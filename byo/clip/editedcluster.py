import logging
import numpy as np
import itertools
import re
from collections import defaultdict

count_functions = {
    'lin' : lambda x : x,
    'uniq' : lambda x : 1,
    'arcsinh' : lambda x : np.arcsinh(x),
    'log' : lambda x : np.log2(x),
}

signature_events = {
    "TC" :  set(["TC","T_"]),
    "GA" :  set(["GA","G_"]),
    "hits" : set(["T_","C_"]),
    "iclip" : set(), # iCLIP is not scored by conversions but by stereotypical 5' end positions.
}

indices = {
    'A' : 0,
    'C' : 1,
    'G' : 2,
    'T' : 3,
    'N' : 0,
}
    
def complement(s):
    comp = { 
        'A' : 'T',
        'C' : 'G',
        'G' : 'C',
        'T' : 'A',
        'N' : 'N',
        '_' : '_' }
    return "".join([comp[x] for x in s])

strand_translate = {'+' : {}, '-' : {}}
for edit in itertools.product('ACGTN_',repeat=2):
    edit = "".join(edit)
    strand_translate['+'][edit] = edit
    strand_translate['-'][edit] = complement(edit)


def requires_read_parsing(f):
    """
    A decorator that wraps functions or properties that 
    require to go through the list of all alignments and 
    parse them.
    """
    def wrapper(self,*args,**kwargs):
        if not self._parsed:
            self.parse_reads()
        return f(self,*args,**kwargs)
    return wrapper
    
from readcluster import ReadCluster
class EditedCluster(ReadCluster):
    """
    Extends the ReadCluster base class. Parses all reads edit and CIGAR strings to
    extract the information of nucleotide matches/mismatches and indels.
    Allows to track signature events (i.e. crosslinking signatures). 
    Builds a matrix representation of the observed nucleotides und depth-of-coverage
    profiles.
    """
    
    def __init__(self,name,reads,strand,chrom,signature={'*' : "TC"},count_func='uniq',**kwargs):
        super(EditedCluster,self).__init__(name,reads,strand,chrom,**kwargs)
        
        self.readcount = count_functions[count_func]
        self.signature = signature

        self.editstat_dict = defaultdict(int)
        self.read_support_dict = defaultdict(int)
        self.sig_support_dict = defaultdict(int)
        
        self.ref = ""
        self.counts = []
        self.n_edits = []
        self.ofs = []
        self._parsed = False

    def parse_reads(self):
        if self._parsed:
            return
        self.nuc_matrix = np.zeros((self.L,4),dtype=np.float32)
        self.signature_counts = np.zeros(self.L,dtype=np.float32)

        signature = {}
        for lib,sig in self.signature.items():
            signature[lib] = signature_events[sig]
        defaultsig = signature.get('*','TC')
        
        # take care of the mismatches
        for pos,edit,w,lib in self.all_edits():
            self.read_support_dict[lib] += w
            edit = strand_translate[self.strand][edit]
            self.editstat_dict[edit] += w

            signature = self.signature.get(lib,defaultsig)
            if edit in signature:
                self.signature_counts[pos] += w
                self.sig_support_dict[lib] += w

            a,b = edit
            if not (a in indices and b in indices):
                # ignore indels for now
                continue

            self.nuc_matrix[pos,indices[a]] -=w
            self.nuc_matrix[pos,indices[b]] +=w
            
        self._parsed = True
        self.counts = np.array(self.counts,dtype=np.float32)

        # now that every read is parsed, we may access self.counts to
        # build the depth-of-coverage profile and collect some statistics
        self.coverage = np.zeros(self.L,dtype=np.float32)
        x0 = self.start
        for start,end,w,n_edits in zip(self.starts,self.ends,self.counts,self.n_edits):
            self.coverage[start-x0:end-x0] += w
            self.editstat_dict['total'] += w
            if n_edits == 0:
                self.editstat_dict['perfect'] += w
            else:
                self.editstat_dict['%d_edit' % n_edits] += w

        # using self.coverage and self.sequence we now add the implicit 
        # nucleotide observations from the matches
        for i,(r,w) in enumerate(zip(self.sequence,self.coverage)):
            self.nuc_matrix[i,indices[r]] += w

        # overwrite conversions with start position counts to make
        # PAR-CLIP code work decently on iCLIP (e.g CCR's)
        if 'iclip' in self.signature.values():
            if self.strand == '-':
                pos = np.array(self.ends,dtype=np.uint32) - self.start - 1
            else:
                pos = np.array(self.starts,dtype=np.uint32) - self.start
               
            self.signature_counts.put(pos,self.counts)
            
    @property
    @requires_read_parsing
    def sequence(self):
        if self.strand == "-":
            return complement(self.ref)
        else:
            return "".join(self.ref)

    @property
    @requires_read_parsing
    def conversions(self):
        return self.signature_counts

    @property
    @requires_read_parsing
    def editstats(self):
        return self.editstat_dict

    @property
    @requires_read_parsing
    def read_support(self):
        return self.read_support_dict

    @property
    @requires_read_parsing
    def sig_support(self):
        return self.sig_support_dict

    def all_edits(self):
        """
        Iterates over all reads in the cluster and yields the edit operations
        (indels, mismatches) at each position in the reference, with the 
        (normalized) multiplicities of the read.
        As a side-effect, this fills the self.counts, self.n_edits, and self.seqs lists,
        which allows to quickly access read multiplicities, total edit distances, and the
        correct reference sequence for subsequent steps.
        """

        self.counts = []
        self.n_edits = []
        self.ref = ["N"] * self.L
        self.ofs = []

        for r in self.reads:
            lib = r.opt('LB')
            # number of edits in this alignment
            n_edits = 0

            # extract mutliplicity of the read and normalize
            # in the desired way (unique only, linear, or arcsinh)
            M = re.search(r"_x(\d+)$",r.qname)
            if M: 
                w = self.readcount(int(M.group(1)))
            else:
                w = 1
            self.counts.append(w)

            # offset within the cluster
            i = r.pos - self.start
            self.ofs.append(i)
            # read sequence as list
            seq = list(r.query)

            # Whew, sometimes I just hate SAM...
            # first we need to go over the CIGAR to get rid of insertions
            j = 0
            cigar_len = 0
            for op,n in r.cigar:
                if op == 0:
                    # match or mismatch
                    j += n
                    cigar_len += n

                elif op == 1:
                    # insertion
                    # remove the inserted nt as the editstr refers 
                    # only to the aligned part of the read
                    n_edits += n
                    for x in range(n):
                        I = seq.pop(j)
                        yield (i+j,"_"+I,w,lib)

                elif op == 2:
                    # deletions are also included in the edit-string
                    # and are otherwise ignored here
                    cigar_len += n

            assert cigar_len == r.alen
            #logging.warning("CIGAR string and read.alen do not match!: '%s'" % str(r))

            # now parse the edit-string
            j = 0
            for ops in re.finditer(r"(\d+)|([A-Z])|(\^[A-Z]+)",r.opt('MD')):
                match,mm,deletion = ops.groups()
                if match:
                    m = int(match)
                    #print "i=%d j=%d m=%d len_counts=%d len_seq=%d" % (i,j,m,len(counts),len(seq))
                    j += m
                elif mm:
                    n_edits += 1
                    #print seq,"MM","i=%d j=%d m=%d len_counts=%d len_seq=%d" % (i,j,m,len(counts),len(seq))
                    yield (i+j,mm+seq[j],w,lib)
                    seq[j] = mm

                    j += 1

                elif deletion:
                    for d in deletion[1:]:
                        yield (i+j,d+"_",w,lib)
                        seq.insert(j,d)
                        j += 1
                        n_edits += 1

            self.n_edits.append(n_edits)
            # store the corrected read sequence as reference sequence
            self.ref[i:i+len(seq)] = seq
