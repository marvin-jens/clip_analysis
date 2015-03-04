import numpy as np
from collections import defaultdict

### Helper functions
def entropy_of_pos(pos,counts):
    data = defaultdict(int)
    for x,c in zip(pos,counts):
        data[x] += c

    p = np.array(data.values(),dtype=np.float32)
    p /= p.sum()

    return - (p * np.log2(p)).sum()

from editedcluster import EditedCluster, requires_read_parsing

def requires_map_parsing(f):
    """
    A decorator that wraps functions or properties that 
    require to go through the list of all alignments and 
    parse their various mapping quality metrics.
    """
    def wrapper(self,*args,**kwargs):
        if not self._map_scored:
            self.parse_map_scores()
        return f(self,*args,**kwargs)
    return wrapper

class ScoredCluster(EditedCluster):

    def __init__(self,name,reads,strand,chrom,**kwargs):
        super(ScoredCluster,self).__init__(name,reads,strand,chrom,**kwargs)

        self.mapq = []
        self.uniqness = []
        self.pp = []
        self.reported_scores = [
            "uniq_reads","n_reads","n_conv_pos","length",
            "max_uniq","avg_uniq","cum_uniq",
            "max_mapq","avg_mapq","cum_mapq",
            "max_pp","avg_pp","cum_pp",
            "conversion_score","signature_density","iclip_score","profile_KL","seq_entropy","entropy_score","edits"
        ]
        self._map_scored = False

        # current behaviour is to switch to iCLIP 
        # as soon as one library is iCLIP
        if 'iclip' in self.signature.values():
            self._scorefuncname = "iclip_score"
        else:
            self._scorefuncname = "conversion_score"

        from byo.clip.statistics import edit_map
        self.edit_keys = [
            # mismatches, deletions
            'AC', 'AG', 'AT', 'A_', 
            'CA', 'CG', 'CT', 'C_',
            'GA', 'GC', 'GT', 'G_',
            'TA', 'TC', 'TG', 'T_',
            '_A', '_C', '_G', '_T'  # insertions
        ]

    @requires_read_parsing
    def parse_map_scores(self):
        for r,w in zip(self.reads,self.counts):
            # record MAPQ mapping quality (0..255)
            self.mapq.append(r.mapq)

            # look at additional mapper output for cues on alignment quality
            tags = dict(r.tags)
            uniqness = 0

            # BWA output
            # get number of best alignments and the hull 
            # (number of alignments at k+1)
            n_best = tags.get('X0',1)
            if 'X1' in tags:
                hull = tags['X1']
                uniqness = (n_best == 1) * 255. / (hull + 1.)
            
            # BOWTIE2 output
            # get the difference between best alignment score and second best
            elif 'AS' in tags:
                AS = tags.get('AS',0)
                XS = tags.get('XS',0)
                uniqness = abs(AS-XS)
            
            self.uniqness.append(uniqness)
            
            # BWA-PSSM output
            # take the posterior probability that this alignment is correct,
            # given the error-model (for example PAR-CLIP conversions)
            pp = tags.get('PP',0.5)
            self.pp.append(pp)
            
        self.uniqness = np.array(self.uniqness)
        self.mapq = np.array(self.mapq)
        self.pp = np.array(self.pp)
        
        self._map_scored = True

    @property
    @requires_read_parsing
    def score(self):
        return getattr(self,self._scorefuncname)

    @property
    @requires_read_parsing
    def conversion_score(self):
        return self.conversions.sum()

    @property
    @requires_read_parsing
    def signature_density(self):
        return self.conversion_score / self.n_reads

    @property
    @requires_read_parsing
    def entropy_5p(self):
        if self.strand == '-':
            return entropy_of_pos(self.ends,self.counts)
        else:
            return entropy_of_pos(self.starts,self.counts)

    @property
    @requires_read_parsing
    def entropy_3p(self):
        if self.strand == '-':
            return entropy_of_pos(self.starts,self.counts)
        else:
            return entropy_of_pos(self.ends,self.counts)

    @property
    @requires_read_parsing
    def iclip_score(self):
        # 1. for a single read, goes to zero the more reads we have
        default = 2. / (self.n_reads + 1.)
        
        # low 5' entropy, even for high numbers of reads is good.
        # high 3' entropy on the other hand is good. The proposed score
        # is 1 if both entropies are the same and goes to 0 if 3' end
        # variability is lower than 5' end, and goes to infinity for
        # 5' entropy close to zero, as read-counts increase.
        return (self.entropy_3p + default) / (self.entropy_5p + default)

    @property
    @requires_read_parsing
    def n_reads(self):
        return self.counts.sum()

    @property
    def uniq_reads(self):
        return len(self.reads)
    
    @property
    @requires_map_parsing
    def cum_uniq(self):
        return (self.counts * self.uniqness).sum()

    @property
    @requires_map_parsing
    def max_uniq(self):
        return self.uniqness.max()

    @property
    @requires_map_parsing
    def avg_uniq(self):
        return (self.counts * self.uniqness).sum() / self.counts.sum()

    @property
    @requires_map_parsing
    def cum_mapq(self):
        return (self.counts * self.mapq).sum()

    @property
    @requires_map_parsing
    def max_mapq(self):
        return self.mapq.max()

    @property
    @requires_map_parsing
    def avg_mapq(self):
        return (self.counts * self.mapq).sum() / self.counts.sum()

    @property
    @requires_map_parsing
    def cum_pp(self):
        return (self.counts * self.pp).sum()

    @property
    @requires_map_parsing
    def max_pp(self):
        return self.pp.max()

    @property
    @requires_map_parsing
    def avg_pp(self):
        return (self.counts * self.pp).sum() / self.counts.sum()

    @property
    @requires_read_parsing
    def profile_KL(self):
        """
        KL-divergence between observed profile and uniform coverage. 
        """
        p = self.coverage / float(self.coverage.sum())
        return np.where(p> 0,p * np.log2(self.L*p),0).sum()

    @property
    @requires_read_parsing
    def seq_entropy(self):
        """
        Shannon entropy of the nucleotide observations (0 for no mismatches)
        """
        p = np.array(self.nuc_matrix,dtype=np.float32) / self.nuc_matrix.sum(axis=1)[:,np.newaxis]
        return np.where(p> 0,-p * np.log2(p),0).sum()

    @property
    def entropy_score(self):
        scale = (1 - (10./(10.+self.n_reads)))*.1
        return self.profile_KL + scale*self.seq_entropy
        
    @property
    def short(self):
        return "%s%s %d:%d " % (self.chrom,self.strand,self.start,self.end)

    @property
    @requires_read_parsing
    def n_conv_pos(self):
        return (self.conversions > 0).sum()
    
    @property
    def edits(self):
        #from byo.clip.statistics import edit_map
        #if self.strand == "+":
            #return ",".join([str(self.editstats[k]) for k in self.edit_keys])
        #else:
            #return ",".join([str(self.editstats[edit_map[k]]) for k in self.edit_keys])

        return ",".join([str(self.editstats[k]) for k in self.edit_keys])
