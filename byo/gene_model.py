# -*- coding: utf-8 -*-
import os,sys
from collections import defaultdict,namedtuple
from bisect import bisect_left,bisect_right
from numpy import array,transpose

from byo.io.gff import gff_importer
from byo.protein import all_orfs,translate
from byo import protein

#from byo.transcript import ProcessedTranscript
#from byo.locus import Locus

import logging
from logging import debug,warning,info,error

class ExonChain(object):
    """
    Implements linked blocks of genomic coordinates (with orientation).
    Features splicing from arbitrary tracks (including sequence) and intersection with
    simple start/end blocks on the genome. The latter returns 3 new ExonChains. As an
    example consider intersection with genomic CDS start and end-coordinates. This will
    return 5'UTR, CDS and 3'UTR as individual ExonChains.
    This class is deliberately kept as light-weight as possible. It has no means of
    generating names/annotation. All of that is handled by the "Transcript" subclass.
    """

    def __init__(self,chrom,sense,exon_starts,exon_ends,system=None):
        if not system:
            from byo.systems import hg18
            system = hg18

        self.system = system
        self.chrom = chrom
        self.sense = sense
        self.exon_starts = array(sorted(exon_starts))
        self.exon_ends = array(sorted(exon_ends))

        self.exon_count = len(self.exon_starts)
        self.start = self.exon_starts[0]
        self.end = self.exon_ends[-1]

        # maybe turn these into properties as well?
        self.exon_bounds = array([self.exon_starts,self.exon_ends]).transpose()
        self.intron_bounds = []
        if self.exon_count > 1:
            self.intron_bounds = array([self.exon_ends[:-1],self.exon_starts[1:]]).transpose()

        self.exon_lengths = self.exon_ends - self.exon_starts
        self.exon_txstarts = array([0,]+ list(self.exon_lengths.cumsum()))

        self.spliced_length = max(self.exon_txstarts[-1],0)
        self.unspliced_length = self.end - self.start

        if sense == "-":
            self.dir = -1
            self.ofs = self.spliced_length #- 1 #oh,oh... :-/
        else:
            self.dir = +1
            self.ofs = 0

        self.tx_start,self.tx_end = [self.start,self.end][::self.dir]
        self.name = "ExonChain (%s:%d-%d%s) %s %s" % (chrom,self.start,self.end,sense,str(exon_starts),str(exon_ends))

    @property
    def intron_chain(self):
        if self.exon_count > 1:
            return ExonChain(self.chrom,self.sense,self.exon_ends[:-1],self.exon_starts[1:],system=self.system)

    def map_to_spliced(self,pos):
        n = bisect_right(self.exon_starts,pos) - 1
        #print "n_left",n
        #print "n_right",bisect_left(self.exon_ends,pos)
        assert (n == bisect_left(self.exon_ends,pos))
        return self.ofs + self.dir * (pos - self.exon_starts[n] + self.exon_txstarts[n])

    def map_block_from_spliced(self,start,end):
        x,y = self.map_from_spliced(start),self.map_from_spliced(end)
        if y > x:
            # semantics change here because start is included, end is excluded
            # in C indexing.
            return y-1,x-1
        else:
            return x,y
        #return min(x,y),max(x,y)

    def map_block_to_spliced(self,start,end):
        x,y = self.map_to_spliced(start),self.map_to_spliced(end)
        if y > x:
            # semantics change here because start is included, end is excluded
            # in C indexing.
            return y-1,x-1
        else:
            return x,y

    def map_from_spliced(self,pos):
        #print pos,self.sense,self.exon_count,self.ofs,self.dir,self.spliced_length
        pos = self.ofs + self.dir * pos
        #assert (0 <= pos < self.spliced_length)
        n = bisect_left(self.exon_txstarts,pos) -1
        #print n,len(self.exon_starts),len(self.exon_txstarts)
        return self.exon_starts[n] + pos - self.exon_txstarts[n]

    def splice(self,track,join = lambda l : "".join(l),get="get_oriented",**kwargs):
        get_func = getattr(track,get)
        return join([get_func(self.chrom,start,end,self.sense,**kwargs) for start,end in self.exon_bounds][::self.dir])

    @property
    def unspliced_sequence(self):
        return self.system.genome.get_oriented(self.chrom,self.start,self.end,self.sense)

    @property
    def spliced_sequence(self):
        return self.splice(self.system.genome)

    @property
    def splice_sites(self):
        """returns tuples of (5'ss, 3'ss)"""
        for bounds in self.intron_bounds[::self.dir]:
            yield bounds[::self.dir]

    @property
    def exons(self):
        for i,(start,end) in enumerate(self.exon_bounds[::self.dir]):
            yield Locus("%s.exon%02d" % (self.name,i+1),self.chrom,self.sense,start,end,"exon")

    @property
    def introns(self):
        for i,(start,end) in enumerate(self.intron_bounds[::self.dir]):
            yield Locus("%s.intron%02d" % (self.name,i+1),self.chrom,self.sense,start,end,"intron")

    @property
    def introns_as_chains(self):
        for i,(start,end) in enumerate(self.intron_bounds[::self.dir]):
            yield "intron", ExonChain(self.chrom,self.sense,[start,],[end,])

    def intersect(self,names,start,end):
        # TODO: Clean this up
        start,end = min(start,end),max(start,end)
        start = max(start,self.start)
        end = min(end,self.end)

        first = bisect_right(self.exon_starts,start) - 1
        last = bisect_right(self.exon_starts,end) -1
        
        f2 = bisect_right(self.exon_ends,start)
        l2 = bisect_right(self.exon_ends,end)
        
        #print first,f2
        #print last,l2

        before_starts = list(self.exon_starts[:f2+1])
        before_ends = list(self.exon_ends[:f2+1])
        
        if not before_starts and not before_ends:
            before_starts = [0]
            before_ends = [0]

        chain_starts = list(self.exon_starts[f2:last+1])
        chain_ends = list(self.exon_ends[f2:last+1])

        if not chain_starts and not chain_ends:
            chain_starts = [0]
            chain_ends = [0]

        #if not chain_starts:
            #print first,last,start,end,self.tx_start,self.tx_end

        after_starts = list(self.exon_starts[last:])
        after_ends = list(self.exon_ends[last:])

        if not after_starts and not after_ends:
            after_starts = [0]
            after_ends = [0]
        
        if first == f2:
            chain_starts[0] = max(start,chain_starts[0])
            before_ends[-1] = min(start,before_ends[-1])
        else:
            # truncate chain, breakpoint between two exons
            before_ends[-1] = before_starts[-1]
            
        if last == l2:
            chain_ends[-1] = min(end,chain_ends[-1])
            after_starts[0] = max(end,after_starts[0])
        else:
            # truncate chain, breakpoint between two exons
            after_starts[0] = after_ends[0]

        #print self.system
        chain = ExonChain(self.chrom,self.sense,chain_starts,chain_ends,system=self.system)
        before = ExonChain(self.chrom,self.sense,before_starts,before_ends,system=self.system)
        after = ExonChain(self.chrom,self.sense,after_starts,after_ends,system=self.system)

        if self.sense == "+":
            return before,chain,after
        else:
            return after,chain,before

    # add some python magic to make things smooth
    def __str__(self):
        exonlist = ",".join(map(str,self.exon_bounds))
        return "%s %s %d %d exons: %s" % (self.chrom,self.sense,self.start,self.end,exonlist)

    def __len__(self):
        """
        zero-length ExonChains will be False in truth value testing.
        so stuff like: "if 5UTR" can work.
        """
        return self.spliced_length

    def __hasitem__(self,pos):
        """check if the coordinate falls into an exon"""
        try:
            x = self.map_to_spliced(pos)
        except AssertionError:
            return False
        else:
            return True

Locus = namedtuple("Locus","name, chrom, sense, start, end, category")

class Transcript(ExonChain):
    def __init__(self,name,chrom,sense,exon_starts,exon_ends,cds,score=0,system=None,description={}):
        # Initialize underlying ExonChain
        super(Transcript,self).__init__(chrom,sense,exon_starts,exon_ends,system=system)
        self.score = score
        self.category = "transcript"
        self.name = name
        self.description = description

        cds_start,cds_end = cds
        if cds_start == cds_end:
            self.UTR5 = None
            self.UTR3 = None
            self.CDS = None
            self.coding = False
            self.tx_class = description.get("tx_class","ncRNA")
        else:
            self.UTR5,self.CDS,self.UTR3 = self.intersect(["UTR5","CDS","UTR3"],cds_start,cds_end)
            self.coding = True
            self.tx_class = description.get("tx_class","coding")
            if self.UTR5: 
                self.UTR5.name = "%s.UTR5" % self.name
                self.UTR5.gene_id = self.gene_id
            if self.CDS: 
                self.CDS.name = "%s.CDS" % self.name
                self.CDS.gene_id = self.gene_id
            if self.UTR3: 
                self.UTR3.name = "%s.UTR3" % self.name
                self.UTR3.gene_id = self.gene_id
        # immitate old-style Transcript API for build_annotation compatibility

    #def __getattr__(self,attr):
        #if attr.

    @property
    def segments(self):
        if self.UTR5:
            UTR = self.UTR5
            yield Locus("%s.5UTR" % self.name,UTR.chrom,UTR.sense,UTR.start,UTR.end,"5UTR")

        if self.UTR3:
            UTR = self.UTR3
            yield Locus("%s.3UTR" % self.name,UTR.chrom,UTR.sense,UTR.start,UTR.end,"3UTR")

        if self.CDS:
            CDS = self.CDS
            yield Locus("%s.CDS" % self.name,CDS.chrom,CDS.sense,CDS.start,CDS.end,"CDS")

        for E in self.exons:
            yield E

        for I in self.introns:
            yield I

    @property
    def features(self):
        start,end = [self.start,self.end][::self.dir]
        yield Locus("%s.TSS" % self.name,self.chrom,self.sense,start,start+1,"TSS")
        yield Locus("%s.PAS" % self.name,self.chrom,self.sense,end,end+1,"PAS")

        for i,(ss5,ss3) in enumerate(self.splice_sites):
            yield Locus("%s.intron%02d.5'SS" % (self.name,i+1),self.chrom,self.sense,ss5,ss5+1,"5'SS")
            yield Locus("%s.intron%02d.3'SS" % (self.name,i+1),self.chrom,self.sense,ss3,ss3+1,"3'SS")

        if self.CDS:
            cds_start,cds_end = [self.CDS.start,self.CDS.end][::self.dir]
            yield Locus("%s.CDSstart" % self.name,self.chrom,self.sense,cds_start,cds_start+1,"CDSstart")
            yield Locus("%s.CDSend" % self.name,self.chrom,self.sense,cds_end,cds_end+1,"CDSend")

    @property
    def gene_id(self):
        gene_id = getattr(self,"name2",None)
        if not gene_id:
            gene_id = self.description.get("gene_id",self.description.get("name2",self.name))
        if not gene_id:
            return self.name
        return gene_id

    @property
    def ORF(self):
        if not self.CDS:
            return ""
        else:
            return protein.translate_aa(self.CDS.spliced_sequence)

    # UCSC table format output
    def __str__(self):
        exonstarts = ",".join([str(s) for s in self.exon_starts])
        exonends = ",".join([str(e) for e in self.exon_ends])
        #exonframes = ",".join([str(f) for f in self.exon_frames])
        # TODO: fix exon-frames
        exonframes = ",".join(["-1" for e in self.exon_ends])

        if self.CDS:
            cds_start = self.CDS.start
            cds_end = self.CDS.end
        else:
            cds_start = self.end
            cds_end = self.end
            
        out = (-1,self.name,self.chrom,self.sense,self.start,self.end,cds_start,cds_end,self.exon_count,exonstarts,exonends,str(self.score),self.gene_id,"unk","unk",exonframes)
        return "\t".join([str(o) for o in out])

ucsc_table_format = "## bin:int \t name:str \t chrom:str \t sense:str \t txStart:int \t txEnd:int \t cdsStart:int \t cdsEnd:int \t exonCount:int \t exonStarts:intlist \t exonEnds:intlist \t id:str \t name2:str \t cdsStartStat:str \t cdsEndStat:str \t exonFrames:intlist"
        
# factory functions to generate transcript models from downloaded gene-models
def transcripts_from_UCSC(fname,system=None,tx_class = None,gene_names = {},fix_chr=True,table_format=ucsc_table_format,**kwargs):
    from byo.io.lazytables import LazyImporter,NamedTupleImporter

    kw = {}
    if tx_class:
        kw = {'tx_class' : tx_class}

    def factory(name="",chrom="",sense="",exonStarts=[],exonEnds=[],cdsStart=-1,cdsEnd=-1,**kwargs):
        k = dict(kw)
        k.update(kwargs)

        if not chrom.startswith('chr') and fix_chr:
            chrom = 'chr'+chrom
        
        # splice in gene ids if given in a separate file (ce6)
        k['name2'] = gene_names.get(name,k.get('name2',k.get('proteinID',name)))
        #print k['name2'],name
        return Transcript(name,chrom,sense,exonStarts,exonEnds,(cdsStart,cdsEnd),description=k,system=system)

    
    #for tx in LazyImporter(path,factory,"name",skip=[0,11,13,14],descr=ucsc_table_format):
        #yield tx

    ucsc_type_hints = dict(txStart="int",txEnd="int",cdsStart="int",cdsEnd="int",exonCount="int",exonStarts="intlist",exonEnds="intlist",exonFrames="intlist")
    return LazyImporter(fname,factory,"name",type_hints=ucsc_type_hints,renames=dict(strand="sense"),parse_comments=True,descr=table_format,**kwargs)
    #return LazyImporter(fname,factory,"name",skip=[0,11,13,14],descr=ucsc_table_format)

def transcripts_from_GTF(fname="/data/BIO2/mjens/HuR/HeLa/RNASeq/cuff/t",ORF_thresh=20):

    for i,args in enumerate(get_exons_for_transcript()):
        if not (i % 1000):
            print i

        tx = ExonChain(*args)
        tx.UTR5 = None
        tx.CDS = None
        tx.UTR3 = None
        #print tx.name
        if ORF_thresh > 0:
            orfs = sorted([ (len(orf),orf,start,end) for (start,end,orf) in all_orfs(tx.spliced_sequence) if len(orf) > ORF_thresh],reverse=True)
            if orfs:
                """
                largest open-reading-frame above threshold becomes CDS
                """
                length,orf,start,end = orfs[0]
                try:
                    cds_start,cds_end = tx.map_block_from_spliced(start,end)
                except AssertionError:
                    print start,end,end-start,tx.spliced_length,tx.tx_start,tx.tx_end,tx.exon_count,tx.sense
                    print tx.map_from_spliced(start)
                    print tx.map_from_spliced(end)

                #cds_start = tx.map_from_spliced(start)
                #cds_end = tx.map_from_spliced(end)

                seg_names = ["UTR5","CDS","UTR3"]
                tx.intersect(seg_names,row.cdsStart,row.cdsEnd,add_as_segments=True)
                #print tx.name,len(orfs),orf

        yield tx

if __name__ == "__main__":
    #for tx in transcripts_from_UCSC():
    for tx in transcripts_from_GTF("/data/BIO2/mjens/HuR/HeLa/RNASeq/cuff/stdout.combined.gtf"):
        #print tx.exon_count
        pass

    #from pprint import pprint
    #E = ExonChain("test","chr1","-",[10,100,300],[30,150,400])
    #pprint( list(E.exons) )
    #pprint( list(E.introns) )

    #CDS = Locus("test.CDS","chr1","-",10,395)
    #print E.spliced_sequence

    #before,chain,after = E.intersect(CDS)
    #print "before",before
    #print "chain",chain
    #print "after",after

    #for tpos in [20,110,300,5,40]:
        #s = E.map_to_spliced(tpos)
        #print "mapped",tpos,s
        #print "inverse",s,E.map_from_spliced(s)

    sys.exit(0)



