from byo.track import Accessor
from logging import debug,warning,error
import re
import numpy


def is_uniq(r):
    try:
        best_hits = r.opt('X0')
        if best_hits < 2:
            return True
    except KeyError:
        # bowtie2 output
        AS = r.opt('AS')
        try:
            XS = int(r.opt('XS'))
        except (KeyError,ValueError):
            XS = 2000

        # if alternative alignment and reported (best) alignment have different score,
        # the mapping is unique
        return abs(XS-AS) != 0
  
class BAMAccessor(Accessor):
    def __init__(self,path,chrom,sense,sense_specific=False,dtype=numpy.uint32,mode="c",filter=None,dim=1,mult=True,unique=True,system='hg19',read_count=False,only_starts=False,**kwargs):
        debug("# BAMAccessor: Loading '%s' coverage for chromosome %s%s from '%s' (sense_specific=%s)" % (str(dtype),chrom,sense,path,sense_specific))
        super(BAMAccessor,self).__init__(path,chrom,sense,sense_specific=sense_specific,system=system,**kwargs)
        self.dtype=dtype
        self.sense_specific = sense_specific
        self.dim = dim
        self.only_starts = only_starts
        self.chrom = chrom
        self.path = path
        self.mult = mult
        self.unique = unique
        self.filter = filter
        import byo.systems
        self.system = getattr(byo.systems,system)

        import pysam
        try:
            self.bam = pysam.Samfile(path,"rb")
        except IOError:
            warning("Could not access '%s'. Switching to dummy mode (only zeros)" % path)
            self.get_data = self.get_dummy
            self.get_sum = self.get_sum_dummy

        if not chrom in self.bam.references:
            warning("chromosome '%s' not in BAM '%s'. Switching to dummy mode (only zeros)" % (chrom,path))
            self.get_data = self.get_dummy
            self.get_sum = self.get_sum_dummy
            
        if read_count:
            self.get_sum = self.get_read_count

        # register for all chroms/strands
        self.covered_strands = "*"
        #self.covered_strands = [chrom+'+' for chrom in self.system.chr_sizes.keys()] + [chrom+'-' for chrom in self.system.chr_sizes.keys()] + [chrom+sense,]
            
        
    def get_sum(self,chrom,start,end,sense,**kwargs):
        #debug("get_sum")
        #try:
            #return self.bam.count(self.chrom,start,end)
        #except ValueError:
            #return 0

        profile = self.get_data(chrom,start,end,sense)
        #print profile
        return profile.sum()
            
    def get_sum_dummy(self,chrom,start,end,sense,**kwargs):
        #debug("get_sum_dummy")
        #try:
        return 0
        #except ValueError:
            #return 0
            
    def get_read_count(self,chrom,start,end,sense):
        
        #print "get_read_count"
        def sense_of_read(read):
            if read.is_reverse:
                return "-"
            else:
                return "+"

        count = 0
        for read in self.bam.fetch(chrom,max(0,start),end):
            if self.unique:
                if not is_uniq(read):
                    continue
            if self.sense_specific and (sense_of_read(read) != sense): continue
            if self.filter and not self.filter(read):
                continue

            w = 1
            if self.mult:
                M = re.search(r"_x(\d+)$",read.qname)
                if M:
                    w = int(M.group(1))

            count += w

        return count
        
    def get_data(self,chrom,start,end,sense,**kwargs):
        #debug("get_data %s:%d-%d" % (chrom,start,end))
        buf = self.get_dummy(chrom,start,end,sense)

        def sense_of_read(read):
            if read.is_reverse:
                return "-"
            else:
                return "+"

        #print "<%s>.get_data(%s,%d,%d,%s)" % (self.path,self.chrom,start,end,sense)
        for read in self.bam.fetch(chrom,max(0,start),end):
            #print read,sense_of_read(read)
            if self.unique:
                if not is_uniq(read):
                    continue
            if self.sense_specific and (sense_of_read(read) != sense): continue
            if self.filter and not self.filter(read):
                continue

            w = 1
            if self.mult:
                M = re.search(r"_x(\d+)$",read.qname)
                if M:
                    w = int(M.group(1))

            s = max(0,read.pos-start)
            e = min(end-start,read.pos-start+read.alen)

            if self.only_starts:
                if sense_of_read(read) == '+':
                    buf[s] += w
                else:
                    buf[e-1] += w
            else:
                buf[s:e] += w

        return buf
        
    def get_dummy(self,chrom,start,end,sense,**kwargs):
        return numpy.zeros(end-start,dtype=self.dtype)
