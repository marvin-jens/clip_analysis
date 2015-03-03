from byo.track import Accessor
from byo import complement, rev_comp
from logging import debug,warning,error

import os
import mmap
import re

class mmap_fasta(object):
    def __init__(self,fname):
        f = file(fname)
        header = f.readline()
        row = f.readline()

        self.ofs = len(header)
        self.lline = len(row)
        self.ldata = len(row.strip())
        self.skip = self.lline-self.ldata
        self.skip_char = row[self.ldata:]
        #print "SKIP",self.skip,self.skip_char
        self.mmap = mmap.mmap(f.fileno(),0,prot=mmap.PROT_READ)

    def __getslice__(self,start,end):
        l_start = start / self.ldata
        l_end = end / self.ldata
        #print "lines",l_start,l_end
        ofs_start = l_start * self.skip + start + self.ofs
        ofs_end = l_end * self.skip + end + self.ofs
        #print "ofs",ofs_start,ofs_end
        
        s = self.mmap[ofs_start:ofs_end].replace(self.skip_char,"")
        L = end-start
        if len(s) == L:
            return s
        else:
            return s+"N"*(L-len(s))
        return 

class indexed_fasta(object):
    def __init__(self,fname):
        self.fname = fname
        self.chrom_stats = {}
        
        # try to load index
        ipath = fname + '.index'
        if os.access(ipath,os.R_OK):
            self.load_index(ipath)
        else:
            self.index()
            self.store_index(ipath)
                
        f = file(fname)
        self.mmap = mmap.mmap(f.fileno(),0,prot=mmap.PROT_READ)

    def index(self):
        debug("# indexed_fasta.index('%s')" % self.fname)

        ofs = 0
        f = file(self.fname)
        chrom = "undef"
        chrom_ofs = 0
        
        for line in f:
            ofs += len(line)
            if line.startswith('>'):
                chrom = line[1:].split()[0].strip()
                chrom_ofs = ofs
            else:
                if not chrom in self.chrom_stats:
                    lline = len(line)
                    ldata = len(line.strip())
                    self.chrom_stats[chrom] = (chrom_ofs,ldata,lline-ldata,line[ldata:])
        f.close()

    def store_index(self,ipath):
        debug("# indexed_fasta.store_index('%s')" % ipath)
        
        # write to tmp-file first and in the end rename in order to have this atomic 
        # otherwise parallel building of the same index may screw it up.
        
        import tempfile
        tmp = tempfile.NamedTemporaryFile(mode="w",dir = os.path.dirname(ipath),delete=False)
        for chrom in sorted(self.chrom_stats.keys()):
            ofs,ldata,skip,skipchar = self.chrom_stats[chrom]
            tmp.write("%s\t%d\t%d\t%d\t%r\n" % (chrom,ofs,ldata,skip,skipchar))
        
        # make sure everything is on disk
        os.fsync(tmp)
        tmp.close()
        
        # make it accessible to everyone
        import stat
        os.chmod(tmp.name, stat.S_IROTH | stat.S_IRGRP | stat.S_IRUSR)
        
        # this is atomic on POSIX as we have created tmp in the same directory, 
        # therefore same filesystem
        os.rename(tmp.name,ipath)
        

    def load_index(self,ipath):
        debug("# indexed_fasta.load_index('%s')" % ipath)
        self.chrom_stats = {}
        for line in file(ipath):
            chrom,ofs,ldata,skip,skipchar = line.rstrip().split('\t')
            self.chrom_stats[chrom] = (int(ofs),int(ldata),int(skip),skipchar[1:-1].decode('string_escape'))
        
    def get_data(self,chrom,start,end,sense):
        if not self.chrom_stats:
            self.index()

        ofs,ldata,skip,skip_char = self.chrom_stats[chrom]
        #print "ldata",ldata
        #print "chromstats",self.chrom_stats[chrom]

        l_start = start / ldata
        l_end = end / ldata
        #print "lines",l_start,l_end
        ofs_start = l_start * skip + start + ofs
        ofs_end = l_end * skip + end + ofs
        #print "ofs",ofs_start,ofs_end,ofs_end - ofs_start
        
        s = self.mmap[ofs_start:ofs_end].replace(skip_char,"")
        L = end-start
        if len(s) < L:
            s+="N"*(L-len(s))
        if sense == "-":
            s = rev_comp(s)
        return s

        
        
if __name__ == "__main__":
    #i = indexed_fasta('/data/BIO2/pcp/systems/hg19/genome/hg19.fna')
    i = indexed_fasta('/data/BIO2/pcp/systems/ce6/genome/ce6.fna')
    #print i.chrom_stats
    print i.get_data('chrIV',0,100,'+')
    #print i.chrom_stats
            
        
        
        
class GenomeAccessor(Accessor):
    def __init__(self,path,chrom,sense,system='hg19',**kwargs):
        #import logging
        #logging.basicConfig(level=logging.DEBUG)

        super(GenomeAccessor,self).__init__(path,chrom,sense,system=system,**kwargs)
        debug("# GenomeAccessor mmap: Loading genomic sequence for chromosome %s from '%s'" % (chrom,path))

        self.system = system
        self.data = None
        import byo.systems
        self.chr_sizes = getattr(byo.systems,system).chr_sizes
        
        # try to access the whole genome, using indexing for fast lookup
        trials = [os.path.join(path,system+'.fna'),os.path.join(path,system+'.fa'),os.path.join(path,chrom+".fa")]
        for fname in trials:
            debug("trying to load '%s'" % fname)
            if os.access(fname,os.R_OK):
                self.data = indexed_fasta(fname)
                break
                
        if not self.data:
            # all fails: return Ns only
            warning("Could not access any of '%s'. Switching to dummy mode (only Ns)" % str(trials))
            self.get_data = self.get_dummy
            self.get_oriented = self.get_dummy
            self.covered_strands = [chrom+'+',chrom+'-']
        else:
            # register for all chroms/strands
            self.covered_strands = [chrom+'+' for chrom in self.data.chrom_stats.keys()] + [chrom+'-' for chrom in self.data.chrom_stats.keys()]

        # TODO: maybe remove this if not needed
        self.get = self.get_oriented

    def load_indexed(self,path):
        ipath = path+'.index'
        if not os.access(ipath,os.R_OK):
            index = self.build_index(path,ipath)
        else:
            index = self.load_index(ipath)

        self.chrom_ofs = index.chrom_ofs
        
    def get_data(self,chrom,start,end,sense):
        #print "indexed_fasta.get_data('%s',%d,%d,%s)" % (chrom,start,end,sense)
        L = self.chr_sizes[chrom]
        if start < 0 or end > L:
            # query violates chromosome bounds. Try to handle gracefully by padding with Ns
            #print "flaky query: start=%d end=%d L=%d" % (start,end,L)
            pad_left = -min(start,0)
            pad_right = -min(L-end,0)
            #print "pad_left=%d pad_right=%d" % (pad_left,pad_right)
            insert = self.data.get_data(chrom,max(0,start),min(end,L),"+")
            #print len(insert),max(0,start),min(end,L),min(end,L)-max(0,start)
            seq = ("N"*pad_left) + insert + ("N"*pad_right)
            #print seq
        else:
            seq = self.data.get_data(chrom,start,end,"+")
            
        if sense == "-":
            seq = complement(seq)

        return seq

    #def get_oriented(self,chrom,start,end,sense):
        #if end < 0:
            #return self.get_dummy(chrom,start,end,sense)
        #elif start < 0:
            #return self.get_dummy(chrom,start,0,sense) + self.get_oriented(chrom,0,end,sense)0
        #else:
            #return self.data.get_data(chrom,start,end,sense)

    def get_dummy(self,chrom,start,end,sense):
        return "N"*int(end-start)



class MSFAccessor(Accessor):
    def __init__(self,path,chrom,sense,system='hg19',offset=0,**kwargs):
        #import logging
        #logging.basicConfig(level=logging.DEBUG)

        super(MSFAccessor,self).__init__(path,chrom,sense,system=system,**kwargs)
        debug("# MSFAccessor mmap: Loading aligned, stitched genomic sequences for chromosome %s from '%s' offset=%d" % (chrom,path,offset))

        self.system = system
        self.data = None
        self.offset = offset
        
        # try to access the whole genome, using indexing for fast lookup
        trials = [os.path.join(path,chrom+'.maf.stitched.cmpl.repeats_lc')]
        for fname in trials:
            debug("trying to load '%s'" % fname)
            if os.access(fname,os.R_OK):
                self.data = indexed_fasta(fname)
                self.species = sorted(self.data.chrom_stats.keys())
                self.extra_data = dict(species=self.species)
                self.covered_strands = [chrom+'+',chrom+'-']
                
                # mangle the "chrom" names to make sure 
                # reference sequences like 'hg19.chr5(+):1-12345' are found as 'hg19'
                for k,v in self.data.chrom_stats.items():
                    parts = k.split(".chr")
                    if len(parts) > 1:
                        self.data.chrom_stats[parts[0]] = v
                break
                
        if not self.data:
            # all fails: return Ns only
            warning("Could not access any of '%s'. Switching to dummy mode (only Ns)" % str(trials))
            self.get_data = self.get_dummy
            self.get_oriented = self.get_dummy
            self.covered_strands = [chrom+'+',chrom+'-']
            self.species = []
        else:
            # register for all chroms/strands
            self.covered_strands = [chrom+'+',chrom+'-']

        # TODO: maybe remove this if not needed
        self.get = self.get_oriented

    def load_indexed(self,path):
        ipath = path+'.index'
        if not os.access(ipath,os.R_OK):
            index = self.build_index(path,ipath)
        else:
            index = self.load_index(ipath)

        self.chrom_ofs = index.chrom_ofs
        
    def get_data(self,chrom,start,end,sense,species=[]):
        start += self.offset
        end += self.offset

        if not species:
            species = self.species
        if start < 0 or end < 0:
            return [self.get_dummy(spc,start,end,sense) for spc in species]
        #UCSC convention: start with 1, end is inclusive
        if sense == "+":
            return [self.data.get_data(spc,start,end,sense).replace('-','N') for spc in species]
        else:
            return [complement(self.data.get_data(spc,start,end,'+').replace('-','N')) for spc in species]

    def get_oriented(self,chrom,start,end,sense,species=[]):
        start += self.offset
        end += self.offset

        if not species:
            species = self.species
        if start < 0 or end < 0:
            return [self.get_dummy(spc,start,end,sense) for spc in species]
        #UCSC convention: start with 1, end is inclusive
        if sense == "+":
            return [self.data.get_data(spc,start,end,sense).replace('-','N') for spc in species]
        else:
            return [rev_comp(self.data.get_data(spc,start,end,'+').replace('-','N')) for spc in species]

    #def get_oriented(self,chrom,start,end,sense,species=[]):
        #if not species:
            #species = self.species

        #if end < 0:
            #return self.get_dummy(chrom,start,end,sense)
        #elif start < 0:
            #return self.get_dummy(chrom,start,end,sense)
            ##NEEDS FIXING FOR LIST OF SEQS
            ##[self.get_dummy(chrom,start,0,sense) + self.get_oriented(spc,0,end,sense) for spc in species]
        #else:
            #return self.get_data(chrom,start,end,sense,species=species)

    def get_dummy(self,chrom,start,end,sense,species=[],**kwargs):
        if not species:
            species = self.species
        
        return ["N"*int(end-start) for spc in species]
