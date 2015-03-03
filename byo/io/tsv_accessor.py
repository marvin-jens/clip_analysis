from byo.track import Accessor
import numpy
import os
from time import time
from logging import debug,warning,error

class SparseFile(Accessor):
    def __init__(self,path,chrom,sense,strand_specific=False,dtype=numpy.int32,default=0,**kwargs):
        self.path = path
        self.sense = sense
        self.chrom = chrom
        self.dtype = dtype
        self.default = default
        
        self.data = {}

        if strand_specific:
            self.fname = chrom+sense+".tsv"
            self.covered_strands = [chrom + sense]
        else:
            self.fname = chrom+".tsv"
            self.covered_strands = [chrom + '+',chrom + '-']

    def load(self):
        t0 = time()
        debug("# SparseFile: loading '%s' array data for chromosome %s%s from '%s'" % (str(self.dtype),self.chrom,self.sense,self.path))

        d = self.data
        t = self.dtype
        for line in file(os.path.join(self.path,self.fname)):
            #if line.startswith("#"): continue
            pos,value = line.split('\t')
            d[int(pos)] = int(value)
            
        debug("# done. loaded %d entries in %.2f seconds." % (len(self.data),time()-t0))
        
    def get_data(self,chrom,start,end,sense):
        if not self.data:
            self.load()
        return numpy.array([self.data.get(pos,self.default) for pos in xrange(start,end)],dtype=self.dtype)\

    def flush(self):
        self.data = {}

        
class TSVAccessor(Accessor):
    supports_write = True

    def __init__(self,path,chrom,sense,strand_specific=True,dtype=numpy.float32,mode="r",dim=1,attribute="score",default=0,**kwargs):
        debug("# TSVAccessor: '%s' array data for chromosome %s%s from '%s' (mode=%s)" % (str(dtype),chrom,sense,path,mode))
        super(TSVAccessor,self).__init__(path,chrom,sense,strand_specific=strand_specific,**kwargs)
        self.dtype=dtype
        self.strand_specific = strand_specific
        self.dim = dim
        self.attribute = attribute
        self.default = default
        self.mode = mode
        self.sense = sense
        
        if strand_specific:
            fname = os.path.join(path,chrom+sense+".tsv")
        else:
            fname = os.path.join(path,chrom+".tsv")

        self.fname = fname
        self.data = {}
        self.dir = int("%s1" % self.sense)
        # load the track for read-only
        if mode == "r" or "a" in mode:
            try:
                from byo.io.lazytables import NamedTupleImporter
                for row in NamedTupleImporter(fname):
                    for pos in xrange(row.start,row.end):
                        self.data[(sense,pos)] = getattr(row,attribute,default)

                debug("# TSVAccessor: Loaded TSV data for %d positions on %s%s from '%s'" % (len(self.data),chrom,sense,fname))
            except IOError:
                warning("Could not access '%s'. Switching to dummy mode (only zeros)" % fname)
        else:
            warning("# Starting new track, potentially overwriting all data in '%s'" % fname)

        # Create a very thin layer around numpy.ndarray that feeds writes to items or slices back
        # to the Accessor via a reference to self.data
        import numpy
        class WriteableView(numpy.ndarray):
            def __setitem__(this,i,x):
                super(WriteableView,this).__setitem__(i,x)
                if x != self.default:
                    self.data[(self.sense,i+this.start)] = x

            def __setslice__(this,i,j,y):
                super(WriteableView,this).__setslice__(i,j,y)                
                for j,x in enumerate(y):
                    if x != self.default:
                        self.data[(self.sense,i+j+this.start)] = x

            def __new__(cls, input_array, start):
                obj = numpy.asarray(input_array).view(cls)
                #obj.dir = dir
                #if dir < 0:
                    #obj.start = start+len(input_array) -1
                #else:
                    #obj.start = start
                obj.start = start                
                
                return obj

            def __array_finalize__(this, obj):
                if obj is None: return
                this.start = getattr(obj, 'start', None)
                this.dir = getattr(obj, 'dir', None)

        self.viewclass = WriteableView

    def get_data(self,chrom,start,end,sense):
        data = numpy.array([self.data.get((sense,pos),self.default) for pos in xrange(start,end)],dtype=self.dtype)
        return self.viewclass(data,start)

    def flush(self):
        if not "w" in self.mode:
            return

        debug("# TSVAccessor.flush(%s)" % self.fname)
        f = file(self.fname,self.mode)
        header = "## start:int\tend:int\t%s:float\n" % self.attribute
        f.write(header)

        lastscore = self.default
        start = -1
        end = -1
        #print self.fname
        for (sense,pos) in sorted(self.data.keys()):
            #print "%s %10d : %4d" % (sense,pos,self.data[(sense,pos)])
            score = self.data[(sense,pos)]
            #print lastscore,start,end
            if score != lastscore or pos > end+1:
                if start >= 0:
                    f.write("%d\t%d\t%s\n" % (start,end+1,str(lastscore)))
                    #print "->","%d\t%d\t%s\n" % (start,end+1,str(lastscore))
                start = pos

            lastscore = score
            end = pos

        f.write("%d\t%d\t%s\n" % (start,end+1,str(lastscore)))
        #print "->","%d\t%d\t%s\n" % (start,end+1,str(lastscore))

        f.close()
        self.data = {}
