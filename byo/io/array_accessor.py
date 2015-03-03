from byo.track import Accessor
from logging import debug,warning,error,info

import os
import numpy

class ZeroSource(object):
    def __init__(self,dtype):
        self.dtype = dtype

    def __getitem__(self,x):
        if isinstance(x,slice):
            return self.__getslice__(x.start,x.stop)
        else:
            return self.dtype(0)

    def __getslice__(self,i,j):
        return numpy.zeros(j - i,dtype=self.dtype)
        
    def flush(self):
        pass
        
class ArrayAccessor(Accessor):
    supports_write = True

    def __init__(self,path,chrom,sense,sense_specific=True,dtype=numpy.float32,mode="r",dim=1,**kwargs):
        debug("# ArrayAccessor mmap: Loading '%s' array data for chromosome %s%s from '%s'" % (str(dtype),chrom,sense,path))
        super(ArrayAccessor,self).__init__(path,chrom,sense,sense_specific=sense_specific,**kwargs)
        # in case we havent been given the actual type, just its name: retrieve it
        if type(dtype) == str:
            dtype = getattr(numpy,dtype)

        self.dtype=dtype
        self.sense_specific = sense_specific
        self.dim = dim

        # preferred mode-codes for how to open an mmap array by assumed intention
        mode_trans = {"r" : "c","rw":"w+","w":"w+"}
        
        if sense_specific:
            fname = os.path.join(path,chrom+sense+".bin")
        else:
            fname = os.path.join(path,chrom+".bin")

        if not "w" in mode:
            # read-only access
            try:
                self.data = numpy.memmap(filename=fname,dtype=dtype,mode=mode_trans[mode])
            except IOError:
                warning("Could not access '%s'. Switching to dummy mode (only zeros)" % fname)
                self.data = ZeroSource(self.dtype)
                #self.get_data = self.get_dummy
        else:
            # create a sparse file
            from byo.io import ensure_path
            ensure_path(path)

            import byo.systems
            system = getattr(byo.systems,kwargs['system'])
            try:
                size = system.chr_sizes[chrom]
            except KeyError:
                warning("Unknown chromosome '%s'. Switching to (read-only) dummy mode (only zeros)" % fname)
                self.data = ZeroSource(self.dtype)
            else:
                info("opening '%s' (size=%d) for write-access" % (fname,size))
                self.data = numpy.memmap(filename=fname,dtype=dtype,mode=mode_trans[mode],shape=size)

    def get_sum(self,chrom,start,end,sense,**kwargs):
        return self.get_data(chrom,start,end,sense).sum()

    def get_data(self,chrom,start,end,sense,**kwargs):
        # TODO properly reshape if dim is != 1
        # or is this done in track?
        return self.data[start:end]

    #def get_dummy(self,start,end,sense):
        #return numpy.zeros(end-start,dtype=self.dtype)

    def flush(self):
        if hasattr(self,"data"):
            self.data.flush()
