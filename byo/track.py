# -*- coding: utf-8 -*-
from logging import getLogger
import os


class Accessor(object):
    supports_write = False

    def __init__(self,path,chrom,sense,sense_specific=False,**kwargs):
        if sense_specific:
            self.covered_strands = [chrom+sense]
        else:
            self.covered_strands = [chrom+'+',chrom+'-']

    def get_data(self,chrom,start,end,sense,**kwargs):
        return []

    def get_oriented(self,chrom,start,end,sense,**kwargs):
        data = self.get_data(chrom,start,end,sense,**kwargs)
        if sense == "-": #  self.sense_specific and
            return data[::-1]
        else:
            return data

    def get_sum(self,chrom,start,end,sense,**kwargs):
        return self.get_data(chrom,start,end,sense,**kwargs).sum()

    def flush(self):
        pass


class DummyCache(object):
    """
    A Dummy dictionary that will claim to have every key
    and return always the same value.
    """
    def __init__(self,retval):
        self.retval = retval

    def __contains__(self,item):
        return True

    def __getitem__(self,name):
        return self.retval

    def __setitem__(self,name,item):
        self.retval = item
        
class Track(object):
    """
    Abstraction of chromosome-wide adressable data like sequences, coverage, scores etc.
    Actual access to the data is delegated to accessor objects which are instantiated on-the-fly for
    each chromosome (strand) upon first access and then cached.
    Use of mmap for the accessors is recommended and implemented for sequences and numpy (C-type)
    arrays.

    See io/track_accessors.py for more examples.
    """

    def __init__(self,path,accessor,sense_specific=True,description="unlabeled track",system="hg18",dim=1,auto_flush=False,mode="r",**kwargs):
        self.path = path
        self.mode = mode
        self.acc_cache = {}
        self.accessor = accessor
        self.kwargs = kwargs
        self.sense_specific = to_bool(sense_specific)
        self.dim = int(dim)
        self.description = description
        self.auto_flush = auto_flush
        self.last_chrom = ""
        self.logger = getLogger("Track('%s')" % path)
        self.system = system

        self.logger.debug("Track(auto_flush=%s)" % (str(auto_flush)))
        kwargs['sense_specific'] = self.sense_specific
        kwargs['mode'] = self.mode
        kwargs['system'] = self.system
        kwargs['description'] = self.description
        kwargs['dim'] = self.dim

        if "w" in mode:
            # make sure the path exists right away so that the accessors 
            # can flush the actual data there!
            from byo.io import ensure_path
            ensure_path(self.path)

    def load(self,chrom,sense):
        # automatically flush buffers whenever a new chromosome is seen. reduces memory-footprint for sorted input data
        if self.auto_flush and chrom != self.last_chrom:
            self.logger.debug("Seen new chromosome %s. Flushing accessor caches." % chrom)
            self.flush_all()
        self.last_chrom = chrom

        ID = chrom+sense
        if not ID in self.acc_cache:
            self.logger.debug("Cache miss for %s%s. creating new accessor" % (chrom,sense))
            acc = self.accessor(self.path,chrom,sense,**(self.kwargs))
            
            # register the accessor for as many chromosomes/strands/contigs as it feels responsible for.
            self.logger.debug("Accessor covers strands '%s'" % str(sorted(set(acc.covered_strands))))
            if acc.covered_strands == '*':
                # hint that this accessors covers the complete genome
                self.logger.debug("disabling auto_flush and registering accessor for everything")
                self.acc_cache = DummyCache(acc)
                self.auto_flush = False
            else:
                for ID in acc.covered_strands:
                    self.acc_cache[ID] = acc

        return self.acc_cache[ID]

    def save(self):
        """
        If a track is opened in "rw" or "w" mode this will save the track-definition config files and flush all accessors.
        Saving of the actual data is performed by the accessors that support writing upon a call to flush.
        """
        if not "w" in self.mode:
            self.logger.warning("save() called on a read-only opened track. Ignored!")
            return

        if not self.accessor.supports_write:
            self.logger.warning("save() called on a track with only read-access supporting accessors. Ignored!")
            return
      
        self.logger.debug("save(): writing '%s'" % self.path)

        def to_str(obj):
            # convert simple data-types to their string representation
            # but classes and more complex types to their names.
            return getattr(obj,"__name__",str(obj))

        kwarg_str = "\n".join(["%s=%s" % (k,to_str(self.kwargs[k])) for k in sorted(self.kwargs.keys()) if k != "mode"])
        file(os.path.join(self.path,"track.rc"),"w+").write(trackrc % dict(accessor=self.accessor.__name__,kwargs=kwarg_str))
        self.flush_all()

    def __del__(self):
        if "w" in self.mode:
            self.save()

    def flush(self,chrom,sense):
        ID = self.get_identifier(chrom,sense)
        if ID in self.acc_cache:
            self.logger.warning("Flushing %s%s" % (chrom,sense))
            del self.acc_cache[ID]

    def flush_all(self):
        for a in self.acc_cache.values():
            a.flush()
        self.acc_cache = {}

    def get(self,chrom,start,end,sense,**kwargs):
        acc = self.load(chrom,sense)
        return acc.get_data(chrom,start,end,sense,**kwargs)

    def get_oriented(self,chrom,start,end,sense,**kwargs):
        acc = self.load(chrom,sense)
        return acc.get_oriented(chrom,start,end,sense,**kwargs)

    def get_sum(self,chrom,start,end,sense):
        #print "TRACK.GETSUM"
        acc = self.load(chrom,sense)
        return acc.get_sum(chrom,start,end,sense)
        
    def get_identifier(self,chrom,sense):
        if self.sense_specific:
            return chrom+sense
        else:
            return chrom

from numpy import float32

trackrc = """
[track]
accessor_type=%(accessor)s

[kwargs]
%(kwargs)s
"""

def to_bool(obj):
    if obj == "False":
        return False
    else:
        return bool(obj)

def load_track(path,default_accessor="ArrayAccessor",**kwargs):
    """
    Factory function that returns a Track instance to access the data stored in <path>.
    It determines the correct ArrayAccessor to use by looking for the config file track.rc
    in the directory.
    """
    from ConfigParser import SafeConfigParser as ConfigParser, NoSectionError
    import os
    import byo.io.track_accessors as track_accessors
    if path.endswith(".bam") and os.path.isfile(path):
        return Track(path,track_accessors.BAMAccessor,description="BAM('%s')" % os.path.basename(path), **kwargs)
        
    track_dict = dict(accessor_type=default_accessor,description=os.path.basename(path))
    try:
        cp = ConfigParser()
        cp.read(os.path.join(path,"track.rc"))
        track_dict.update(dict(cp.items("track")))
        kwargs.update(dict(cp.items("kwargs")))
        kwargs.update(track_dict)
        
    #except IOError,NoSectionError:
    # TODO: Do not just silence everything. Might shadow important issues
    except:
        pass
    #from pprint import pprint
    #pprint(kwargs)
    accessor = getattr(track_accessors,track_dict["accessor_type"])
    return Track(path,accessor,**kwargs)

# convenience function for often-needed track type
def array_track(path,dtype=float32,**kwargs):
    return load_track(path,dtype=dtype,description=path,**kwargs)


if __name__ == "__main__":
    import logging
    logging.basicConfig(level=logging.DEBUG)
    print "lets do some track-testing"
    from byo.io.track_accessors import ArrayAccessor,TSVAccessor
    import numpy

    #t = Track('test_track',ArrayAccessor,mode="w",dtype=numpy.uint32,sense_specific=False)
    #t.get_oriented('chr1',200,300,"-")[:] = numpy.arange(100)

    #t.save()
    
    #t = load_track('test_track')
    #print t.get_oriented('chr1',200,300,"-")
    
    t = Track('test_track_tsv',TSVAccessor,mode="w+",dtype=numpy.uint32,sense_specific=False)
    t.get_oriented('chr1',200,300,"-")[:] = numpy.arange(100)

    t.save()
    
    t = load_track('test_track_tsv')
    
    print t.get_oriented('chr1',200,300,"+")

    
    
    