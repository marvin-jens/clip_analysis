#!/usr/bin/env python
# -*- coding: utf-8 -*-
from byo.io.tables import Importer,Exporter
from byo.io.gff import gff_importer,dict_from_attrstr
from byo.io.lazytables import NamedTupleImporter
from collections import namedtuple


import os,sys,logging
from logging import debug,warning,error,info

from byo.gene_model import *

class LookupCache(object):
    def __init__(self,accessor,**kwargs):
        self.filtered_index = {}
        self.accessor = accessor
        self.name = "dummy"

    @staticmethod
    def process(x):
        return x.rstrip().split("$")

    def __filtered_and_cached(self,key):
        if not key in self.filtered_index:
            self.filtered_index[key] = self.process(self.accessor.index[key])

        return self.filtered_index[key]

    def get_data(self,start,end,sense):
        return [self.__filtered_and_cached(i) for i in self.accessor.data[start:end]]
      

from byo.io.array_accessor import ArrayAccessor
from numpy import uint32, array


class AnnotationAccessor(ArrayAccessor):
    """
    Annotation Tracks consist of a genome-wide integer array with indices to a table with string-lists referencing all annotated features.
    optionially the lookup_table can be processed by a custom filter, which keeps filtering for features of interest out of the inner loops.
    """

    def __init__(self,path,chrom,sense,sense_specific=True,annotation_filter=None,filter_class=LookupCache,dtype=uint32,empty="",**kwargs):

        debug("# AnnotationAccessor mmap: Loading '%s' annotation data for chromosome %s" % (str(dtype),chrom+sense))
        super(AnnotationAccessor,self).__init__(path,chrom,sense,dtype=dtype,sense_specific=sense_specific,dim=getattr(annotation_filter,"dim",1))

        self.annotation_filter = filter_class(self,**kwargs)
        if annotation_filter:
            self.annotation_filter.process = annotation_filter

        self.empty = self.annotation_filter.process(empty)            
        self.load_index(self.index_name(path,chrom,sense),empty)


    def index_name(self,path,chrom,sense):
        return os.path.join(path,chrom+sense+"_index.table")

    def load_index(self,fname,empty):
        debug("# AnnotationAccessor: loading index from '%s'" % fname)
        try:
            self.index = [empty] + [line.rstrip() for line in file(fname)]
        except IOError:
            warning("Could not access '%s'. Switching to dummy mode (only empty)" % fname)
            self.index = [empty]

        debug("# AnnotationAccessor: loaded %d feature combinations" % len(self.index))

    def get_data(self,chrom,start,end,sense):
        return self.annotation_filter.get_data(start,end,sense)

    def get_dummy(self,chrom,start,end,sense):
        return self.annotation_filter.get_data(start,end,sense)
        #print "AnnAccessor.get_dummy"
        #return [self.empty] * (end-start)

        

def transcript_source(source):
    for tx in source:
        start,end = [tx.start,tx.end][::tx.dir]

        # offsets to get the sense informatino +/- 1 madness right. (hopefully)
        a = int(-.5 * (1 - tx.dir))
        b = int( .5 * (1 + tx.dir))
        c = int(-.5 * (1 + tx.dir))
        d = int( .5 * (1 - tx.dir))

        # promoter regions
        yield (tx.chrom,tx.sense,min(start - 500*tx.dir,start),max(start - 500*tx.dir,start),"regulation/promoter_500/", tx.name)
        yield (tx.chrom,tx.sense,min(start - 1000*tx.dir,start),max(start - 1000*tx.dir,start),"regulation/promoter_1000/", tx.name)

        # downstream (possible 3'UTR extension) regions
        yield (tx.chrom,tx.sense,min(end + 500*tx.dir,end),max(end + 500*tx.dir,end),"transcript/downstream_500/", tx.name)
        yield (tx.chrom,tx.sense,min(end + 1000*tx.dir,end),max(end + 1000*tx.dir,end),"transcript/downstream_1000/", tx.name)
        
        # the transcript itself, transcription start and poly adenylation sites
        yield (tx.chrom,tx.sense,tx.start,tx.end,"transcript/%s/" % tx.tx_class,tx.name)
        
        yield (tx.chrom,tx.sense,start+a,start+b,"transcript/processing/TSS/",tx.name)
        yield (tx.chrom,tx.sense,end+a,end+b,"transcript/processing/polyA_site/",tx.name)

        # a gene locus, the transcript region
        yield (tx.chrom,tx.sense,tx.start,tx.end,"gene/", tx.gene_id)
        
        # all exons
        for i,(e_start,e_end) in enumerate(tx.exon_bounds[::tx.dir]):
            yield (tx.chrom,tx.sense,e_start,e_end,"transcript/processing/exon/","%s.exon%02d" % (tx.name,i+1))
            
        # all introns and associated splice-sites
        for i,(ss5,ss3) in enumerate(tx.splice_sites):
            i_start,i_end = min(ss5,ss3),max(ss5,ss3)
            yield (tx.chrom,tx.sense,i_start,i_end,"transcript/processing/intron/","%s.intron%02d" % (tx.name,i+1))
            
            yield (tx.chrom,tx.sense,ss5+a,ss5+b,"transcript/processing/SS5/","%s.intron%02d" % (tx.name,i+1))
            yield (tx.chrom,tx.sense,ss3+c,ss3+d,"transcript/processing/SS3/","%s.intron%02d" % (tx.name,i+1))
            

        if tx.CDS:
            cds_start,cds_end = [tx.CDS.start,tx.CDS.end][::tx.dir]
            
            yield (tx.chrom,tx.sense,min(cds_start,cds_end),max(cds_start,cds_end),"transcript/translation/CDS/",tx.name)
            yield (tx.chrom,tx.sense,cds_start+a,cds_start+b,"transcript/translation/start/",tx.name)
            yield (tx.chrom,tx.sense,cds_end+a,cds_end+b,"transcript/translation/stop/",tx.name)

            if tx.UTR5:
                yield (tx.chrom,tx.sense,tx.UTR5.start,tx.UTR5.end,"transcript/translation/UTR5/",tx.name)
            
            if tx.UTR3:
                yield (tx.chrom,tx.sense,tx.UTR3.start,tx.UTR3.end,"transcript/translation/UTR3/",tx.name)

def gff_source(fname,prefix="transcript/noncoding",**kwargs):
    d,f = os.path.split(fname)
    parts = f.split(".")
    if len(parts) > 2:
        prefix = os.path.join("/".join(parts[:-2]))
        
    for gff in gff_importer(fname,**kwargs):
        attrs = dict_from_attrstr(gff.attr_str)
        name = attrs.get("Name",attrs.get("ID",attrs.get("name",attrs.get("transcript_id",gff.attr_str))))
        path = os.path.join(prefix,gff.type)
        chrom = gff.chrom
        
        yield (chrom,gff.sense,gff.start-1,gff.end,path,name)

        

def pileup_loci(loci_list):
    """
    returns all distinct combinations of overlapping (or single) loci with start 
    and end-coordinates of the block.
    Works by sorting loci by start, then end positions. Then, iterates over all loci, keeping a 
    sorted list of end-coordinates to detect which elements to kick out and which to keep.
    Should be near optimal, probably ~O(N log(N)), dominated by the initial sorting.
    """
    from bisect import bisect_left
    
    # will keep 'short' elements on top
    loci = []
    ends = []
    last_start = 0
    
    # in-place sorting!
    loci_list.sort()

    for l in loci_list:
        chrom,sense,start,end,path,name = l
        if start > last_start:
            # find everything that ends before this new block and yield the combinations
            while ends and (ends[0] < start):
                if ends[0] - last_start > 0:
                    yield loci,last_start,ends[0]
                # advance the "cursor" until the end of the smallest element
                last_start = ends[0]
                # and kick it out. (Could be more than one element of same size)
                while ends and ends[0] == last_start:
                    ends.pop(0)
                    loci.pop(0)

            # maybe there is still some space, yield until the new locus begins
            if last_start < start:
                if loci: yield loci,last_start,start

        # now we are at the start position of the new block
        last_start = start

        # drop everything that ended before this position
        while ends and ends[0] <= last_start:
            ends.pop(0)
            loci.pop(0)

        # find appropriate position by end-position
        pos = bisect_left(ends,end)
        loci.insert(pos,l)
        ends.insert(pos,end)

    # we're through, yield the remaining blocks away
    while ends:
        if ends[0] - last_start > 0:
            yield loci,last_start,ends[0]

        last_start = ends[0]
        loci.pop(0)
        ends.pop(0)


def compress_combination(loci):
    by_path = defaultdict(set)
    for chrom,sense,start,end,path,name in loci:
        by_path[path].add(name)

    return [os.path.join(p,",".join(sorted(list(by_path[p])))) for p in sorted(by_path.keys())]


#class AnnotationLocus(object):
    #def __init__(tx,locus,annotation_id):
        ## TODO: transitional hack to support deprecated Transcript model from transcript.py
        #if hasattr(locus,"seq_id"):
            #tx.seq_id = locus.seq_id
        #else:
            #tx.seq_id = locus.chrom
        #tx.start = int(locus.start)
        #tx.end = int(locus.end)
        #tx.sense = locus.sense
        #tx.annotation_id = annotation_id

#def add_prefix(prefix,source):
    #T = namedtuple("annotation_build","seq_id start end sense annotation_id")

    #for item in source:
        #annotation_id = prefix + item.annotation_id
        #yield T(item.seq_id,item.start,item.end,item.sense,annotation_id)

#def spliced_transcript_source(source,prefix="refGene"):
    #for item in source:
        #yield AnnotationLocus(item,"%s.%s#%s" % (prefix,item.category,item.name))
        #if hasattr(item,"segments"):
            #for seg in item.segments:
                #yield AnnotationLocus(seg,"%s.%s#%s" % (prefix,seg.category,seg.name))
        #if hasattr(item,"features"):
            #for feat in item.features:
                #yield AnnotationLocus(feat,"%s.%s#%s" % (prefix,feat.category,feat.name))

#def rnaGene_source(transcripts):
    #for tx in transcripts:
        #yield AnnotationLocus(tx,"rnaGene.%s#%s" % (tx.type,tx.name))

#def evoFold_source(loci):
    #for l in loci:
        #yield AnnotationLocus(l,"evoFold#%s" % l.name)

##def gff_source(fname,prefix,split_attr=True):
    ##T = namedtuple("annotation_build","seq_id start end sense annotation_id")
    ##for gff in gff_importer(fname):
        ##ID = gff.attr_str.strip()
        ##if split_attr:
            ##s = gff.attr_str.replace("transcript_id ","ID=")
            ##for attr in s.split(';'):
                ##if attr.strip() and "=" in attr:
                    ##name,value = attr.split("=")
                    ##name,value = name.strip(),value.strip()
                    ##if name == "ID" or name == "name":
                        ##ID = value[1:-1]
        ##if not prefix:
            ##prefix = gff.source
        ##yield T(gff.chrom,gff.start-1,gff.end,gff.sense,"%s#%s.%s" % (prefix,ID,gff.type))

#def PARCLIP(path):
    #T = namedtuple("annotation_build","seq_id start end sense annotation_id")

    #for site in NamedTupleImporter(
            #path,
            #cols=["seq_id","start","end","score","sense","annotation_id"],
            #casts=[str,int,int,float,str,str],skip=[1,2,7]
        #):

        #ID = "SID_" + site.annotation_id.split('" (')[0][18:]
        #annotation_id = "%s|score=%.3f" % (ID,site.score)

        #yield T(site.seq_id,site.start,site.end,site.sense,annotation_id)

#def cached(source,src_file):
    #cache_dir = "annotation_build_cache"
    #cache_file = os.path.join(cache_dir,"cached_"+os.path.basename(src_file))

    #if not os.access(cache_dir,os.F_OK):
        #os.makedirs(cache_dir)

    #if not os.access(cache_file,os.R_OK) or os.stat(cache_file).st_ctime < os.stat(src_file).st_mtime:
        #debug("precaching loci from %s" % src_file)
        #Exporter(source(src_file),cache_file,cols=["seq_id","start","end","sense","annotation_id"]).save()

    #return NamedTupleImporter(cache_file)


if __name__ == "__main__":
    from byo.systems import hg19 as hg19
    loci = []
    for i,a in enumerate(transcript_source(hg19.get_refGenes())):
        loci.append(a)
        if i > 1000:
            break

    from pprint import pprint
    for loci,start,end in pileup_loci(loci):
        print start,end
        #pprint(loci)
        pprint(compress_combination(loci))
        
