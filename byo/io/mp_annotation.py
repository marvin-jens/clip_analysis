#!/usr/bin/env python
# -*- coding: utf-8 -*-
from logging import getLogger,debug,warning,info,error
from multiprocessing import Pool,Manager
import os

def pileup_loci(loci_list):
    """
    returns all distinct combinations of overlapping (or single) loci with start 
    and end-coordinates of the block.
    Works by sorting loci by start, then end positions. Then, iterates over all loci, keeping a 
    sorted list of end-coordinates to detect which elements to kick out and which to keep.
    Should be near optimal, probably ~O(N log(N)), dominated by the initial sorting.
    """
    from collections import namedtuple
    L = namedtuple("locus","start,end,annotation_id")

    from bisect import bisect_left
    
    # will keep 'short' elements on top
    loci = []
    ends = []
    last_start = 0

    for l in sorted(loci_list):
        l = L(*l)
        if l.start > last_start:
            # find everything that ends before this new block and yield the combinations
            while ends and (ends[0] < l.start):
                yield loci,last_start,ends[0]
                # advance the "cursor" until the end of the smallest element
                last_start = ends[0]
                # and kick it out. (Could be more than one element of same size)
                while ends and ends[0] == last_start:
                    ends.pop(0)
                    loci.pop(0)

            # maybe there is still some space, yield until the new locus begins
            if last_start < l.start:
                if loci: yield loci,last_start,l.start

        # now we are at the start position of the new block
        last_start = l.start

        # find appropriate position by end-position
        pos = bisect_left(ends,l.end)
        loci.insert(pos,l)
        ends.insert(pos,l.end)

    # we're through, yield the remaining blocks away
    while ends:
        yield loci,last_start,ends[0]
        loci.pop(0)
        ends.pop(0)


def do_chrom((track_path,chr_key,size,loci_list)):
    if not loci_list:
        return

    log = getLogger("mp_annotation.%s" % chr_key)

    import numpy
    buf_name = os.path.join(track_path,chr_key+".bin")
    log.debug("# saving buffer to %s" % buf_name)
    buf = numpy.memmap(buf_name,dtype=numpy.uint32,mode="w+",shape=(size,))

    ind_name = os.path.join(track_path,chr_key+"_index.table")
    log.debug("# saving index table to %s" % ind_name)
    index = file(ind_name,"w")

    N = 0
    ann_dict = {"":0}
    
    log.debug("# scanning along the chromosome for unique feature combinations")
    for loci,start,end in pileup_loci(loci_list):
        ann_key = "$".join(sorted([l.annotation_id for l in loci]))

        if not ann_key in ann_dict:
            # this must be a new feature combination
            N += 1
            ann_dict[ann_key] = N
            index.write("%s\n" % ann_key)

        # assign an increasing number to each unique combination
        buf[start:end] = ann_dict[ann_key]

        if (N % 1000) == 0:
            log.debug("#%.1f %% scanned: found %d keys until pos %d/%d. " % (end*100./size,N,end,size))

    log.debug("# freeing buffers")
    buf = None
    index = None        


def both_strands(chrom_and_size):
    for chrom,size in chrom_and_size:
        yield chrom+"+",size
        yield chrom+"-",size
    
def build_mmap_annotation_mp(track_path,loci_sources,genome=None,processes=4):
    pool = Pool(processes=processes)
    manager = Manager()

    # prepare filesystem
    if not os.access(track_path,os.F_OK):
        os.makedirs(track_path)


    debug("# presorting loci")
    from collections import defaultdict
    chr_loci = defaultdict(dict)
    
    for src in loci_sources:
        debug("# from source %s" % src)
        n = 0
        for l in src:
            key = (l.start,l.end,l.annotation_id)
            chr_loci[l.seq_id+l.sense][key] = l
            n += 1
            if hasattr(l,"segments"):
                for s in l.segments:
                    key = (s.start,s.end,s.annotation_id)
                    chr_loci[s.seq_id+s.sense][key] = s
                    n += 1

            if hasattr(l,"features"):
                for f in l.features:
                    key = (f.start,f.end,f.annotation_id)
                    chr_loci[f.seq_id+f.sense][key] = f
                    n += 1

        debug("# %d loci imported" % n)
        debug("# %d strands" % len(chr_loci))

    debug("# building annotation maps for the following strands: %s" % chr_loci.keys())
    result = pool.map(do_chrom,[(track_path,chr_key,size,chr_loci[chr_key].keys()) for chr_key,size in both_strands(genome.chr_sizes.items())])
    pool.close()
    pool.join()

    return result

if __name__ == "__main__":
    from collections import namedtuple
    L = namedtuple("locus","start,end,name")
    
    test_loci = [
        L(0,5,"first"),
        L(0,10,"second_half_ov_w_first"),
        L(0,100,"third_long_contains_first_and_second_and_fourth_overlaps_fifths"),
        L(20,50,"fourth_overlaps_third"),
        L(20,50,"B_fourth_overlaps_third"),
        L(90,110,"fifth_overlaps_thirds_end"),
        L(200,1000,"tx"),
        L(200,201,"tss"),
        L(1000,1001,"pas"),
    ]
    for loci,start,end in pileup_loci(test_loci):
        print start,"-------------",end
        for l in loci:
            print l
        