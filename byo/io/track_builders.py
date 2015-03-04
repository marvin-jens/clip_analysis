#!/usr/bin/env python
# -*- coding: utf-8 -*-
from byo.io.tables import Importer,BioData
from numpy import *
import sys,os
from logging import debug,warning,info,error
from mp_annotation import build_mmap_annotation_mp

def filter_unique(src):
    last_id = None
    last_r = None
    count = 0
    for r in src:
        if r.read_id != last_id:
            if count == 1:
                yield last_r
            last_r = r
            count = 0
            last_id = r.read_id
        count +=1

    if count == 1:
        yield last_r

def both_strands(chrom_and_size):
    for chrom,size in chrom_and_size:
        yield chrom+"+",size
        yield chrom+"-",size


def build_score_track(source,track_path,sense_specific=False,genome=None,dtype=float32):
    from byo.io.wiggle import wig_reader

    # prepare filesystem
    if not os.access(track_path,os.F_OK):
        os.makedirs(track_path)

    # allocate chromosome maps
    debug("# preparing chromosome buffers...")
    chr_cov = {}
    if sense_specific:
        for chrom,size in both_strands(genome.chr_sizes.items()):
            chr_cov[chrom] = zeros(size,dtype=dtype)
    else:
        for chrom,size in genome.chr_sizes.items():
            chr_cov[chrom] = zeros(size,dtype=dtype)

    debug("# collecting data")
    for item in source:
        if not item.chrom in chr_cov:
            warning("# ignoring wiggle data for unknown chromosome %s" % item.chrom)
            continue

        for pos in xrange(item.start,item.end):
            chr_cov[item.chrom][pos] = item.score

    # store track data to disk
    debug("# saving to disk...")
    for name,cov in chr_cov.items():
        fname = os.path.join(track_path,name+".bin")
        debug("# %s" % fname)
        cov.tofile(fname)


def build_score_track_from_wig(wig_path,track_path,genome=None,dtype=float):
    from byo.io.wiggle import wig_reader
    from glob import glob

    # prepare filesystem
    if not os.access(track_path,os.F_OK):
        os.makedirs(track_path)

    # allocate chromosome maps
    debug("# preparing chromosome buffers...")
    chr_cov = {}
    for chrom,size in genome.chr_sizes.items():
        chr_cov[chrom] = zeros(size,dtype=dtype)

    # prepare parsing of wig files
    for wig_file in glob(wig_path):
        for chrom,block in wig_reader(wig_file):
            if not chrom in chr_cov:
                warning("# ignoring wiggle data for unknown chromosome %s" % chrom)
                continue

            for pos,score in block:
                chr_cov[chrom][pos] = score

    # store track data to disk
    debug("# saving to disk...")
    for name,cov in chr_cov.items():
        fname = os.path.join(track_path,name+".bin")
        debug("# %s" % fname)
        cov.tofile(fname)


def build_coverage_track(read_src,track_path,unique=True,genome=None):
    # read_src must provide objects with read_id, chr, start, end and sense attributes

    # prepare filesystem
    if not os.access(track_path,os.F_OK):
        os.makedirs(track_path)

    # allocate chromosome maps
    debug("# preparing chromosome buffers...")
    chr_cov = {}
    for chrom,size in genome.chr_sizes.items():
        chr_cov[chrom+"+"] = zeros(size,dtype=uint32)
        chr_cov[chrom+"-"] = zeros(size,dtype=uint32)

    if unique:
        read_src = filter_unique(read_src)

    # count into chromosomal coverage map
    debug("# parsing reads...")
    for i,b in enumerate(read_src):
        name = b.read_id.split("|")[0].strip()

        b.count = 1
        edit = 0
        if "_x" in name:
            nameparts = name.split("_")
            count,edit = nameparts[-2:]
            b.count = int(count[1:])

        b.start -= 1 # using a zero-based, half-open index-convention

        chrom = b.chrom+b.sense
        if not chrom in chr_cov:
            warning("# ignoring read mapping to unknown chromosome %s" % chrom)
            continue

        c = chr_cov[chrom]
        try:
            for j in arange(b.start,b.end):
                c[j] += b.count
        except:
            print "index mess:",b.start,b.end,b.chrom,b.sense,"chr_size",chr_sizes[b.chr]

    # store track data to disk
    debug("# saving to disk...")
    for name,cov in chr_cov.items():
        fname = os.path.join(track_path,name+".bin")
        debug("# %s" % fname)
        cov.tofile(fname)

def build_mmap_annotation(track_path,loci_sources,genome=None):
    # prepare filesystem
    if not os.access(track_path,os.F_OK):
        os.makedirs(track_path)

    from collections import defaultdict
    from byo.lowlevel.annotation_core import AnnotatedSeq
    from struct import pack

    debug("# presorting loci")
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

    def do_chrom((chr_key,size)):
        debug("# allocating %d elements raw buffer for %s..." % (size,chr_key))
        #buf = zeros(size,dtype=uint32)
        buf = memmap(filename,dtype=uint32,mode="w+",shape=(size,))

        debug("# allocating %d elements annotation buffer for %s..." % (size,chr_key))
        chr_ann = AnnotatedSeq(size)

        loci_list = chr_loci[chr_key].keys()
        debug("# populating annotation buffer with %d loci" % len(loci_list) )
        chr_ann.add_annotation_list(loci_list)

        # scan along the chromosome and hash every distinct list of features
        anno_dict = {"" : 0}
        anno_list = [""]

        buf_name = os.path.join(track_path,chr_key+".bin")
        debug("# saving buffer to %s" % buf_name)

        ind_name = os.path.join(track_path,chr_key+"_index.table")
        debug("# saving index table to %s" % ind_name)
        index = file(ind_name,"w")

        debug("# scanning along the chromosome for unique feature combinations")
        for x in xrange(size):
            ann = chr_ann.whats_at(x)
            ann_key = "$".join(ann)

            if not ann_key in anno_dict:
                N = len(anno_list)
                anno_dict[ann_key] = N
                anno_list.append(ann_key)
                index.write("%s\n" % ann_key)

                if (N % 100) == 0:
                    debug("#%.1f %% of %s scanned: found %d keys until pos %d/%d. " % (x*100./size,N,chr_key,x,size))
            key = anno_dict[ann_key]
            # keep sparseness of files
            if key:
                buf[x] = key

        debug("# saving buffer to %s" % buf_name)
        buf.tofile(buf_name)

        debug("# freeing buffers")
        buf = None
        index = None
        chr_ann.free()
        chr_ann = None

    for chr_key,size in chr_key,size in both_strands(genome.chr_sizes.items()):
        do_chrom(chr_key,size)
        import gc
        gc.collect()
        debug("#done %s" % chr_key)

