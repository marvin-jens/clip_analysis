#!/usr/bin/env python
# -*- coding: utf-8 -*-
from byo.annotation import *
from byo.gene_model import transcripts_from_UCSC
from byo.io.track_builders import build_mmap_annotation, build_mmap_annotation_mp
import numpy
import gzip

import os,sys,logging,re
logging.basicConfig(level=logging.DEBUG,format='%(asctime)10s '+logging.BASIC_FORMAT)

from optparse import *

usage = """

  %prog [options] <features.gff> <gene_models.ucsc> <... >

Read features from supplied gffs/gene models from ucsc tables and build a genome-wide annotation database for fast scanning.
"""

parser = OptionParser(usage=usage)
parser.add_option("-S","--system",dest="system",type=str,default="hg19",help="model system (hg18,ce6,...)")
parser.add_option("-p","--processes",dest="processes",default=None,type=int,help="number of processes to use. (default=cpu_count())")
parser.add_option("-d","--destination",dest="destination",default="system_annotation",help="destination (default: annotation_track)")
parser.add_option("-G","--gene_names",dest="gene_names",type=str,default="",help="additional gene name translation table")
parser.add_option("-F","--dontfixchr",dest="fix_chr",default=True,action="store_false",help="don't add 'chr' to chrom identifiers that lack it (e.g for Contigs)")
options,args = parser.parse_args()

import byo.systems
system = getattr(byo.systems,options.system)

if not system.chr_sizes:
    logging.error("Uninitialized model system '%s'. Run make as described in the output of bootstrap.py, aborting!" % options.system)
    sys.exit(1)

gene_names = {}
if not options.gene_names:
    gn_path = os.path.join(system.root,"annotation","gene_names")
else:
    gn_path = options.gene_names

if os.path.exists(gn_path):
    for line in file(gn_path):
        name,gene = line.rstrip().split('\t')
        gene_names[name] = gene

def do_chrom((chr_key,size,fname)):
    log = logging.getLogger("mp_annotation.%s" % chr_key)
    log.debug("# reading loci from %s" % fname)

    loci_list = []
    for line in gzip.open(fname,"rb"):
        chrom,sense,start,end,path,name = line.rstrip().split('\t')
        loci_list.append((chrom,sense,int(start),int(end),path,name))

    if not loci_list:
        log.debug("# no loci to process. exiting" % fname)
        return
    else:
        log.debug("# preparing to build annotation map for %d loci on '%s'" % (len(loci_list),chr_key))

    buf_name = os.path.join(options.destination,chr_key+".bin")
    log.debug("# saving buffer to %s" % buf_name)
    buf = numpy.memmap(buf_name,dtype=numpy.uint32,mode="w+",shape=(size,))

    ind_name = os.path.join(options.destination,chr_key+"_index.table")
    log.debug("# saving index table to %s" % ind_name)
    index = file(ind_name,"w")

    N = 0
    ann_dict = {"":0}
    
    log.debug("# scanning along the chromosome for unique feature combinations")
    for loci,start,end in pileup_loci(loci_list):
        ann_key = "$".join(compress_combination(loci))
        
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


# prepare filesystem
if not os.access(options.destination,os.F_OK):
    os.makedirs(options.destination)


from multiprocessing import Pool
pool = Pool(processes=options.processes)

debug("# presorting loci")
from collections import defaultdict
chr_loci = defaultdict(list)

for x in args:
    #take the filename except the extension as prefix
    if not os.access(x,os.R_OK):
        logging.warning("cannot access '%s'" % x)
        continue

    fname = os.path.basename(x)
    
    if x.lower().endswith(".gff"):
        #prefix = "/".join(os.path.basename(x).split(".")[:-1])
        src = gff_source(x,fix_chr=options.fix_chr)

    elif x.lower().endswith(".ucsc"):

        if fname.startswith('transcript'):
            tx_class = fname.split(".")[1]
        else:
            tx_class = None

        src = transcript_source(transcripts_from_UCSC(x,system=options.system,tx_class=tx_class,gene_names=gene_names,fix_chr=options.fix_chr))
    else:
        logging.warning("unknown annotation type '%s'. skipping." % fname)
        continue

    debug("# from source %s" % fname)
    n = 0
    for l in src:
        chrom,sense,start,end,path,name = l
        chr_loci[chrom+sense].append(l)
        n += 1

    debug("# %d features imported" % n)
    debug("# %d strands" % len(chr_loci))


debug("before storing annotation, using '%d' bytes for chr_loci dictionary" % (sys.getsizeof(chr_loci)) )
for strand in sorted(chr_loci.keys()):
    loci_list = chr_loci[strand]
    if not strand[:-1] in system.chr_sizes:
        warning("skipping unknown chromosome/strand '%s'" % strand)
        del chr_loci[strand]
        continue

    if not loci_list:
        debug("skipping '%s' lacking any annotation")
        del chr_loci[strand]
        continue
        
    fname = os.path.join(options.destination,"%s.gz" % strand)
    debug("storing precompiled annotation in '%s'" % fname)
    buf = gzip.open(fname,"wb")
    for l in loci_list:
        buf.write("%s\t%s\t%d\t%d\t%s\t%s\n" % l)

    buf.close()
    chr_loci[strand] = fname
    
import sys
debug("after storing annotation, now using '%d' bytes for chr_loci dictionary" % (sys.getsizeof(chr_loci)) )

strands = sorted(chr_loci.keys(),key = lambda x : len(chr_loci[x]),reverse=True)
sizes = [system.chr_sizes[s[:-1]] for s in strands]

import gc
gc.collect()

if chr_loci:
    debug("putting chromosome strands into pool.map()")
    result = pool.map(do_chrom,[(chr_key,size,chr_loci[chr_key]) for chr_key,size in zip(strands,sizes)])
    #result = map(do_chrom,[(chr_key,size,chr_loci[chr_key]) for chr_key,size in zip(strands,sizes)])
    debug("pool.close()")
    pool.close()
    pool.join()

