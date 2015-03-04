#!/usr/bin/env python
from byo import rev_comp,complement
from byo.clip.bamsource import BAM_ClusterGenerator
from byo.clip.scoredcluster import ScoredCluster

import os
import sys
import re
import numpy
numpy.seterr(all='ignore')
from numpy import *
from optparse import *

usage = """
  %prog [options] <sorted_alignments.bam> [another_library.bam] [...]

reads alignment data and produces clusters with various quality measures.
Optionally also produces coverage + conversion tracks for the rest of the toolchain.
If multiple input files are given --mode should be a comma-separated list of conversion 
signatures like -m TC,GA,GA,TC. Default is TC for everything.
When evaluating multiple libraries at once, the --min_support and --require_conversions 
switch allows for consensus filtering.
"""

parser = OptionParser(usage=usage)
parser.add_option("-l","--logfile",dest="logfile",default="parclip.log",help="logfile to use")
parser.add_option("-S","--system",dest="system",type=str,default="hg19",help="model system database (hg19|ce6|mm9)")
parser.add_option("-t","--track",dest="track",default="",help="path to write tracks into")
parser.add_option("-F","--fp_flag",dest="fp_flag",default="markov",help="the flag representing false-positive (default='markov'). Only used for gathering statistics.")
parser.add_option("-d","--discard",dest="discard",default="",help="discard clusters overlapping with this track")
parser.add_option("-M","--mask",dest="mask",default=[],action="append",help="can be used multiple times to pass genomic intervals (samtools 0-based syntax chrX:start-end) that should be skipped for cluster-building (e.g. ribosomal DNA)")
parser.add_option("-m","--mode",dest="mode",default="TC",help="conversion events to look for")
parser.add_option("-s","--stats",dest="stats",default="",help="filename to dump conversion statistics to")
parser.add_option("-C","--cluster-stats",dest="cluster_stats",default="",help="filename for storing detailed cluster scores")
parser.add_option("-c","--conv",dest="minconv",type=int,default=1,help="minimal number of conversion events to be reported (default=1)")
parser.add_option("","--require_conversions",dest="require_conversions",default=False,action="store_true",help="require conversions in other libraries to keep a cluster")
parser.add_option("","--min_support",dest="min_support",type=int,default=1,help="minimum number of libraries with read support")
parser.add_option("","--min_reads",dest="min_reads",type=int,default=2,help="minimum number of reads in a reported cluster (default=2)")
parser.add_option("","--ccr",dest="ccr",type=int,default=0,help="extract strongest crosslink centered regions with this flank-size and report on stderr")
parser.add_option("-R","--region",dest="region",default="",help="only scan this (samtools style, 0-based end exclusive) region, e.g. chr1:1234567-1234667")
parser.add_option("-G","--groups",dest="groups",default="",help="comma-separated list of group identifiers to assign to input BAMs instead of filename. Allows to combine technical replicates in consensus calling. Number of groups must match number of BAMs and comes in the same order.")
parser.add_option("-u","--min_uniq",dest="min_uniq",type=float,default=1,help="minimal alignment score diff. between best and second to consider as 'unique'. default=1")
parser.add_option("-r","--read_count",dest="count_func",default='uniq',choices=['uniq','lin','arcsinh'],help="read count normalization: uniq: count each unique (distinct) read only once (default), lin: linear, arcsinh: linear until ~10, then log-like")


options,args = parser.parse_args()

import logging
logging.basicConfig(level=logging.INFO,filename=options.logfile,format='%(asctime)s %(levelname)-8s %(name)s %(message)s')
logger = logging.getLogger("clip.py")

if len(args) < 1:
    parser.error("insufficient arguments")
    sys.exit(1)

import importlib
system = importlib.import_module("byo.systems.{options.system}".format(**locals()))

from byo.track import Track
from byo.io.track_accessors import ArrayAccessor,TSVAccessor

if options.track:
    cov_track = Track(os.path.join(options.track,"coverage"),TSVAccessor,mode="w",dtype=numpy.uint32,auto_flush=True)
    conv_track = Track(os.path.join(options.track,"conversions"),TSVAccessor,mode="w",auto_flush=True)

def update_tracks(cluster):
    cov_track.get(cluster.chrom,cluster.start,cluster.end,cluster.strand)[:] = cluster.coverage
    conv_track.get(cluster.chrom,cluster.start,cluster.end,cluster.strand)[:] = cluster.conversions


if options.discard:
    def nonempty_channel(line):
        return bool(line.strip())

    discard = system.get_annotation_channel(path=options.discard,channel=nonempty_channel,auto_flush=True)

masked = []
if options.mask:
    for mask in options.mask:
        chrom,start,end = mask.replace("-",":").split(":")
        masked.append((chrom,int(start),int(end)))

from byo.clip.statistics import Stats
all_stats = Stats(path=options.stats)
kept_stats = Stats()

if options.cluster_stats:
    raw_stats = file(options.cluster_stats,'w')

def debug_output(cluster):
    if cluster.strand == "+":
        return
    
    def conv_str(conv):
        l = []
        for c in conv:
            if c:
                l.append(".")
            else:
                l.append(' ')
        return "".join(l)
    
    ref = system.genome.get(cluster.chrom,cluster.start,cluster.end,cluster.strand)
    print cluster
    print ">>>>"
    print ref
    print cluster.sequence
    print conv_str(cluster.conversions)
    print "-----"

from collections import defaultdict
N = defaultdict(int)

try:
    # make or parse a list of group identifiers for the BAM files
    lib_groups = options.groups.split(",")
    if '' in lib_groups: lib_groups.remove('')
    if not lib_groups:
        lib_groups = [os.path.basename(bam) for bam in args]
        
    # and link them to the respective signature to scan for within these
    signature_dict = dict(zip(lib_groups,options.mode.split(",")))

    # build a generator for ScoredCluster objects from BAM coverage
    cluster_generator = BAM_ClusterGenerator(
        bams=args,groups=lib_groups,region=options.region,cluster_type=ScoredCluster,
        # these are passed down as kwargs to the constructor of each ScoredCluster
        signature=signature_dict,count_func=options.count_func)

    # go over each cluster as it is generated, and report it if minimal requirements are met
    for cluster in cluster_generator:
        if cluster.chrom.startswith('rRNA'): continue

        if options.track:
            update_tracks(cluster)

        N["total"] += 1 
        all_stats.update_stats("total",cluster.editstats,cluster.strand)

        decoy_hit = bool(re.search(options.fp_flag,cluster.chrom))
        if decoy_hit:
            N["decoy_hits"] += 1
            all_stats.update_stats("decoy",cluster.editstats,cluster.strand)
        else:
            N["reference_hits"] += 1
            all_stats.update_stats("reference",cluster.editstats,cluster.strand)

        # now we can start to discard clusters without skewing the 'total' statistics
        if cluster.uniq_reads < options.min_reads:
            N["discarded_few_reads"] += 1
            continue

        total_conv = cluster.conversions.sum()
        if total_conv < options.minconv:
            N["discarded_few_conversions"] += 1
            #print "few_conversions"
            continue

        if cluster.max_uniq < options.min_uniq:
            N["discarded_non_uniq"] += 1
            #print "skip"
            continue
            
        read_support = len(cluster.read_support.keys())
        if read_support < options.min_support:
            logger.debug("cluster '%s' has insufficient multilibrary read coverage: %d. Discarded." % (cluster.short,read_support))
            N["discarded_insuff_read_support"] += 1
            continue

        conv_support = len(cluster.sig_support.keys())
        if options.require_conversions and conv_support < options.min_support:
            logger.debug("cluster '%s' has insufficient multilibrary conversions: %d. Discarded." % (cluster.short,conv_support))
            N["discarded_insuff_conv_support"] += 1
            continue
            
        if options.discard:
            discard_hits = array(discard.get(cluster.chrom,cluster.start,cluster.end,cluster.strand)).sum()
            
            if discard_hits:
                logger.debug("cluster '%s' overlaps %d elements from the discard track. Discarded." % (cluster.short,discard_hits))
                N["discarded_discard_track_overlap"] += 1
                continue

        if masked:
            hit = False
            for chrom,start,end in masked:
                if chrom != cluster.chrom:
                    continue
                if cluster.end < start:
                    continue
                if cluster.start > end:
                    continue
                hit = (chrom,start,end)

            if hit:
                logger.debug("warning '%s' overlaps masked region '%s'. Discarded." % (cluster.short,str(hit)))
                N["discarded_masked_region"] += 1
                continue

        N["kept"] += 1
        kept_stats.update_stats("total",cluster.editstats,cluster.strand)
        
        if decoy_hit:
            # update statistics of kept and reported clusters
            N["decoy_kept"] += 1
            kept_stats.update_stats("decoy",cluster.editstats,cluster.strand)
        else:
            N["reference_kept"] += 1
            kept_stats.update_stats("reference",cluster.editstats,cluster.strand)
        
        conv_peak = cluster.conversions.argmax() + cluster.start
        total_edits = cluster.editstats['edits']
        total_entropy = cluster.entropy_score
        
        cid = cluster.name #"CID_%06d" % N["kept"]
        attr_str = 'name="%s"; n_reads="%s"; n_uniq="%s"; ccr="%d";' % \
        (cid,cluster.n_reads,cluster.uniq_reads,conv_peak+1)

        out = [cluster.chrom,"clip.py","cluster",cluster.start+1,cluster.end,cluster.score,cluster.strand,total_entropy,attr_str]
        print "\t".join([str(o) for o in out])

        if options.ccr:
            ccr = [cluster.chrom,"clip.py","CCR",conv_peak+1-options.ccr,conv_peak+1+options.ccr,cluster.score,cluster.strand,total_entropy,attr_str]
            sys.stderr.write("%s\n" % "\t".join([str(o) for o in ccr]))

        if options.cluster_stats:
            if N["kept"] == 1:
                header = ["cid","decoy"] + cluster.reported_scores
                raw_stats.write("#%s\n" % "\t".join([str(h) for h in header]))

            raw = [cid,str(decoy_hit)] + [str(getattr(cluster,s)) for s in cluster.reported_scores]
            raw_stats.write("\t".join([str(r) for r in raw])+'\n')
        
except KeyboardInterrupt:
    pass

if options.stats:
    all_stats.save_dict(title="cluster building statistics",stats=N)
    all_stats.save_editstats(title="conversion statistics")
    all_stats.save_editstats(title="kept cluster conversion statistics",editstats = kept_stats.editstats)
    all_stats.flush()

