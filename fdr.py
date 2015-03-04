#!/usr/bin/env python
import os,sys
from time import time
import itertools
from numpy import *
from optparse import *


seterr(all='ignore')
usage = """
  %prog [options] <cluster_stats.tsv> <all_clusters.gff> 

Analyze cluster statistics and estimate false discovery rates based on reference/decoy hit distributions.
"""

parser = OptionParser(usage=usage)
parser.add_option("-l","--logfile",dest="logfile",default="fdr.log",help="logfile to use")
parser.add_option("-s","--statfile",dest="statfile",default="",help="write statistics to file")
parser.add_option("-c","--cutoff",dest="max_fdr",default=.05,type=float,help="maximally tolerated FDR (default=.05)")
parser.add_option("-d","--decoy",dest="decoy_scale",default=1.,type=float,help="efficiency of decoy compared to reference (typically < 1.0)")
parser.add_option("-2","",dest="combinations",default=False,action="store_true",help="also explore score combinations (takes a long time, off by default)")
parser.add_option("-p","--pdf",dest="pdf",default="",help="render PDFs to this path")

options,args = parser.parse_args()

if len(args) < 1:
    parser.error("insufficient arguments")
    sys.exit(1)

if os.path.exists(options.logfile):
    os.remove(options.logfile)

import logging

logging.basicConfig(level=logging.DEBUG,filename=options.logfile)
logger = logging.getLogger()
false_positives = []
true_positives = []

from collections import defaultdict
stats_by_map_qual = defaultdict(lambda : defaultdict(list))
stats_by_id = {}
raw_data = []

from time import time
class ClusterSet(object):
    def __init__(self,stat_source):
        """
        A ClusterSet object is created from the table data in 
        *.cluster_stats.tsv output of parclip.py
        It offers methods to apply scoring functions to the raw quality
        metrics computed by parclip.py and furthermore to select cutoffs
        on these scores or raw metrics to deplete false positives 
        (FP=decoy hits) more strongly than true positives (TP=reference 
        hits) and thus satisfy a false discovery rate (FDR) limit.
        """
        t0 = time()

        # load all cluster stats into a list
        self.all_data = list(stat_source)

        # create a dictionary with all the columns (cid,map_qual,sense, ...)
        self.data = defaultdict(list)
        if not self.all_data:
            self.N = 0
            self.tp_index = set()
            self.fp_index = set()
            self.unfiltered_fdr = 0
        else:
            columns = self.all_data[0]._fields
            for c in self.all_data:
                for k,v in c._asdict().items():
                    self.data[k].append(v)

            # convert columns to arrays for speed and convenience
            for k,l in self.data.items():
                self.data[k] = array(l)

            self.N = len(self.data['cid'])
            flags = self.data['sense']
            
            self.tp_index = (flags == 'reference').nonzero()[0]
            self.fp_index = (flags == 'decoy').nonzero()[0]

            self.unfiltered_fdr = self.fdr()
        
        logger.debug("loaded %d cluster stats in %.3f sec. tp=%d fp=%d" % (self.N,time()-t0,len(self.tp_index),len(self.fp_index)))
        
        #bp,fdr,kept,loss
        self.kept_by_method = {('unfiltered',):(arange(self.N),None,self.unfiltered_fdr,self.N,0), ('none',):([],-1,1,1)}

    def fdr(self,indices=set()):
        """
        estimate the false discovery rate for (a subset) of a cluster set
        based on annotation as reference aligning or decoy hit.
        """
        if not self.N:
            return 0.5

        if not len(indices):
            indices = arange(self.N)
            
        flags = take(self.data['sense'],indices)
        return self._fdr((flags == 'reference').sum(),(flags == 'decoy').sum())

    @staticmethod
    def _fdr(real,decoy):
        #print "_FDR",real,decoy
        fp = decoy / options.decoy_scale
        
        # assume as many fp on reference than on (scaled) decoy
        tp = max(0,real - fp)

        # pseudo counts: in absence of data, est. FDR -> 0.5
        return (fp + 1.) / (tp + fp + 2.)
        
    def score(self,name,func):
        """
        evaluate an arbitrary scoring function on the cluster quality data and
        keep the results as new column.
        """
        
        t0 = time()
        self.data[name] = array([func(c) for c in self.all_data])
        logger.debug("computed scoring function '%s' in %.3f sec." % (name,time()-t0))
        
    def find_cutoff(self,column,max_fdr=options.max_fdr):
        """
        search a cutoff on the given column (score) which would deplete
        FP more than TP and thus satisfy max_fdr.
        """
        if (column,) in self.kept_by_method:
            return self.kept_by_method[(column,)][1:]

        t0 = time()
        from bisect import bisect_left
        logger.debug("searching cutoff on '%s' to satisfy FDR <= %.1f %%" % (column,max_fdr*100))

        data = self.data.get(column)
        breakpoints = sorted(set(data))
    
        #print data,self.tp_index,self.fp_index
        data_tp = data.take(self.tp_index)
        data_fp = data.take(self.fp_index)
        
        #print len(data_tp),len(data_fp)
        # indices in sorted order (ascending)
        I = data.argsort()
        I_tp = data_tp.argsort()
        I_fp = data_fp.argsort()
        
        D = data.take(I)
        TP = data_tp.take(I_tp)
        FP = data_fp.take(I_fp)
        
        total = self.N
        tp0 = len(data_tp)
        fp0 = len(data_fp)
        
        for bp in breakpoints:
            #above = (data >= bp).nonzero()[0]
            above = bisect_left(D,bp)
            tp_above = bisect_left(TP,bp)
            fp_above = bisect_left(FP,bp)
            
            kept = total - above
            tp = tp0 - tp_above
            fp = fp0 - fp_above
            
            fdr = self._fdr(tp,fp)
            #fdr = self.fdr(above)
            
            #print "bp=%f above=%d data[above]=%f tp_above=%d fp_above=%d" % (bp,above,D[above],tp_above,fp_above)
            #print "kept=%d tp=%d fp=%d fdr=%f" % (kept,tp,fp,fdr)
            
            
            if fdr < max_fdr:
                loss = 1. - kept/float(self.N)
                logger.debug("KEEP %d with cutoff %s >= %f with FDR=%.2f%% loss=%.2f%% in %.3f sec." % (kept,column,bp,fdr*100,loss*100,time()-t0))
                # as this is monotonous in 1d, we can safely return now, solutions can only get worse.
                self.kept_by_method[(column,)] = (I[above:],bp,fdr,kept,loss)
                return bp,fdr,kept,loss
            

        logger.debug("exhausted all %d cutoffs on '%s' without satisfying FDR limit in %.3f sec." % (len(breakpoints),column,time()-t0))


    def find_cutoff_combination(self,columns,max_fdr=options.max_fdr,max_cutoffs=[]):
        """
        search a cutoff combination on the given columns (scores) which
        deplete FP more than TP and thus satisfy max_fdr.
        """
        if tuple(columns) in self.kept_by_method:
            return self.kept_by_method[tuple(columns)][1:]
        
        t0 = time()
        columns_str = ",".join(columns)
        from itertools import izip
        from bisect import bisect_left

        # the actual values in the columns to be searched for cutoffs
        data = array([self.data.get(col) for col in columns])
        
        # only the values of real/decoy hits
        data_tp = array([d.take(self.tp_index) for d in data])
        data_fp = array([d.take(self.fp_index) for d in data])
            
        #print data_tp.shape,data_fp.shape

        # indices in sorted order (ascending)
        I = array([d.argsort() for d in data])
        I_tp = array([d_tp.argsort() for d_tp in data_tp])
        I_fp = array([d_fp.argsort() for d_fp in data_fp])

        #print "I_tp,I_fp",I_tp.shape,I_fp.shape
        
        # the values in sorted order to bisect the breakpoints
        D = array([d.take(i) for d,i in izip(data,I)])
        TP = array([d.take(i) for d,i in izip(data_tp,I_tp)])
        FP = array([d.take(i) for d,i in izip(data_fp,I_fp)])

        #print TP.shape,FP.shape
        
        total = self.N
        tp0 = len(self.tp_index) #array([len(d_tp) for d_tp in data_tp])
        fp0 = len(self.fp_index) #array([len(d_fp) for d_fp in data_fp])

        if len(max_cutoffs):
            #print "## columns=%s max_cutoffs = %s" % (columns_str,str(max_cutoffs))
            breakpoints = []
            for d,maxcut in izip(data,max_cutoffs):
                breakpoints.append(sorted(set([x for x in d if x <= maxcut])))
            #print "breakpoints below cutoff",breakpoints
        else:
            breakpoints = [sorted(set(d)) for d in bp_candidates]
        
        N_comb = product(array([len(bp) for bp in breakpoints]))
        data = data.transpose()
       
        def drop_below(I_all,_I,above):
            I_left = I_all.copy()
            #print "I_left, _I, above:",len(I_all),_I.shape,above
            #print sorted(I_left)
            for i,start in izip(_I,above):
                # dicard the indices that belong to data BELOW the threshold
                #print "dropping",sorted(set(i[:start]))
                I_left -= set(i[:start])
            #print "left",len(I_left)
            return I_left

        all_ind = set(arange(self.N))
        all_ind_tp = set(arange(tp0))
        all_ind_fp = set(arange(fp0))

        logger.debug("COMB: searching best cutoff on '%s' (%d combinations) to satisfy FDR <= %.1f %% (setup in %.3f sec.)" % (columns_str,N_comb,max_fdr*100,time()-t0))

        best_kept = 0
        best = None


        for bp in itertools.product(*breakpoints):
            above = [bisect_left(d,b) for d,b in izip(D,bp)]
            #print above
            tp_above = [bisect_left(t,b) for t,b in izip(TP,bp)]
            #print tp_above
            fp_above = [bisect_left(f,b) for f,b in izip(FP,bp)]
            #print fp_above    

            #print "I_tp,I_fp",I_tp.shape,I_fp.shape
            #print "before drop",len(all_ind),len(all_ind_tp),len(all_ind_fp)
            #print "all"
            I_left = drop_below(all_ind,I,above)
            #print "tp"
            I_tp_left = drop_below(all_ind_tp,I_tp,tp_above)
            #print "fp"
            I_fp_left = drop_below(all_ind_fp,I_fp,fp_above)
                            
            kept = len(I_left)
            tp = len(I_tp_left)
            fp = len(I_fp_left)
                
            fdr = self._fdr(tp,fp)

            #print "bp=%s" % (",".join([str(b) for b in bp]))
            #print "kept=%d tp=%d fp=%d fdr=%f" % (kept,tp,fp,fdr)
                
            if fdr < max_fdr:
                loss = 1. - kept/float(self.N)
                if kept > best_kept:
                    best = (kept,bp,fdr,loss,I_left)
                    best_kept = kept

        if best:
            kept,bp,fdr,loss,I_left = best
            logger.debug("COMB: KEEP %d with best cutoff combination (%s) >= (%s) with FDR=%.2f%% loss=%.2f%% in %.3f sec." % (kept,columns_str,",".join([str(b) for b in bp]),fdr*100,loss*100,time()-t0))
            self.kept_by_method[tuple(columns)] = (array(sorted(I_left)),bp,fdr,kept,loss)
            return bp,fdr,kept,loss
        else:
            logger.debug("COMB: exhausted all %d cutoffs on '%s' without satisfying FDR limit in %.3f sec." % (N_comb,columns_str,time()-t0))



    def find_cutoff_pair(self,columns,max_fdr=options.max_fdr,max_cutoffs=[],min_kept=0):
        """
        search a cutoff combination on the two given columns (scores)
        which deplete FP more than TP and thus satisfy max_fdr.
        """
        
        if tuple(columns) in self.kept_by_method:
            return self.kept_by_method[tuple(columns)][1:]
        
        t0 = time()
        columns_str = ",".join(columns)
        from itertools import izip
        from bisect import bisect_left

        # the actual values in the columns to be searched for cutoffs
        data = array([self.data.get(col) for col in columns])
        
        # only the values of real/decoy hits
        data_tp = array([d.take(self.tp_index) for d in data])
        data_fp = array([d.take(self.fp_index) for d in data])
            
        #print data_tp.shape,data_fp.shape

        # indices in sorted order (ascending)
        I = array([d.argsort() for d in data])
        I_tp = array([d_tp.argsort() for d_tp in data_tp])
        I_fp = array([d_fp.argsort() for d_fp in data_fp])

        #print "I_tp,I_fp",I_tp.shape,I_fp.shape
        
        # the values in sorted order to bisect the breakpoints
        D = array([d.take(i) for d,i in izip(data,I)])
        TP = array([d.take(i) for d,i in izip(data_tp,I_tp)])
        FP = array([d.take(i) for d,i in izip(data_fp,I_fp)])

        #print TP.shape,FP.shape
        
        total = self.N
        tp0 = len(self.tp_index) #array([len(d_tp) for d_tp in data_tp])
        fp0 = len(self.fp_index) #array([len(d_fp) for d_fp in data_fp])

        if len(max_cutoffs):
            #print "## columns=%s max_cutoffs = %s" % (columns_str,str(max_cutoffs))
            breakpoints = []
            for d,maxcut in izip(data,max_cutoffs):
                breakpoints.append(sorted(set([x for x in d if x <= maxcut])))
            #print "breakpoints below cutoff",breakpoints
        else:
            breakpoints = [sorted(set(d)) for d in bp_candidates]
        
        N_comb = product(array([len(bp) for bp in breakpoints]))
        data = data.transpose()
       
        all_ind = set(arange(self.N))
        all_ind_tp = set(arange(tp0))
        all_ind_fp = set(arange(fp0))

        logger.debug("COMB2: searching best cutoff on '%s' (%d combinations) to satisfy FDR <= %.1f %% (setup in %.3f sec.)" % (columns_str,N_comb,max_fdr*100,time()-t0))

        best_kept = 0
        best = None

        all_ind0 = all_ind.copy()
        all_ind0_tp = all_ind_tp.copy()
        all_ind0_fp = all_ind_fp.copy()

        above0_last = 0
        above0_last_tp = 0
        above0_last_fp = 0

        for b0 in breakpoints[0]:
            above0 = bisect_left(D[0],b0)
            if above0 == above0_last:
                continue
            above0_tp = bisect_left(TP[0],b0)
            above0_fp = bisect_left(FP[0],b0)

            newly_lost0 = I[0][above0_last:above0]
            newly_lost0_tp = I_tp[0][above0_last_tp:above0_tp]
            newly_lost0_fp = I_fp[0][above0_last_fp:above0_fp]

            all_ind0 -= set(newly_lost0)
            all_ind0_tp -= set(newly_lost0_tp)
            all_ind0_fp -= set(newly_lost0_fp)
            
            kept0 = len(all_ind0)
            if kept0 < min_kept:
                logger.debug("Aborting search in this column because we are already below min_kept")
                break
           
            above0_last = above0
            above0_last_tp = above0_tp
            above0_last_fp = above0_fp

            all_ind1 = all_ind0.copy()
            all_ind1_tp = all_ind0_tp.copy()
            all_ind1_fp = all_ind0_fp.copy()

            above1_last = 0
            above1_last_tp = 0
            above1_last_fp = 0

            for b1 in breakpoints[1]:
                bp = (b0,b1)
                above1 = bisect_left(D[1],b1)
                if above1 == above1_last:
                    continue

                above1_tp = bisect_left(TP[1],b1)
                above1_fp = bisect_left(FP[1],b1)

                newly_lost1 = I[1][above1_last:above1]
                newly_lost1_tp = I_tp[1][above1_last_tp:above1_tp]
                newly_lost1_fp = I_fp[1][above1_last_fp:above1_fp]

                all_ind1 -= set(newly_lost1)
                all_ind1_tp -= set(newly_lost1_tp)
                all_ind1_fp -= set(newly_lost1_fp)
                
                above1_last = above1
                above1_last_tp = above1_tp
                above1_last_fp = above1_fp
                            
                kept = len(all_ind1)
                if kept < min_kept:
                    logger.debug("Aborting search in this row because we are already below min_kept")
                    break

                tp = len(all_ind1_tp)
                fp = len(all_ind1_fp)
                    
                fdr = self._fdr(tp,fp)

                #print "bp=%s" % (",".join([str(b) for b in bp]))
                #print "kept=%d tp=%d fp=%d fdr=%f" % (kept,tp,fp,fdr)
                    
                if fdr < max_fdr:
                    loss = 1. - kept/float(self.N)
                    if kept > best_kept:
                        best = (kept,bp,fdr,loss,all_ind1.copy())
                        best_kept = kept

        if best:
            kept,bp,fdr,loss,I_left = best
            logger.debug("COMB2: KEEP %d with best cutoff combination (%s) >= (%s) with FDR=%.2f%% loss=%.2f%% in %.3f sec." % (kept,columns_str,",".join([str(b) for b in bp]),fdr*100,loss*100,time()-t0))
            self.kept_by_method[tuple(columns)] = (array(sorted(I_left)),bp,fdr,kept,loss)
            return bp,fdr,kept,loss
        else:
            logger.debug("COMB2: exhausted all %d cutoffs on '%s' without satisfying FDR limit in %.3f sec." % (N_comb,columns_str,time()-t0))



    def filter_by_method(self,source,method):
        kept_indices,bp,fdr,kept,loss = self.kept_by_method[method]
        kept_indices = set(kept_indices)
        #print method,kept_indices
        N = defaultdict(int)

        for i,gff in enumerate(source):
            N["total"] += 1
            gff = list(gff)
            gff[8] = gff[8].rstrip() + ' filtered_by="%s"; est_FDR="%.3f";' % (",".join(method),5)

            out = "\t".join([str(g) for g in gff])
            if i in kept_indices:
                if i in self.fp_index:
                    N['decoy'] += 1
                    sys.stderr.write(out+'\n')
                else:
                    N['kept'] += 1
                    print out
            else:
                N["low_score"] += 1
        
        self.filter_stats = N

scoring_functions = {
    'signature' : lambda s : s.n_signature,
    'sig_density' : lambda s : float(s.n_signature)/float(s.n_reads),
    'sig_density_pseudo' : lambda s : float(0.04+s.n_signature+1)/float(s.n_reads+2+0.04),
    'length' : lambda s : s.length,
    'entropy' : lambda s : s.total_entropy,
    '4su' : lambda s : s.edits[13]+s.edits[15],
    '6sg' : lambda s : s.edits[8]+s.edits[11],
    '4su_density' : lambda s : (s.edits[13]+s.edits[15])/float(s.n_reads),
    '6sg_density' : lambda s : (s.edits[8]+s.edits[11])/float(s.n_reads),
    'indels' : lambda s : s.edits[3]+s.edits[7]+s.edits[11]+s.edits[15]+s.edits[16]+s.edits[17]+s.edits[18]+s.edits[19],
    'deletions' : lambda s : s.edits[3]+s.edits[7]+s.edits[11]+s.edits[15],
    'all_edits' : lambda s: array(s.edits)[:20].sum(),
    'n_best' : lambda s : s._asdict().get('n_best',0),
    'n_reads' : lambda s : s._asdict().get('n_reads',0),
    'n_uniq' : lambda s : s._asdict().get('n_uniq',0),
    'uniqness' : lambda s : s._asdict().get('uniqness',0),
    'cum_uniq' : lambda s : s._asdict().get('cum_uniq',0),
    'avg_mapq' : lambda s : s._asdict().get('avg_mapq',0),
    'best_mapq' : lambda s : s._asdict().get('best_mapq',0),
    'mean_pp' : lambda s : float(s._asdict().get('mean_pp',0)),
    'median_pp' : lambda s : float(s._asdict().get('median_pp',0)),
    'max_pp' : lambda s : float(s._asdict().get('max_pp',0)),
}

T0 = time()
type_hints = dict(cid='str',sense='str',map_qual='str',edits='intlist')
from byo.io.lazytables import NamedTupleImporter

CS = ClusterSet(NamedTupleImporter(args[0],type_hints=type_hints,parse_comments=True,default_cast='float'))
unfiltered_fdr = CS.unfiltered_fdr 
logger.info("unfiltered FDR=%.2f %%" % (unfiltered_fdr*100))
            
if unfiltered_fdr <= options.max_fdr:
    logger.info("no filtering needed!")
    kept,method,fdr,loss,bp = CS.N,("unfiltered",),unfiltered_fdr,0,0
else:
    candidates = []
    maxcutoffs = {}
    min_kept = 0
    for name,func in scoring_functions.items():
        CS.score(name,func)
        res = CS.find_cutoff(name)
        if res:
            bp,fdr,kept,loss = res
            maxcutoffs[name] = bp
            min_kept = max(min_kept,kept)
            candidates.append((kept,1,(name,),fdr,loss,bp))

    for base_func in ["n_reads","n_uniq","signature"]:
        for supp_func in scoring_functions.keys():
            if supp_func == base_func:
                continue

            cutoff_bound = [
                maxcutoffs.get(base_func,CS.data.get(base_func).max()),
                maxcutoffs.get(supp_func,CS.data.get(supp_func).max())
            ]
            res = CS.find_cutoff_pair([base_func,supp_func],max_cutoffs=array(cutoff_bound),min_kept=min_kept)
            if res:
                bp,fdr,kept,loss = res
                min_kept = max(min_kept,kept)
                candidates.append((kept,.5,(base_func,supp_func),fdr,loss,bp))

        if options.combinations:
            for name1,name2 in itertools.combinations(sorted(scoring_functions.keys()),2):
                cutoff_bound = [maxcutoffs.get(name1,CS.data.get(name1).max()),maxcutoffs.get(name2,CS.data.get(name2).max())]
                res = CS.find_cutoff_pair([name1,name2],max_cutoffs=array(cutoff_bound))
                if res:
                    bp,fdr,kept,loss = res
                    candidates.append((kept,.5,(name1,name2),fdr,loss,bp))
            
    if not candidates:
        logger.warning("no method was able to satisfy FDR limit, bailing out without clusters")
        kept,method,fdr,loss,bp = 0,("none",),unfiltered_fdr,1.,0
    else:
        ### rank all filterings by how many clusters they keep
        ranked = sorted(candidates,reverse=True)
        logger.info("FDR filtering completed in %d sec. with the following candidate scoring functions:" % (time() -T0))
        for kept,x,method,fdr,loss,bp in ranked:
            logger.info("%30s >= %30s : keeps %d clusters (loss=%.3f %%) at FDR <= %.2f %%" % (str(method),str(bp),kept,loss*100,fdr*100))

        ### get best solution
        kept,x,method,fdr,loss,bp = ranked[0]
        logger.info("selected method: %s" % str(method))
        

from byo.io.gff import gff_importer,dict_from_attrstr
N = defaultdict(int)
CS.filter_by_method(gff_importer(args[1]),method)

if options.statfile:
    N = CS.filter_stats

    from byo.parclip.statistics import Stats
    stats = Stats(path=options.statfile)
    d = dict(raw_clusters=N['total'],kept=N['kept'],fp_remaining=N['decoy'],FDR=fdr,raw_FDR=unfiltered_fdr,loss=100*loss)
    d["FDR_%s" % ",".join(method)] = fdr

    stats.save_dict(title="cluster diag",stats=d)
    stats.flush()
        
