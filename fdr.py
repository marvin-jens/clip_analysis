#!/usr/bin/env python
import os,sys
from time import time
import itertools
from numpy import *
from optparse import *


seterr(all='ignore')
usage = """
  %prog [options] <cluster_stats.tsv> <all_clusters.gff> 

Analyze cluster statistics and estimate false discovery rates based on sense/antisense distributions.
"""

parser = OptionParser(usage=usage)
parser.add_option("-l","--logfile",dest="logfile",default="fdr.log",help="logfile to use")
parser.add_option("-s","--statfile",dest="statfile",default="",help="write statistics to file")
parser.add_option("-c","--cutoff",dest="max_fdr",default=.05,type=float,help="maximally tolerated FDR (default=.05)")
parser.add_option("-d","--decoy",dest="decoy_scale",default=1.,type=float,help="efficiency of decoy compared to true positives (1 for random genome, .5 for antisense)")
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
        (FP=decoy hits) more strongly than true positives (TP=transcript 
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
            
            self.tp_index = (flags == 'transcript').nonzero()[0]
            self.fp_index = (flags == 'antisense').nonzero()[0]

            self.unfiltered_fdr = self.fdr()
        
        logger.debug("loaded %d cluster stats in %.3f sec. tp=%d fp=%d" % (self.N,time()-t0,len(self.tp_index),len(self.fp_index)))
        
        #bp,fdr,kept,loss
        self.kept_by_method = {('unfiltered',):(arange(self.N),None,self.unfiltered_fdr,self.N,0), ('none',):([],-1,1,1)}

    def fdr(self,indices=set()):
        """
        estimate the false discovery rate for (a subset) of a cluster set
        based on annotation as transcript aligning or decoy hit.
        """
        if not self.N:
            return 0.5

        if not len(indices):
            indices = arange(self.N)
            
        flags = take(self.data['sense'],indices)
        
        ## assume all antisense are fp and correct for decoy efficiency
        #fp = (flags == 'antisense').sum() / options.decoy_scale
        
        ## assume as many fp on sense than on antisense
        #tp = max(0,(flags == 'transcript').sum() - fp)
        #fp *=2

        ##f = (fp + 1.) / (tp + fp + 2.)
        #f = (fp + 1.) / (tp + fp + 2.)
        return self._fdr((flags == 'transcript').sum(),(flags == 'antisense').sum())

    @staticmethod
    def _fdr(real,decoy):
        #print "_FDR",real,decoy
        fp = decoy / options.decoy_scale
        
        # assume as many fp on real than on (scaled) decoy
        tp = max(0,real - fp)
        fp *=2

        #f = (fp + 1.) / (tp + fp + 2.)
        f = (fp + 1.) / (tp + fp + 2.)
        
        return f
        
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
                    N['antisense'] += 1
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
    'mapqual' : lambda s : {'UNIQ40_K0' : 10,'UNIQ10_K0': 9,'UNIQ1_K0': 8,'UNIQ40_K1' : 7,'UNIQ10_K1' : 6,'UNIQ1_K1' : 5, 'UNIQ40_K2' : 4,'UNIQ10_K2' : 3,'UNIQ1_K2' : 2, 'DUBIOUS' : 1}[s._asdict().get('map_qual')]
    #'fail' : lambda s : 0,
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
    d = dict(raw_clusters=N['total'],kept=N['kept'],fp_remaining=N['antisense'],FDR=fdr,raw_FDR=unfiltered_fdr,loss=100*loss)
    d["FDR_%s" % ",".join(method)] = fdr

    stats.save_dict(title="cluster diag",stats=d)
    stats.flush()
    # here comes the really ugly part. Make the tables with all the statistics.
    
    #f = file(options.statfile,"a")
    #f.write(">FDR (false discovery rate filtering)\n# mapping quality class\tNO FILTER\t%s\n## --figsize=15,10 --labelspace=.2\n" % "\t".join(top_scorings))
    #f.write("ALL\t%.3f" % (raw_fdr*100))

    #for left,pri,name,tp_left,fp_left,ot_left,avg_fdr,raw_fdr in results[:top]:
        #logger.info("-> %d (out of %d) clusters left after filtering by '%s', reducing FDR from %.3f %% to %.3f %%" % (left,N['total'],name,raw_fdr*100,avg_fdr*100))
        #f.write("\t%.3f" % (avg_fdr*100))

    #f.write("\n")

    #for map_qual in ['NOHULL_K0','NOHULL_K1','NOHULL_K2','HULL_K0','HULL_K1','HULL_K2']:
        #if not map_qual in stats_by_map_qual:
            #continue

        #raw_fdr = score_performance[top_scorings[0]][map_qual][1][0]*100
        #f.write("%s\t%.3f\t%s" % (map_qual,raw_fdr,"\t".join(["%.3f" % (score_performance[n][map_qual][10]*100) for n in top_scorings])))
        #f.write("\n")

    #f.write(">Number of clusters\n# mapping quality class\tNO FILTER\t%s\n## --figsize=15,10 --labelspace=.2\n" % "\t".join(top_scorings))
    #f.write("ALL\t%.3f" % N['total'])
    #for left,pri,name,tp_left,fp_left,ot_left,avg_fdr,raw_fdr in results[:top]:
        #f.write("\t%.3f" % left)

    #f.write("\n")

    #for map_qual in ['NOHULL_K0','NOHULL_K1','NOHULL_K2','HULL_K0','HULL_K1','HULL_K2']:
        #if not map_qual in stats_by_map_qual:
            #continue

        #total_at_qual = int(np.array(score_performance[top_scorings[0]][map_qual][4:7]).sum())
        #f.write("%s\t%.3f\t%s" % (map_qual,total_at_qual,"\t".join([str(np.array(score_performance[n][map_qual][7:10]).sum()) for n in top_scorings])))
        #f.write("\n")
    
  
  

        
sys.exit(0)
import bisect

def compute_fdr_loss(tp_sig,fp_sig):
    
    fp = len(fp_sig)
    tp = len(tp_sig)

    breakpoints = sorted([0,]+list(set(fp_sig) | set(tp_sig)))

    fdr = []
    loss = []
    for b in breakpoints:
        f = fp - bisect.bisect_left(fp_sig,b)
        l = bisect.bisect_left(tp_sig,b)
        t = tp - l

        if tp > 0:\
            LLL = float(l)/tp
        else:
            LLL = 1.
        logger.debug("b=%.3f fp=%d / %d tp=%d / %d -> fdr=%.3f loss=%.3f" % (b,f,fp,t,tp,float(f+1)/(f+t+2),LLL))

        fdr.append(float(f+1)/(f+t+2))
        loss.append(LLL)

    return breakpoints,np.array(fdr),np.array(loss)


score_performance = defaultdict(dict)

colors = {
    'NOHULL_K0' : 'k',
    'NOHULL_K1' : 'g',
    'NOHULL_K2' : 'y',
    'HULL_K0' : 'c',
    'HULL_K1' : 'm',
    'HULL_K2' : 'r',
}

if options.pdf:
    import matplotlib
    matplotlib.use('pdf')
    from pylab import figure,title,plot,axhline,axvline,legend,ylim,text,twinx,ylabel,xlabel

#M = []
#M_names = [s.cid for s in raw_data]
#M_flags = array([s.sense for s in raw_data])

#for name in sorted(scoring_functions.keys()):
    #func = scoring_functions[name]

    #M.append([func(s) for s in raw_data])

#M = transpose(array(M))
#import byo.parclip.pca as pca
#pca.Center(M)



#P = pca.PCA(M)
#P.npc = 2
#res = transpose(P.pc())

#def print_component(C):
    #c_abs = abs(C)
    #print c_abs
    #ind = c_abs.argsort()[::-1]
    #print ind
    
    #bla = zip(C,sorted(scoring_functions.keys()))
    #print bla
    
    #for c,name in [bla[i] for i in ind]:
        #print name,"%.4f" % c
        
#print res.shape
#print "first component"
#print_component(P.Vt[0])

#print "second component"
#print_component(P.Vt[1])

#x,y = res

#from pylab import *

#i_tx = (M_flags == 'transcript').nonzero()[0]
#i_decoy = (M_flags == 'antisense').nonzero()[0] 
#print M_flags
#print i_tx
#scatter(x,y,color='gray')
#scatter(x.take(i_tx),y.take(i_tx),s=5,color='blue')
#scatter(x.take(i_decoy),y.take(i_decoy),s=5,color='red')

#show()

#sys.exit(1)

for name,func in scoring_functions.items():
    if name == 'fail':
        continue
    #print "-> scoring function",name

    # PLOTTING
    if options.pdf:
        figure(figsize=(15,10))
        ylim(0,100)    
        ylabel("loss [%]")
        xlabel("cluster score: '%s'" % name)
        ax2 = twinx()
        ylabel("FDR [%]")
        
        maxcut = 0
    
    works = False
    
    for map_qual in sorted(stats_by_map_qual.keys()):
        stats = stats_by_map_qual[map_qual]
        
        fp = len(stats['antisense'])
        tp = len(stats['transcript'])
        ot = len(stats['overlap_tx']) + len(stats['intergenic'])
    
        fp_sig = sorted([func(s) for s in stats['antisense']])
        tp_sig = sorted([func(s) for s in stats['transcript']])
        ot_sig = sorted([func(s) for s in stats['overlap_tx']]+[func(s) for s in stats['intergenic']])
    
        breakpoints,fdr,loss = compute_fdr_loss(tp_sig,fp_sig)
        logger.debug("map_qual %s (fp=%d tp=%d) -> unfiltered FDR=%.2f %%" % (map_qual,fp,tp,fdr[0]*100))
        
        cutoff_cand = (fdr < options.max_fdr).nonzero()[0]
        if len(cutoff_cand):
            cutoff_i = cutoff_cand[0]
            cutoff = breakpoints[cutoff_i]

            t_left = tp - bisect.bisect_left(tp_sig,cutoff)
            f_left = fp - bisect.bisect_left(fp_sig,cutoff)
            o_left = ot - bisect.bisect_left(ot_sig,cutoff)

            select_fdr = fdr[cutoff_i]
            logger.debug("%s FDR(%.3f) = %.2f %%  (< %.3f%%) for %d fp and %d tp (loss=%.2f %%) %d other (%d total)" % (name,cutoff,fdr[cutoff_i] * 100,options.max_fdr,f_left,t_left,loss[cutoff_i]*100,o_left,ot))

            #tp_left += t_left
            works = True
            
        else:
            logger.warning("method '%s' is unable to satisfy FDR constraint in '%s'" % (name,map_qual))
            t_left = 0
            f_left = 0
            o_left = 0
            cutoff = breakpoints[len(breakpoints)-1]+1
            select_fdr = 0

        score_performance[name][map_qual] = (breakpoints,fdr,loss,cutoff,tp,fp,ot,t_left,f_left,o_left,select_fdr)
        #print cutoff_i,cutoff
        
        if options.pdf:
            ax2.plot(breakpoints,100*fdr,colors[map_qual]+'-', drawstyle='steps-mid',label="%s FDR" % map_qual)#,linestyle='step')
            plot(breakpoints,25*loss,colors[map_qual]+'--', drawstyle='steps-post',label="%s loss" % map_qual)#,linestyle='step')
            #print array(100*loss,dtype=int)
            
            if cutoff < breakpoints[-1]:
                axvline(cutoff,color=colors[map_qual],linewidth=.1)
                maxcut = min(max(cutoff*2,maxcut),breakpoints[-1])
        
            if select_fdr == 0:
                maxcut = breakpoints[-1]

    if options.pdf:
        base = os.path.basename(args[0]).split(".cluster_stats.tsv")[0]        
        title(base)
        oc = options.max_fdr*100
        ax2.axhline(oc,color="gray",linestyle='dashed',linewidth=.1)
        #text("FDR-cutoff",-1,oc)
        
        legend(loc='upper-left')
        ax2.set_ylim((0,25))
        xlim(0,maxcut)
        #show()
        #if works:
        from byo.io import ensure_path
        savefig(os.path.join(ensure_path(options.pdf),"FDR_%s_%s.pdf" % (name,base)))
        # END PLOTTING

results = []
#sys.exit(1)
# use as additional cues when sorting results. Only used to break ties in a non-confusing way!
priorities = {
    'signature' : 2,
    'signature_density' : 2,
    '4su_density' : 1,
    '6sg_density' : 1,
    '4su' : 0,
    '6sg' : 0,
}

for name,d in score_performance.items():
    tp_left = 0
    fp_left = 0
    ot_left = 0
    avg_fdr = 0.0
    raw_fdr = 0.0
    total = 0
    total_tp_fp = 0
    
    for map_qual,(breakpoints,fdr,loss,cutoff,tp,fp,ot,t_left,f_left,o_left,select_fdr) in d.items():
        tp_left += t_left
        fp_left += f_left
        ot_left += o_left
        
        total += fp+tp+ot
        total_tp_fp += fp+tp

        avg_fdr += select_fdr * (t_left+f_left)
        raw_fdr += fdr[0] * (tp+fp)

    left = tp_left+ot_left

    if (tp_left+fp_left):
        avg_fdr /= (tp_left+fp_left)
    else:
        avg_fdr = 0
        
    raw_fdr /= total_tp_fp
    results.append((left,priorities.get(name,-1),name,tp_left,fp_left,ot_left,avg_fdr,raw_fdr))

    
results = sorted(results,reverse=True)
if not results:
    results = [(0,0,"fail",0,0,0,1.,1.)]
    total = .01
    
top = 3
top_scorings = [r[2] for r in results[:top]]

# now the actual filtering
best_scorer = top_scorings[0]
best_func = scoring_functions[best_scorer]

fp_above_cutoff = 0
from byo.io.gff import gff_importer,dict_from_attrstr
N = defaultdict(int)
for gff in gff_importer(args[1]):
    #print gff
    N["total"] += 1
    attrs = dict_from_attrstr(gff.attr_str)
    cid = attrs['name']

    if not cid in stats_by_id:
        logger.warning("trying to filter the wrong cluster set? unknown id '%s'" % cid)
        N["unknown ID"] += 1
        continue
    
    map_qual = attrs['map_qual']
    breakpoints,fdr,loss,cutoff,tp,fp,ot,t_left,f_left,o_left,select_fdr = score_performance[best_scorer][map_qual]
    
    score = best_func(stats_by_id[cid])

    i = bisect.bisect_left(breakpoints,score)
    if i < len(fdr):
        est_fdr = fdr[i]
    else:
        est_fdr = fdr[-1]
    #if est_fdr > options.max_fdr:
        #N["FDR_too_large"] += 1
    if score < cutoff:
        N["low_score"] += 1
        continue
    
    
    est_fdr *= 100
    gff = list(gff)
    gff[7] = est_fdr
    gff[8] = gff[8].rstrip() + ' filtered_by="%s"; est_FDR="%.3f";' % (best_scorer,est_fdr)

    out = "\t".join([str(g) for g in gff])
    if stats_by_id[cid].sense == "antisense":
        N['antisense'] += 1
        sys.stderr.write(out+'\n')
        fp_above_cutoff += 1
    else:
        N['kept'] += 1
        print out


