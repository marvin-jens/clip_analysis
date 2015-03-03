import numpy as np
from collections import defaultdict

### Helper functions
def entropy_of_pos(pos,counts):
    data = defaultdict(int)
    for x,c in zip(pos,counts):
        data[x] += c

    p = np.array(data.values(),dtype=np.float32)
    p /= p.sum()

    return - (p * np.log2(p)).sum()


from editedcluster import EditedCluster

class ScoredCluster(EditedCluster):

    def __init__(self,name,reads,strand,chrom,**kwargs):
        super(ScoredCluster,self).__init__(name,reads,strand,chrom,**kwargs)

        self.mapq = []
        self.uniqness = []
        self.pp = []
        self.reported_scores = [
            "uniq_reads","n_reads","n_conv_pos","length",
            "max_uniq","avg_uniq","cum_uniq",
            "max_mapq","avg_mapq","cum_mapq",
            "max_pp","avg_pp","cum_pp",
            "score","signature_density","iclip_score","profile_KL","seq_entropy","entropy_score"
        ]
        self._map_scored = False

        # current behaviour is to switch to iCLIP 
        # as soon as one library is iCLIP
        if 'iclip' in self.signature.values():
            self._scorefuncname = "iclip_score"
        else:
            self._scorefuncname = "conversion_score"
                
    def parse_map_scores(self):
        if not self._parsed: self.parse_reads()
        for r,w in zip(self.reads,self.counts):
            # record MAPQ mapping quality (0..255)
            self.mapq.append(r.mapq)

            # look at additional mapper output for cues on alignment quality
            tags = dict(r.tags)
            uniqness = 0

            # BWA output
            # get number of best alignments and the hull 
            # (number of alignments at k+1)
            n_best = tags.get('X0',1)
            if 'X1' in tags:
                hull = tags['X1']
                uniqness = (n_best == 1) * 255. / (hull + 1.)
            
            # BOWTIE2 output
            # get the difference between best alignment score and second best
            elif 'AS' in tags:
                AS = tags.get('AS',0)
                XS = tags.get('XS',0)
                uniqness = abs(AS-XS)
            
            self.uniqness.append(uniqness)
            
            # BWA-PSSM output
            # take the posterior probability that this alignment is correct,
            # given the error-model (for example PAR-CLIP conversions)
            pp = tags.get('PP',0.5)
            self.pp.append(pp)
            
        self.uniqness = np.array(self.uniqness)
        self.mapq = np.array(self.mapq)
        self.pp = np.array(self.pp)
        
        self._map_scored = True

    @property
    def score(self):
        return getattr(self,self._scorefuncname)

    @property
    def conversion_score(self):
        if not self._parsed: self.parse_reads()
        return self.conversions.sum()

    @property
    def signature_density(self):
        return self.conversion_score / self.n_reads

    @property
    def entropy_5p(self):
        if not self._parsed: self.parse_reads()
        if self.strand == '-':
            return entropy_of_pos(self.ends,self.counts)
        else:
            return entropy_of_pos(self.starts,self.counts)
    @property
    def entropy_3p(self):
        if not self._parsed: self.parse_reads()
        if self.strand == '-':
            return entropy_of_pos(self.starts,self.counts)
        else:
            return entropy_of_pos(self.ends,self.counts)

    @property
    def iclip_score(self):
        # 1. for a single read, goes to zero the more reads we have
        default = 2. / (self.n_reads + 1.)
        
        # low 5' entropy, even for high numbers of reads is good.
        # high 3' entropy on the other hand is good. The proposed score
        # is 1 if both entropies are the same and goes to 0 if 3' end
        # variability is lower than 5' end, and goes to infinity for
        # 5' entropy close to zero, as read-counts increase.
        return (self.entropy_3p + default) / (self.entropy_5p + default)

    @property
    def n_reads(self):
        if not self._parsed: self.parse_reads()
        return self.counts.sum()
    
    @property
    def uniq_reads(self):
        return len(self.reads)
    
    @property
    def cum_uniq(self):
        if not self._map_scored: self.parse_map_scores()
        return (self.counts * self.uniqness).sum()
    @property
    def max_uniq(self):
        if not self._map_scored: self.parse_map_scores()
        return self.uniqness.max()
    @property
    def avg_uniq(self):
        if not self._map_scored: self.parse_map_scores()
        return (self.counts * self.uniqness).sum() / self.counts.sum()



    @property
    def cum_mapq(self):
        if not self._map_scored: self.parse_map_scores()
        return (self.counts * self.mapq).sum()
    @property
    def max_mapq(self):
        if not self._map_scored: self.parse_map_scores()
        return self.mapq.max()
    @property
    def avg_mapq(self):
        if not self._map_scored: self.parse_map_scores()
        return (self.counts * self.mapq).sum() / self.counts.sum()



    @property
    def cum_pp(self):
        if not self._map_scored: self.parse_map_scores()
        return (self.counts * self.pp).sum()
    @property
    def max_pp(self):
        if not self._map_scored: self.parse_map_scores()
        return self.pp.max()
    @property
    def avg_pp(self):
        if not self._map_scored: self.parse_map_scores()
        return (self.counts * self.pp).sum() / self.counts.sum()

        
    @property
    def profile_KL(self):
        """
        KL-divergence between observed profile and uniform coverage. 
        """
        if not self._parsed: self.parse_reads()
        
        p = self.coverage / float(self.coverage.sum())
        return np.where(p> 0,p * np.log2(self.L*p),0).sum()
    @property
    def seq_entropy(self):
        """
        Shannon entropy of the nucleotide observations (0 for no mismatches)
        """
        if not self._parsed: self.parse_reads()
        
        p = np.array(self.nuc_matrix,dtype=np.float32) / self.nuc_matrix.sum(axis=1)[:,np.newaxis]
        return np.where(p> 0,-p * np.log2(p),0).sum()

    @property
    def entropy_score(self):
        scale = (1 - (10./(10.+self.n_reads)))*.1
        return self.profile_KL + scale*self.seq_entropy
        
    @property
    def short(self):
        return "%s%s %d:%d " % (self.chrom,self.strand,self.start,self.end)

    @property
    def n_conv_pos(self):
        if not self._parsed: self.parse_reads()
        return (self.conversions > 0).sum()
