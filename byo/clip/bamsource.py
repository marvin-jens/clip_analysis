from readcluster import ReadCluster

class BAM_ClusterGenerator(object):
    """
    This class allows to iterate over alignments from multiple, sorted BAM files.
    Reads that contiguously cover the reference sequence are collected into
    clusters and yielded by the instances __iter__ function.
    Use it as a generator:
    
    for cluster in BAM_ClusterGenerator(bams=["rep1.bam","rep2.bam","rep3.bam"],stranded=True):
        print "found a cluster of %d reads on %s strand of chrom %s" % (len(cluster.reads),cluster.strand,cluster.chrom)
        print cluster

    each read gets an addtional tag 'LB', referring to the bam's basename 
    ('rep1','rep2','rep3' in the above example). Alternatively, this can be replaced by a 'group' a 
    bam belongs to. Example: 
    
        BAM_ClusterGenerator(bams=["sampleA_rep1","sampleA_rep2","sampleB_rep1","sampleB_rep2"],groups=["A","A","B","B"])
        
    Here, the LB tags will either be "A" or "B", and the replicates are treated as if they were one file.
    By default BAM_CLusterGenerator will yield objects of type ReadCluster (see clustertypes.py). Use the kwargs
    cluster_type=MyClusterClass to change this. Any additional kwarg will be passed on to each new cluster upon 
    instantiation.
    """
    def __init__(self,bams=[],groups=[],stranded=True,region="",cluster_type=ReadCluster,**kwargs):
        import pysam
        import logging
        import itertools
        import os

        self.logger = logging.getLogger("BAM_ClusterGenerator")
        self.region = region
        self.stranded = stranded
        self._plus = "+"
        if not stranded:
            self._plus = "*"

        self.bams = [pysam.Samfile(name,"rb") for name in bams]

        if groups:
            if len(groups) != len(bams):
                reason = "number of supplied library groups (%d) and BAM files (%d) does not match!" % (len(groups),len(bams))
                self.logger.error(reason)
                raise ValueError(reason)
            
            self.names = groups
        else:
            self.names = [os.path.basename(bam) for bam in bams]
            
        self.N_clusters = 0
        self.cluster_type = cluster_type
        self.kwargs = kwargs

    def fetch(self):
        """
        Iterates jointly over multiple sorted bams, tagging reads with their respective library-names
        """
        if self.region:
            iters = [bam.fetch(region=self.region).__iter__() for bam in self.bams]
        else:
            iters = [bam.fetch().__iter__() for bam in self.bams]

        curr = [i.next() for i in iters]
        names = self.names
        
        def head():
            pos = sorted([(read.tid,read.pos,i) for i,read in enumerate(curr)])
            return pos[0][2]

        N = 0
        from time import time
        t0 = time()
        
        while curr:
            i = head()
            read = curr[i]
            read.tags = read.tags + [('LB',names[i])]

            yield read
            
            N += 1
            if N and not (N % 100000):
                dt = time() - t0
                self.logger.info("processing %.2fk reads/s" % (N/(dt*1000.)))
                
            try:
                curr[i] = iters[i].next()
            except StopIteration:
                curr.pop(i)
                iters.pop(i)
                names.pop(i)

    def make_cluster(self,reads,strand,tid):
        self.N_clusters += 1
        return self.cluster_type("CID_%06d" % self.N_clusters,reads,strand,self.bams[0].getrname(tid),**self.kwargs)
        
    def __iter__(self):
        """
        Uses fetch() to collect reads that contiguously cover the reference. 
        Yields lists of such reads that form a cluster, if desired strands
        are treated separately.
        """

        lastref = -1
        lastpos_plus = 0
        lastpos_minus = 0
        reads_plus = []
        reads_minus = []

        # iterate over all reads
        for read in self.fetch():
            # continous coverage or new chrom marks a new cluster
            new_chrom = read.tid != lastref
            
            if self.stranded and read.is_reverse:
                if read.pos > lastpos_minus or new_chrom:
                    if reads_minus:
                        yield self.make_cluster(reads_minus,"-",lastref)
                    reads_minus = []

                    if new_chrom:
                        if reads_plus:
                            yield self.make_cluster(reads_plus,self._plus,lastref)
                        reads_plus = []
                        lastpos_plus = 0
                    
                lastpos_minus = max((1-new_chrom)*lastpos_minus,read.aend)
                reads_minus.append(read)
            else:
                if read.pos > lastpos_plus or new_chrom:
                    if reads_plus:
                        yield self.make_cluster(reads_plus,self._plus,lastref)
                    reads_plus = []
                    if new_chrom:
                        if reads_minus:
                            yield self.make_cluster(reads_minus,"-",lastref)
                        reads_minus = []
                        lastpos_minus = 0
                    
                lastpos_plus = max((1-new_chrom)*lastpos_plus,read.aend)
                reads_plus.append(read)

            lastref = read.tid

        if reads_plus:
            yield self.make_cluster(reads_plus,self._plus,lastref)

        if reads_minus:
            yield self.make_cluster(reads_minus,"-",lastref)

