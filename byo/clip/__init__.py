if __name__ == "__main__":
    from byo.systems import hg19 as system
    def debug_output(cluster):
        #if cluster.strand == "+":
            #return
        
        def conv_str(conv):
            l = []
            for c in conv:
                if c:
                    l.append(".")
                else:
                    l.append(' ')
            return "".join(l)
        
        #ref = system.genome.get(cluster.chrom,cluster.start,cluster.end,cluster.strand)
        print cluster
        print ">>>>"
        #print ref
        print cluster.sequence
        print conv_str(cluster.conversions)
        print "n_reads    ", cluster.n_reads
        print "n_uniq     ", cluster.uniq_reads
        print "profile_KL ", cluster.profile_KL
        print "seq_entropy", cluster.seq_entropy
        print "5'- entropy", cluster.entropy_5p
        print "3'- entropy", cluster.entropy_3p
        print "iclip_score", cluster.iclip_score
        #print "read copies", cluster.counts
        #print cluster.coverage
        #print cluster.nuc_matrix
        #print "editstats dict", cluster.editstats
        
        
        print "-----"

    from bamsource import BAM_ClusterGenerator
    from scoredcluster import ScoredCluster
    #from multiprocessing import Pool
    
    print "testing CLIP/read cluster code"
    import logging
    logging.basicConfig(level=logging.DEBUG)
    
    def handle_cluster(cluster):
        if len(cluster.reads) < 2:
            return

        debug_output(cluster)

    #P = Pool()
    #for cluster in BAM_ClusterGenerator(bams=["rep1.bam","rep2.bam","rep3.bam"],stranded=True):
    #for cluster in BAM_ClusterGenerator(bams=["../../lebedeva2011/mapping/filtered.hg19.elavl1_6sg_SILAC.bam"],stranded=True,cluster_type=ScoredCluster,signature={'*' : "GA"}):
        #handle_cluster(cluster)

    for cluster in BAM_ClusterGenerator(bams=["/home/mjens/pcp/pariclip/mapping/filtered.hg19.nSR100_iCLIP_rep2.bam"],stranded=True,cluster_type=ScoredCluster,signature={'*' : "iclip"},count_func='lin'):
        handle_cluster(cluster)
        
    #P.map(handle_cluster,BAM_ClusterGenerator(bams=["../../lebedeva2011/mapping/filtered.hg19.elavl1_6sg_SILAC.bam"],stranded=True,cluster_type=ScoredCluster,signature={'*' : "GA"}))

