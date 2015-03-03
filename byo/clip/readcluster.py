class ReadCluster(object):
    """
    The baseclass of all read clusters. Represents a pile of reads
    that contiguously cover the reference. Does nothing much except 
    determining its start and end coordinates. __str__ produces BED6
    which allows for print() to do something useful already.
    """
    def __init__(self,name,reads,strand,chrom,**kwargs):
        self.name = name
        self.reads = reads
        self.strand = strand
        self.chrom = chrom
        self.kwargs = kwargs
        
        self.starts = [r.pos for r in reads]
        self.ends = [r.pos + r.alen for r in reads]

        # reads are sorted by start-position
        self.start = self.starts[0]
        self.end = max(self.ends)
        self.L = self.end-self.start
        self.length = self.L

    @property
    def score(self):
        return len(self.reads)

    def __str__(self):
        """
        gives a BED6 compatible representation of the read cluster
        """
        return "\t".join([self.chrom,str(self.start),str(self.end),self.name,str(self.score),self.strand])

