from codons import codon_table

def translate(seq,frame=0):
    started = False
    for i in xrange(frame,len(seq)-frame-4,3):
        cod = seq[i:i+3].upper()
        a = codon_table.get(cod,"")

        if started and not a:
            # translating and encountered stop codon
            break

        if not started and a == 'M':
            # scanning for start codon
            started = True

        if started and a:
            # within the ORF
            yield a

def translate_aa(seq,frame=0):
    return "".join(translate(seq,frame))

def all_orfs(seq):
    for frame in range(3):
        started = False
        start = 0
        orf = []

        for i in xrange(frame,len(seq)-frame-4,3):
            cod = seq[i:i+3].upper()
            a = codon_table.get(cod,"")

            if started and not a:
                # translating and encountered stop codon
                #print "* terminating",orf
                yield (start,i+3,"".join(orf))
                started = False
                orf = []

            if (not started) and a == 'M':
                # scanning for start codon
                started = True
                start = i
                orf = []
                #print "* starting at %d" % i,orf

            if started and a:
                # within the ORF
                orf.append(a)
                #print "* elongating",orf

def find_ORF(seq,len_thresh=30):
    """
    find all open-reading frames, sort by length (if above threshold) and return the longest one.
    """
    orfs = sorted([ (len(orf),orf,start,end) for (start,end,orf) in all_orfs(seq) if len(orf) >= len_thresh],reverse=True)
    if not orfs:
        return "",-1,-1
    else:
        L,ORF,start,end = orfs[0]
        return ORF,start,end
