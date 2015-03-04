# -*- coding: utf-8 -*-

def k_mers(k,alphabet=['a','c','g','t',]):
    if k == 1:
        prefix = [""]
    else:
        prefix = k_mers(k-1,alphabet)

    for pre in prefix:
        for a in alphabet:
            yield pre+a

COMPLEMENT = {
    'a' : 't',
    't' : 'a',
    'c' : 'g',
    'g' : 'c',
    'k' : 'm',
    'm' : 'k',
    'r' : 'y',
    'y' : 'r',
    's' : 's',
    'w' : 'w',
    'b' : 'v',
    'v' : 'b',
    'h' : 'd',
    'd' : 'h',
    'n' : 'n',
    'A' : 'T',
    'T' : 'A',
    'C' : 'G',
    'G' : 'C',
    'K' : 'M',
    'M' : 'K',
    'R' : 'Y',
    'Y' : 'R',
    'S' : 'S',
    'W' : 'W',
    'B' : 'V',
    'V' : 'B',
    'H' : 'D',
    'D' : 'H',
    'N' : 'N',
}

def complement(s):
    return "".join([COMPLEMENT[x] for x in s])

def rev_comp(seq):
    return complement(seq)[::-1]

if __name__=='__main__':
    seq = 'attaCGtTTTTGCCGCTTAaaaaaaaaa'
    print seq
    print complement(seq)
    print rev_comp(seq)

