#!/usr/bin/env python
# -*- coding: utf-8 -*-
from sequence_data.io import fasta_chunks
from collections import defaultdict
from sequence_data import k_mers
import sys,os

from optparse import *

from numpy import *
from numpy.random import random,seed

usage = """
  %prog [options] <frequency_table>

"""

parser = OptionParser(usage=usage)
parser.add_option("-l","--length",dest="length",type=int,default=1000,help="how much sequence to generate")
parser.add_option("-s","--seed",dest="seed",type=int,default=None,help="seed value for the pseudo random number generator")
parser.add_option("-H","--header",dest="header",default=">chrmarkov",help="header line (default='>chrmarkov')")
options,args = parser.parse_args()

seed(options.seed)

def load_freqs(fname):
    elements = []
    freqs = []
    
    for line in file(fname):
        if line.startswith("#"):
            continue

        elem,count = line.split("\t")[:2]
        elements.append(list(elem))
        freqs.append(float(count))
        k = len(elem)

    elements = k_mers(2)
    freqs = array(freqs,dtype=float)
    
    l = int(len(freqs)**(1./k))
    freqs = freqs.reshape((l,l))    

    return k,elements,freqs

def draw(elements,cum_p):
    mask = ((cum_p - random()) > 0)
    i = mask.nonzero()[0][0]
    return elements[i]

def generate(L,k,elements,freqs):
    freqs /= freqs.sum(axis=1)[:,newaxis]
    
    cum_p = freqs.cumsum(axis=1)


    p = freqs.sum(axis=1)
    
    nucs = list(k_mers(1))
    nuc_nums = range(len(nucs))
    last = draw(nuc_nums,p.cumsum())
    yield nucs[last]

    for i in xrange(L-1):
        next = draw(nuc_nums,cum_p[last])
        last = next
        yield nucs[next]

def fold(src):
    line = []
    for nt in src:
        line.append(nt)
        if len(line) >= 80:
            print "".join(line)
            line = []
    if line:
        print "".join(line)

if options.header:
    print options.header
    
fold(generate(options.length,*load_freqs(args[0])))
