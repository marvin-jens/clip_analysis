#!/usr/bin/env python
import os
import sys
import logging
logging.basicConfig(level=logging.DEBUG,format="%(asctime)s - %(levelname)s - %(message)s")

from optparse import *

usage = """
usage: %prog [options] [input.gff] > annotated.gff 

Reads GFF from file or stdin and adds annotations from the model-system database specified via -s <system> to the attr field.
"""

parser = OptionParser(usage=usage)
parser.add_option("-S","--system",dest="system",default="hg19",help="model system (default=hg19)")
parser.add_option("-a","--ann_path",dest="ann_path",default="",help="custom path to compiled annotation data")
parser.add_option("-b","--bounds",dest="bounds",default="off",choices=["left","right","both","off"],help="annotate start and end coordinate only (strand corrected). Can be 'off' (default), 'left', 'right' or 'both (creating attributes LEFT_.. and RIGHT_..")
#parser.add_option("-H","--hierarchical",dest="hierarchical",default=False,action="store_true",help="annotate only one feature that comes first in the hierarchy")
options,args = parser.parse_args()

from byo.annotation.annotator import Annotator
annotation = Annotator(options.system,options.ann_path)

if args:
    _input = file(args[0])
else:
    _input = sys.stdin

from byo.io.gff import gff_importer
try:
    annotation.reset()
    for gff in gff_importer(_input):
        if options.bounds != 'off':
            left_str,left_classes = annotation.get_str(gff.chrom,gff.start-2,gff.start+1,gff.sense)
            right_str,right_classes = annotation.get_str(gff.chrom,gff.end-2,gff.end+1,gff.sense)
            
            if gff.sense == "-":
                left_str,right_str = right_str,left_str
                
            if options.bounds == 'left':
                ann_str = left_str
            elif options.bounds == 'right':
                ann_str = right_str
            else:
                left_str = "LEFT_"+ left_str.replace(" "," LEFT_")
                right_str = "RIGHT_"+ right_str.replace(" "," RIGHT_")
                
                ann_str = left_str + " " + right_str
        else:
            ann_str,ann_classes = annotation.get_str(gff.chrom,gff.start-1,gff.end,gff.sense)
            
            
        attr_str = "%s %s" % (gff.attr_str,ann_str)
        out = [gff.chrom,gff.source,gff.type,gff.start,gff.end,gff.score,gff.sense,gff.frame,attr_str]
        print "\t".join([str(o) for o in out])
except KeyboardInterrupt:
    pass
