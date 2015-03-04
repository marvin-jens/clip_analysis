# -*- coding: utf-8 -*-
from byo.track import Track
from byo.io.track_accessors import Accessor, ArrayAccessor, GenomeAccessor, AnnotationAccessor
from byo.io.lazytables import NamedTupleImporter as Importer
import byo.config

import os,re
from logging import debug,warning,error,info
from numpy import uint32, float32

root = byo.config.systems['hg19']

genome = Track(os.path.join(root,"genome"),GenomeAccessor,system='hg19')

def get_annotation_track(path=os.path.join(root,"system_annotation"),accessor=AnnotationAccessor,**kwargs):
    return Track(path,accessor,**kwargs)

def get_refGenes(path = os.path.join(root,'annotation','wgEncodeGencodeBasicV17.ucsc'))):
    from byo.gene_model import transcripts_from_UCSC
    import byo
    return transcripts_from_UCSC(path,system=byo.hg19)

chr_sizes = {}
try:
    for c in Importer(os.path.join(root,"chrom.sizes"),skip=[2],descr="## chrom:str \t size:int \t fileName:str"):
        chr_sizes[c.chrom] = c.size
except IOError:
    loaded = False
else:
    loaded = True
