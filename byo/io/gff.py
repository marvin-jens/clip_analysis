# -*- coding: utf-8 -*-

from byo.io.lazytables import NamedTupleImporter as Importer
def gff_importer(path,fix_chr=True,**kwargs):
    cols = ['chrom','source','type','start','end','score','sense','frame','attr_str']
    casts = [str,str,str,int,int,str,str,str,str]

    # some people save 3 bytes in C.elegans and call the chromosomes I,II,II,IV,X
    if not fix_chr:
        mangle_chrom = lambda x : x
    else:
        def mangle_chrom(data):
            if not data[0].startswith('chr'):
                data[0] = 'chr'+data[0]
            return data
    
    for row in Importer(path,cols=cols,casts=casts,mangling=mangle_chrom,parse_headers=False,**kwargs):
        yield row

def format_gff(chrom,source,type,start,end,score,sense,frame,attr_str):
    out = [chrom,source,type,start,end,score,sense,frame,attr_str]
    return "\t".join([str(o) for o in out])

import re
from collections import defaultdict
def dict_from_attrstr(attrs):
    attr_dict = defaultdict(lambda : None)
    for m in re.finditer(r'(\w+)([ ="]+)(\S+?)(["\;\s]+|$)',attrs):
        #print m.groups()
        key,bla1,value,bla2 = m.groups()
        attr_dict[key] = value

    return attr_dict
        
if __name__ == "__main__":
    tests = [
        """name="CID_012656"; length="21"; ent="1.780"; enrichment="-inf"; cov="63750"; edits="3364.0"; conv="1760"; ccr="48652938";""",
        """ID=Transcript:F58E1.5.T2;transcribed=true;normscore=1;overlapping_wormbase_transcript=F58E1.5;prediction_status=Partially_confirmed""",
        """ID=MI0022705;Alias=MI0022705;Name=hsa-mir-6859-1"""
    ]
    for t in tests:
        print t
        print dict_from_attrstr(t)
