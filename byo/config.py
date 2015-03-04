# -*- coding: utf-8 -*-
# global default parameters independent of model-system data
import os
import logging

os.environ['PCP'] = os.path.abspath(os.path.join(os.path.dirname(__file__),"../"))

systems = {}
system_root = "${PCP}"

try:
    from byo.localconfig import *
except ImportError:
    logging.warning("You have not configured the pipeline by creating 'byo/localconfig.py'! Please read the manual before proceeding!")
    import sys
    sys.exit(1)


# recursively crawl dictionaries and expand variables like ${HOME} 
def _expand_recursive(d):
    for k,v in d.items():
        if not k.startswith("__"):
            if type(v) == dict:
                _expand_recursive(v)
            if type(v) == str:
                d[k] = os.path.expandvars(v)
                #print k,v

_expand_recursive(globals())
