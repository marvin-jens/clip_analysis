# -*- coding: utf-8 -*-
from byo.io.tables import makelist
from collections import namedtuple
from pprint import pprint
import itertools
import logging

class LazyImporter(object):
    def __init__(self,src,dst,key_attr,**kwargs):
        self.dst = dst
        self.tuple_data = {}
        self.object_data = {}
        self.ordered_keys = []
        self.key_attr = key_attr

        self.kw = kwargs
        self.load(src)

    def load(self,src,**kwargs):
        if hasattr(self.key_attr,"__call__"):
            getter = self.key_attr
        else:
            getter = lambda item : getattr(item,self.key_attr)
            
        self.kw.update(kwargs)
        for item in NamedTupleImporter(src,**self.kw):

            key = getter(item)
            self.tuple_data[key] = item
            self.ordered_keys.append(key)

    def __getitem__(self,key):
        if key in self.object_data:
            return self.object_data[key]

        obj = self.dst(**self.tuple_data[key]._asdict())
        self.object_data[key] = obj

        return obj

    def __iter__(self):
        for key in self.ordered_keys:
            yield self[key]

    def __contains__(self,key):
        return key in self.tuple_data

    def __len__(self):
        return len(self.tuple_data)

class NamedTupleImporter(object):
    """
    """

    casttypes = {
        'int' : int,
        'intlist' : lambda x : makelist(x,int),
        'floatlist' : lambda x : makelist(x,float),
        'strlist' : lambda x : makelist(x,str),
        'float' : float,
        'str' : lambda x : str(x).strip()
    }

    def __init__(self,src,dst=None,cols = [],casts = [],skip = [],type_hints={},renames={},parse_comments=False,parse_headers=True,default_cast='str',comment_char="#",sep="\t",descr="",skip_broken=True,keep_line=False,debug=False,mangling=lambda x : x):
        self.src_name = src
        if type(src) == file:
            self.src = src
        else:
            self.src = file(src)

        self.mangling = mangling
        self.debug = debug
        self.logger = logging.getLogger("NamedTupleImporter('%s')" % src)
        self.dst_type = dst
        self.cols = cols
        self.casts = casts
        self.skip = sorted(skip,reverse=True)
        self.type_hints = type_hints
        self.renames = renames
        self.comment_char = "#"
        self.parse_comments = parse_comments
        self.parse_headers = parse_headers
        self.skip_broken = skip_broken
        self.keep_line = keep_line
        if keep_line:
            self.casts.append(str)
            self.cols.append('line')
        
        self.default_cast = default_cast
        self.sep = sep

        #self.casttypes = casttypes
        self.casttypes['intlist'].__name__ = "intlist"
        self.casttypes['floatlist'].__name__ = "floatlist"
        self.casttypes['strlist'].__name__ = "strlist"

        # alternatively you can pass a shorthand table description
        if descr:
            self.load_descr(descr)
        #print "#src: ",src

    def load_descr(self,l):
        try:
            # try parsing the new table description line
            n_cols = []
            n_casts = []
            for colstr in l.replace("#","").split("\t"):
                #print "#",colstr
                key,typestr = colstr.strip().split(":")
                n_cols.append(key)
                n_casts.append(self.casttypes[typestr])

            skip = [s for s in self.skip if s < len(n_cols)]

            for i in skip:
                n_cols.pop(i)
                n_casts.pop(i)
                    
        except ValueError:
            self.logger.warning("ignoring malformed table description line '%s' in '%s'" % (l,self.src_name))
        else:
            # everything went well. Accept the new definition
            self.header_line = l
            self.cols = n_cols
            self.casts = n_casts
            self.skip = skip
        
        if self.keep_line:
            self.casts.append(str)
            self.cols.append('line')
            
    def __str__(self):
        return "io.lazytables.Importer('%s')" % self.src_name

    def __iter__(self):
        types = self.casttypes
        cols = self.cols
        casts = self.casts
        comment = self.comment_char
        sep = self.sep

        import re
        
        self.tuple_type = namedtuple('namedtupleimporter_predefined',",".join(cols))

        for l in self.src:
            if l.startswith("##") and self.parse_headers:
                    self.load_descr(l)
                    cols = self.cols
                    casts = self.casts
                    
                    self.tuple_type = namedtuple('namedtupleimporter_autodetect',",".join(cols))

            elif l.startswith(comment):
                if not self.parse_comments:
                    continue
                #self.logger.warning("loading from comments '%s'" % l)

                self.cols = [self.renames.get(c,c) for c in re.split("\s+",l[1:].strip())]
                self.casts = [self.casttypes[self.type_hints.get(c,self.default_cast)] for c in self.cols]
                if self.keep_line:
                    self.casts.append(str)
                    self.cols.append('line')
                
                cols = self.cols
                casts = self.casts
                self.skip = []
                                  
                self.tuple_type = namedtuple('namedtupleimporter_autodetect',",".join(cols))
                #if self.debug:
                    #pprint(zip(cols,casts))
            
            else:
                l = l.rstrip()
                if not l:
                    continue

                vals = l.split(sep)
                #print vals,len(vals)
                for i in self.skip:
                    vals.pop(i)
                #print vals,self.skip
                try:
                    if self.debug:
                        pprint(list(itertools.izip_longest(cols,vals,casts,fillvalue="l")))
                    if self.keep_line:
                        #print vals
                        data = [c(v) for c,v in itertools.izip_longest(casts,vals,fillvalue=str)][:len(casts)]
                        data[-1] = l
                    else:
                        data = [c(v) for c,v in itertools.izip_longest(casts,vals,fillvalue="")]
                    yield self.tuple_type(*self.mangling(data))

                except ValueError as E:
                    self.logger.warning("cannot parse line '%s'" % l)

                    if not self.skip_broken: 
                        raise E
        
    def read(self):
        return list(self.__iter__())
