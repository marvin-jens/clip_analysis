# -*- coding: utf-8 -*-
# generate self-contained flat files that can be reimported without hassle
import os,sys
from datetime import datetime
from logging import debug,warning,error

class Unsupported(Exception):
    pass

def check_supported(item):
    if not hasattr(item,"collect_data"):
        raise Unsupported()

def makelist(x,T):
    x = x.strip()
    x = x.replace(";",",")

    if not x:
        return []

    if x.endswith(','):
        x = x[:-1]

    return [T(v.strip()) for v in x.split(',')]


class BioData(object):
    """
    Implements autodecorating instances. 
    Basically a dictionary with dot-syntax data access and a thin interface
    for file-serialization and everyday handling.
    """
    def __init__(self,data_cols=None,*argc,**kwargs):
#         super(BioData,self).__init__(self,*argc,**kwargs)
#         self._types = {}

        if not data_cols:
            self._cols = kwargs.keys()
        else:
            self._cols = data_cols

        for k,v in kwargs.items():
            self.__setattr__(k,v)

        self.post_init(*argc,**kwargs)

        self._types = property(lambda self : dict([(k,type(getattr(self,k))) for k in self._cols]))

    def post_init(self,*argc,**kwargs):
        """overload this to do whatever additional parsing/decorating needs to be done"""
        pass

    def collect_data(self,cols=None):
        if not cols: cols = self._cols
        values = [self.__getattribute__(k) for k in cols]
        return dict(zip(cols,values))

    def __repr__(self):
        return "BioData: %s" % repr(self.collect_data())

class Exporter(object):
    """
    Takes an iterable sequence of instances supporting the BioData interface
    and generates a flat-file including format comment-line.
    Importer will be able to construct equivalent BioData instances from this 
    flat file.

    (Basically a *very* simple but comfortable pickle with human-readable 
    flat-files)
    """
    def __init__(self,src,dst,comment=None,cols=None):
        self.src = src
        self.dst = file(dst,"w",buffering=1)
        self.comment = comment

        # get a representative item
        x = iter(src).next()
        debug("exporting %s objects from %s" % (str(type(x)),src))

        if not cols:
            self.cols = x._cols
        else:
            self.cols = cols

        def translate(x):
            if x == 'int64':
                return 'int'
            else:
                return x

        self.colstrs = ["%s:%s" % (k,translate(type(getattr(x,k)).__name__)) for k in self.cols]

    def save(self):
        self.write_header()
        for item in self.src:
            values = [str(getattr(item,k)) for k in self.cols]
            valstr = "\t".join(values)
            self.dst.write("%s\n" % valstr)

        self.dst.close()

    def save_fasta(self,id_attr,seq_attr):
        for item in self.src:
            self.dst.write("> %s\n" % getattr(item,id_attr))
            self.dst.write("%s\n" % getattr(item,seq_attr))
        self.dst.close()

    def write_header(self):
        colstr = "\t".join(self.colstrs)
        host = os.uname()[1]
        date = datetime.now()
        
        self.dst.write("# run by Marvin Jens on %s at %s\n" % (host,date))
        if self.comment:
            self.dst.write("# %s\n" % self.comment)
        self.dst.write("# %s\n" % " ".join(sys.argv))
        self.dst.write("## %s\n" % colstr)

    def store(self,item):
#         check_supported(item)
        #data = item.collect_data()
        values = [getattr(item,k) for k in self.cols]
        valstr = "\t".join(values)
        self.dst.write("%s\n" % valstr)
        #self.dst.flush()

class Importer(object):
    """
    Factory for BioData instances. Reads a flat file with tabular layout, 
    optionally parses the column information from comments formatted as this:
        
    ## name1:type1 \t name2:type2 ...

    Alternatively accepts a list of column names and cast functions (e.g. types) 
    to do the same job for really flat files without the format comment-line.
    
    Also supports skipping of some columns (by index) to save some code and 
    memory.
    
    Supports iterator style usage for slick processing of large data-sets (In 
    this case the intermediate BioData instances should probably be destroyed
    before they consume too much memory).
    """
    
    casttypes = {
        'int' : int,
        'intlist' : lambda x : makelist(x,int),
        'floatlist' : lambda x : makelist(x,float),
        'strlist' : lambda x : makelist(x,str),
        'float' : float,
        'str' : lambda x : str(x).strip()
    }

    def __init__(self,src,dst=BioData,cols = [],casts = [],skip = [],skip_incomplete=False,skip_lines=0,guess_columns=False,skip_errors=True):
        self.src_name = src
        self.src = file(src)
        self.dst_type = dst
        self.cols = cols
        self.casts = casts
        self.skip = skip
        self.skip_incomplete = skip_incomplete
        self.skip_errors = skip_errors

        #self.casttypes = casttypes
        self.casttypes['intlist'].__name__ = "intlist"
        self.casttypes['floatlist'].__name__ = "floatlist"
        self.casttypes['strlist'].__name__ = "strlist"

        #print self.casttypes['intlist']
        if guess_columns:
            self.guess(guess_columns)
            skip_lines = max(1,skip_lines)

        if skip_lines:
            [self.src.readline() for i in range(skip_lines)]

        #print "#src: ",src

    def __str__(self):
        return "io.tables.Importer('%s')" % self.src_name

    def guess(self,cols):
        """
        attempt to guess column names and types from first two lines (like MaxQuant .txt files)
        """
        def skip(l):
            for j,i in enumerate(self.skip):
#                 print l,i-j
                l.pop(i-j)
            return l

        def pythonize(s):
            s = s.replace(" ","_")
            s = s.replace("(","")
            s = s.replace(")","")
            s = s.replace("[","")
            s = s.replace("]","")
            s = s.replace("/","_")
            s = s.replace(".","")
            s = s.replace("-","_")
            s = s.replace("+","plus")
            s = s.replace("%","percent")
            s = s.replace("__","_")
            s = s.replace("__","_")
            return s

        f = file(self.src_name)
        colnames = skip(f.readline().split('\t'))
        colnames = [pythonize(c.strip()) for c in colnames]

        for i in range(cols): f.readline()

        examples = skip(f.readline().split('\t'))
        casts = [str] * len(colnames)
        cast_candidates = [self.casttypes[c] for c in ['int','float','intlist','floatlist',]]
        for i,ex in enumerate(examples):
            for c in cast_candidates:
                try:
                    c(ex.strip())
                except(ValueError):
                    pass
                else:
                    casts[i] = c
                    break

        self.cols = colnames
        self.casts = casts

    def __iter__(self):
        types = self.casttypes
        cols = self.cols
        casts = self.casts

        def skip(l):
            for j,i in enumerate(self.skip):
#                 print l,i-j
                l.pop(i-j)
            return l

        for l in self.src:
            if not l.strip():
                # skip blank lines
                continue
            if l.startswith("##"):
                self.header_line = l
                cols = []
                casts = []
                for colstr in l.replace("#","").split("\t"):
                    #print "#",colstr
                    key,typestr = colstr.strip().split(":")
                    cols.append(key)
                    casts.append(types[typestr])

                cols = skip(cols)
                casts = skip(casts)

            data = {}
            invalid = {
                int : -1,
                float : -1,
                types['intlist'] : [],
                types['floatlist'] : [],
            }
            if not l.startswith("#"):
                #item = self.dst_type()
                vals = l.replace("#","").split("\t")
                s_vals = skip(vals)
                #print "#vlas",vals
                if self.skip_incomplete and len(vals) != len(cols):
                    warning("# warning, incomplete line! %d %d" %(len(s_vals),len(cols)))
                    continue #incomplete line
                err = False
                for i,(k,v,c) in enumerate(zip(cols,s_vals,casts)):
                    try:
                        data[k] = c(v)
                    except ValueError,V:
                        warning("# %d %s %s panic %s" % (i,k,v,V))
                        data[k] = invalid[c]
                        err = True

            if data and not (err and self.skip_errors):
                x = self.dst_type(data_cols=cols,**data)
                x.input_line = l
                yield x

    def read(self):
        return list(self.__iter__())


def hammer_exporter():
    def producer():
        #print "lalala"
        for i in xrange(int(1E6)):
            #print "*"
            b = BioData(n=i,n2 = i**2,l = ['a']*(i % 10), c = tuple())
            yield b

    Exporter(producer(),"test.table").save()
    

if __name__ == "__main__":
    #test_collection = [ BioData(a = i,b=float(i/10.),c="la"*i) for i in range(10)]
    #Exporter(test_collection,"test.table").save()
    #for i in Importer("test.table"):
        #print i
        
    #new_collection = list(Importer("test.table"))
    hammer_exporter()
    #from cProfile import run
    #run('hammer_exporter()')