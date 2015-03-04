# -*- coding: utf-8 -*-
import os
def ensure_path(path,silent=False):
    path = os.path.abspath(os.path.expandvars(path))
    if not silent and (not os.path.isdir(path) and path):
        try:
            os.makedirs(path)
        except OSError:
            # path already exists, don't care
            pass
    return path

def stat_chunks(src):
    print src
    from numpy import array
    import re
    class Item(object):
        def __init__(self,title,col_labels,row_labels,data,local_options):
            self.title = title
            #self.name = name
            self.col_labels = col_labels
            self.row_labels = row_labels
            self.data = array(data)
            self.local_options = local_options

            from string import maketrans
            table = maketrans(" +-:","_pm_")
            self.safe_title = title.translate(table)
            for rl,d in zip(row_labels,data):
                rl = rl.translate(table)
                if len(d) == 1:
                    d = d[0]
                setattr(self,rl,d)
                
        def __repr__(self):
            attr_str = ", ".join(["%s=%s" % (k,str(v)) for k,v in self.__dict__.items() if not k.startswith("__")])
            return "[Item: %s {%s}]" % (self.title,attr_str)

            #print self.title, self.safe_title

    title = ""
    col_labels = ["NO LABELS"]
    data = []
    row_labels = []
    local_options = []

    for line in src:
        if line.startswith(">"):
            if data:
                yield Item(title,col_labels,row_labels,data,local_options)
                data = []
                row_labels = []
                local_options = []

            title = line[1:].strip()
            continue

        if line.startswith("##"):
            local_options = re.split("\s+",line[2:].rstrip())
            continue

        if line.startswith("#"):
            col_labels = [l.strip() for l in line[1:].rstrip().split("\t")]
            continue

        if line.strip():
            parts = line.rstrip().split("\t")
            row_labels.append(parts[0])

            data.append([float(p) for p in parts[1:]])

    if data:
        yield Item(title,col_labels,row_labels,data,local_options)
        data = []
        row_labels = []
        local_options = []
    

def qfa_chunks(lines):
    """
    Iterates over FASTQ text lines from a file like object and yields
    NamedTuple instances of the FASTQ with attributes 
    qfa.name, qfa.seq, qfa.qual
    """
    from collections import namedtuple
    QFA = namedtuple("qfa_tuple","name,seq,qual")
    
    I = lines.__iter__()

    try:
        while I:
            name = I.next().rstrip()[1:]
            seq = I.next().rstrip()
            plus = I.next().rstrip()[1:]
            qual = I.next().rstrip()
        
            yield QFA(name,seq,qual)
    except StopIteration:
        pass

def fasta_chunks(lines,strip=True,fuse=True):
    chunk = ""
    data = []

    for l in lines:
        if l.startswith("#"): continue
        if l.startswith(">"):
            if data and chunk:
                #print chunk
                yield chunk,"".join(data)

                if strip:
                    data = []
                else:
                    data = [l]

            chunk = l[1:].strip()
        else:
            if fuse:
                data.append(l.rstrip())
            else:
                data.append(l)

    if data and chunk:
        yield chunk,"".join(data)

from byo import k_mers
from numpy import array,uint8

nums = {}
for i,b in enumerate(k_mers(1)):
    nums[b] = i

def numeric(seq):
    seq = seq.lower().replace("n","a")
    return array([nums[s] for s in seq],dtype=uint8)

def fasta_as_nums(fname):
    #logger.debug("scan_fasta('%s')" % fname)
    #bases = {}
        #bases[i] = b

    for fa_id,seq in fasta_chunks(file(fname)):
        #if 'n' in seq: continue
        yield fa_id,numeric(seq)

   #logger.debug("read %d sequences of length %d" % (N,L))


import re
def wigfix_chunks(lines):
    import numpy
    regexp = r"^fixedStep\s+chrom=(?P<chrom>\w+)\s+start=(?P<start>\d+)\s+step=(?P<step>\d+).*"

    chrom = "none"
    start = 0
    step = 1
    linebuf = []

    def emit_block(linebuf,step,start):
        scores = numpy.array(linebuf,dtype=numpy.float32).repeat(step)
        return chrom,scores,start-1,start-1+len(scores)

    for l in lines:
        l = l.strip()
        if not l or l.startswith('track'):
            continue

        if l.startswith("f"):
            info = re.match(regexp,l)

            if linebuf:
                yield emit_block(linebuf,step,start)
                linebuf = []

            d = info.groupdict()
            chrom = d['chrom']
            start = int(d['start'])
            step = int(d['step'])

        else:
            linebuf.append(l)

    if linebuf:
        yield emit_block(linebuf,step,start)
        linebuf = []

