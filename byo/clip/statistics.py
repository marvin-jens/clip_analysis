
from collections import defaultdict


class Stats(object):
    def __init__(self,path="",stats=None):
        if not stats:
            stats = defaultdict(lambda : defaultdict(int))
        self.editstats = stats
        self.path = path
        self.store = None

    def update_stats(self,flags,stats,sense):
        fstats = self.editstats[flags]
        if sense == "-":
            for k,v in stats.items():
                if not "N" in k: fstats[edit_map.get(k,k)] += v
        else:
            for k,v in stats.items():
                if not "N" in k: fstats[k] += v

        self.editstats[flags] = fstats

    def save_dict(self,title="count statistics",stats={},labels=["key","value"]):
        block = ["%s\t%s" % (k,stats[k]) for k in sorted(stats.keys())]

        self.prep_storage()            
        self.store.write(stat_block % dict(title=title,labels="\t".join(labels),block="\n".join(block)))
        

    def prep_storage(self):
        from byo.io import ensure_path
        import os
        if not self.store:
            ensure_path(os.path.dirname(self.path))
            self.store = file(self.path,"w")
        
    def save_editstats(self,title="conversion statistics",editstats=None):

        if not editstats:
            editstats = self.editstats
        keys = ["intergenic","transcript","antisense"]
        #if len(editstats.keys()) < 2:
            #keys = ["all"]

        events = ['total','uniq','perfect','edits','AC','AG','AT','CA','CG','CT','GA','GC','GT','TA','TC','TG','_A','_C','_G','_T','A_','C_','G_','T_']
            
        labels = []
        data = defaultdict(list)

        for key in keys:
            #if not key in editstats: continue
            stats = editstats[key]
            if key == "antisense":
                sign = -1
            else:
                sign = 1
            for k in events:
                data[k].append(sign * stats[k])

            labels.append(key)

        block = ["%s\t%s" % (translation.get(k,k),"\t".join([str(v) for v in data[k]])) for k in events]

        self.prep_storage()
        self.store.write(stat_block % dict(title=title,labels="\t".join(labels),block="\n".join(block)))

    def flush(self):
        if self.store:
            self.store.close()

edit_map = {
    'AC':'TG',
    'AG':'TC',
    'AT':'TA',
    'CA':'GT',
    'CG':'GC',
    'CT':'GA',
    'GA':'CT',
    'GC':'CG',
    'GT':'CA',
    'TA':'AT',
    'TC':'AG',
    'TG':'AC',
    '_A':'_T',
    '_C':'_G',
    '_G':'_C',
    '_T':'_A',
    'A_':'T_',
    'C_':'G_',
    'G_':'C_',
    'T_':'A_',
    'total':'total',
    'uniq':'uniq',
    'perfect' : 'perfect',
    'edits' : 'edits',
}

translation = {
    'total':'total',
    'perfect':'perfect',
    'edits':'edits',
    'uniq' : 'distinct',
    'AC':"A:C",
    'AG':"A:G",
    'AT':"A:T",
    'CA':"C:A",
    'CG':"C:G",
    'CT':"C:T",
    'GA':"G:A",
    'GC':"G:C",
    'GT':"G:T",
    'TA':"T:A",
    'TC':"T:C",
    'TG':"T:G",
    '_A':"A_ins",
    '_C':"C_ins",
    '_G':"G_ins",
    '_T':"T_ins",
    'A_':"A_del",
    'C_':"C_del",
    'G_':"G_del",
    'T_':"T_del",
}

stat_block = """
>%(title)s
#read alignments and edits\t%(labels)s
## --rotate=90 --figsize=14,8 --color_cycle=gray,blue,red --spacing=.5
%(block)s
"""
        