import re,os
from collections import defaultdict

#from byo.annotation import LookupCache

#class AnnotationFilter(LookupCache)

class Annotator(object):
    def __init__(self,system_name,ann_path="",**kwargs):
        from byo import systems
        system = getattr(systems,system_name)
        if ann_path:
            self.annotation = system.get_annotation_track(path=ann_path,annotation_filter=self.channel,**kwargs)
        else:
            self.annotation = system.get_annotation_track(annotation_filter=self.channel,**kwargs)

        self.reset()

    def reset(self):
        self.path_counts = defaultdict(float)

    @staticmethod
    def channel(data):
        categories = set()
        features = set()
        elements = set()
        
        for a in data.split("$"):
            if not a:
                continue
           
            path,name = os.path.split(a)
            categories.add(path)
            root = path.split("/")[0]
            for n in name.split(","):
                features.add("%s/%s" % (root,n))
                elements.add("%s/%s" % (path,n))

        if categories:    
            return dict(categories = categories,features=features,elements=elements)

        return {}

    @staticmethod
    def compact(s):
        return re.sub(r"(transcript/processing/)|(transcript/translation/)","",s)

    @staticmethod
    def compact_feat(s):
        return re.sub(r"transcript/","",s)

    def get(self,chrom,start,end,sense):
        cat_counts = defaultdict(int)
        as_cat_counts = defaultdict(int)
        feat_counts = defaultdict(int)
        as_feat_counts = defaultdict(int)
        
        elem_counts = defaultdict(int)
        as_elem_counts = defaultdict(int)
        
        for a in self.annotation.get(chrom,start,end,sense):
            if a:
                for cat in a['categories']:
                    cat_counts[cat] += 1
                for f in a['features']:
                    feat_counts[f] += 1
                for e in a['elements']:
                    elem_counts[e] += 1

        antisense = {'+' : '-', '-':'+'}
        for a in self.annotation.get(chrom,start,end,antisense[sense]):
            if a:
                for cat in a['categories']:
                    as_cat_counts[cat] += 1
                    #cat_counts["antisense/"+cat] += 1
                for f in a['features']:
                    as_feat_counts[f] += 1
                for e in a['elements']:
                    as_elem_counts[e] += 1

        for cat in as_cat_counts.keys():
            if cat.startswith("transcript"):
                cat_counts["antisense"] += 1

        classes = sorted(cat_counts.keys())
        as_classes = sorted(as_cat_counts.keys())
        features = sorted(feat_counts.keys())
        as_features = sorted(as_feat_counts.keys())
        
        if not classes:
            classes = ["intergenic"]
            cat_counts = dict(intergenic=1)

        #self.add_counts(cat_counts,length=end-start)
        #self.add_counts(as_cat_counts,prefix="antisense",length=end-start)
        self.add_counts(elem_counts,length=end-start)
        self.add_counts(as_elem_counts,prefix="antisense",length=end-start)
         
            
        return classes,features,as_classes,as_features

    def hierarchical_category(self,chrom,start,end,sense):

        categories = [
            ('structural','mt-transcript',('chrM',None),'#0cadbf',10),
            ('structural','RMRP',("RMRP",None),'#7d8bbf',6),
            ('noncoding.regulatory','miRNA',("miRNA",None),'#F25C05',4),
            ('structural','snoRNA',("(snoRNA)|(SNORA)",None),'#ab8705',5),
            ('structural','snRNA',("(snRNA)|(RNU12)",None),'#93de05',6),
            ('structural','scRNA',("scRNA",None),'#93816e',6),
            ('structural','vtRNA',("(vtRNA)|(VTRNA)",None),'#9372c1',6),
            ('structural','srpRNA',("srpRNA",None),'#93a8de',6),
            ('structural','RNAse-P',("RPPH",None),'#2c58ab',6),
            ('structural','(pseudo-)tRNA',("tRNA",None),'#93a605',7),
            ('structural','(pseudo-)rRNA',("(rRNA)|(RNA5SP)|(RNA5-8SP)|(RNA18SP)|(RNA28SP)",None),'#3e5916',8),
            ('structural','other_structural',("Y_RNA",None),'#00b792',9),
            ('coding','CDS',(r"transcript/coding(.*?)exon(.*?)CDS",None),'#d19743',1),
            ('coding',"5'UTR",(r"transcript/coding(.*?)exon(.*?)UTR5",None),'#ffca0e',0),
            ('coding',"3'UTR",(r"transcript/coding(.*?)exon(.*?)UTR3",None),'#ab362e',2),
            ('coding','intronic',("intron","exon"),'#ffedb7',3),   
            ('coding','other_coding',(None,None),'#ffc03a',3.1),
            #('tx_island',("bed2gff.py",None),'#F5EFFA',13),
            ('misc','antisense',("(antisense)",None), '#808080',11.6),
            ('noncoding','lincRNA/long ncRNA',("(transcript/(noncoding)|(ncRNA))|(lincRNA)","(miRNA)|(snoRNA)|(snRNA)|(tRNA)|(rRNA)|(bed2gff.py)|(transcript/coding)"),'#ff8080',10),
            ('coding',"3'extension",("downstream",None),'#7d087f',2.5),
            ('misc','promoter/enhancer',("(promoter)|(enhancer)",None),'#f5ff00',-1),
            ('misc','repeat',("repeat","transcript"),'#bfbfbf',12),
            ('misc','intergenic',("(intergenic)","transcript"),'#FFFFFF',11.5),
            ('misc','other',(None,None),'#e2e2e2',11),
        ]
        
        classes,features,as_classes,as_features =  self.get(chrom,start,end,sense)
        #print classes
        
        cls = ",".join(classes + [chrom])
        matches = []
        #print line
        for group,name,(yes,no),color,rank in categories:
            #print cls
            #print "Y",yes
            #print "N",no
            
            if not yes:
                continue

            if re.search(yes,cls):
                if not no or (not re.search(no,cls)):
                    return name
                    #matches.append(name)
                    #break

        return "None (this should not happen)"
    
        

    def get_str(self,chrom,start,end,sense):
        
        classes,features,as_classes,as_features = self.get(chrom,start,end,sense)
        #print features
        #self.print_counts()

        attrs = [
            'CLASSES="%s";' % self.compact(",".join(classes)),
            'AS_CLASSES="%s";' % self.compact(",".join(as_classes)),
            'SENSE="%s";' % self.compact_feat(",".join(features)),
            'ANTISENSE="%s";' % self.compact_feat(",".join(as_features))
        ]
        return ' '.join(attrs),classes

    def add_counts(self,count_dict,prefix="",normalize=True,length=1):
        paths = count_dict.keys()
        
        #weights = {}
        #decomposed = defaultdict(list)
        
        w = 1
        todo = set()
        for p in paths:
            p = os.path.join(prefix,p)
            for m in re.finditer("/",p):
                #self.path_counts[p[:m.start()]] += w
                todo.add(p[:m.start()])
            todo.add(p)

        for p in todo:
            self.path_counts[p] += w

            

    def print_counts(self):
        for count,path in sorted([(c,p) for p,c in self.path_counts.items()],reverse=True):
            print "%6d '%s'" % (count,path)
        #Z = 1
        #if normalize:
            #for count in count_dict.
        #for path,count in count_dict.items():
        #pass