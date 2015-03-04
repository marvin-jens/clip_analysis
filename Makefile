.DELETE_ON_ERROR:

all: .references .bwa_indices .compiled_annotations
	    echo "done!"

.references: reference/hg19/hg19.fa
	    touch .references

.bwa_indices: reference/hg19/.bwa_hg19 
	    touch .bwa_indices

.compiled_annotations: annotation/hg19/.compiled_annotation
	    touch .compiled_annotations

reference/hg19/hg19.fa:
	    echo "downloading human genome hg19 from UCSC..."
	    -mkdir -p reference/hg19 2> /dev/null
	    wget -c http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz -O - | gunzip | tar -xf - -O > reference/hg19/hg19.fa

reference/hg19/.bwa_hg19: reference/hg19/hg19.fa
	    echo "building BWA index..."
	    -bwa 2>&1 | grep Version > reference/hg19/bwa_hg19.jk
	    cd reference/hg19; bwa index -a bwtsw hg19.fa >> bwa_hg19.jk 2> bwa_hg19.jk2
	    touch reference/hg19/.bwa_hg19

reference/hg19/chrom.sizes:
	    -mkdir reference 2> /dev/null
	    wget -c http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes  -O - > reference/hg19/chrom.sizes

annotation/hg19/wgEncodeGencodeBasicV17.ucsc:
	    echo "downloading GENCODE17 annotation for hg19 from UCSC.."
	    -mkdir annotation 2> /dev/null
	    wget 'http://genome.ucsc.edu/cgi-bin/hgTables?clade=mammal&org=Human&db=hg19&hgta_group=genes&hgta_track=wgEncodeGencodeBasicV17&hgta_table=wgEncodeGencodeBasicV17&hgta_regionType=genome&hgta_outputType=primaryTable' --post-data='hgta_doTopSubmit=get_output' -O - > annotation/hg19/wgEncodeGencodeBasicV17.ucsc

annotation/hg19/.compiled_annotation: reference/hg19/chrom.sizes annotation/hg19/wgEncodeGencodeBasicV17.ucsc
	    echo "compiling annotation for fast lookup..."
	    -mkdir -p annotation/hg19/compiled 2> /dev/null
	    ./build_annotation_track.py -S hg19 annotation/hg19/*.ucsc -d annotation/hg19/compiled
	    touch annotation/hg19/.compiled_annotation

