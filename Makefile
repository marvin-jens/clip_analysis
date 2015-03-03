.DELETE_ON_ERROR:

all: .reference .bwa_indices .compiled_annotation
	    touch bootstrap_complete 

.reference:
	    -mkdir reference 2> /dev/null
	    wget -c http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz -O - | gunzip | tar -xf - -O > reference/hg19.fa
	    wget -c http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes  -O - > reference/hg19.chrom.sizes
	    touch .reference

annotation/wgEncodeGencodeBasicV17.ucsc:
	    -mkdir annotation 2> /dev/null
	    echo  wgEncodeGencodeBasicV17
	    wget 'http://genome.ucsc.edu/cgi-bin/hgTables?clade=mammal&org=Human&db=hg19&hgta_group=genes&hgta_track=wgEncodeGencodeBasicV17&hgta_table=wgEncodeGencodeBasicV17&hgta_regionType=genome&hgta_outputType=primaryTable' --post-data='hgta_doTopSubmit=get_output' -O - > annotation/wgEncodeGencodeBasicV17.ucsc


.compiled_annotation: chrom.sizes annotation/wgEncodeGencodeBasicV17.ucsc  annotation/repeat.repmask.ucsc_repmask_track.gff  annotation/refGene.ucsc annotation/transcript.noncoding.mirbase20.gff annotation/transcript.noncoding.ucsc_trna_track.gff annotation/transcript.noncoding.snoRNA.ucsc_srna_track.gff annotation/transcript.noncoding.lincRNA.ucsc_lincrna_track.ucsc annotation/transcript.noncoding.ucsc_evofold_track.gff 
	    build_annotation_track.py -S hg19 annotation/*.ucsc annotation/*.gff
	    touch .compiled_annotation


.bwa_genome: .reference 
	    -bwa 2>&1 | grep Version > reference/bwa_hg19.jk
	    cd reference; bwa index -a bwtsw hg19.fa >> bwa_hg19.jk 2> bwa_hg19.jk2
	    touch .bwa_genome

.bwa_indices: .bwa_genome 
	    touch .bwa_indices

bwa aln -t 8 -n 1 -l 100 -k 1 ~/pcp/systems/hg19/genome/hg19.fna \
  collapsed.elavl1_4su_test_data.fastq.gz 2> bwa_aln.log | \
bwa samse ~/pcp/systems/hg19/genome/hg19.fna - \
  collapsed.elavl1_4su_test_data.fastq.gz 2> bwa_samse.log | \
samtools view -hbuS -F 4 - | samtools sort - hg19.elavl1_4su_test_data