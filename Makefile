.DELETE_ON_ERROR:

## The default make target is to execute the analysis of the test data 
## as outlined in the MiMB chapter
chapter: fdr_filtered.elavl1_4su_test_data.gff
	@echo done! `wc -l fdr_filtered.elavl1_4su_test_data.gff` filtered clusters.

fdr_filtered.elavl1_4su_test_data.gff: elavl1_4su_test_data.gff
	./fdr.py elavl1_4su_test_data.tsv elavl1_4su_test_data.gff \
	--cutoff 0.05 --decoy 0.1 --logfile fdr.log \
	> fdr_filtered.elavl1_4su_test_data.gff 2> filtered_out.elavl1_4su_test_data.gff

elavl1_4su_test_data.gff: hg19.elavl1_4su_test_data.bam.bai
	./clip.py -S hg19 hg19.elavl1_4su_test_data.bam \
	-l clip.log --fp_flag markov --cluster-stats elavl1_4su_test_data.tsv \
	> elavl1_4su_test_data.gff

hg19.elavl1_4su_test_data.bam.bai: collapsed.elavl1_4su_test_data.fastq.gz hg19_decoy
	bwa aln -t 8 -n 1 -l 100 -k 1 reference/hg19/hg19_w_decoy.fa \
	 collapsed.elavl1_4su_test_data.fastq.gz 2> bwa_aln.log | \
	bwa samse reference/hg19/hg19_w_decoy.fa - \
	 collapsed.elavl1_4su_test_data.fastq.gz 2> bwa_samse.log | \
	samtools view -hbuS -F 4 - | samtools sort - hg19.elavl1_4su_test_data
	samtools index hg19.elavl1_4su_test_data.bam

collapsed.elavl1_4su_test_data.fastq.gz: trimmed.elavl1_4su_test_data.fastq.gz
	zcat trimmed.elavl1_4su_test_data.fastq.gz | \
	./collapse_reads.pl elavl1_4su 2> collapsing.log | \
	gzip > collapsed.elavl1_4su_test_data.fastq.gz

trimmed.elavl1_4su_test_data.fastq.gz: elavl1_4su_test_data.fastq.gz
	echo '>adapter\nTCGTATGCCGTCTTCTGCTTGT' > adapter.fa
	flexbar -a adapter.fa -ae RIGHT --adapter-threshold=2 \
	-r elavl1_4su_test_data.fastq.gz \
	-t trimmed.elavl1_4su_test_data \
	-f sanger --pre-trim-phred=3 \
	--min-read-length=15 -u 2 -z GZ --threads=8


## Human genome assembly hg19
hg19_reference: reference/hg19/hg19.fa reference/hg19/chrom.sizes reference/hg19/.bwa_hg19
	touch hg19_reference

reference/hg19/hg19.fa:
	echo "downloading human genome hg19 from UCSC..."
	-mkdir -p reference/hg19 2> /dev/null
	wget -c http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz -O - | gunzip | tar -xf - -O > reference/hg19/hg19.fa

reference/hg19/chrom.sizes:
	-mkdir -p reference/hg19 2> /dev/null
	wget -c http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes  -O - > reference/hg19/chrom.sizes

reference/hg19/.bwa_hg19: reference/hg19/hg19.fa
	echo "building BWA index..."
	-bwa 2>&1 | grep Version > reference/hg19/bwa_hg19.jk
	cd reference/hg19; bwa index -a bwtsw hg19.fa >> bwa_hg19.jk 2> bwa_hg19.jk2
	touch reference/hg19/.bwa_hg19


## Adding ~600Mb of random decoy sequence (decoy scale 0.1)
hg19_decoy: reference/hg19/hg19_w_decoy.fa reference/hg19/.bwa_hg19_w_decoy
	touch hg19_decoy

reference/hg19/hg19_w_decoy.fa: reference/hg19/hg19.fa
	./markov.py -l 500000000 2mer_freqs/hg19.2mer > rnd_hg19.fa
	./markov.py -H ">chrmarkov_utr3" -l 10000000 2mer_freqs/utr3.2mer > rnd_3utr.fa
	./markov.py -H ">chrmarkov_intron" -l 100000000 2mer_freqs/introns.2mer > rnd_introns.fa
	./markov.py -H ">chrmarkov_CDS" -l 10000000 2mer_freqs/utr3.2mer > rnd_cds.fa
	cat rnd_hg19.fa rnd_3utr.fa rnd_introns.fa rnd_cds.fa reference/hg19/hg19.fa > reference/hg19/hg19_w_decoy.fa

reference/hg19/.bwa_hg19_w_decoy: reference/hg19/hg19_w_decoy.fa
	echo "building BWA index..."
	-bwa 2>&1 | grep Version > reference/hg19/bwa_hg19_w_decoy.jk
	cd reference/hg19; bwa index -a bwtsw hg19_w_decoy.fa >> bwa_hg19_w_decoy.jk 2> bwa_hg19_w_decoy.jk2
	touch reference/hg19/.bwa_hg19_w_decoy

## Downloading transcript annotation (and optionally compiling it into a custom format for fast lookup)
hg19_annotation: annotation/hg19/wgEncodeGencodeBasicV17.ucsc
	touch hg19_annotation

annotation/hg19/wgEncodeGencodeBasicV17.ucsc:
	echo "downloading GENCODE17 annotation for hg19 from UCSC.."
	-mkdir annotation 2> /dev/null
	wget 'http://genome.ucsc.edu/cgi-bin/hgTables?clade=mammal&org=Human&db=hg19&hgta_group=genes&hgta_track=wgEncodeGencodeBasicV17&hgta_table=wgEncodeGencodeBasicV17&hgta_regionType=genome&hgta_outputType=primaryTable' --post-data='hgta_doTopSubmit=get_output' -O - > annotation/hg19/wgEncodeGencodeBasicV17.ucsc

hg19_compiled_annotation: hg19_annotation reference/hg19/chrom.sizes
	echo "compiling annotation for fast lookup..."
	-mkdir -p annotation/hg19/compiled 2> /dev/null
	./build_annotation_track.py -S hg19 annotation/hg19/*.ucsc -d annotation/hg19/compiled
	touch hg19_compiled_annotation


