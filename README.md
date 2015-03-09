# A pipeline for PAR-CLIP data analysis
This repository contains the python scripts to analyze (PAR-)CLIP data, as referenced in [2015 Methods in Molecular Biology (MiMB)](http://www.springer.com/series/7651). 

## Acknowledgements
The procedures for PAR-CLIP data analysis implemented in these scripts have been developed over the past 5 years in the Rajewsky lab and have received input and work from many people. Among them are Jonas Maaskola, Sebastian Mackowiak, Minnie Zhuo Fang,  Andranik Ivanov, and Nikolaus Rajewsky. They have also been adapted and developed in close collaboration with the wet lab, where PAR-CLIP experiments have been performed by Svetlana Lebedeva, Agnieszka Rybak-Wolf, Yasuhiro Murakawa (Landthaler lab), Kerstin Baethge (Landthaler lab), Anna-Carina Jungkamp and Stefanie Grosswendt. The author further acknowledges many fruitful discussions with Markus Landthaler.

## Authorship and legal information
The scripts, in the form provided in this archive, were implemented by Dr. Marvin Jens, currently working as a post-doc in Nikolaus Rajewsky's lab at the [Max-Delbrueck-Center for Molecular Biology](https://www.mdc-berlin.de/1151037/en/research/research_teams/systems_biology_of_gene_regulatory_elements), with the exception of `collapse_reads.pl`, provided by Sebastian Mackowiak. The code is hereby released under the terms of the [GNU General Public License](LICENSE). If you use this code for a publication, please cite the MiMB chapter. If there is a newer publication referring to a potentially improved or extended method, this README will be updated. The scripts are well-tested, but come without any guarantee! It is your own responsibility to sanity-check your results and use our code with consideration. Apart of that, please be welcome to use it, play around with it, and contact me if you have improvements, suggestions, or bug-reports.

## Installation
The first step is to clone this git repository:
```
git clone https://github.com/marvin-jens/clip_analysis
```
The code is entirely written in python and known to run with python 2.7 on Ubuntu 14.04 LTS. The `numpy` and `pysam` python libraries (available either through your distributions package manager, or via the Python Package Index) are required. On Ubuntu, the following commands should do the trick:

```
sudo apt-get install python-pip python-numpy
pip install pysam
```

Further 3rd party requirements are the Flexible Adapter Remover `flexbar`, the Burrows-Wheeler aligner `bwa`,`samtools`, and (optionally) `bwa pssm`. It is recommended to download the latest versions from the corresponding project websites:

* [FLEXBAR](http://sourceforge.net/projects/flexbar/)
* [BWA](http://bio-bwa.sourceforge.net/)
* [SAMTOOLS](http://www.htslib.org/)
* [BWA-PSSM](http://bwa-pssm.binf.ku.dk/) (optional)

To run these tools efficiently (esp. `bwa index`) it is recommended to have at least 8GB of RAM in your workstation, and many cores are helpful for adapter removal and mapping.

## Reference sequences and annotations
There is a `Makefile` included, which streamlines downloading reference sequences and transcript annotations from UCSC. You can also use it to build mapping decoys and the corresponding indices, and run the code from the MiMB chapter. It currently supports the following targets:
* hg19_reference
* hg19_decoy
* hg19_annotation
* hg19_compiled_annotation (Optional. Takes ~12Gb of RAM to build)

If you run `make chapter` or just `make`, it will download and build the pre-requisites for the analysis of the test data and then run the commands from the book.

## Compiled annotation [OPTIONAL]
If you make the *hg19_compiled_annotation* target, this will break down each annotated transcript into a set of functional features (exons, introns, UTRs, splice-sites, etc.), sort them on their genomic coordinates, and build a big, sparse, binary file for each strand in the reference (chromosome+strand) for fast random lookup. This is an entirely custom annotation infrastructure that circumvents any use of databases. If this procedure works for you, the `gff_annotate.py` script can be used to quickly add annotation information to your PAR-CLIP read clusters (or any other GFF).

## Other tools
There is no single 'best' way to analyze your data. It always depends on the exact question you are trying to answer and on the details of the project. There are other excellent tools available to analyze PAR-CLIP data, that may be more suited to your needs or easier to work with. Especially noteworthy are [PARalyzer](https://ohlerlab.mdc-berlin.de/software/PARalyzer_85/), which features a kernel density estimation based segmentation of clusters into binding sites, and [CLIPZ](http://www.clipz.unibas.ch/), a web-browser accessible environment for data analysis and exploration that does not require the use of commandline tools.
