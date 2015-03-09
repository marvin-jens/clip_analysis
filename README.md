# A pipeline for PAR-CLIP data analysis
This repository contains the python scripts to analyze (PAR-)CLIP data, as referenced in [2015 Methods in Molecular Biology (MiMB)](http://www.springer.com/series/7651). 

## Acknowledgements
The procedures for PAR-CLIP data analysis implemented in these scripts have been developed over the past 5 years in the Rajewsky lab and have received input and work from many people. Among them are Jonas Maaskola, Sebastian Mackowiak, Minnie Zhuo Fang,  Andranik Ivanov, and Nikolaus Rajewsky. They have also been adapted and developed in close collaboration with the wet lab, where PAR-CLIP experiments have been performed by Svetlana Lebedeva, Agnieszka Rybak-Wolf, Yasuhiro Murakawa (Landthaler lab), Anna-Carina Jungkamp and Stefanie Grosswendt. The author further acknowledges many fruitful discussions with Markus Landthaler.

## Authorship and legal information
The scripts, in the form provided in this archive, are the work of Dr. Marvin Jens, currently working as a post-doc in Nikolaus Rajewsky's lab at the [Max-Delbrueck-Center for Molecular Biology](https://www.mdc-berlin.de/1151037/en/research/research_teams/systems_biology_of_gene_regulatory_elements), with the exception of `collapse_reads.pl`, provided by Sebastian Mackowiak. The code is hereby released under the terms of the [GNU General Public License](LICENSE). If you use this code for a publication, please cite the MiMB chapter. If there is a newer publication referring to a potentially improved or extended method, this README will be updated. The scripts are well-tested, but come without any guarantee! It is your own responsibility to sanity-check your results and use our code with consideration. Apart of that, please be welcome to use it, play around with it, and contact me if you have improvements, suggestions, or bug-reports.

## Installation
The code is purely python and known to run with python 2.7 on Ubuntu 14.04 LTS. It depends on the `numpy` and `pysam` python libraries (available either through your distributions package manager, or via the Python Package Index). On Ubuntu, the following commands should do the trick:

```bash
sudo apt-get install python-pip python-numpy
pip install pysam
```

Further, 3rd party requirements are the Flexible adapter remover `flexbar`, the Burrows-Wheeler aligner `bwa`,`samtools`, and (optionally) `bwa pssm`. It is recommended to download the latest versions from the corresponding project websites:

* [FLEXBAR](http://sourceforge.net/projects/flexbar/)
* [bwa](http://bio-bwa.sourceforge.net/)
* [samtools](http://www.htslib.org/)
* [bwa pssm](http://bwa-pssm.binf.ku.dk/)

It is recommended to have at least 8GB of RAM in your workstation for running these tools efficiently (esp. `bwa index`), and many cores are helpful, too.

## Reference sequences and annotations
There is a `Makefile` included, which streamlines downloading reference sequences and transcript annotations from UCSC. You can also use it to build mapping decoys and the corresponding indices, and run the code from the MiMB chapter. It currently supports the following targets:
* `hg19_reference, hg19_`



