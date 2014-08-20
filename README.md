#ContaVect
=========

## Motivation
Contavect was developped to quantify and caracterize DNA contaminants from gene therapy vector production after NGS sequencing. This automated pipeline can however be used for wider pourpose requiring to identify map NGS datasets consisting of a mix of DNA sequences on multiple references. It combine several features such as reference homologies masking, fastq filtering/adapter trimming, short read alignments, SAM file splitting and generating human readable output.

##Principle


### Dependencies:

The programm was developed under Linux Mint 16/17 and require a python 2.7 environment.
The following dependencies are required for proper program excecution:

* [bwa 0.7.0+](https://github.com/lh3/bwa)
* [samtools 0.1.17+](https://github.com/samtools/samtools)
* blast+ (not required is not no reference masking is to be performed) 

In addtion 2 third party python packages are also needed 

* [Biopython](https://github.com/biopython/biopython)
* [pysam](https://github.com/pysam-developers/pysam)

## Get IsFinder

1. Clone the repository
``` bash
$ git clone https://github.com/a-slide/ContaVect/ my_folder/
```

2. Make the main script excecutable
``` bash
$ sudo chmod u+x IsFinder.py
```

3.Compile the ssw aligner (and add the dynamic library it to the path)

If you wish to perform a step of adapter trimming before mapping you need to complile the dynamic library ssw.so to be able to use the Smith watermann algorithm forked from mengyao's [Complete-Striped-Smith-Waterman-Library](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library)

To use the dynamic library libssw.so you may need to modify the LD_LIBRARY_PATH environment
variable to include the library directory (export LD_LIBRARY_PATH=$PWD) or for definitive
inclusion of the lib edit /etc/ld.so.conf and add the path or the directory containing the
library and update the cache by using /sbin/ldconfig as root

## Usage

Prepare the configuration file to include your files and settings as indicated in the template Conf.txt file provided with the source files

``` bash
$ ./ContaVect.py Conf.txt......
```

## Pipeline development logbook
TO BE INCLUDED
* [Logbook]()

## Authors and Contact

Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)




