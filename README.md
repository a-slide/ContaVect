#ContaVect

## Motivation
Contavect was developped to quantify and caracterize DNA contaminants from gene therapy vector production after NGS sequencing. This automated pipeline can however be used for wider pourpose requiring to identify map NGS datasets consisting of a mix of DNA sequences on multiple references. It combine several features such as reference homologies masking, fastq filtering/adapter trimming, short read alignments, SAM file splitting and generating human readable output.

##Principle

Contavect a python pipeline composed of several modules linked together to analyse NGS Data. Here is a description of the overall workflow principle :

* Each reference fasta file is parsed to identify all sequences within it and a Reference object is initialised to save the reference characteristics, the name and the output required.
* Facultative: Homologies between references can be masked iteratively, starting by the last reference which is masked by all the others then to the penultimate masked by all others except the last and and so forth until there is only 1 reference remaining. This is done using blastn from blast+ package
* Facultative: Fastq can be filtered by mean quality and adapters can be trimmed using an homemade fully integrated fastq filter parallel procecing module written in python and C.
* If needed an index for bwa will be generated from the modified reference files or from the original one after being merged together in a temporary directory.
* Fastq sequences are then aligned against the bwa merged reference genome index and a temporary sam file is generated
* Aligned reads from the sam file are splited and attributed to the reference Object for which a hit was found. or to one of the following garbage reads categories: unmapped, lowMapq, secondary.
* Each reference will then generates the output required in the configuration file (Bam, sam, bedgraph, bed, covgraph and/or variant report).
* Finally a distribution report and a log file are generated 

For more information, a comprehensive developper documentation can be generated from ContaVect.dox using [Doxygen](https://github.com/doxygen/doxygen) with [doxypy](https://github.com/0xCAFEBABE/doxypy).

## Dependencies:

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

2. Enter the root of the program folder and make the main script excecutable
``` bash
$ sudo chmod u+x ContaVect.py
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
$ ./ContaVect.py Conf.txt 
```
No command line option available, everything is in the Configuration file

## Pipeline development logbook
* [Logbook](http://nbviewer.ipython.org/github/a-slide/ContaVect/blob/master/doc/Logbook.ipynb)

## Authors and Contact

Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)




