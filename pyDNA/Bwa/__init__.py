# -*- coding: utf-8 -*-
"""
@package Bwa
@brief **Wrapper for BWA mem**
Please see the BWA user manual document for further details
[MANUAL](http://bio-bwa.sourceforge.net/bwa.shtml)
Basically the top level function Mem.align processes sequences as follow:

* If a bwa index is provided it will attempt to validate it by creating and IndexWrapper.ExistingIndex object.
* If the validation fails of if no index was provided a new index will be created by using IndexWrapper.ExistingIndex from a reference fasta file
(or a list a fasta files that will be combined in a single reference)
* A instance of MemWrapper.Aligner is then created by passing the Index object as an argument.
* A single or a pair of fastq files are then aligned against the reference through MemWrapper.Aligner
* Finally, results are piped into a sam file

@copyright [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""
__all__ = ["Mem", "IndexWrapper", "MemWrapper"]
