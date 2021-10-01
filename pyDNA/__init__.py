# -*- coding: utf-8 -*-
"""
@package    pyDNA
@brief      **Generic collection of python Utilities and Wrapper for DNA/NGS data manipulation for python 2.7**
For now it includes BWA, Blast and Smith–Waterman wrappers as well as higher level functionnalities:

* FastqFT : Filter fastq file based on quality and adapter trimming
* RefMasker : Align a fasta reference against several fastq queries and mask homologies in the reference file
* pySamTools :  Manipulate aligned reads manipulation though pysam
* Utilities : Library of simple generic functions to manipulate paths and file, interact with command line interpreter and so on

@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

__all__ = ["Utilities", "Blast", "Bwa", "FastqFT", "pySamTools", "Ssw", "RefMasker", "Ungzip"]
__version__ = 0.1
