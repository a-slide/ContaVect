# -*- coding: utf-8 -*-
"""
@package Blast
@brief **Wrapper for Blastn from NCBI Blast**
Please see the Blast+ user manual document for further details
[MANUAL](http://www.ncbi.nlm.nih.gov/books/NBK1762/)
Basically the top level function Blastn.align processes sequences as follow:

* If a blastn database is provided it will attempt to validate it by creating and MakeblastdbWrapper.ExistingDB object.
* If the validation fails of if no database was provided a new database will be created by using MakeblastdbWrapper.NewDB from a reference fasta file
* A instance of BlastnWrapper.Aligner is then created by passing the DB object as an argument.
* A list of fasta file (each can contains many sequences) is then submitted one by one to the balst database though  BlastnWrapper.Aligner
* For each query fasta a list of BlastHit objects containing informations for blast hits found will be returned
* Hit lists are combined into a single flat list and returned at the end of Blastn.align execution.

@copyright [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""
__all__ = ["Blastn", "BlastnWrapper", "MakeblastdbWrapper", "BlastHit"]
