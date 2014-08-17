"""
@package Bwa
@brief **Wrapper for BWA mem**
Please see the BWA user manual document for further details
[MANUAL](http://bio-bwa.sourceforge.net/bwa.shtml)
* To use the wrapper, it is first needed to generate a bwa index either by using NewIndex if no
index was already created or ExistingIndex if an index is already available.
* A instance of Aligner from MemWrapper can then be created by passing the index object as an argument.
* The Aligner method align can finally be used as many times as desired with different queries.
It will returned each time a path to a sam file

Alternatively the Mem module can be called to facilitate the usage of the wrapper.

@copyright [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""
__all__ = ["Mem", "IndexWrapper", "MemWrapper"]
