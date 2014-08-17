"""
@package Blast
@brief **Wrapper for Blastn from NCBI Blast**
Please see the Blast+ user manual document for further details
[MANUAL](http://www.ncbi.nlm.nih.gov/books/NBK1762/)
* To use the wrapper, it is first needed to generate a blastn database either by using NewDB if no
database was already created or ExistingDB if a database is already available.
* A instance of Blastn can then be created by passing the DB object as an argument.
* The blastn method align can finally be used as many times as desired with different queries.
It will returned each time a list of BlastHit objects containg informations for each blast found.

@copyright [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""
__all__ = ["Blastn", "BlastnWrapper", "MakeblastdbWrapper", "BlastHit"]
