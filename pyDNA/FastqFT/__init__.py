# -*- coding: utf-8 -*-
"""
@package    FastqFT
@brief	    **Module for fastq files adapters trimming and quality filtering**
for NGS paired end data. Basically the top level function FastqFilter process fastq as follow:

* Read a fastq file and add sequence pairs in a queue
* Process fastq pairs from the queue in parrallel threads through a Quality Filter Object and/or an Adapter Trimmer Object. If the sequences were not fitered out add them to a second queue
* Collect filterer pairs from the queue and write them back in a new fastq file.

@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

__all__ = ["PairwiseAligner", "AdapterTrimmer", "QualityFilter","FastqFilter"]
