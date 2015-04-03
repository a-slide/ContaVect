# -*- coding: utf-8 -*-

"""
@package    Quade
@brief      Contain the template of the empty configuration file for Quade
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

def write_example_conf():

    with open ("ContaVect_conf_file.txt", 'wb') as fp:
        fp.write ("""
###################################################################################################
#                                CONTAVECT CONFIGURATION FILE                                     #
###################################################################################################
# Values can by customized with users values, but the file template must remain unchanged,
# otherwise the program will not be able to load default values.
# File path should be indicated as absolute path preferably and should not contain blank spaces
# Values identified with '**' in the descriptor are not recommended to be modified

###################################################################################################
# REFERENCE DEFINITION

# These sections corresponds to the definition of References sequence objects. It is possible to
# include as many independent reference as required by duplicating a entire Ref section. If you
# performed a Reference Masking step before or if you provided an existing BWA index, indicate each
# individual reference used. All references have to be unique and to be organized as follow :
#   * [referenceX] = Reference identifier section, where X is the reference id starting from 1 for
#     the first reference and incrementing by 1 for each additional reference
#   * name = A name to define the reference (STRING)
#   * fasta = A path to a fasta file that can contain multiple sequences, can be gziped (STRING)
#   * output_type = A list of files to be generated for this reference among the following choices:
#     bam, sam, bedgraph, bed and covgraph separated by blank spaces (See Readme for more
#     details about these output format) (LIST OF STRING)

[reference1]
name : AAV
fasta : ../
output_type : bam    sam    bedgraph    bed    covgraph

[reference2]
name : Backbone
fasta : ../
output_type : bam    sam    bedgraph

###################################################################################################
# SAMPLE DEFINITIONS

# These sections corresponds to the definition of Samples to be analyzed against the references.
# Similar to Reference Sections, sample sections can be duplicated to include as many independant
# sample as required. All sample name and index have to be unique and to be organized as follow :
#   * [sampleX] = Sample identifier section, where X is the sample id number starting from 1 for the
#     first sample and incrementing by 1 for each additional sample
#   * name = Unique identifier that will be used to prefix the read files (STRING)
#   * R1_path = Valid path to the fastq(.gz) file containing the forward reads of the pair (STRING)
#   * R2_path = Valid path to the fastq(.gz) file containing the reverse reads of the pair (STRING)

[sample1]
name : S1
R1_path : ../
R2_path : ../

[sample2]
name : S2
R1_path : ../
R2_path : ../

###################################################################################################
[Bwa_Alignment]

# This section group options to parameter bwa index and alignment with bwa mem

# Path of a bwa index for the reference indicated above. If not available just leave parameter
# blank and a new index will be generated (STRING)
bwa_index :

# Command line options that will be used by bwa mem aligner and bwa index (STRING) **
bwa_mem_opt : -M
bwa_index_opt :

# Path to bwa mem and bwa index executable (not required if added to system path) (STRING) **
bwa_aligner  : bwa mem
bwa_indexer : bwa index

###################################################################################################
[Output]

# Minimal MAPQ score and size of alignment of a read to be considered valid (POSITIVE INTEGER)
min_mapq : 30
min_size : 25

# Write Garbage reads (unmapped, low mapq and secondary) in BAM and or SAM files (BOOLEAN)
unmapped_bam = True
unmapped_sam = False
""")
