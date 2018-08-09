# -*- coding: utf-8 -*-

"""
@package    ContaVect
@brief      Contain the template of the empty configuration file for ContaVect
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@authors     Adrien Leger - 2014 & Mathieu Bolteau - 2018
* [Github](https://github.com/a-slide)
* [Github](https://github.com/mablt)
* INSERM 1089 - University of Nantes
"""

def write_conf_file (number_reference):
	"""
	Create an empty configuration file for ContaVect. By default, the file contains one 
	reference. 
	"""

	# Initialize variables
	counter = 1
	reference_to_add = str()

	# Write parameters for reference informations
	while counter<=number_reference-1:
		reference_to_add += """
[Ref{}]
name : 
fasta :
output :\n\n""".format(str(counter+1))
		counter+=1

	# Write the empty configuration file
	with open ("Conf.txt", 'w') as f:
		f.write ("""
#########################################################################################
#                                                                                       #
#                          CONTAVECT CONFIGURATION FILE                                 #
#                                                                                       #
#########################################################################################                                                    

#########################################################################################
#                                 INSTRUCTIONS      									#
#########################################################################################

#   Values can by customized with users values, but the file template must remain
#   unchanged, otherwise the program will not be able to load default values.
#   Please follow the following recommendations :
#   
#   - File path should be indicated as absolute path preferably and should not contain
#     blank spaces  
#   - Values identified with '**' in the descriptor are not recommended to be modified              
#                                                                                       


#########################################################################################
#                                 PARAMETERS      										#
#########################################################################################

#########################
#    GENERAL SECTION    #
#########################

#   General parameters of the program 
[General]

#   Path where to create a folder in which all files generated during the program will
#   be stored. Can be leaved empty if you want to work in the current directory (string)
outdir : 

#   String that will be used to prefix all data generated by the program. Cannot
#   contain any '\\' character. If leaved blank 'out' will be used as a prefix (string)
outprefix : 

###############################
#    REFERENCE DEFINITIONS    #
###############################

#   These sections corresponds to the definition of References sequence objects. They
#   are defined by a name, a path to a fasta file that can contain multiple sequences, 
#   and finally a list of output files for this reference. To be properly interpreted
#   Items in the list needs to be choose among :'bam', 'sam', 'bedgraph', 'bed',
#   'covgraph'.

#   BAM and SAM correspond to the standard BAM/SAM format. They will be sorted and with
#   a bam index. BAM is much more compact than SAM but not human readeable.

#   BEDGRAPH and BED are data recapitulating the coverage of reads over the sequence.
#   bed report the coverage of all positions in the reference while bedgraph is more
#   concise since uncovered regions are not reported and contigous positions of same
#   coverage are reported only once. for large sequences bedgraph only should be used.

#   COVGRAPH is a graphical representation of the coverage as a svg file and should be
#   used only for small sized reference (> 50 000 pb). If selected the program will
#   generate a graphics per sequence mapped in the reference fasta file.

#   It is possible to include as many independant reference as required by duplicating
#   a entire Ref section and updating the id number (Ref1, Ref2, Ref3 ... ). The order
#   of references is CRITICAL if you want to perform a reference masking since it will
#   start by the last reference masked by all the others then the penultimate masked by
#   all others except the last and and so on until there is only one reference remaining

[Ref1]
#   Name of the reference (string)
name : 

#   Path to the reference fasta file that can contain multiple sequences.
#   Can be gzipped (string)
fasta : 

#   List of output required for the reference separated by any blank space. See above
#   for valid values (list of strings)
output : 

"""
+reference_to_add+
"""
#####################
#    FASTQ FILES    #
#####################

#   This program is dedicated to paired en data only and will not work if a pair of fastq
#   is not provided. Fastq can be pretreated for quality or adapter trimming or it can
#   be done directly within the pipeline with an homemade efficient module.
[Fastq]

# Path to the file containing the forward reads (string)
R1 : 
# Path to the file containing the reverse reads (string)
R2 : 

#################################################
#    REFERENCE MASKING OPTIONS (FACULTATIVE)    #
#################################################

#   This section group options to parameters a facultative step of homologies masking
#   in ref as explain in the program documentation and above. It is based on the retrieval
#   of similarities using blastn from blast+ followed by an hard masking of hits before
#   rewritting of a new reference fasta file 
[Ref_Masking]

#   Flag to activate or deactivate the reference masking. If False the following options
of the section won't be parsed (Boolean) 
ref_masking : 

# ** Command line options that will be used by blastn for alignment (string) 
blastn_opt : -evalue 0.1
blastn_threads : 2

# ** Command line options that will be used by mkblastdb for database creation (string)
mkblastdb_opt : 

# ** Path to the blastn executable (not required is added to system path) (string)
blastn : blastn

# ** Path to the mkblastdb executable (not required is added to system path) (string)
mkblastdb : mkblastdb

###############################################
#    FASTQ FILTERING OPTIONS (FACULTATIVE)    #
###############################################

#   This section group options to parameters a facultative step of fastq filtering
#   in ref as explain in the program documentation and above. 2 types o filtering can be  
#   made together or individually : A filtering of bad quality reads and a trimming of
#   primers homologies found in reads. 
[Fastq_Filtering]

#   Quality format associated with reads autorized values are : 'fastq-sanger' 'fastq-solexa' or 
#   'fastq-illumina'. See Biopython QualityIO documentation for more details. (string)
input_qual : fastq-sanger

#   If True a quality filtering base on the general read quality will be performed (Boolean)
quality_filtering : 

#   Minimal mean quality of a read to be kept for alignment (Integer)
min_qual : 25

#   If True a step of adaptor trimming will be performed (Boolean)
adapter_trimming : 

#   List of adapter sequences separated by a blank space  (list of strings)
adapters : 

#   Try to find matched for complementary sequence of all adapters (Boolean)
find_rc : 

# ** Minimal fraction of the read that should remain after adapter trimming (Positive Float)
min_read_len : 0.6

# ** Minimal fraction of the adapter that have to be aligned on reads to be consider as a valid match (Positive Float)
min_match_len : 0.8

# ** Minimal score per base of the adapter to be consider as a valid match (Positive Float)
min_match_score : 1.4

# ** Bonus score for SSW alignment in case of a match (Positive Integer)
ssw_match : 2

# ** Malus score for SSW alignment in case of a match (Positive Integer)
ssw_mismatch : 2

# ** Malus score for SSW alignment in case of gap opening (Positive Integer)
ssw_gap_open : 3

# ** Malus score for SSW alignment in case of gap extension(Positive Integer)
ssw_gap_extend : 1

#############################
#   BWA ALIGNMENT OPTIONS   #
#############################

#   This section group options to parameters bwa index and alignment with mem 
[Bwa_Alignment]

#   Path of a bwa index for the reference indicated above. If not available just leave parameter
#   blank and a new index will be generated (string) 
bwa_index : 

# ** Command line options that will be used by bwa mem for alignment (string) 
bwa_mem_opt : -M
bwa_threads : 2

# ** Command line options that will be used by bwa incdex for reference indexation (string) 
bwa_index_opt : 

# ** Path to bwa mem executable (not required is added to system path) (string)
bwa_aligner  : bwa mem

# ** Path to bwa index executable (not required is added to system path) (string)
bwa_indexer : bwa index


#######################
#   OUTPUT OPTIONS    #
#######################

#   This section group Miscelianous options for configuring program output
[Output]

#   This value correspond to the minimal quality mapping of a read to be considered valid (Positive Integer)
min_mapq : 30

#   This value correspond to the minimal size of mapping (Positive integer)
min_size : 25

#   If True Garbage reads (unmapped, low mapq and secondary) will be written in bam files (Boolean)
unmapped_bam = 

#   If True Garbage reads (unmapped, low mapq and secondary) will be written in sam files (Boolean)
unmapped_sam = 

#   Minimal depth of coverage to report the region in bedgraph, bed and covgraph (Positive Integer)
cov_min_depth : 4

#   Minimal depth of coverage to report a variant position (Positive Integer)
var_min_depth : 500

#   Minimal fraquency of a base at a given position to be considered as frequent (Positive Float)
var_min_freq : 0.2""")
