#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
@package    Contavect
@brief      **Main Class of Contavect program**. Contains The Main class and top level instructions
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

# TODO Raise flag if a facultative package is not properly imported
	
# IMPORTS
try:
	# Standard library packages import
	from os import path, remove # Mandatory package
	from time import time # Mandatory package
	import ConfigParser # Mandatory package
	from sys import argv # Mandatory package
	import csv # Mandatory package
	import optparse #Mandatory package

	# Third party packages
	import pysam # Mandatory package
	import Bio # Mandatory package

	# Local Package import
	from pyDNA.Utilities import mkdir, file_basename, file_name, expand_file, rm_blank, is_gziped # Mandatory package
	from pyDNA.Blast import Blastn # if not imported = not ref masking
	from pyDNA.RefMasker import mask # if not imported = not ref masking
	from pyDNA.FastqFT.FastqFilter import FastqFilter # if not imported = not fasta filter
	from pyDNA.FastqFT.QualityFilter import QualityFilter # if not imported = not fasta filter
	from pyDNA.FastqFT.AdapterTrimmer import AdapterTrimmer # if not imported = not fasta filter
	from pyDNA.Ssw import ssw_wrap # if not imported = not fasta filter
	from pyDNA.Bwa import Mem # Mandatory package
	from pyDNA.pySamTools import Bam, Coverage, Variant # if not imported = not requested output
	from ContaVect_src.Reference import Reference, Sequence # Mandatory package
	from Conf_file import write_example_conf #if not  imported = not creation of configuration file

except ImportError as E:
	print (E)
	print ("Please verify your dependencies. See Readme for more informations\n")
	exit()

#~~~~~~~MAIN FUNCTION~~~~~~~#


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Main(object):
	"""
	Main class of the program. In a first time the class is initialize with values parsed from the
	configuration file. Then, the pipeline of analysis is launch through the 'run' method
	"""
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

	#~~~~~~~CLASS FIELDS~~~~~~~#

	VERSION = "ContaVect 0.2.1"
	USAGE = "Usage: %prog -c Conf.txt [-i -h]"

	#~~~~~~~CLASS METHODS~~~~~~~#

	@classmethod
	def class_init (self):
		"""
        init class method for instantiation from command line. Parse arguments parse CL arguments
        """

        # Define parser usage, options
		parser = optparse.OptionParser(usage = self.USAGE, version = self.VERSION)
		parser.add_option('-i', dest="number_ref", type="int", 
			help="Generate an example configuration file with the number of references necessary and exit [Mandatory]")
		parser.add_option('-c', dest="conf_file",
			help="Path to the configuration file [Mandatory]")

		# Parse arguments
		options, args = parser.parse_args()

		return Main(options.number_ref, options.conf_file)

	#~~~~~~~FONDAMENTAL METHODS~~~~~~~#

	def __init__ (self, number_ref=None, conf_file=None):
		"""
		Parse command line argument and configuration file, then verify a limited number of
		critical values.
		"""
		
		# Create a example conf file if needed with the number of reference necessary
		if number_ref:
			print "Create an example configuration file with",number_ref,"references in the current folder"
			write_example_conf(number_ref)
			exit()
		
		# TODO manage import flags define default values and verify file existence
		# CL ARGUMENTS AND CONF FILE PARSING

		self.conf = ConfigParser.RawConfigParser(allow_no_value=True)
		self.conf.read(conf_file)
		if not self.conf.sections():
			print ("Empty or invalid configuration file. See Readme for more informations\n")
			exit()
		try:
			# Mandatory paramaters
			self.outdir = rm_blank(self.conf.get("General", "outdir"), replace="_")
			if not self.outdir:
				self.outdir = "./"
			if self.outdir[-1] != "/":
				self.outdir += "/"
			self.outprefix = rm_blank(self.conf.get("General", "outprefix"), replace="_")
			if not self.outdir:
				self.outdir = "out"

			self.ref_masking = self.conf.getboolean("Ref_Masking", "ref_masking")
			self.R1 = rm_blank(self.conf.get("Fastq", "R1"), replace="\ ")
			self.R2 = rm_blank(self.conf.get("Fastq", "R2"), replace="\ ")
			self.input_qual = self.conf.get("Fastq_Filtering", "input_qual")
			self.quality_filtering = self.conf.getboolean("Fastq_Filtering", "quality_filtering")
			self.adapter_trimming = self.conf.getboolean("Fastq_Filtering", "adapter_trimming")
			self.bwa_index = rm_blank(self.conf.get("Bwa_Alignment", "bwa_index"), replace="\ ")
			self.bwa_mem_opt = self.conf.get("Bwa_Alignment", "bwa_mem_opt")
			self.bwa_threads = self.conf.get("Bwa_Alignment", "bwa_threads")
			self.bwa_index_opt = self.conf.get("Bwa_Alignment", "bwa_index_opt")
			self.bwa_aligner = self.conf.get("Bwa_Alignment", "bwa_aligner")
			self.bwa_indexer = self.conf.get("Bwa_Alignment", "bwa_indexer")
			self.min_mapq = self.conf.getint("Output", "min_mapq")
			self.min_size = self.conf.getint("Output", "min_size")
			self.unmapped_bam = self.conf.getboolean("Output", "unmapped_bam")
			self.unmapped_sam = self.conf.getboolean("Output", "unmapped_sam")
			self.cov_min_depth = self.conf.getint("Output", "cov_min_depth")
			self.var_min_depth = self.conf.getint("Output", "var_min_depth")
			self.var_min_freq = self.conf.getfloat("Output", "var_min_freq")

			# Conditional paramaters
			if self.ref_masking:
				self.blastn_opt = self.conf.get("Ref_Masking", "blastn_opt")
				self.blastn_threads = self.conf.get("Ref_Masking", "blastn_threads")
				self.mkblastdb_opt = self.conf.get("Ref_Masking", "mkblastdb_opt")
				self.blastn = self.conf.get("Ref_Masking", "blastn")
				self.mkblastdb = self.conf.get("Ref_Masking", "mkblastdb")
			if self.quality_filtering:
				self.min_qual = self.conf.getint("Fastq_Filtering", "min_qual")
			if self.adapter_trimming:
				self.adapters = self.conf.get("Fastq_Filtering", "adapters").split()
				self.find_rc = self.conf.getboolean("Fastq_Filtering", "find_rc")
				self.min_read_len = self.conf.getfloat("Fastq_Filtering", "min_read_len")
				self.min_match_len = self.conf.getfloat("Fastq_Filtering", "min_match_len")
				self.min_match_score = self.conf.getfloat("Fastq_Filtering", "min_match_score")
				self.ssw_match = self.conf.getint("Fastq_Filtering", "ssw_match")
				self.ssw_mismatch = self.conf.getint("Fastq_Filtering", "ssw_mismatch")
				self.ssw_gap_open = self.conf.getint("Fastq_Filtering", "ssw_gap_open")
				self.ssw_gap_extend = self.conf.getint("Fastq_Filtering", "ssw_gap_extend")

			# More complicated import in a list of dictionnary for references informations
			self.raw_ref_list =[]

			for i in range(1,100):
				ref_id = "Ref"+str(i)
				if not self.conf.has_section(ref_id):
					break
				ref = { 'name'   : rm_blank(self.conf.get(ref_id, "name"), replace="_"),
						'fasta'  : rm_blank(self.conf.get(ref_id, "fasta"), replace="\ "),
						'output' : self.conf.get(ref_id, "output").split(),}
				self.raw_ref_list.append(ref)

		except ConfigParser.NoOptionError as E:
			print (E)
			print ("An option is missing in the configuration file")
			print ("Please report to the descriptions in the configuration file\n")
			exit()
		except ConfigParser.NoSectionError as E:
			print (E)
			print ("An section is missing in the configuration file")
			print ("Please report to the descriptions in the configuration file\n")
			exit()
		except ValueError as E:
			print (E)
			print ("One of the value in the configuration file is not in the correct format")
			print ("Please report to the descriptions in the configuration file\n")
			exit()

	def __repr__(self):
		msg = "MAIN CLASS\n"
		msg+= "\tParameters list\n"
		for i, j in self.__dict__.items():
			msg+="\t{}\t{}\n".format(i, j)
		return (msg)

	def __str__(self):
		return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

	def get(self, key):
		return self.__dict__[key]

	def set(self, key, value):
		self.__dict__[key] = value

	#~~~~~~~PUBLIC METHODS~~~~~~~#

	def __call__(self):
		"""
		Launch the complete pipeline of analyse:

		* Reference importation/parsing
		* Facultative step of reference masking to remove homologies between reference sequences
		* Facultative step of Fastq quality Filtering/ adapter trimming
		* Facultative step of reference indexing for bwa from merged references
		* Short read alignment with bwa mem
		* Spliting of sam to attribute reads to each original references (or unmmapped)
		* Output per reference bam, sam, bedgraph, bed, covgraph, variant call
		* Output distribution table and graph
		"""
		stime = time()
		self.outdir = mkdir(path.abspath(self.outdir))

		print ("\n##### PARSE REFERENCES #####\n")
		# Create CV_Reference.Reference object for each reference easily accessible through
		# Reference class methods
		
		if self.ref_masking or not self.bwa_index:
			self.ref_dir = mkdir(path.join(self.outdir, "references/"))
			self.index_dir = mkdir(path.join(self.outdir, "bwa_index/"))
			self._extract_ref(expand=True)
		else:
			self.ref_dir = ""
			self.index_dir = ""
			self._extract_ref(expand=False)
		
		# Reference Masking
		if self.ref_masking:
			print ("\n##### REFERENCE HOMOLOGIES MASKING #####\n")
			self.db_dir = mkdir(path.join(self.outdir, "blast_db/"))
			ref_list = self._iterative_masker()
			# Erase existing index value if ref masking was performed
			bwa_index = None

		# Fastq Filtering
		if self.quality_filtering or self.adapter_trimming:
			print ("\n##### FASTQ FILTERING #####\n")
			self.fastq_dir = mkdir(path.join(self.outdir, "fastq/"))
			self.R1, self.R2 = self._fastq_filter()

		# BWA alignment
		print ("\n##### READ REFERENCES AND ALIGN WITH BWA #####\n")
		# An index will be generated if no index was provided
		self.result_dir = mkdir(path.join(self.outdir, "results/"))
		
		self.sam = Mem.align (
			self.R1, self.R2,
			index = self.bwa_index,
			ref = Reference.allFasta(),
			align_opt = self.bwa_mem_opt,
			index_opt = self.bwa_index_opt,
			aligner = self.bwa_aligner,
			align_threads = self.bwa_threads,
			indexer = self.bwa_indexer,
			align_outdir = self.result_dir,
			index_outdir = self.index_dir,
			align_outname = self.outprefix+".sam",
			index_outname = self.outprefix+".idx")

		print ("\n##### FILTER ALIGNED READS AND ASSIGN A REFERENCE #####\n")
		# Split the output sam file according to each reference
		self._sam_spliter ()

		print ("\n##### GENERATE OUTPUT FOR EACH REFERENCE #####\n")
		# Deal with garbage read dictionnary
		self._garbage_output()
		# Ask references to generate the output they were configured to
		Reference.mk_output_global(self.result_dir+self.outprefix)
		# Create a distribution table
		self._distribution_output()
		self._make_report()
		
		print ("\n##### DONE #####\n")
		print ("Total execution time = {}s".format(round(time()-stime, 2)))

	##~~~~~~~PRIVATE METHODS~~~~~~~#

	def _extract_ref(self, expand=True):
		"""
		Import and expand fasta references and associated flags in a Reference object
		expand the file if Gziped to avoid multiple compression/decompression during execution
		if require for next operations
		"""
		for ref in self.raw_ref_list:
			# Expand fasta if needed
			if expand:
				ref_fasta = expand_file(infile=ref['fasta'], outdir=self.ref_dir)
			else:
				ref_fasta = ref['fasta']
			
			# Create a Reference object
			Ref = Reference(
				name = ref['name'],
				ref_fasta = ref_fasta,
				bam_maker = Bam.BamMaker(
					make_bam = 'bam' in ref['output'],
					make_sam = 'sam' in ref['output']),
				cov_maker = Coverage.CoverageMaker(
					min_depth=self.cov_min_depth,
					make_bedgraph = 'bedgraph' in ref['output'],
					make_bed = 'bed' in ref['output'],
					make_covgraph = 'covgraph' in ref['output']),
				var_maker = Variant.VariantMaker(
					min_depth=self.var_min_depth,
					min_freq=self.var_min_freq,
					make_freqvar = 'variant' in ref['output']))
			
			## Test if all seq in ref are longer than 3000 for compatibility with bwa 
			#for seq in Ref.seq_dict.values():
				#if seq.length < 3000:
					#import_and_pad (
		
			print (repr(Ref))

	def _iterative_masker (self): #### TODO The fuction directly manipulate reference field= change that
		"""
		Mask references homologies iteratively, starting by the last reference which is masked by
		all the others then to the penultimate masked by all others except the last and and so
		forth until there is only 1 reference remaining
		"""
		# Iterate over index in Reference.instances staring by the last one until the 2nd one
		for i in range(Reference.countInstances()-1, 0, -1):

			# Extract subject and query_list from ref_list
			subject = Reference.Instances[i]
			query_list = Reference.Instances[0:i]
			print ("\n# PROCESSING REFERENCE {} #\n".format(subject.name))

			# Perform a blast of query list against subject
			hit_list = Blastn.align (
				query_list = [ref.ref_fasta for ref in query_list],
				subject_fasta = subject.ref_fasta,
				align_opt = self.blastn_opt,
				num_threads = self.blastn_threads,
				db_opt = self.mkblastdb_opt,
				db_outdir = self.db_dir,
				db_outname = subject.name)

			# Masking hits in suject fasta if hits in hit_list
			subject.ref_fasta = mask (
				subject_fasta= subject.ref_fasta,
				hit_list = hit_list,
				ref_outdir = self.ref_dir,
				ref_outname = "masked_{}.fa".format(subject.name),
				compress_ouput = False)

	def _fastq_filter (self):
		"""
		Filter fastq with FastqFilterPP
		"""
		# Define a quality filter object

		if self.quality_filtering:
			self.qFilter = QualityFilter (self.min_qual)
		else:
			self.qFilter = None

		# Define a adapter trimmer object
		if self.adapter_trimming:
			
			self.ssw_aligner = ssw_wrap.Aligner(
				match = self.ssw_match,
				mismatch = self.ssw_mismatch,
				gap_open = self.ssw_gap_open,
				gap_extend = self.ssw_gap_extend)
			
			self.trimmer = AdapterTrimmer(
				Aligner = self.ssw_aligner,
				adapters = self.adapters,
				find_rc = self.find_rc,
				min_read_len = self.min_read_len,
				min_match_len = self.min_match_len,
				min_match_score = self.min_match_score)
		else:
			self.trimmer = None

		# Filter fastq for quality and adapter with FastqFilter
		self.fFilter = FastqFilter (
			self.R1, self.R2,
			quality_filter = self.qFilter,
			adapter_trimmer = self.trimmer,
			outdir = self.fastq_dir,
			input_qual = self.input_qual,
			compress_output=False)

		# Print a simple result
		print ("Pairs processed: {}\t Pairs passed : {}\t in {} s".format(
			self.fFilter.getCTypeVal('total'),
			self.fFilter.getCTypeVal('total_pass'),
			self.fFilter.get('exec_time')))

		# Write a detailed report in a logfile
		output = "{}{}_FastqFilter_report.txt".format(self.fastq_dir, self.outprefix)
		with open (output, "wb") as outfile:
			outfile.write(repr(self.fFilter))

		return self.fFilter.getTrimmed()

	def _sam_spliter (self):
		"""
		"""
		with pysam.Samfile(self.sam, "r" ) as samfile:
			self.bam_header = samfile.header

			# Give the header of the sam file to all Reference.Instances to respect the same order
			# references in sorted bam files
			Reference.set_global("bam_header", self.bam_header)

			# Create a dict to collect unmapped and low quality reads
			Secondary = Sequence (name = 'Secondary', length = 0)
			Unmapped = Sequence (name = 'Unmapped', length = 0)
			LowMapq = Sequence (name = 'LowMapq', length = 0)
			self.garbage_read = [Secondary, Unmapped, LowMapq]

			for read in samfile:
				# Always remove secondary alignments
				if read.is_secondary:
					Secondary.add_read(read)
				# Filter Unmapped reads
				elif read.tid == -1:
					Unmapped.add_read(read)
				# Filter Low MAPQ reads
				elif read.mapq < self.min_mapq:
					LowMapq.add_read(read)
				# Filter short map ##### FOR FUTURE CREATE A SEPARATE CATEGORY
				elif len(read.query_alignment_sequence) < self.min_size:
					Unmapped.add_read(read)    
				# Finally if all is fine attribute the read to a Reference
				else:
					Reference.addRead(samfile.getrname(read.tid), read)
		
		# Removing the original sam file which is no longer needed
		remove(self.sam)
		self.sam = None

	def _garbage_output (self):
		"""
		Output bam /sam for garbage reads
		"""

		# Define a generic Bam.BamMaker resulting in unsorted bam/sam for all garbage reads
		bam_maker = Bam.BamMaker(
			sort=False,
			make_index=False,
			make_bam = self.unmapped_bam,
			make_sam = self.unmapped_sam)

		for seq in self.garbage_read:
			print "Processing garbage reads :{}\tReads aligned :{} ".format(seq.name, seq.nread)
			bam_maker.make(
				header = self.bam_header,
				read_col = seq.read_list,
				outpath = self.result_dir+self.outprefix,
				ref_name = seq.name)

	def _distribution_output (self):
		"""
		"""
		output = "{}{}_Reference_distribution.csv".format(self.result_dir, self.outprefix)
		with open(output, 'wb') as csvfile:
			writer = csv.writer(csvfile, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
			# Table for all reference
			writer.writerow(["Ref name","length","nread","RPKB"])
			for ref in Reference.getInstances():
				writer.writerow([ref.name, len(ref), ref.nread, float(ref.nread)/len(ref)*1000])
			# Add a line for garbage reads excluding the secondary alignments
			nread = sum([seq.nread for seq in self.garbage_read[1:]])
			writer.writerow(["Unmaped_and LowMapq","NA",nread,"NA"])

		output = "{}{}_Sequence_distribution.csv".format(self.result_dir, self.outprefix)
		with open(output, 'wb') as csvfile:
			writer = csv.writer(csvfile, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
			# Table decomposing Sequence per Reference
			writer.writerow(["Seq name","length","nread","RPKB"])
			for ref in Reference.getInstances():
				for seq in ref.seq_dict.values():
					writer.writerow([seq.name, len(seq), seq.nread, float(seq.nread)/len(seq)*1000])
			# Add a lines for garbage reads including the secondary alignments
			for seq in self.garbage_read:
				writer.writerow([seq.name, "NA", seq.nread, "NA"])
				
	def _make_report (self):
		"""
		"""
		output = "{}{}_parameters.txt".format(self.result_dir, self.outprefix)
		with open(output, 'wb') as outfile:
		
			# References options
			outfile.write("################## REFERENCES ##################\n\n")
			outfile.write(Reference.reprInstances()) 
			
			if self.ref_masking:
				outfile.write("Reference homologies were masked with RefMasker\n")
				outfile.write("blastn options : {}\n".format(self.blastn_opt))
				outfile.write("makeblastdb options : {}\n".format(self.mkblastdb_opt))
				
			else:
				outfile.write("No Reference homologies masking done\n")
			
			# Fastq options
			outfile.write("\n################## FASTQ FILES ##################\n\n")
			outfile.write("R1 : {}\n".format(self.R1))
			outfile.write("R2 : {}\n\n".format(self.R2))
			
			if self.quality_filtering or self.adapter_trimming:
				outfile.write(repr(self.fFilter)+"\n")
				if self.quality_filtering:
					outfile.write(repr (self.qFilter)+"\n")
				if self.adapter_trimming:
					outfile.write(repr (self.ssw_aligner)+"\n")
					outfile.write(repr (self.trimmer)+"\n")
			else:
				outfile.write("\nNo Fastq Filtering done\n")
				
			# bwa alignment options
			outfile.write("\n################## BWA ALIGNMENT ##################\n\n")
			outfile.write("index file : {}\n".format(self.bwa_index))
			outfile.write("bwa index options: {}\n".format(self.bwa_index_opt))
			outfile.write("bwa mem option: {}\n".format(self.bwa_mem_opt))
			outfile.write("bwa threads : {}\n".format(self.bwa_threads))
			
			# Output Options
			outfile.write("\n################## OUTPUT ##################\n\n")
			outfile.write("Minimal MAPQ score : {}\n".format(self.min_mapq))
			outfile.write("Write garbage reads to sam: {}\n".format(str(self.unmapped_sam)))
			outfile.write("Write garbage reads to bam: {}\n".format(str(self.unmapped_bam)))
			outfile.write("Minimal depth for Coverage output : {}\n".format(self.cov_min_depth))
			outfile.write("Minimal depth for Variant output : {}\n".format(self.var_min_depth))
			outfile.write("Minimal Variant frequency : {}\n".format(self.var_min_freq))

#~~~~~~~TOP LEVEL INSTRUCTIONS~~~~~~~#
if __name__ == '__main__':

	main = Main.class_init()
	main()
