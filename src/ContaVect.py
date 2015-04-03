#!/usr/bin/env python2.7
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

# IMPORTS
try:
    # Standard library packages import
    from os import path, remove
    from time import time
    from datetime import datetime
    import ConfigParser
    import optparse
    import csv

    # Third party packages
    import pysam # Mandatory package
    import Bio # Mandatory package

    # Local Package import
    from pyDNA.FileUtils import is_readable_file, rm_blank
    from Sample import Sample
    from Reference import Reference
    from Sequence import Sequence

    #from pyDNA.Utilities import mkdir, file_basename, file_name, expand_file, rm_blank, is_gziped # Mandatory package
    #from pyDNA.Blast import Blastn # if not imported = not ref masking
    #from pyDNA.RefMasker import mask # if not imported = not ref masking
    #from pyDNA.FastqFT.FastqFilter import FastqFilter # if not imported = not fasta filter
    #from pyDNA.FastqFT.QualityFilter import QualityFilter # if not imported = not fasta filter
    #from pyDNA.FastqFT.AdapterTrimmer import AdapterTrimmer # if not imported = not fasta filter
    #from pyDNA.Ssw import ssw_wrap # if not imported = not fasta filter
    #from pyDNA.Bwa import Mem # Mandatory package
    #from pyDNA.pySamTools import Bam, Coverage, Variant # if not imported = not requested output
    #from ContaVect_src.Reference import Reference, Sequence # Mandatory package

except ImportError as E:
    print (E)
    print ("Please verify your dependencies. See Readme for more informations\n")
    sys.exit()

#~~~~~~~MAIN FUNCTION~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class ContaVect(object):
    """
    Main class of the program. In a first time the class is initialize with values parsed from the
    configuration file. Then, the pipeline of analysis is launch through the 'run' method
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~CLASS FIELDS~~~~~~~#

    VERSION = "ContaVect 0.3"
    USAGE = "Usage: %prog -c Conf.txt [-i -h]"

    #~~~~~~~CLASS METHODS~~~~~~~#

    @classmethod
    def class_init (self):
        """
        init class method for instantiation from command line. Parse arguments parse CL arguments
        """

        # Define parser usage, options
        optparser = optparse.OptionParser(usage = self.USAGE, version = self.VERSION)
        optparser.add_option('-c', dest="conf_file",
            help= "Path to the configuration file [Mandatory]")
        optparser.add_option('-i', dest="init_conf", action='store_true',
            help= "Generate an example configuration file and exit [Facultative]")

        # Parse arguments
        options, args = optparser.parse_args()

        return ContaVect(options.conf_file, options.init_conf)

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self, conf_file=None, init_conf=None):
        """
        Initialization function, parse options from configuration file and verify their values.
        All self.variables are initialized explicitly in init.
        """

        # Create a example conf file if needed
        if init_conf:
            print("Create an example configuration file in the current folder")
            write_example_conf()
            sys.exit(0)

        print("Initialize Quade\n")
        # Parse the configuration file and verify the values of variables
        try:

            # Verify if conf file was given and is valid
            assert conf_file, "A path to the configuration file is mandatory"
            is_readable_file(conf_file)
            self.conf = conf_file

            # Define a configuration file parser object and load the configuration file
            cp = ConfigParser.RawConfigParser(allow_no_value=True)
            cp.read(self.conf)

            print("\tParse Reference sequences")
            # Iterate only on sections starting by "reference", create Reference objects
            # And store them in a list
            self.reference_list = []
            for reference in [i for i in cp.sections() if i.startswith("reference")]:
                # Create Reference objects
                self.reference_list.append (
                    Reference (
                        name = rm_blank(cp.get(reference, "name"), replace ='_'),
                        fasta = rm_blank(cp.get(reference, "fasta"), replace ='\ '),
                        output_type = cp.get(reference, "output_type").split()))

            print("\tParse Samples")
            # Iterate only on sections starting by "sample", create Sample objects
            # And store them in a list
            self.sample_list = []
            for sample in [i for i in cp.sections() if i.startswith("sample")]:
                # Create Sample objects
                self.sample_list.append (Sample (
                    name = rm_blank(cp.get(sample, "name"), replace ='_'),
                    R1_path = rm_blank(cp.get(sample, "R1_path"), replace ='\ '),
                    R2_path = rm_blank(cp.get(sample, "R2_path"), replace ='\ '))

            print("\tParse BWA options")
            # Bwa options
            self.bwa_mem_opt = cp.get("Bwa_Alignment", "bwa_mem_opt")
            self.bwa_index_opt = cp.get("Bwa_Alignment", "bwa_index_opt")
            self.bwa_aligner = cp.get("Bwa_Alignment", "bwa_aligner")
            self.bwa_indexer = cp.get("Bwa_Alignment", "bwa_indexer")
            self.bwa_index = rm_blank(cp.get("Bwa_Alignment", "bwa_index"), replace ='\ ')

            if not bwa_index:
                bwa_index = make_index ([reference.fasta for reference in self.reference_list]) ################################## A VOIR SI INDEXATION ICI OU DANS CALL

            print("\tParse output options")
            # Output options
            self.min_mapq = cp.getint("Output", "min_mapq")
            self.min_size = cp.getint("Output", "min_size")
            self.unmapped_bam = cp.getboolean("Output", "unmapped_bam")
            self.unmapped_sam = cp.getboolean("Output", "unmapped_sam")

            assert self.min_mapq
            assert self.min_size ###########################################################################"

        # Handle the many possible errors occurring during conf file parsing or variable test
        except (ConfigParser.NoOptionError, ConfigParser.NoSectionError) as E:
            print ("Option or section missing. Report to the template configuration file\n" + E.message)
            sys.exit(1)
        except (ValueError, AssertionError) as E:
            print ("One of the value in the configuration file is not correct\n" + E.message)
            sys.exit(1)
        except (IOError) as E:
            print ("One of the file is incorrect or unreadable\n" + E.message)
            sys.exit(1)

        print ("\tAll configuration file parameters are valid\n")

    def __str__(self):
        msg = "CONTAVECT CLASS\n\tParameters list\n"
        # list all values in object dict in alphabetical order
        keylist = [key for key in self.__dict__.keys()]
        keylist.sort()
        for key in keylist:
            msg+="\t{}\t{}\n".format(key, self.__dict__[key])
        return (msg)

    def __repr__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def __call__(self):
        """
        Launch the complete pipeline of analyse:
        * BWA indexation if needed
        * Short read alignment with bwa mem
        * Spliting of sam to attribute reads to each original references (or unmmapped)
        * Output per reference bam, sam, bedgraph, bed, covgraph, variant call
        * Output distribution table and graph
        """
        stime = time()

        for sample in self.sample_list:

            outdir = mkdir(sample.name)

            print ("\n##### PARSE REFERENCES #####\n")

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

    contavect = ContaVect.class_init()
    contavect()
