#!/usr/bin/env python
# -*- coding: utf-8 -*-

# IMPORTS
try:
    # Standard library packages import
    from os import path
    from time import time
    import ConfigParser
    from sys import argv
    import csv

    # Third party packages
    import pysam
    import Bio

    # Local Package import
    from Utilities import mkdir, file_basename, file_name, expand_file, rm_blank
    from Blast import Blastn
    import RefMasker
    from FastqFT.FastqFilterPP import FastqFilterPP
    from FastqFT.QualityFilter import QualityFilter
    from FastqFT.AdapterTrimmer import AdapterTrimmer
    from Ssw import ssw_wrap
    from Bwa import Mem
    from CV_Reference import Reference, Sequence
    from pySamTools import Bam, Coverage2, Variant

except ImportError as E:
    print (E)
    print ("Please verify your dependencies. See Readme for more informations\n")
    exit()


#~~~~~~~MAIN FUNCTION~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Main(object):
    """
    @class  Reference
    @brief  Object oriented class containing informations of reference
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__ (self):
        """
        """

        # CL ARGUMENTS AND CONF FILE PARSING

        if len(argv) != 2:
            print ("Please provide the path to the Configuration file as an unique argument\n")
            exit()

        self.conf = ConfigParser.RawConfigParser(allow_no_value=True)
        self.conf.read(argv[1])
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
            self.bwa_index_opt = self.conf.get("Bwa_Alignment", "bwa_index_opt")
            self.bwa_aligner = self.conf.get("Bwa_Alignment", "bwa_aligner")
            self.bwa_indexer = self.conf.get("Bwa_Alignment", "bwa_indexer")
            self.min_mapq = self.conf.getint("Output", "min_mapq")
            self.unmapped_bam = self.conf.getboolean("Output", "unmapped_bam")
            self.unmapped_sam = self.conf.getboolean("Output", "unmapped_sam")
            self.cov_min_depth = self.conf.getint("Output", "cov_min_depth")
            self.var_min_depth = self.conf.getint("Output", "var_min_depth")
            self.var_min_freq = self.conf.getfloat("Output", "var_min_freq")

            # Conditional paramaters
            if self.ref_masking:
                self.blastn_opt = self.conf.get("Ref_Masking", "blastn_opt")
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
    # TODO verify if files are valid and numeric value not negative

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

    def run(self):
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
        self.ref_dir = mkdir(path.join(self.outdir, "references/"))
        self.db_dir = mkdir(path.join(self.outdir, "blast_db/"))
        self.fastq_dir = mkdir(path.join(self.outdir, "fastq/"))
        self.align_dir = mkdir(path.join(self.outdir, "bwa_align/"))
        self.index_dir = mkdir(path.join(self.outdir, "bwa_index/"))
        self.result_dir = mkdir(path.join(self.outdir, "results/"))

        print ("\n##### PARSE REFERENCES #####\n")
        # Create CV_Reference.Reference object for each reference easily accessible through
        # Reference class methods
        self._extract_ref()

        # Reference Masking
        if self.ref_masking:
            print ("\n##### REFERENCE HOMOLOGIES MASKING #####\n")
            ref_list = self._iterative_masker()
            # Erase existing index value if ref masking was performed
            bwa_index = None

        # Fastq Filtering
        if self.quality_filtering or self.adapter_trimming:
            print ("\n##### FASTQ FILTERING #####\n")
            self._fastq_filter()

        # BWA alignment
        print ("\n##### READ REFERENCES AND ALIGN WITH BWA #####\n")
        # An index will be generated if no index was provided
        self.sam = Mem.align (
            self.R1, self.R2,
            index = self.bwa_index,
            ref = Reference.allFasta(),
            align_opt = self.bwa_mem_opt,
            index_opt = self.bwa_index_opt,
            aligner = self.bwa_aligner,
            indexer = self.bwa_indexer,
            align_outdir = self.align_dir,
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


        print ("\n##### DONE #####\n")
        print ("Total execution time = {}s".format(round(time()-stime, 2)))

    ##~~~~~~~PRIVATE METHODS~~~~~~~#

    def _extract_ref(self):
        """
        Import and expand fasta references and associated flags in a Reference object
        expand the file if Gziped to avoid multiple compression/decompression during execution
        """

        for ref in self.raw_ref_list:
            Ref = Reference(
                name = ref['name'],
                ref_fasta = expand_file(fp=ref['fasta'], outdir=self.ref_dir, copy_ungz=True),
                bam_maker = Bam.BamMaker(
                    make_bam = 'bam' in ref['output'],
                    make_sam = 'sam' in ref['output']),
                cov_maker = Coverage2.CoverageMaker(
                    min_depth=self.cov_min_depth,
                    make_bedgraph = 'bedgraph' in ref['output'],
                    make_bed = 'bed' in ref['output'],
                    make_covgraph = 'covgraph' in ref['output']),
                var_maker = Variant.VariantMaker(
                    min_depth=self.var_min_depth,
                    min_freq=self.var_min_freq,
                    make_freqvar = 'variant' in ref['output']))
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
                db_opt = self.mkblastdb_opt,
                db_outdir = self.db_dir,
                db_outname = subject.name)

            # Masking hits in suject fasta if hits in hit_list
            subject.fasta_path = RefMasker.mask (
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
            qFilter = QualityFilter (self.min_qual)
        else:
            qFilter = None

        # Define a adapter trimmer object
        if self.adapter_trimming:
            trimmer = AdapterTrimmer(
                Aligner = ssw_wrap.Aligner(
                    match = self.ssw_match,
                    mismatch = self.ssw_mismatch,
                    gap_open = self.ssw_gap_open,
                    gap_extend = self.ssw_gap_extend),
                adapters = self.adapters,
                find_rc = self.find_rc,
                min_read_len = self.min_read_len,
                min_match_len = self.min_match_len,
                min_match_score = self.min_match_score)
        else:
            trimmer = None

        # Filter fastq for quality and adapter with FastqFilterPP
        fFilter = FastqFilterPP (
            self.R1, self.R2,
            quality_filter = qFilter,
            adapter_trimmer = trimmer,
            outdir = self.fastq_dir,
            input_qual = self.input_qual,
            compress_output=False)

        print (repr(fFilter))
        return fFilter.getTrimmed()

    def _sam_spliter (self):
        """
        """
        samfile = pysam.Samfile(self.sam, "r" )
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
            # Finally if all is fine attribute the read to a Reference
            else:
                Reference.addRead(samfile.getrname(read.tid), read)

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
        output = "{}{}_Ref_distribution.csv".format(self.result_dir, self.outprefix)

        with open(output, 'wb') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t', quoting=csv.QUOTE_MINIMAL)

            # Table for all reference
            writer.writerow(["ALL REFERENCES "])
            writer.writerow(["Ref name","length","nread","RPKB"])
            for ref in Reference.getInstances():
                writer.writerow([ref.name, len(ref), ref.nread, float(ref.nread)/len(ref)*1000])
            # Add a line for garbage reads
            nread = sum([seq.nread for seq in self.garbage_read])
            writer.writerow(["Garbage_Reads","NA",nread,"NA"])

            # Table decomposing Sequence per Reference
            for ref in Reference.getInstances():
                writer.writerow([""])
                writer.writerow(["REFERENCE "+ref.name])
                writer.writerow(["Seq name","length","nread","RPKB"])
                for seq in ref.seq_dict.values():
                    writer.writerow([seq.name, len(seq), seq.nread, float(seq.nread)/len(seq)*1000])

            # Add a lines for garbage reads
            writer.writerow([""])
            writer.writerow(["REFERENCE Garbage reads"])
            writer.writerow(["Seq name","length","nread","RPKB"])
            for seq in self.garbage_read:
                writer.writerow([seq.name, "NA", seq.nread, "NA"])

#~~~~~~~TOP LEVEL INSTRUCTIONS~~~~~~~#
if __name__ == '__main__':

    main = Main()
    main.run()
