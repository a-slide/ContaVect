#!/usr/bin/env python
# -*- coding: utf-8 -*-

#~~~~~~~GLOBAL IMPORTS~~~~~~~#

# Standard library packages import
from os import path
#from tempfile import mkdtemp
#from shutil import rmtree
from time import time, sleep
import ConfigParser
from sys import argv

# Third party packages
import pysam

# Local Package import
from Utilities import mkdir, file_basename, file_name, expand_file
from Blast import Blastn
import RefMasker
from FastqFT.FastqFilterPP import FastqFilterPP
from FastqFT.QualityFilter import QualityFilter
from FastqFT.AdapterTrimmer import AdapterTrimmer
from Ssw import ssw_wrap
from FastaReader2 import FastaReader2 as FastaReader
from Bwa import Mem
from CV_Reference import Reference

#~~~~~~~MAIN FUNCTION~~~~~~~#

def main(conf):
    """
    """
    stime = time()

    # Create a main directory for file created during the programm
    if conf.get("General", "outdir"):
        main_dir = path.abspath(conf.get("General", "outdir"))+"/"
        mkdir(main_dir)
    else:
        main_dir = path.abspath ("./")+"/"

    # Get boolean value for pipeline orientation
    ref_masking = conf.getboolean("Ref_Masking", "ref_masking")
    quality_filtering = conf.getboolean("Fastq_Filtering", "quality_filtering")
    adapter_trimming = conf.getboolean("Fastq_Filtering", "adapter_trimming")
    bwa_index = conf.get("Bwa Alignment", "index")

    print ("\n##### EXPAND AND PARSE REFERENCES #####\n")
    # Expand the reference sequence to avoid multiple decompression during Program execution

    ref_dir = path.join(main_dir+"references/")
    mkdir(ref_dir)
    extract_ref(conf, ref_dir) # Reference object are instancied and accessible through the class methods

    # Reference Masking
    if ref_masking:
        print ("\n##### REFERENCE HOMOLOGIES MASKING #####\n")
        db_dir = path.join(main_dir, "blast_db/")
        mkdir (db_dir)
        ref_list = iterative_masker(conf, db_dir, ref_dir)
        # Erase existing index value if ref masking was performed
        bwa_index = None

    # Fastq Filtering
    if quality_filtering or adapter_trimming:
        print ("\n##### FASTQ FILTERING #####\n")
        fastq_dir = path.join(main_dir+"fastq/")
        mkdir(fastq_dir)
        R1, R2 = fastq_filter(conf, fastq_dir, quality_filtering, adapter_trimming)
    else:
        R1, R2 = conf.get("Fastq", "R1"), conf.get("Fastq", "R2")

    # BWA alignment
    print ("\n##### READ REFERENCES AND ALIGN WITH BWA #####\n")
    # Create a folder for aligned read and index file storage
    align_dir = path.join(main_dir+"bwa_align/")
    index_dir = path.join(main_dir+"bwa_index/")
    mkdir(align_dir)
    mkdir(index_dir)

    # An index will be generated if no index was provided
    sam = Mem.align (
        R1, R2,
        index = bwa_index,
        ref = [ref.fasta_path for ref in Reference.Instances],
        align_opt = conf.get("Bwa Alignment", "mem_opt"),
        index_opt = conf.get("Bwa Alignment", "index_opt"),
        aligner = conf.get("Bwa Alignment", "bwa_mem"),
        indexer = conf.get("Bwa Alignment", "bwa_index"),
        align_outdir= align_dir,
        index_outdir = index_dir,
        align_outname = conf.get("General", "outprefix")+".sam",
        index_outname = conf.get("General", "outprefix")+".idx")

    # Split the output sam file according to each reference
    sam_spliter (conf, sam)

    # Ask references to generate the output they were configured to
    result_dir = path.join(main_dir+"results/")
    mkdir(result_dir)
    for ref in Reference.Instances:
        ref.mk_output(outpath=result_dir+conf.get("General", "outprefix"))

    print ("\n##### DONE #####\n")
    print ("Total execution time = {}s".format(round(time()-stime, 2)))

#~~~~~~~HELPER FUNCTIONS~~~~~~~#

def extract_ref(conf, refdir="./"):
    """
    Import and expand fasta references and associated flags in a Reference object
    """
    i = 1
    while True:
        # Iterate over Ref in the config file until no more reference is found
        ref_id = "Ref"+str(i)

        if not conf.has_section(ref_id):
            break
        try:
            fasta = expand_file(conf.get(ref_id, "fasta"), outdir=refdir, copy_ungz=True)
            sam = conf.getboolean (ref_id, "sam")
            covgraph = conf.getboolean (ref_id, "covgraph")
            bedgraph = conf.getboolean (ref_id, "bedgraph")
            vcf = conf.getboolean (ref_id, "pileup")
        except ValueError:
            sleep (2)
            print ("ERROR INVALID VALUES, SKIPPING THE REFERENCE NUMBER "+str(i))
        else:
            a = Reference(fasta, sam, covgraph, bedgraph, vcf)
            print (repr(a))
            i+=1

def iterative_masker (conf, db_dir="./", ref_dir="./"):
    """
    """
    # Iterate over index in Reference.instances staring by the last one until the 2nd one
    for i in range(Reference.countInstances()-1, 0, -1):
        # Extract subject and query_list from ref_list
        subject = Reference.Instances[i]
        query_list = Reference.Instances[0:i]
        print ("\n# PROCESSING REFERENCE {} #\n".format(subject.name))

        # Perform a blast of query list against subject
        hit_list = Blastn.align (
            [ref.fasta_path for ref in query_list],
            subject_fasta=subject.fasta_path,
            align_opt=conf.get("Ref_Masking", "blastn_opt"),
            db_opt=conf.get("Ref_Masking", "mkblastdb_opt"),
            db_outdir = db_dir,
            db_outname = subject.name)

        # Masking hits in suject fasta if hits in hit_list
        subject.fasta_path = RefMasker.mask (
            subject.fasta_path,
            hit_list,
            ref_outdir=ref_dir,
            ref_outname="masked_{}.fa".format(subject.name),
            compress_ouput=False )

def fastq_filter (conf, outdir="./", quality_filtering=False, adapter_trimming=False):
    """
    """
    # Define a quality filter object
    qFilter = None
    if quality_filtering:
        qFilter = QualityFilter (conf.getfloat("Fastq_Filtering", "min_qual"))

    # Define a adapter trimmer object
    trimmer = None
    if adapter_trimming:
        sswAligner = ssw_wrap.Aligner(
            match=conf.getint("Fastq_Filtering", "ssw_match"),
            mismatch=conf.getint("Fastq_Filtering", "ssw_mismatch"),
            gap_open=conf.getint("Fastq_Filtering", "ssw_gap_open"),
            gap_extend=conf.getint("Fastq_Filtering", "ssw_gap_extend"))
        trimmer = AdapterTrimmer(
            aligner=sswAligner,
            adapter_path=conf.get("Fastq_Filtering", "adapter_path"),
            min_read_len=conf.getfloat("Fastq_Filtering", "min_read_len"),
            min_match_len=conf.getfloat("Fastq_Filtering", "min_match_len"),
            min_match_score=conf.getfloat("Fastq_Filtering", "min_match_score"))

    # Filter fastq for quality and adapter with FastqFilterPP
    fFilter = FastqFilterPP (
        R1=conf.get("Fastq", "R1"),
        R2=conf.get("Fastq", "R2"),
        quality_filter=qFilter,
        adapter_trimmer=trimmer,
        outdir=outdir,
        input_qual=conf.get("Fastq_Filtering", "input_qual"),
        compress_output=False)

    print (repr(fFilter))
    return fFilter.getTrimmed()


def sam_spliter (conf, sam):
    """
    """
    samfile = pysam.Samfile( sam, "r" )

    #TODO define an unmapped reference
    unmapped = []

    for read in samfile:
        if not read.is_secondary:
            if read.mapq < 30:
                unmapped.append(read)
            else:
                try:
                    Reference.addRead(samfile.getrname(read.tid), read)
                except ValueError:
                    unmapped.append(read)

    print ("Number of unmapped reads : {}\n".format(len (unmapped)))


#~~~~~~~TOP LEVEL INSTRUCTIONS~~~~~~~#
if __name__ == '__main__':

    conf = ConfigParser.RawConfigParser(allow_no_value=True)
    conf.read(argv[1])
    main(conf)

