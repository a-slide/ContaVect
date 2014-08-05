#~~~~~~~GLOBAL IMPORTS~~~~~~~#

# Standard library packages import
from os import path
from tempfile import mkdtemp
from shutil import rmtree
from time import time

# Local Package import
from Utilities import mkdir, file_basename, file_name, expand_filelist
from Blast import Blastn
import RefMasker
from FastqFT.FastqFilterPP import FastqFilterPP
from FastqFT.QualityFilter import QualityFilter
from FastqFT.AdapterTrimmer import AdapterTrimmer
from Ssw import ssw_wrap
from FastaReader2 import FastaReader2 as FastaReader
from Bwa import Mem

#~~~~~~~MAIN FUNCTION~~~~~~~#

def maindef(ref_list, R1, R2, outdir="./", outprefix="out"):

    stime = time()

    # Reference Masking
    ref_masking = True
    blastn_opt=None
    mkblastdb_opt=None

    # Fastq Filtering
    quality_filtering = True
    min_qual = 30
    input_qual = "fastq-sanger"
    adapter_trimming = True
    adapter_path = "/media/analyse/Pharmaco_AAV/documents_pharmaco_AAV/Programming/Python/ContaVect/test/adapter.fa"
    # ... Other parameters not available for now

    # Alignment
    bwa_index = None
    mem_opt = None
    index_opt = None

    # Create a main directory for file created during the programm
    main_dir = path.abspath(outdir)+"/"
    mkdir(main_dir)

    print ("\n##### EXPAND AND REGROUP REFERENCES #####\n")
    # Expand the reference sequence to avoid multiple decompression during Program execution
    ref_dir = path.join(main_dir+"references/")
    mkdir(ref_dir)
    ref_list = expand_filelist (ref_list, outdir=ref_dir, copy_ungz=True)
    #~print "\nReference list:\n" + "\n".join([ref for ref in ref_list]) + "\n"

    # Reference Masking
    if ref_masking:
        print ("\n##### REFERENCE HOMOLOGIES MASKING #####\n")
        db_dir = path.join(main_dir, "blast_db/")
        mkdir (db_dir)
        ref_list = iterative_masker(ref_list, db_dir, ref_dir, blastn_opt, mkblastdb_opt)

    print "\nReference list:\n" + "\n".join([ref for ref in ref_list]) + "\n"

    # Fastq Filtering
    if quality_filtering or adapter_trimming:
        print ("\n##### FASTQ FILTERING #####\n")
        fastq_dir = path.join(main_dir+"fastq/")
        mkdir(fastq_dir)
        R1, R2 = fastq_filter(R1, R2, fastq_dir, quality_filtering, min_qual, input_qual, adapter_trimming, adapter_path, maxqueue)

    print "R1 = " + R1
    print "R2 = " + R2

    # BWA alignment
    print ("\n##### READ REFERENCES AND ALIGN WITH BWA #####\n")

    # Create a folder for aligned read storage
    align_dir = path.join(main_dir+"bwa_align/")
    mkdir(align_dir)

    if ref_masking or not bwa_index:
        # Read fasta ref to map seq names and merge then together in a single fasta file
        reader = FastaReader(ref_list, write_merge=True, output=ref_dir+"merged.fa", write_summary=True)

        # Create a folder for index storage
        index_dir = path.join(main_dir+"bwa_index/")
        mkdir(index_dir)

        # Index and Align with Bwa Mem
        Mem.align (R1,R2, ref_fasta=reader.getMergeRef(),
                   align_opt=mem_opt, align_outdir=align_dir, align_outname=outprefix+".sam",
                   index_opt=index_opt, index_outdir=index_dir, index_outname = "merged.fa")

    else:
        # Read fasta ref to map seq names
        reader = FastaReader(ref_list, write_summary=True)

        # Align with Bwa Mem against
        Mem.align (R1,R2, ref_index=bwa_index,
                   align_opt=mem_opt, align_outdir=align_dir, align_outname=outprefix+".sam")

    print ("\n##### DONE #####\n")
    print ("Total execution time = {}s".format(round(time()-stime, 2)))

#~~~~~~~HELPER FUNCTIONS~~~~~~~#

def iterative_masker (
    ref_list,
    db_dir="./",
    ref_dir="./",
    align_opt=None,
    db_opt=None):
    new_ref_list = []

    while (len(ref_list)>=2):
        # Extract subject and query_list from ref_list
        subject = ref_list [-1]
        query_list = ref_list [0:-1]
        subject_name = file_basename(subject)
        print ("\n# PROCESSING REFERENCE {} #\n".format(subject_name ))

        # Perform a blast of query list against subject
        hit_list = Blastn.align (
            query_list,
            subject_fasta=subject,
            evalue=0.5,
            align_opt=None,
            db_opt=None,
            db_outdir = db_dir,
            db_outname = subject_name)

        # Masking hits in suject fasta
        new_ref_list.append (RefMasker.mask (
            subject,
            hit_list,
            ref_outdir=ref_dir,
            ref_outname="masked_{}.fa".format(subject_name),
            compress_ouput=False ))

        # Remove the last last element of the list
        ref_list.pop()

    # Adding the last reference that was not poped out of the list
    new_ref_list.extend(ref_list)
    # Reverse the list to be in the original order
    new_ref_list.reverse()
    return new_ref_list


def fastq_filter (
    R1, R2,
    outdir="./",
    quality_filtering=False,
    min_qual=0,
    input_qual="fastq-sanger",
    fastq_trimming=False,
    adapter_path=""):

    # Define a quality filter object
    if quality_filtering:
        qFilter = QualityFilter (min_qual)
    else:
        qFilter = None

    # Define a adapter trimmer object
    if fastq_trimming:
        sswAligner = ssw_wrap.Aligner() # Instanciate with a default score matrix
        trimmer = AdapterTrimmer(sswAligner, adapter_path) # Instanciate with a default score filters
    else:
        trimmer = None

    # Filter fastq for quality and adapter with FastqFilterPP
    fFilter = FastqFilterPP (R1, R2,
        quality_filter=qFilter,
        adapter_trimmer=trimmer,
        outdir=outdir,
        input_qual=input_qual,
        maxqueue=maxqueue,
        compress_output=False)
    
    print (repr(fFilter))
    return (fFilter.R1_out, fFilter.R2_out)
