from os import path
from tempfile import mkdtemp
from shutil import rmtree

from Utilities import mkdir, file_basename, expand_filelist
from Blast import Blastn
import RefMasker

from FastqFT.FastqFilterPP import FastqFilterPP
from FastqFT.QualityFilter import QualityFilter
from FastqFT.AdapterTrimmer import AdapterTrimmer
from Ssw import ssw_wrap

#~~~~~~~MAIN FUNCTION~~~~~~~#

def maindef(ref_list, R1, R2, outdir="./"):

    # Reference Masking
    ref_masking = True
    align_opt=None
    db_opt=None

    # Fastq Filtering
    quality_filtering = True
    min_qual = 30
    input_qual = "fastq-sanger"
    fastq_trimming = False
    adapter_path = "/home/adrien/Programming/Python/DNA-tools/test/adapter.fa"
    # ... Other parameters not available for now

    # Alignment
    bwa_index=None


    # Create a main directory for file created during the programm
    main_dir = path.abspath(outdir)+"/"
    mkdir(main_dir)

    print ("\n##### EXPAND AND REGROUP REFERENCES #####\n")
    # Expand the reference sequence to avoid multiple decompression during Program execution
    working_dir = path.join(main_dir+"references/")
    mkdir(working_dir)
    ref_list = expand_filelist (ref_list, outdir=working_dir, copy_ungz=True)
    #~print "\nReference list:\n" + "\n".join([ref for ref in ref_list]) + "\n"

    # Reference Masking
    if ref_masking:
        print ("\n##### REFERENCE HOMOLOGIES MASKING #####\n")
        #~working_dir = path.join(main_dir+"references")
        #~mkdir(mask_dir)
        working_dir = main_dir
        ref_list = iterative_masker(ref_list, working_dir, align_opt, db_opt)

    print "\nReference list:\n" + "\n".join([ref for ref in ref_list]) + "\n"

    # Fastq Filtering
    if quality_filtering or fastq_trimming:
        print ("\n##### FASTQ FILTERING #####\n")
        working_dir = path.join(main_dir+"filtered_fastq/")
        mkdir(working_dir)
        R1, R2 = fastq_filter(R1, R2, working_dir, quality_filtering, min_qual, input_qual, fastq_trimming, adapter_path)

    print "R1 = " + R1
    print "R2 = " + R2

    # BWA alignment
    print ("\n##### PREPARING FILES FOR ALIGNMENT WITH BWA #####\n")
    if ref_masking or not bwa_index:
        pass

    else:
        pass
        # Generate the list of ref names

    #print ("\n##### ALIGNMENT WITH BWA #####\n")
    #sam = bwa...

    #sam.split

#~~~~~~~HELPER FUNCTIONS~~~~~~~#

def iterative_masker (
    ref_list,
    main_dir="./",
    align_opt=None,
    db_opt=None):
    new_ref_list = []

    while (len(ref_list)>=2):
        # Extract subject and query_list from ref_list
        subject = ref_list [-1]
        query_list = ref_list [0:-1]
        subject_name = file_basename(subject)
        print ("\n# PROCESSING REFERENCE {} #\n".format(subject_name ))

        # Create folder for secure file writting
        db_outdir = path.join(main_dir, "blast_db/")
        ref_outdir = path.join(main_dir, "references/")
        mkdir (db_outdir)
        mkdir (ref_outdir)

        # Perform a blast of query list against subject
        hit_list = Blastn.align (
            query_list,
            subject_fasta=subject,
            evalue=0.5,
            align_opt=None,
            db_opt=None,
            db_outdir = db_outdir,
            db_outname = subject_name)

        # Masking hits in suject fasta
        new_ref_list.append (RefMasker.mask (
            subject,
            hit_list,
            ref_outdir=ref_outdir,
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
        trimmer = AdapterTrimmer(sswaligner, adapter_path) # Instanciate with a default score filters
    else:
        trimmer = None

    # Filter fastq for quality and adapter with FastqFilterPP
    fFilter = FastqFilterPP (R1, R2, qFilter, trimmer, outdir, input_qual)
    print (repr(fFilter))
    return (fFilter.R1_out, fFilter.R2_out)
