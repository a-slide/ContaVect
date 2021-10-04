#~~~~~~~GLOBAL IMPORTS~~~~~~~#

# Standard library packages import
from os import remove, path

# Local library packages import
#from SamSpliter import SamSpliter
from .MemWrapper import Aligner
from .IndexWrapper import NewIndex, ExistingIndex
from pyDNA.Utilities import mkdir

#~~~~~~~MAIN METHODS~~~~~~~#

def align  (R1,
            R2='',
            index = '',
            ref = '',
            aligner = "bwa mem",
            align_opt="",
            align_threads = 1,
            align_outdir= "./bwa_align/",
            align_outname= "out.sam",
            indexer = "bwa index",
            index_opt="",
            index_outdir = "./bwa_index/",
            index_outname = "out"):
    """
    Main function of the package allowing to validate an existing index or to create a new one,
    then perform a alignment of single or paired fastq sequences against the index. Finally a sam
    file is returned for further analysis. If an valid existing index was given all index option
    and ref_fasta are not required.
    @param R1 Path to the file containing fastq sequences (can be gzipped)
    @param R2 Facultative path to the file containing paired fastq sequence (can be gzipped)
    @param index Index files basename if available
    @param ref Path of the fasta file containing the reference sequence (can be gzipped)
    This parameter can also be a list of fasta file (gzipped or not) in this case all references
    will be merged into a single fasta reference
    @param aligner Path ot the bwa mem executable. Not required if bwa if added to your path
    @param align_opt Bwa mem command line options as a string
    @param align_outdir Directory where to store the sam file
    @param align_outname Name of the output sam file
    @param indexer Path ot the bwa index executable. Not required if bwa if added to your path
    @param index_opt Bwa index command line options as a string
    @param index_outdir Directory where to store the index files
    @param index_outname Basename of the index file
    @return Path of the output sam file
    """
    # Try to import an existing index
    try:
        if not index:
            raise Exception("No index provided")

        print("Existing index provided")
        idx = ExistingIndex(index)

    # If no index or if an error occured during validation of the existing index = create a new one
    except Exception as E:
        print (E)

        # Verify the presence of the reference fasta file
        if not ref:
            raise Exception("Invalid or no fasta file provided. Cannot create an index")

        print("Generating index...")
        mkdir(index_outdir)
        index_path = path.join(index_outdir, index_outname)
        idx = NewIndex(ref, index_path, index_opt, indexer)

    # Create a Aligner object
    mem = Aligner(idx, align_opt, aligner, align_threads)
    #~print (repr(mem))
    mkdir(align_outdir)

    # Align the reference index with R1 fastq (and R2)
    align_path = path.join(align_outdir, align_outname)
    return (mem.align(R1, R2, align_path))
