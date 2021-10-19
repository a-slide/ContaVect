#~~~~~~~GLOBAL IMPORTS~~~~~~~#
# Standard library packages import
from os import path, remove, rmdir
from multiprocessing import cpu_count
from time import time

# Local library packages
from pyDNA.Utilities import run_command, file_basename, make_cmd_str

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Aligner(object):
    """
    @class  Aligner
    @brief  Perform a alignement of a query fastq file against a bwa index. Results are written in
    a sam file whose path is returned at the end of the alignment
    BWA 0.7.5+ needs to be install and eventually added to the path.
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __repr__(self):
        msg = "BWA MEM WRAPPER\n"
        msg += "Bwa mem path : {}\n".format(self.aligner)
        msg += "Options : {}\n".format(self.align_opt)
        msg += repr(self.Index)
        return msg

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def __init__ (self, Index, align_opt="", aligner = "bwa mem", bwa_threads = 1):
        """
        Initialize the object and index the reference genome if necessary
        @param Index Bwa index object NewIndex or ExistingIndex
        @param align_opt Bwa mem command line options as a string
        @param bwa_mem Path ot the bwa mem executable. Not required if bwa if added to your path
        """
        # Creating object variables
        self.aligner = aligner
        self.Index = Index
        self.bwa_threads = bwa_threads
        # if bwa_threads == 0 use all cores on the node
        if self.bwa_threads == 0:
            self.bwa_threads = cpu_count()
        # By default the option t is setted to the max number of available threads
        self.align_opt = "{} -t {}".format(align_opt, bwa_threads)

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def align(self, R1_path, R2_path="", out_path="./out.sam"):
        """
        Align query fastq against a subject database and return a list of BlastHit object
        @param R1_path Path to the file containing fastq sequences
        @param R2_path Facultative path to the file containing paired fastq sequence
        @param out_path Path to the output sam file
        @return A list of BlastHit objects if at least one hit was found
        @exception (SystemError,OSerror) May be returned by run_command in case of invalid command line.
        """
        # Build the command line
        cmd = "{} {} {} {} {} {} ".format(
            self.aligner,
            self.align_opt,
            self.Index.index_path,
            R1_path,
            R2_path,
            "> "+out_path)

        # Execute bwa mem (Can raise a SystemError) and verify if stdout is not None
        print(("Align against {} index with bwa mem".format(file_basename (self.Index.index_path))))
        print(cmd)
        stderr_list = run_command(cmd, stdin=None, ret_stderr=True, ret_stdout=False).decode().split("\n")
    
        # In bwa stderr return a report of alignment just print the most important one
        print((stderr_list[0]))
        print(("\n".join(stderr_list[-4:])))
        return out_path
