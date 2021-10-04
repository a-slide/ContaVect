#~~~~~~~GLOBAL IMPORTS~~~~~~~#
# Standard library packages import
import gzip
from os import close, remove, path
from multiprocessing import cpu_count
from time import time
from tempfile import mkstemp

# Local library packages
from pyDNA.Utilities import run_command, file_basename, fgunzip
from .BlastHit import BlastHit

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Aligner(object):
    """
    @class  Aligner
    @brief  Perform de blastn of a DNA query against a blast database. If hits are found, a list of
    BlastHit objects is returned.
    Blast+ 2.8+ needs to be install and eventually added to the path.
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __repr__(self):
        msg = "BLASTN WRAPPER\n"
        msg += "Blastn path : {}\n".format(self.blastn)
        msg += "Options : {}\n".format(self.blastn_opt)
        msg += repr(self.Blastdb)
        return msg

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def __init__ (self, Blastdb, blastn_opt="", blastn="blastn", num_threads=1):
        """
        Initialize the object and index the reference genome if necessary
        @param Blastdb Blast database object NewDB or ExistingDB
        @param blastn_opt Blastn command line options as a string
        @param blastn Path ot the bwa executable. If bwa if already added to your system path
        do not change the default value
        """
        # Creating object variables
        self.blastn = blastn
        self.Blastdb = Blastdb
        self.num_threads = num_threads

        # init an option dict and attribute defaut options
        # if num_threads == 0 use all cores on the node
        if self.num_threads == 0:
            self.num_threads = cpu_count()
        self.blastn_opt = "{} -num_threads {} -task {} -outfmt {} -dust {} -db {}".format(
            blastn_opt, self.num_threads, "blastn", 6, "no", self.Blastdb.db_path)

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def align (self, query):
        """
        Blast query against a subject database and return a list of BlastHit object
        @param  query Path to a fasta file containing the query sequences
        @param  evalue  Cutoff used in blast to select valid hits
        @return A list of BlastHit objects if at least one hit was found
        @exception (SystemError,OSerror) May be returned by run_command in case of invalid command line
        """
        # Build the command line string
        query_name = file_basename(query)

        # If the fasta file is compressed = extract the file in a temporary file
        if query[-2:].lower() == "gz":
            print ("Extracting the compressed fasta in a temporary file")
            fd, tmp_path = mkstemp()
            fgunzip (in_path=query, out_path=tmp_path)
            blastn_opt = self.blastn_opt + " -query {}".format(tmp_path)
            hits_list = self._align(query_name, blastn_opt)
            close(fd)
            remove(tmp_path)
            return hits_list
        # Else just proceed by using the fasta reference
        else:
            blastn_opt = self.blastn_opt + " -query {}".format(query)
            return self._align(query_name, blastn_opt)

    def _align (self, query_name, blastn_opt):

        print(("Blast {} against {} database with blastn".format(query_name, file_basename (self.Blastdb.db_path))), end=' ')

        # Build the command line string
        cmd = "{} {}".format(self.blastn, blastn_opt)

        # Run the command line without stdin and asking only stdout
        blast_lines = run_command(cmd, stdin=None, ret_stderr=False, ret_stdout=True).splitlines()

        for line in blast_lines:
            # Parse each result lines and create a BlastHit object
            h = line.split()
            BlastHit(h[0], h[1] , h[2], h[3], h[4], h[5], h[6], h[7], h[8], h[9], h[10], h[11])

        # Sumarize the hit count in the different references
        print(("\t{} hits found".format(BlastHit.count_total())))

        # Get the list of hits from BlastHit class and reset the class list.
        hits_list = BlastHit.get()
        BlastHit.reset_list()
        return hits_list
