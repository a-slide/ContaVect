#~~~~~~~GLOBAL IMPORTS~~~~~~~#
# Standard library packages import
from os import path, remove
from time import time

# Local library packages
from pyDNA.Utilities import run_command, file_basename, make_cmd_str, merge_files

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class NewIndex(object):
    """
    @class NewIndex
    @brief Wrapper for bwa index. Create a reference index from a fasta file
    BWA 0.7.5+ needs to be install and eventually added to the path.
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __repr__(self):
        msg = "BWA INDEX WRAPPER (NEW INDEX)\n"
        msg += "bwa index path : {}\n".format(self.indexer)
        msg += "Blastn database path : {}\n".format(self.index_path)
        msg += "Options : {}\n".format(self.index_opt)
        return msg

    def __str__(self):
        return "\n<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def __init__ (self, ref, index_path="./out.idx", index_opt="", bwa_index = "bwa index"):
        """
        Initialize the object and index the reference genome if necessary
        @param ref Path of the fasta file containing the reference sequence (can be gzipped)
        This parameter can also be a list of fasta file (gzipped or not) in this case all references
        will be merged into a single fasta reference
        @param index_path Outname for the bwa index files basename.
        @param index_opt bwa index dictionnary of option arguments such as "-t 5". The option flag
        have to be the key (without "-") and the the option value in the dictionnary value. If no
        value is requested after the option flag "None" had to be asigned to the value field.
        @param bwa index Path ot the bwa index executable. Not required if bwa if added to your path
        """
        # Creating object variables
        self.indexer = bwa_index
        self.index_path = index_path
        self.index_opt = "{} -p {}".format (index_opt, self.index_path)

        try:
            # If only one ref = use it directly to make an index
            if isinstance(ref, str):
                self.ref = ref
                self._make_index()

            # If list at one element = same thing
            elif isinstance(ref, list) and len(ref) == 1:
                self.ref = ref[0]
                self._make_index()

            # If severel references, merged them, make index and remove the merged reference file
            elif isinstance(ref, list):
                print("Merge references files for indexation")
                self.ref = merge_files(ref, outpath="./out.fa", compress_output=False)
                self._make_index()
                #remove ("./out.fa")
            else:
                raise TypeError("ref variable is neither a list nor a string")

        except Exception as E:
            self._remove_index_files()
            raise Exception (E.message+"Impossible to generate a valid index from the reference sequence")

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _make_index(self):
        """
        Create a bwa index from ref using bwa index
        """
        # Build the command line
        cmd = "{} {} {}".format(self.indexer, self.index_opt, self.ref)

        print("Creating a BWA index")

        # Run the command line without stdin and asking both stdout and stderr
        start_time = time()
        stderr = run_command(cmd, stdin=None, ret_stderr=True, ret_stdout=False)

        print (stderr)
        print(("Index created in {}s\n".format(round(time()-start_time, 3))))

    def _remove_index_files(self):
        """
        Remove index files in case of exception during indexing
        """
        print("Remove index files")

        # Removing Index file
        for ext in ["amb", "ann", "bwt", "pac", "sa"]:
            f = "{}.{}".format(self.index_path, ext)
            if path.isfile (f):
                remove (f)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class ExistingIndex(object):
    """
    @class ExistingIndex
    @brief Import an existing bwa index + verify the existence of files
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    def __repr__(self):
        msg = "BWA INDEX WRAPPER (EXISTING INDEX)\n"
        msg += "Bwa Index Path : {}\n".format(self.index_path)
        return msg

    def __str__(self):
        return "\n<Instance of {} from {}>\n".format(self.__class__.__name__, self.__module__)

    def __init__ (self, index_path):
        """
        @param index_path The index path is the name of any of the index files up to but not
        including the final "amb", "ann", "bwt", "pac" and "sa"
        """
        # Creating object variables
        self.index_path = index_path

        print ("Checking index files")
        # Checking if all index files needed by bwa are
        for ext in ["amb", "ann", "bwt", "pac", "sa"]:
            f = "{}.{}".format(self.index_path, ext)

            if not path.isfile (f):
                raise Exception ("Invalid database : {} does not exist".format(f))
            if path.getsize(f) == 0:
                raise Exception ("Invalid database : {} is empty".format(f))

        print ("All index files are valid")
