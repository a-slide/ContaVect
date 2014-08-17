"""
@package    Utilities
@brief      Contains several usefull functions to interact with OS environement and to parse files
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

#~~~~~~~COMMAND LINE UTILITIES~~~~~~~#

def run_command(cmd, stdin=None, ret_stderr=False, ret_stdout=True):
    """
    Run a command line in the default shell and return the standard output
    @param  cmd A command line string formated as a string
    @param  stdinput    Facultative parameters to redirect an object to the standard input
    @param  ret_stderr  If True the standard error output will be returned
    @param  ret_stdout  If True the standard output will be returned
    @note If ret_stderr and ret_stdout are True a tuple will be returned and if both are False
    None will be returned
    @return If no standard error return the standard output as a string
    @exception  OSError Raise if a message is return on the standard error output
    @exception  (ValueError,OSError) May be raise by Popen
    """
    # Function specific imports
    from subprocess import Popen, PIPE

    # Execute the command line in the default shell
    if stdin:
        proc = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        stdout, stderr = proc.communicate(input=stdin)
    else:
        proc = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        stdout, stderr = proc.communicate()

    if proc.returncode == 1:
        msg = "An error occured during execution of following command :\n"
        msg += "COMMAND : {}\n".format(cmd)
        msg += "STDERR : {}\n".format(stderr)
        raise Exception (msg)

    # Else return data according to user choices is returned
    if ret_stdout and ret_stderr:
        return stdout, stderr
    elif ret_stdout:
        return stdout
    elif ret_stderr:
        return stderr
    else:
        return None

def make_cmd_str(prog_name, opt_dict={}, opt_list=[]):
    """
    Create a Unix like command line string from a
    @param prog_name Name (if added to the system path) or path of the programm
    @param opt_dict Dictionnary of option arguments such as "-t 5". The option flag have to
    be the key (without "-") and the the option value in the dictionnary value. If no value is
    requested after the option flag "None" had to be asigned to the value field.
    @param opt_list List of simple command line arguments
    @exemple make_cmd_str("bwa", {"b":None, t":6, "i":"../idx/seq.fa"}, ["../read1", "../read2"])
    """

    # Start the string by the name of the program
    cmd = "{} ".format(prog_name)

    # Add options arguments from opt_dict
    if opt_dict:
        for key, value in opt_dict.items():
            if value:
                cmd += "-{} {} ".format(key, value)
            else:
                cmd += "-{} ".format(key)

    # Add arguments from opt_list
    if opt_list:
        for value in opt_list:
            cmd += "{} ".format(value)

    return cmd

#~~~~~~~FILE MANIPULATION~~~~~~~#

def copyFile(src, dest):
    """
    Copy a single file to a destination file or folder (with error handling/reporting)
    @param src Source file path
    @param dest Path of the folder where to copy the source file
    """
    import shutil

    try:
        shutil.copy(src, dest)
    # eg. src and dest are the same file
    except shutil.Error as e:
        print('Error: %s' % e)
    # eg. source or destination doesn't exist
    except IOError as e:
        print('Error: %s' % e.strerror)


def fgzip(in_path, out_path=None):
    """
    @param in_path Path of the input uncompressed file
    @param out_path Path of the output compressed file (facultative)
    @exception  OSError Can be raise by open
    """
    # Function specific imports
    import gzip
    from os import remove, path

    # Generate a automatic name if none is given
    if not out_path:
        out_path = in_path +".gz"

    # Try to initialize handle for
    try:
        in_handle = open(in_path, "rb")
        out_handle = gzip.open(out_path, "wb")
        # Write input file in output file
        print ("Compressing {}".format(in_path))
        out_handle.write (in_handle.read())
        # Close both files
        in_handle.close()
        out_handle.close()
        return path.abspath(out_path)

    except IOError as E:
        print(E)
        if path.isfile (out_path):
            try:
                remove (out_path)
            except OSError:
                print "Can't remove {}".format(out_path)

def fgunzip(in_path, out_path=None):
    """
    @param in_path Path of the input compressed file
    @param out_path Path of the output uncompressed file (facultative)
    @exception  OSError Can be raise by open
    """
    # Function specific imports
    import gzip
    from os import remove, path

    # Generate a automatic name without .gz extension if none is given
    if not out_path:
        out_path = in_path[0:-3]

    try:
        # Try to initialize handle for
        in_handle = gzip.GzipFile(in_path, 'rb')
        out_handle = open(out_path, "wb")
        # Write input file in output file
        print ("Uncompressing {}".format(in_path))
        out_handle.write (in_handle.read())
        # Close both files
        out_handle.close()
        in_handle.close()
        return path.abspath(out_path)

    except IOError as E:
        print(E)
        if path.isfile (out_path):
            try:
                remove (out_path)
            except OSError:
                print "Can't remove {}".format(out_path)

def expand_file (fp, outdir="./", copy_ungz=False):
    """
    Iterate over a list of files and expand files in outdir if the files are gzipped
    Else the file won't be modified and won't be moved from it's current directory
    @param fp Path to a file eventually gzipped
    @param outdir Path of the directory in which to uncompress or copy the files
    @param copy_ungz Copy uncompressed file in the outdir without modification
    @return A path to an uncompressed file
    """
    # Function specific imports
    from os import path, link
    from shutil import copy

    assert path.isfile(fp), "{} is not a valid file".format(fp)

    # Extract if gziped
    if fp[-2:].lower() == "gz":
        out_path = path.join (outdir, file_name(fp)[:-3])
        fgunzip (fp, out_path)

    # Copy to the outdir if requested
    elif copy_ungz :
        out_path = path.join (outdir, file_name(fp))
        # try to create a hard link else just recopy
        try:
            link(fp, out_path)
        except Exception:
            print ('link failed')
            copyFile(fp, outdir)
    else:
        out_path = file_name(fp)

    return path.abspath(out_path)

def expand_filelist (file_list, outdir="./", copy_ungz=False):
    """
    Iterate over a list of files and expand files in outdir if the files are gzipped
    Else the file won't be modified and won't be moved from it's current directory
    @param file_list List of path to files eventually gzipped
    @param outdir Path of the directory in which to uncompress the gziped files
    @param copy_ungz Copy uncompressed file in the outdir without modification
    @return A list of path to uncompressed files including the name of original files that do not
    needed to be uncompressed.
    """
    # Function specific imports
    from os import path, link
    from shutil import copy

    new_file_list = []

    for fp in file_list:
        new_file_list.append(expand_file(fp, outdir, copy_ungz))
    return new_file_list

def mkdir(fp):
    """
    Create a directory at the indicated path\n
    Reproduce the ability of UNIX "mkdir -p" command
    (ie if the path already exits no exception will be raised).
    @param  fp path name where the folder should be created
    @exception  OSError Can be raise by os.mkdir
    """
    # Function specific imports
    from os import mkdir, path

    if path.exists(fp) and path.isdir(fp):
        #print ("'{}' already exist in the current directory".format(fp))
        return fp
    else:
        #print ("Creating '{}' in the current directory".format(fp))
        mkdir(fp)
        return fp

def merge_files (inpath_list, outpath="out", compress_output=True, bufsize = 100000):
    """
    Merge a list of text file (gzip or not) in a single file taht can be compress or not
    @param input_list List of files to merge
    @param outpath Destination file
    @param compress_output Gzip the output file. Slower if true
    @param bufline Size of the output file write buffer in line (positive integer)
    @return path of the output merged file
    """
    # Standard library import
    import gzip
    from sys import stdout
    from os import path
    from time import time

    stime = time()
    # Creating and storing a file for writting output
    outpath = path.abspath(outpath)+".gz" if compress_output else path.abspath(outpath)
    openout = gzip.open if compress_output else open

    with openout(outpath, "wb") as out_handle:
        # Iterate over files in the input list
        for inpath in inpath_list:

            # Open according to the compression
            openin = gzip.open if file_extension(inpath) == "gz" else open
            with openin (inpath, "rb") as in_handle:
                stdout.write("\t+ {}  ".format(file_name(inpath)))
                stdout.flush()

                # Init a line counter and a text buffer
                lineno = 0
                linebuf = ""

                # Store line in the buffer until the line size is full then flush in out_handle
                for line in in_handle:
                    lineno += 1
                    linebuf += line
                    if lineno % bufsize == 0:
                        out_handle.write(linebuf)
                        linebuf = ""
                    if lineno % 1000000 == 0:
                        stdout.write("*")
                        stdout.flush()

                # Flush the remaining lines in the buffer
                stdout.write("*\n")
                stdout.flush()
                out_handle.write(linebuf)

    print ("{} files merged in {}s\n".format (len(inpath_list), round(time()-stime,3)))
    return outpath

def file_basename (path):
    """
    @param path Filepath as a string
    @return The basename of a file without folder location and extension
    """
    return path.rpartition('/')[2].partition('.')[0]

def file_extension (path):
    """
    @param path Filepath as a string
    @return The extension of a file in lowercase
    """
    return path.rpartition(".")[2].lower()

def file_name (path):
    """
    @param path Filepath as a string
    @return The complete name of a file with the extension but without folder location
    """
    return path.rpartition("/")[2]

def dir_name (path):
    """
    @param path Filepath as a string
    @return The complete path where is located the file without the file name
    """
    return path.rpartition("/")[0].rpartition("/")[2]

def rm_blank (name, replace=""):
    """
    @param name String with blanck spaces
    @param replace character of replacement for blanks (Default None)
    @return String with blanks removed and replace. Blanks at extremities are always removed
    and nor replaced
    """
    return replace.join(name.split())

#~~~~~~~SEQUENCE UTILITIES~~~~~~~#

def import_seq(filename, col_type="dict", seq_type="fasta"):
    """
    Import sequences from a fasta files in a list of biopython SeqRecord
    @param filename Valid path to a fasta file. Can contains several sequences and can be gzipped
    @param col_type Type of the collection where SeqReccord entries will be added ("list" or "dict").
    @param seq_type Type of the sequence file to parse (see Biopython seqIO for supported format)
    @return A list or a dictionnary containing all seqReccord objects from the fastq file
    @exception IOError  Raise if the path in invalid or unreadeable
    """
    # Require the Third party package Biopython
    from Bio import SeqIO
    import gzip

    # Try to open the file first gz compressed and uncompressed
    try:

        # Verify if the type of the input sequence is valid
        seq_type = seq_type.lower()
        allowed_seq = ["fasta", "genbank", "gb", "fastq-illumina", "fastq-solexa" , "fastq",
        "fastq-sanger", "embl ", "abi ", "seqxml", "sff", "uniprot-xml"]
        assert seq_type in allowed_seq , "The input file format have to be in the following list : "+ ", ".join(allowed_seq)

        # Verify if the type of the output collection is valid
        col_type = col_type.lower()
        allowed_types = ["dict", "list"]
        assert col_type  in allowed_types, "The output collection type have to be in the following list : "+ ", ".join(allowed_types)

        # Open gzipped or uncompressed file
        if file_extension(filename) == "gz":
            #print("\tUncompressing and extracting data")
            handle = gzip.open(filename, "r")
        else:
            #print("\tExtracting data")
            handle = open(filename, "r")

        # Create the collection
        if col_type == "list":
            seq_col = list(SeqIO.parse(handle, seq_type))
        else:
            seq_col = SeqIO.to_dict(SeqIO.parse(handle, seq_type))

        # Close file, verify if the collection is filled and returned it
        handle.close()
        assert seq_col, 'The collection contains no SeqRecord after file parsing. Exit'
        return seq_col

    except IOError as E:
        print(E)
        exit()

    except AssertionError as E:
        print (E)
        exit()

def count_seq (filename, seq_type="fasta"):
    """
    Count the number of sequences in a fastq or a fastq file
    @param filename Path to a valid readeable file
    @param file_type Should be either fastq or fastq. Default fasta
    """
    #Standard library import
    import gzip

    # Verify if the file is fasta or fastq type
    assert seq_type in ["fasta", "fastq"], "The file has to be either fastq or fasta format"

    # Open the file
    if file_extension(filename) == "gz":
        fp = gzip.open(filename, "rb")
    else:
        fp = open(filename, "rb")

    # line no counter
    nline = 0

    # FASTA Find a start line seq character ">" an increment the counter each time
    if seq_type ==  "fasta":
        for line in fp:
            if line[0] == ">":
                nline+=1
        fp.close()
        return nline

    # FASTQ No motif to find, but 4 lines correspond to 1 sequence
    else:
        for line in fp:
            nline+=1
        fp.close()
        return nline/4

def DNA_reverse_comp (sequence, AmbiguousBase=True):
    """
    Generate the reverese complementary sequence of a given DNA sequence
    @param
    """
    if AmbiguousBase:
        compl = {'A':'T','T':'A','G':'C','C':'G','Y':'R','R':'Y','S':'S','W':'W','K':'M','M':'K',
               'B':'V','V':'B','D':'H','H':'D','N':'N'}
    else:
        compl = {'A':'T','T':'A','G':'C','C':'G','N':'N'}

    compl_sequence=""
    for base in sequence:
        try:
            compl_sequence += compl[base.upper()]
        except KeyError:
            compl_sequence += 'N'

    return compl_sequence[::-1]


#~~~~~~~GRAPHICAL UTILIES~~~~~~~#

def fill_between_graph (X, Y, basename="out", img_type="png", title=None, xlabel=None, ylabel=None,
    baseline=0, xsize=15, ysize=10, dpi=100, fill_color='green'):
    """
    Trace a generic fill between graph with matplotlib pyplot
    @param X List of values for x axis
    @param Y List of values for y axis
    @param title Title of graph (facultative)
    @param xlabel Label for x axis (facultative)
    @param ylabel Label for y axis (facultative)
    @param basename Output basename of the image file (Default "out")
    @param img_type Type of the image file (Default "png")
    @param baseline lower value of the colorated area (Default 0)
    @param xsize Width of the graphics (Default 15)
    @param ysize Heigth of the graphics (Default 10)
    @param dpi Resolution of the graphics (Default 100)
    @param fill_color Color of the filled area (Default 'green')
    """

    # Require the Third party package matplotlib
    from matplotlib import pyplot as plt

    # Create a figure object and adding details
    fig = plt.figure(figsize=(xsize, ysize), dpi=dpi)
    if title:
        plt.title(title)
    if xlabel:
        plt.xlabel(xlabel)
    if ylabel:
        plt.ylabel(ylabel)

    # Plot an area representing the coverage depth
    plt.fill_between(X, Y, baseline, facecolor=fill_color, alpha=0.5)

    # Tweak spacing to prevent clipping of ylabel
    plt.subplots_adjust(left=0.15)

    # Export figure to file
    try:
        fig.savefig(basename+"."+img_type, format = img_type)
    except ValueError as E:
        print (E)
        print ("Saving file as png")
        fig.savefig(basename+".png", format = "png")
