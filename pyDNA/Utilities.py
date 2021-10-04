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
        for key, value in list(opt_dict.items()):
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

def is_readable_file (file_path):
    """
    Verify the readability of a file or list of file
    """
    from os import access, R_OK
    
    if isinstance (file_path, str):
        file_path = [file_path]
    
    if isinstance (file_path, list):
        for fp in file_path :
            if not access(fp, R_OK):
                raise ValueError ("{} is not a valid file".format(fp))
                
    else:
        raise TypeError ("File path shoud be a path or a list of paths")
        

def is_gziped (fp):
    """
    @param fp path to a files eventually gzipped
    @return True is ended by gz extension else false
    """
    return fp[-2:].lower() == "gz"

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
        print(('Error: %s' % e))
    # eg. source or destination doesn't exist
    except IOError as e:
        print(('Error: %s' % e.strerror))


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
        print(("Compressing {}".format(in_path)))
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
                print("Can't remove {}".format(out_path))

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
        print(("Uncompressing {}".format(in_path)))
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
                print("Can't remove {}".format(out_path))

def expand_file (infile, outdir="./"):
    """
    expand file in outdir if the file are gzipped
    Else the file won't be modified and won't be moved from it's current directory
    @param infile Path to a file eventually gzipped
    @param outdir Path of the directory in which to uncompress or copy the files
    @return A path to an uncompressed file
    """
    # Function specific imports
    from os import path

    assert path.isfile(infile), "{} is not a valid file".format(infile)

    # Extract if gziped
    if is_gziped(infile):
        return fgunzip (in_path=infile, out_path=path.join(outdir,file_name(infile)[:-3]))
    
    # else just return the original file path
    else:
        return infile

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

    with openout(outpath, "w") as out_handle:
        # Iterate over files in the input list
        for inpath in inpath_list:

            # Open according to the compression
            openin = gzip.open if is_gziped(inpath) else open
            with openin (inpath, "r") as in_handle:
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

    print(("{} files merged in {}s\n".format (len(inpath_list), round(time()-stime,3))))
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
        if is_gziped(filename):
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
    if is_gziped(filename):
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
    @param sequence DNA sequence
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


def gb_to_bed(gb_file, track_description="User Supplied Track", features_type = []):
    """
    write a bed file containing annotations from a gb file 
    @param gb_file path to a gb file containing a single molecule
    @param track_description Line of text descitpion for the bed track (60 char max)
    @param feature_type Restrict to the list of feature indicated 
    """
    # Specific imports
    from Bio import SeqIO
    from unicodedata import normalize
    
    # open in and out files
    with open (gb_file, "rb") as gb:
        outf = file_basename(gb_file)+".bed"
        with open(outf, 'w') as bed:
            
            # parse record
            record = SeqIO.read( gb, "genbank")
            print(("{} features to be parsed\n".format(len(record.features))))
            
            # write bed header
            bed.write("""track name={} description="{}" visibility=2 colorByStrand="255,0,0 0,0,255"\n""".format(
                record.id,
                track_description))
                
            for feature in record.features:
                
                # if specific deatures are required
                if features_type and feature.type not in features_type:                
                    continue
                
                start = feature.location.start.position
                stop = feature.location.end.position
                strand = '-' if feature.strand < 0 else '+'
                
                # try several alternative key for feature name else skip the feature 
                if 'gene' in feature.qualifiers:
                    name = feature.qualifiers['gene'][0]
                elif 'locus_tag' in feature.qualifiers:
                    name = feature.qualifiers['locus_tag'][0]
                elif 'label' in feature.qualifiers:
                    name = feature.qualifiers['label'][0]
                else:
                    continue
                    
                # Normalize name to asci
                name = str(name.decode('ascii', 'ignore'))
                bed.write("{0}\t{1}\t{2}\t{3}\t1000\t{4}\t\n".format(
                    record.id,
                    start,
                    stop,
                    name,
                    strand))
    return outf

def fetch_count_read (alignment_file, seq_name, start, end):
    """
    Count the number of read that are at least partly overlapping a specified chromosomic region
    @param alignment_file Path to a sam or a bam file
    @param seq_name Name of the sequence where read are to be aligned on
    @param start Start genomic coordinates of the area of alignment
    @param end End End genomic coordinates of the area of alignment
    """
    # Specific imports
    from pysam import AlignmentFile
    
    # Init a generator on the sam or bam file with pysam
    if alignment_file[-3:].lower() == "bam":
        al = AlignmentFile(alignment_file, "rb")
        
    elif alignment_file[-3:].lower() == "sam":
        al = AlignmentFile(alignment_file, "r")
    
    else:
        raise Exception("Wrong file format (sam or bam)") 
    
    # Count read aligned at least partly on the specified region
    n = 0
    for i in al.fetch(seq_name, start, end):
        n += 1
        
    al.close()
    
    return n
    
    
def fetch_all_bam (bam_pattern, seq_name, coord_list, outname):
    """
    for all bam files matching the pattern, count the number of read overlapping list of coordinates
    and generate a file report
    @param bam_pattern Pattern to match in bam file to be included in the analysis 
    @param seq_name Name of the sequence where read are to be aligned on
    @param coord_list list of coordinate start and end where to find overlapping reads
    @param outname Name of the output csv file repport
    """
    import os, csv, glob
    import pysam
    
    # Find bam matching the patern and sorting the list alphabetically
    bam_list = list(glob.iglob("*"+bam_pattern+"*"+".bam"))
    bam_list.sort()
    print ("Files analysed :")
    print (bam_list)
    
    # Generate a table header
    header = ["Coordinates"]
    for start, end in coord_list:
        header.append("{}:{}".format(start, end))
    
    # Generate a list to store hit found
    all_hits = []
    all_hits.append(header)

    # Fetch file for coordinates in coord list 
    for bam in bam_list:
        print(("\nAnalysing file : "+bam))
        
        if not os.path.isfile(bam+".bai"):
            print ("\tGenerate a Bam index")
            pysam.index(bam)
        
        hits = [bam]
        for start, end in coord_list:
            hits.append(fetch_count_read(bam, seq_name, start, end))
        
        all_hits.append(hits)
    
    # Finally write a new table
    print ("\nWrite results in a csv table") 
    with open(outname, 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for i in all_hits:
            writer.writerow(i)

    print("Done")

    return

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
