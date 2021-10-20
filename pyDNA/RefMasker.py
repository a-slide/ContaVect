# TODO : TRANSFORM INTO A CLASS AND CREATE A REPORT OF REGION TRIMMED

#~~~~~~~GLOBAL IMPORTS~~~~~~~#

# Standard library packages import
from os import remove, path
import gzip
from time import time
from sys import stdout

# Third party package import
from Bio import SeqIO
from Bio.Seq import MutableSeq

# Local library packages import
from pyDNA.Utilities import import_seq, file_basename, mkdir
from .Blast import Blastn

#~~~~~~~MAIN METHODS~~~~~~~#

def mask (  subject_fasta,
            hit_list,
            ref_outdir="./references/",
            ref_outname="masked_ref.fa",
            compress_ouput=True ):
    """
    Import a reference fasta sequence, Mask positions indicated by hits from a hit_list and write
    the modified fasta sequence in a new file.
    @param subject_fasta Fasta sequence of the subject to edit (can be gzipped)
    @param hit_list List of hit objects. Hits need at least 3 fields named s_id, s_start and s_end
    coresponding to the name of the sequence matched, and the hit start/end (0 based).
    @param ref_outdir Directory where the masked reference will be created
    @param ref_outname Name of the masked reference
    @param compress_ouput If true the output will be gzipped
    @return A path to the modified sequence if the hit list was valid.
    """

    # Test if object the first object of hit_list have the require s_id, s_start and s_end fields
    try:
        a = hit_list[0].s_id
        a = hit_list[0].s_start
        a = hit_list[0].s_end

    except IndexError:
        print ("No hit found, The subject fasta file will not be edited")
        return subject_fasta
    except AttributeError as E:
        print ("The list provided does not contain suitable hit object, The subject fasta file will not be edited")
        return subject_fasta
    
    # Initialize output folder
    mkdir(ref_outdir)

    # Initialize input fasta file
    if subject_fasta[-2:].lower() == "gz":
        in_handle = gzip.open(subject_fasta, "r")
    else:
        in_handle = open(subject_fasta, "r")

    # Initialize output fasta file
    if compress_ouput:
        ref_path = path.join (ref_outdir, ref_outname+".gz")
        out_handle = gzip.open(ref_path, 'w')
    else:
        ref_path = path.join (ref_outdir, ref_outname)
        out_handle = open(ref_path, 'w')

    # Generate a list of ref that will need to be modified
    id_list = list({hit.s_id:0 for hit in hit_list}.keys())
    id_list_decode = [contig.decode() for contig in id_list]

    # Iterate over record in the subject fasta file
    print(("Masking hit positions and writting a new reference for {} ".format(ref_outname)))
    i=j=0
    start_time = time()
    for record in SeqIO.parse(in_handle, "fasta"):
        # Progress Marker
        stdout.write("*")
        stdout.flush()

        # Check if the record is in the list of record to modify
        if record.id in id_list_decode:
            i+=1
            # print("Hit found in {}. Editing the sequence".format(record.id))
            # Casting Seq type to MutableSeq Type to allow string editing
            record.seq = MutableSeq(record.seq)

            # For each hit in the list of hit found
            for hit in hit_list:
                # print("Record id: {}".format(record.id))
                # print("Hit id: {}".format(hit.s_id))
                try:
                    hit.s_id = hit.s_id.decode()
                except:
                    pass
                if record.id == hit.s_id:
                    # For all position between start and end coordinates modify the base by N
                    for position in range (hit.s_start, hit.s_end):
                        record.seq[position]= 'n'
        else:
            j+=1
            # print("No hit found in {}".format(record.id))

        # Finally write the sequence modified or not
        out_handle.write(record.format("fasta"))
    print("")
    # Report informations
    print(("{} sequence(s) from {} modified in {}s".format(i,ref_outname, round(time()-start_time),2)))

    # Close files and return the masked ref path
    in_handle.close()
    out_handle.close()
    return ref_path
