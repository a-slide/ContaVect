# -*- coding: utf-8 -*-

"""
@package    ContaVect
@brief
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

# Standard library packages import
from gzip import open as gopen
from sys import stdout
from collections import OrderedDict

# Third party packages import
from Bio import SeqIO

# Local Package import
from pyDNA.FileUtils import is_gziped, is_readable_file
from Sequence import Sequence

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Reference(object):
    """
    @class  Reference
    @brief  Object oriented class containing informations of reference
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~CLASS FIELDS~~~~~~~#

    Instances = [] # Class field used for instance tracking
    id_count = 1

    #~~~~~~~CLASS METHODS~~~~~~~#

    @ classmethod
    def next_id (self):
        cur_id = self.id_count
        self.id_count +=1
        return cur_id

    @ classmethod
    def countInstances (self):
        return len(self.Instances)

    @ classmethod
    def refLen (self):
        length=0
        for ref in self.Instances:
            for seq in ref.seq_dict.values():
                length+=seq.length
        return length

    @ classmethod
    def allSeqList (self):
        if self.Instances:
            all_seq_list=[]
            for ref in self.Instances:
                all_seq_list.extend([name for name in ref.seq_dict.keys()])
            return all_seq_list
        else:
            return []

    @ classmethod
    def allName (self):
        return [ref.name for ref in self.Instances]

    @ classmethod
    def allFasta (self):
        return [ref.fasta for ref in self.Instances]

    @ classmethod
    def getInstances (self):
        return self.Instances

    @ classmethod
    def printInstances (self):
        for ref in self.Instances:
            print (repr(ref))

    @ classmethod
    def reprInstances (self):
        msg = ""
        for ref in self.Instances:
            msg+= repr(ref)+"\n"
        return msg

    @ classmethod
    def resetInstances (self):
        print "Clearing Reference instances list"
        self.Instances = []
        self.id_count = 0

    @ classmethod
    def addRead (self, seq, read):
        """
        class method attibuting a read to a sequence by searching this sequence name in
        all Reference instances.
        """
        for ref in self.Instances:
            if seq in ref.seq_dict:
                ref.seq_dict[seq].add_read(read)
                ref.nread+=1
                return

        raise ValueError, "Seq name not found in references"

    @ classmethod
    def set_global (self, key, value):
        for ref in self.Instances:
            ref.set(key, value)

    @ classmethod
    def mk_output_global (self, outpath="./out"):
        for ref in self.Instances:
            ref.mk_output (outpath)

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    ###### Reference will create a sorted bam and a bai by default
    ###### For the other output type it will create specific objects


    def __init__(self, name, fasta, output_type):
        """
        @param name Name of the reference
        @param fasta Path to the fasta reference needed to determine sequence name associated
        with this reference
        @param output_type type of output file to generate
        """

        print ("Creating reference object {}".format(name))

        # test fasta readability
        is_readable_file(fasta)

        # Create self variables
        self.name = name
        self.id = self.next_id()
        self.fasta = fasta
        self.output_type = output_type #############################################################################

        # Define additional variables
        self.seq_dict = {}
        self.bam_header = ""
        self.nread = 0

        # Parse the fasta reference
        self.seq_dict = self._fasta_reader()

        # Add the instance to the class instance tracking list
        self.Instances.append(self)

    def __str__(self):
        msg = "REFERENCE {}".format(self.id)
        msg+= "\tName: {}\n".format(self.name)
        msg+= "\tFasta_path: {}\n".format(self.fasta)
        msg+= "\tTotal reference length: {}\n".format(len(self))
        msg+= "\tNumber of sequences: {}\n".format(len(self.seq_dict))
        #msg+= "\tSequence list:\n"
        #for seq in self.seq_dict.values():
        #    msg+= "\t* {}".format(repr(seq))

        for i in [self.bam_maker, self.cov_maker, self.var_maker]:
            msg+= repr(i)

        if self.nread:
            msg+= "\tTotal read mapped: {}\n".format(self.nread)
        return (msg)

    def __repr__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def __len__(self):
        return sum([len(seq) for seq in self.seq_dict.values()])

    def get(self, key):
        return self.__dict__[key]

    def set(self, key, value):
        self.__dict__[key] = value

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def mk_output (self, outpath="./out"):
        """
        Create output files according to user specifications
        """
        print "Processing reference :{}\tReads aligned :{} ".format(self.name, self.nread)
        print "\tPreparing data..."
        # Generate a simple dictionary associating seq name and read_list
        read_dict = {name: seq.read_list for name, seq in self.seq_dict.items()}

        # Call Behavior methods
        self.bam_maker.make(self.bam_header, read_dict, outpath, self.name)
        self.cov_maker.make(self.bam_maker.bam, self.bam_maker.bai, outpath, self.name)
        self.var_maker.make(self.bam_maker.bam, self.bam_maker.bai, outpath, self.name)

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _fasta_reader(self):
        """
        Read seq names a fasta file and verify the absence of duplicates
        Create a dictionary CV_Reference.Sequence object indexed by sequence name
        """

        # Init a file pointer
        try:
            fp = gopen(self.fasta,"rb") if is_gziped(self.fasta) else open(self.fasta,"rb")

        except Exception as E:
            fp.close()
            raise Exception (E.message+"Can not create a list of sequence from{}".format(self.name))

        # Create a dict that remenbered the order of items added
        seq_dict=OrderedDict()
        for seq in SeqIO.parse(fp, "fasta"):
            # verify the absence of the sequence name in the current list and in other ref seq_dict
            if seq.id in self.allSeqList():
                raise Exception ("{} is duplicated\n".format(seq.id))
            else:
                seq_dict[seq.id] = Sequence(name=seq.id, length=len(seq))
            stdout.write("*")
            stdout.flush()

        print ("")
        fp.close()

        return seq_dict
