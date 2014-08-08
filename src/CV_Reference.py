#~~~~~~~GLOBAL IMPORTS~~~~~~~#

# Standard library packages import
import gzip
from sys import stdout
from os import path

# Third party packages import
from Bio import SeqIO

# Local Package import
from Utilities import file_basename

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Reference(object):
    """
    @class  Reference
    @brief  Object oriented class containing informations of reference
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~CLASS FIELDS~~~~~~~#

    Instances = [] # Class field used for instance tracking
    id_count = 0

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
            for seq in ref.seq_list:
                length+=seq.length
        return length

    @ classmethod
    def allSeqList (self):
        if self.Instances:
            all_seq_list=[]
            for ref in self.Instances:
                all_seq_list.extend([seq.name for seq in ref.seq_list])
            return all_seq_list
        else:
            return []

    @ classmethod
    def allName (self):
        return [ref.name for ref in self.Instances]

    @ classmethod
    def allFasta (self):
        return [ref.fasta_path for ref in self.Instances]

    @ classmethod
    def getInstances (self):
        return self.Instances

    @ classmethod
    def printInstances (self):
        for ref in self.Instances:
            print (repr(ref))

    @ classmethod
    def resetInstances (self):
        print "Clearing Reference instances list"
        self.Instances = []
        self.id_count = 0
    
    @ classmethod
    def addRead (self, seqname, read):
        for ref in self.Instances:
            for seq in ref.seq_list:
                if seqname == seq.name:
                    ref.read_list.append(read)
                    ref.nread+=1
                    seq.nread+=1
                    return 1
        return 0
    

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self,
            fasta_path,
            mk_sam=False,
            mk_covgraph=False,
            mk_bedgraph=False,
            mk_vcf=False):
        """
        """
        # Init object variables 
        self.name = file_basename(fasta_path)
        print ("Creating reference object {}".format(self.name))
        self.id = self.next_id()
        self.fasta_path = fasta_path
        self.seq_list = []
        self.read_list = []
        self.nread = 0
        self.mk_sam = mk_sam
        self.mk_covgraph = mk_covgraph
        self.mk_bedgraph = mk_bedgraph
        self.mk_vcf = mk_vcf
                
        # parse the fasta reference
        self.seq_list = self._fasta_reader()
        # Add the instance to the class instance tracking list
        self.Instances.append(self)

    def __repr__(self):
        msg = "REF {}".format(self.id)
        msg+= "\tName: {}\n".format(self.name)
        msg+= "\tFasta_path: {}\n".format(self.fasta_path)
        msg+= "\tSequence list:\n"
        for seq in self.seq_list:
            msg+= "\t{}\n".format(repr(seq))
        msg+= "\tRequired output:"
        if self.mk_sam:
            msg+= "  Sam file"
        if self.mk_covgraph:
            msg += "  Coverage graph"
        if self.mk_bedgraph:
            msg += "  Bedgraph file"
        if self.mk_vcf:
            msg += "  VCF file"
        msg+="\n"
        return (msg)

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def _fasta_reader(self):
        """
        Read seq names a fasta file and verify the absence of duplicates
        """
        try:
            if self.fasta_path[-2:].lower() == "gz":
                fp = gzip.open(self.fasta_path,"rb")
            else:
                fp = open(self.fasta_path,"rb")

        except Exception as E:
            fp.close()
            raise Exception (E.message+"Can not create a list of sequence from{}".format(self.name))

        seq_list=[]
        for seq in SeqIO.parse(fp, "fasta"):
            # verify the absence of the sequence name in the current list and in other ref seq_list
            if seq.id in self.allSeqList():
                raise Exception ("{} is duplicated\n".format(seq.id))
            else:
                seq_list.append(Sequence(name=seq.id, length=len(seq)))

            stdout.write("*")
            stdout.flush()

        print ("")
        fp.close()

        return seq_list

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Sequence(object):
    """
    @class  Sequence
    @brief  Object oriented class containing informations of a sequence from a reference file
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    def __init__(self, name, length):
        self.name = name
        self.length = length
        self.nread = 0
    
    def __repr__(self):
        msg = "{}({}bp) : {} read(s)".format(self.name, self.length, self.nread)
        return (msg)

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)
