# TODO : Store lenght of each ref and seq somewhere for future RPKM 

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
    def allSeqList (self):
        if self.Instances:
            all_seq_list=[]
            for ref in self.Instances:
                all_seq_list.extend(ref.seq_list)
            return all_seq_list
        else:
            return []

    @ classmethod
    def seqOrigin (self, seqname):
        for ref in self.Instances:
            if seqname in ref.seq_list:
                return ref.name
        return None
    
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
        
    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self,
            fasta_path,
            mk_sam=False,
            mk_covgraph=False,
            mk_bedgraph=False,
            mk_vcf=False):
        """
        """
        self.name = file_basename(fasta_path)
        print ("Creating reference object {}".format(self.name))
        self.id = self.next_id()
        self.fasta_path = path.abspath(fasta_path)
        self.seq_list = self._fasta_reader()
        #self.bam = []
        self.mk_sam = mk_sam
        self.mk_covgraph = mk_covgraph
        self.mk_bedgraph = mk_bedgraph
        self.mk_vcf = mk_vcf
        
        # Add the instance to the class instance tracking list
        self.Instances.append(self)

    def __repr__(self):
        msg = "REF {}".format(self.id)
        msg+= "\tName: {}\n".format(self.name)
        msg+= "\tFasta_path: {}\n".format(self.fasta_path)
        msg+= "\tSequence list:  {}\n".format("  ".join(self.seq_list))
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
            if seq.id in seq_list:
                raise Exception ("{} is duplicateda in the current reference\n".format(seq.id, self.name))
            if seq.id in self.allSeqList():
                raise Exception ("{} is duplicated in another reference\n".format(seq.id))
            else:
                seq_list.append(seq.id)

            stdout.write("*")
            stdout.flush()
        
        print ("")
        fp.close()
        
        return seq_list


