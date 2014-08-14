#~~~~~~~GLOBAL IMPORTS~~~~~~~#

# Standard library packages import
import gzip
from sys import stdout

# Third party packages import
from Bio import SeqIO

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
        return [ref.ref_fasta for ref in self.Instances]

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
    def set(self, key, value):
        for ref in reference.Instances:
            ref.set(key, value)
    
    @ classmethod
    def mk_output (self, outpath="./out"):
        for ref in reference.Instances:
            ref.mk_output (outpath)

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self, name, ref_fasta, Bam, Coverage, Variant):
        """
        @param name Name of the reference
        @param ref_fasta Path to the fasta reference needed to determine sequence name associated
        with this reference
        @param Bam BamMaker object to create bam sam and bai
        @param Covgraph Coverage object to create a coverage graph, bedgraph and bed
        @param Variant Variant object to create frequent variant report file
        """

        # Store object variables
        print ("Creating reference object {}".format(name))
        self.name = name
        self.id = self.next_id()
        self.ref_fasta = ref_fasta
        self.Bam = Bam
        self.Coverage = Coverage
        self.Variant = Variant

        # Define additional variables
        self.seq_dict = {}
        self.bam_header = ""
        self.nread = 0

        # Parse the fasta reference
        self.seq_dict = self._fasta_reader()

        # Add the instance to the class instance tracking list
        self.Instances.append(self)

    def __repr__(self):
        msg = "REFERENCE {}".format(self.id)
        msg+= "\tName: {}\n".format(self.name)
        msg+= "\tFasta_path: {}\n".format(self.ref_fasta)
        msg+= "\tTotal reference length: {}\n".format(sum([seq.length for seq in self.seq_dict.values()]))
        msg+= "\tSequence list:\n"
        for seq in self.seq_dict.values():
            msg+= "\t* {}".format(repr(seq))
        
        for i in [self.Bam, self.Coverage, self.Variant]:
            msg+= repr(i)
        
        if self.nread:
            msg+= "\tTotal read mapped: {}\n".format(self.nread)
        return (msg)

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)
        
    def get(self, key):
        return self.__dict__[key]

    def set(self, key, value):
        self.__dict__[key] = value

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def mk_output (self, bam_header, outpath="./out"):
        """
        Create output files according to user specifications
        """
        print "Processing reference :{}\tReads aligned :{} ".format(self.name, self.nread)
        print "\tPreparing data..."
        # Generate a simple dictionary associating seq name and read_list
        read_dict = {name: seq.read_list for name, seq in self.seq_dict.items()}
        # Generate a simple dictionary associating seq name and coverage_list
        cov_dict = {name: seq.mk_coverage() for name, seq in self.seq_dict.items()}

        # Call Behaviour methods
        self.Bam.make(bam_header, read_dict, outpath+self.name)
        self.CovGraph.make(cov_dict, outpath, self.name)
        self.BedGraph.make(cov_dict, outpath, self.name)
        self.PileUp.make(Bam.bam, outpath, self.name)

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _fasta_reader(self):
        """
        Read seq names a fasta file and verify the absence of duplicates
        Create a dictionary CV_Reference.Sequence object indexed by sequence name
        """
        
        # Init a file pointer
        try:
            if self.ref_fasta[-2:].lower() == "gz":
                fp = gzip.open(self.ref_fasta,"rb")
            else:
                fp = open(self.ref_fasta,"rb")

        except Exception as E:
            fp.close()
            raise Exception (E.message+"Can not create a list of sequence from{}".format(self.name))
        
        seq_dict={}
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
        self.read_list = []

    def __repr__(self):
        msg = "{} ({} bp)".format(self.name, self.length)
        if self.nread:
            return (msg + "{} read(s)\n".format(self.nread))
        else:
            return (msg + "\n")

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)
        
    def get(self, key):
        return self.__dict__[key]

    def set(self, key, value):
        self.__dict__[key] = value

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def add_read (self, read):
        """
        Add a read to read_list and update the counter
        """
        self.read_list.append(read)
        self.nread+=1

    def sort_read (self):
        """
        sort read in read_list acording to their leftmost position
        """
        self.read_list.sort(key = lambda x: x.pos)

    def mk_coverage (self):
        """
        Create a coverage depth over the length of the sequence
        """
        # Init a null coverage list
        coverage = [0 for i in range(self.length)]

        # Fill the coverage with read values
        for read in self.read_list :
            for i in range (read.pos, read.pos+read.alen):
                coverage[i]+=1

        return coverage
