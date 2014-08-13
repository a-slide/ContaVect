#~~~~~~~GLOBAL IMPORTS~~~~~~~#

# Standard library packages import
import gzip
from sys import stdout
from os import path
from random import randint
import csv

# Third party packages import
from Bio import SeqIO
import pysam
from matplotlib import pyplot

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
        """
        class method attibuting a read to a sequence by searching this sequence name in
        all Reference instances.
        """
        for ref in self.Instances:
            for seq in ref.seq_list:
                if seqname == seq.name:
                    seq.read_list.append(read)
                    ref.nread+=1
                    seq.nread+=1
                    return
        raise ValueError, "Seq name not found in references"

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self,
            fasta_path,
            mk_sam=False,
            mk_covgraph=False,
            mk_bedgraph=False,
            mk_pileup=False):
        """
        """
        # Init object variables
        self.name = file_basename(fasta_path)
        print ("Creating reference object {}".format(self.name))
        self.id = self.next_id()
        self.fasta_path = fasta_path
        self.seq_list = []
        self.nread = 0

        # Store flags and path of output files
        self.mk_sam = mk_sam
        self.sam =""
        self.mk_covgraph = mk_covgraph
        self.covgraph = ""
        self.mk_bedgraph = mk_bedgraph
        self.bedgraph = ""
        self.mk_pileup = mk_pileup
        pileup = ""

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
        msg+= "\tTotal reference length: {}\n".format(sum([seq.length for seq in self.seq_list]))
        msg+= "\tTotal read mapped: {}\n".format(sum([seq.nread for seq in self.seq_list]))
        msg+= "\tRequired output:"
        if self.mk_sam:
            msg+= "  Sam file"
        if self.mk_covgraph:
            msg += "  Coverage graph"
        if self.mk_bedgraph:
            msg += "  Bedgraph file"
        if self.mk_pileup:
            msg += "  Pileup file"
        msg+="\n"
        return (msg)

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def mk_output (self, outpath="./out"):
        """
        Create output files according to user specifications
        """
        print "Generate output for reference : {}".format(self.name)
        # Sort read in list of each sequence
        print ("\tSort sequence by coordinates")
        for seq in self.seq_list:
            seq._sort_read()

        # Minimal output
        print ("\tCreate a bam file...")
        self.bam = self._mk_bam_sam(bam=True, outpath=outpath)
        print ("\tCreate a bam index...")
        pysam.index(self.bam)

        # Facultative output
        if self.mk_sam:
            print ("\tCreate a sam file...")
            self.sam = self._mk_bam_sam(bam=False, outpath=outpath)
        if self.mk_covgraph:
            print ("\tCreate a coverage graph...")
            self.covgraph = self._mk_covgraph(outpath=outpath)
        if self.mk_bedgraph:
            print ("\tCreate a bedGraph file...")
            self.sam = self._mk_bedgraph(outpath=outpath, min_depth=3)
        if self.mk_pileup:
            print ("\tCreate a pileUp file...")
            self.sam = self._mk_pileup(outpath=outpath, min_depth=500, min_freq=0.01)

    #~~~~~~~PRIVATE METHODS~~~~~~~#

################################################## TO DO = get methods out of the class (if possible)

    def _fasta_reader(self):
        """
        Read seq names a fasta file and verify the absence of duplicates
        Integrated in
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


    def _mk_bam_sam(self, bam=True, outpath="./out"):
        """
        """
        # Prepare the header for pysam
        header = {}
        header['HD'] =  {'VN': '1.5', 'SO' : 'coordinate'}
        header['SQ'] = [{'LN': seq.length, 'SN': seq.name} for seq in self.seq_list]

        # Init a file and write the header in sam or bam format
        if bam:
            outname = "{}_{}.bam".format(outpath, self.name)
            outfile = pysam.Samfile(outname, "wb", header=header)
        else:
            outname = "{}_{}.sam".format(outpath, self.name)
            outfile = pysam.Samfile(outname, "wh", header=header)

        # Write reads
        for id, seq in enumerate(self.seq_list):
            for read in seq.read_list:
                a = pysam.AlignedRead()
                a.qname = read.qname
                a.seq = read.seq
                a.flag = read.flag
                a.rname = id
                a.pos = read.pos
                a.mapq = read.mapq
                a.cigar = read.cigar
                a.mrnm = -1
                a.mpos = -1
                a.isize = 0
                a.qual = read.qual
                a.tags = read.tags
                outfile.write(a)
        outfile.close()
        return outname


    def _mk_covgraph (self, outpath="./out"):
        """
        """
        # If needed generate the coverage over seq
        for seq in self.seq_list:
            if not seq.coverage:
                seq._coverage()

            # Create a figure object and adding details
            fig = pyplot.figure(figsize=(50, 10), dpi=200)
            pyplot.title("Coverage of reads over {}".format (self.name))
            pyplot.ylabel('Count')
            pyplot.xlabel('Position')

            # List of numbers for x axis positions
            x = [i+1 for i in range (seq.length)]

            # Plot an area representing the coverage depth
            pyplot.fill(x,seq.coverage, facecolor='green', alpha=0.5)

            # Tweak spacing to prevent clipping of ylabel
            pyplot.subplots_adjust(left=0.15)

            # Export figure to file
            outname = "{}_{}_{}.svg".format(outpath, self.name, seq.name)
            fig.savefig(outname)
            return outname


    def _mk_bedgraph (self, outpath="./out", min_depth=0):
        """
        """
        outname =  "{}_{}.bedgraph".format(outpath, self.name)

        with open (outname, "wb") as outfile:

            # Write bedGraph header
            outfile.write ("track type={} name={} color={},{},{} visibility={}\n".format(
                "bedGraph", self.name, randint(0,255), randint(0,255), randint(0,255),"full"))

            # If needed generate the coverage over seq
            for seq in self.seq_list:
                if not seq.coverage:
                    seq._coverage()

                start = -1
                depth_prec = 0

                # Start to iterate over the coverage list
                for position, depth in enumerate (seq.coverage):
                    if depth >= min_depth:
                        if depth != depth_prec:
                            if start != -1:
                                outfile.write("{}\t{}\t{}\t{}\n".format(
                                    seq.name, start, position-1, depth_prec))
                            start = position
                            depth_prec = depth

                    else:
                        if start != -1:
                            outfile.write("{}\t{}\t{}\t{}\n".format(
                            seq.name, start, position-1, depth_prec))
                        start = -1
                        depth_prec = 0
        return outname

    def _mk_pileup (self, outpath="./out", min_depth=1000, min_freq=0.01):
        """
        """
        # Create a list and for results collecting
        out_list =[]

        # Open a handle on the bamfile and create an output file
        bamfile = pysam.Samfile( self.bam, "rb")
        # Open a PileUpColumn generator with an high limit of depth to avoid cuts
        # Analyse only if the sequencing depth is sufficient
        for PileUpCol in bamfile.pileup(max_depth=1000000):
            if PileUpCol.n > min_depth:

                # Fill a dictionary for positions of reads overlapping each positions.
                posDic = {'A':0,'C':0,'T':0,'G':0,'Del':0,'N':0}
                for read in PileUpCol.pileups:
                    if read.is_del:
                        posDic['Del']+=1
                    else:
                        try:
                            posDic[read.alignment.seq[read.qpos].upper()]+=1
                        except IndexError:
                            posDic['N']+=1

                # Define a threshold above which variations are not considered
                threshold = int(PileUpCol.n*min_freq)

                # If more than one frequent DNA base or indel was found at this position
                if sum([1 for base, count in posDic.items() if count >= threshold]) >= 2:
                    out_list.append([
                        self.name,
                        bamfile.getrname(PileUpCol.tid),
                        PileUpCol.pos,
                        PileUpCol.n,
                        posDic['A'],
                        posDic['T'],
                        posDic['C'],
                        posDic['G'],
                        posDic['N'],
                        posDic['Del'],
                        round(posDic['A']/float(PileUpCol.n), 3),
                        round(posDic['T']/float(PileUpCol.n), 3),
                        round(posDic['C']/float(PileUpCol.n), 3),
                        round(posDic['G']/float(PileUpCol.n), 3),
                        round(posDic['N']/float(PileUpCol.n), 3),
                        round(posDic['Del']/float(PileUpCol.n), 3)])

        # Create a file to write out the list if results were found
        if out_list:
            outname =  "{}_{}_pileup.csv".format(outpath, self.name)
            with open(outname, 'wb') as csvfile:
                writer = csv.writer(csvfile, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
                writer.writerow(["Ref","Seq","Pos","Total","CountA","CountT","CountC","CountG",
                    "CountN", "CountDel","FreqA","FreqT","FreqC","FreqG","FreqN","FreqDel"])
                for line in out_list:
                    writer.writerow(line)
            return outname

        else:
            print("\t  No frequent variation found")
            return ""

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
        self.coverage = []

    def __repr__(self):
        msg = "{}({}bp) : {} read(s)\n".format(self.name, self.length, self.nread)
        return (msg)

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)


    def _sort_read (self):
        """
        sort read in read_list acording to their leftmost position
        """
        self.read_list.sort(key = lambda x: x.pos)

    def _coverage (self):
        """
        Create a coverage depth over the length of the
        """
        # Init a null coverage list
        self.coverage = [0 for i in range(self.length)]

        # Fill the coverage with read values
        for read in self.read_list :
            for i in range (read.pos, read.pos+read.alen):
                self.coverage[i]+=1
