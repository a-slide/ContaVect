#~~~~~~~GLOBAL IMPORTS~~~~~~~#

# Third party packages import
import pysam

# Local Package import
from pyDNA.Utilities import fill_between_graph

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class CoverageMaker (object):
    """
    Generate a Bedgraph Bed and or Coverage graphics from a sorted bam file.
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__ (self, min_depth=0, make_bedgraph=True ,make_bed=False, make_covgraph=False):
        """
        Create a CoverageMaker object
        @param min_depth Minimal depth to report. Elsewhere it will be considered to be null
        @param make_bedgraph Allow the object to create bedgraph file = region below the min depth
        are not written and contiguous same coverage areas are grouped in a single entry.
        @param make_bed Allow the object to create bedgraph file = One entry is created for each
        bases of the sequence even if bellow the minimal depth. In this case the depth will be 0
        @param make_covgraph Allow the object to create coverage graphics (1 per sequence mapped in
        the bam file) using Matplotlib.
        """
        # Creating object variables
        self.min_depth = min_depth

        self.make_bedgraph = make_bedgraph
        self.bedgraph = ""

        self.make_bed = make_bed
        self.bed = ""

        self.make_covgraph = make_covgraph
        self.covgraph_list = []

    def __repr__(self):
        msg = "\tCOVERAGE MAKER\n"

        if not self.make_bed and not self.make_bedgraph and not self.make_covgraph:
            msg += "\t\tNo output requested\n"
            return msg

        msg+= "\t\tMinimal depth : {}\n".format(self.min_depth)
        msg+= "\t\tOutput requested :"
        if self.make_bedgraph:
            msg+= "\tBedGraph"
        if self.make_bed:
            msg+= "\tBed"
        if self.make_covgraph:
            msg+= "\tCovGraph"
        msg+= "\n"
        return msg

    def __str__(self):
        return "\n<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def get(self, key):
        return self.__dict__[key]

    def set(self, key, value):
        self.__dict__[key] = value

    #~~~~~~~PUBLIC METHODS~~~~~~~#
    
    def __call__ (self, bam_path, bai_path, outpath="./out", ref_name = "ref"):
        """
        @param bam_path Path to the bam file that will be used for coverage extraction
        @param bai_path Path to the bam index file. I won't be used but it's existence is
        required by pileup
        @param outpath Basename of the path where to output files
        @param ref_name Name of the reference genome containing the sequence listed in the bam file
        """
        self.make(bam_path, bai_path, outpath, ref_name)
    
    def make (self, bam_path, bai_path, outpath="./out", ref_name = "ref"):
        """
        @param bam_path Path to the bam file that will be used for coverage extraction
        @param bai_path Path to the bam index file. I won't be used but it's existence is
        required by pileup
        @param outpath Basename of the path where to output files
        @param ref_name Name of the reference genome containing the sequence listed in the bam file
        """
        if not self.make_bed and not self.make_bedgraph and not self.make_covgraph:
            return

        if self.make_bedgraph:
            print ("\tCreate a bedGraph File...")
            self.bedgraph = self._make_bedgraph(bam_path, outpath, ref_name)

        if self.make_bed:
            print ("\tCreate a bed File...")
            self.bed = self._make_bed(bam_path, outpath, ref_name)

        if self.make_covgraph:
            print ("\tCreate coverage graphics...")
            self.covgraph_list = self._make_covgraph(bam_path, outpath, ref_name)

    ##~~~~~~~PRIVATE METHODS~~~~~~~#

    def _make_bedgraph (self, bam_path, outpath="./out", ref_name = "ref"):
        """
        Write in a file a bedgraph report with region below the min depth not written and
        contiguous base with a same coverage are grouped in a single entry
        """

        # Define a path where to write the bedgraph
        bedgraph = "{}_{}.bedgraph".format(outpath, ref_name)

        # Open a handle on the to read the bam file and another to write the bedgraph
        with pysam.Samfile(bam_path, "rb") as bamfile:
            with open (bedgraph, "wb") as outfile:
            
            # TODO NAME SHOULD INCLUDE USER PREFIX 
            # Write bedGraph header
                outfile.write ("track type=bedGraph name={} color=0,0,0\n".format(ref_name).encode())

                # Iterate over all seq in bam header
                for seq_name in bamfile.references:

                    # Open a PileUpColumn generator with an high limit of depth to avoid cuts
                    start = depth_prec = pos_prec = -1
                    for PileUpCol in bamfile.pileup(seq_name, max_depth=10000000):

                        # If bellow the threshold write only if the last position was above
                        if PileUpCol.n < self.min_depth:
                            if start != -1:
                                outfile.write("{}\t{}\t{}\t{}\n".format(seq_name, start, pos_prec+1, depth_prec).encode())
                            start = -1
                            depth_prec = -1

                        # If above the threshold but missing bases
                        elif PileUpCol.pos != pos_prec+1:
                            if start != -1:
                                outfile.write("{}\t{}\t{}\t{}\n".format(seq_name, start, pos_prec+1, depth_prec).encode())
                            start = PileUpCol.pos
                            depth_prec = PileUpCol.n

                        # And finally if there was a change in depth.
                        elif PileUpCol.n != depth_prec:
                            if start != -1:
                                outfile.write("{}\t{}\t{}\t{}\n".format(seq_name, start, pos_prec+1, depth_prec).encode())
                            start = PileUpCol.pos
                            depth_prec = PileUpCol.n

                        # Update the position tracker
                        pos_prec = PileUpCol.pos

                    # Write the last entry for the sequence if needed
                    if start != -1:
                        outfile.write("{}\t{}\t{}\t{}\n".format(seq_name, start, pos_prec+1, depth_prec).encode())
        return bedgraph


    def _make_bed (self, bam_path, outpath="./out", ref_name = "ref"):
        """
        Write in a file a bed report with entries created for each bases of the sequence even if
        bellow the minimal depth. In this case the depth will be 0. A sequence larger than
        1 000 000 will no be included due to large memory needed for this representation.
        Bedgraph should be use instead
        """

        # Define a path where to write the bedgraph
        bed = "{}_{}.bed".format(outpath, ref_name)

        # Open a handle on the to read the bam file and another to write the bedgraph
        with pysam.Samfile(bam_path, "rb") as bam:
            with open (bed, "wb") as outfile:

            # Write bed header
                outfile.write("track type=bed name={} color=0,0,0\n".format(ref_name).encode())

                # Iterate over all seq in bam header
                for seq_dict in bam.header['SQ']:
                    seq_name = seq_dict['SN']
                    seq_len = seq_dict['LN']

                    # If the sequence if larger than 1M bp it will no be included in the bed
                    if seq_len >= 1000000:
                        continue

                    # Open a PileUpColumn generator with an high limit of depth to avoid cuts
                    pos_prec = -1
                    write = False
                    for PileUpCol in bam.pileup(seq_name, max_depth=10000000):

                        # Write the missing positions since pileup do not report uncovered regions
                        if PileUpCol.pos != pos_prec+1:
                            for i in range(pos_prec+1, PileUpCol.pos):
                                outfile.write("{}\t{}\t{}\t{}\n".format(seq_name, i, i+1, 0).encode())
                                write = True

                        # Write the current position
                        if PileUpCol.n < self.min_depth:
                            outfile.write("{}\t{}\t{}\t{}\n".format(seq_name, PileUpCol.pos, PileUpCol.pos+1, 0).encode())
                            write = True
                        else:
                            outfile.write("{}\t{}\t{}\t{}\n".format(seq_name, PileUpCol.pos, PileUpCol.pos+1, PileUpCol.n).encode())
                            write = True

                        # Update the position tracker
                        pos_prec = PileUpCol.pos

                    # Write the last entries for the sequence if needed
                    if write and pos_prec < seq_len-1:
                        for i in range(pos_prec+1, seq_len-1):
                            outfile.write("{}\t{}\t{}\t{}\n".format(seq_name, i, i+1, 0).encode())
        return bed

    def _make_covgraph (self, bam_path, outpath="./out", ref_name = "ref"):
        """
        Create a coverage graphics directly using matplotlib. Value below the min depth will be
        represented as 0. A sequence larger than 50 000 will no be included due to large memory
        needed for this representation and lisibility of graphics. Bedgraph should be use instead
        """
        covgraph_list = []

        # Open a handle on the to read the bam file
        with pysam.Samfile(bam_path, "rb") as bam:

            # Iterate over all seq in bam header
            for seq_dict in bam.header['SQ']:
                seq_name = seq_dict['SN']
                seq_len = seq_dict['LN']

                # If the sequence if larger than 50 000 bp the graphics will not be created
                if seq_len >= 50000:
                    continue

                # Open a PileUpColumn generator with an high limit of depth to avoid cuts
                pos_prec = -1
                coverage = []
                for PileUpCol in bam.pileup(seq_name, max_depth=10000000):

                    # Write the missing positions since pileup do not report uncovered regions
                    if PileUpCol.pos != pos_prec+1:
                        for i in range(pos_prec+1, PileUpCol.pos):
                            coverage.append(0)

                    # Write the current position
                    if PileUpCol.n < self.min_depth:
                        coverage.append(0)
                    else:
                        coverage.append(PileUpCol.n)

                    # Update the position tracker
                    pos_prec = PileUpCol.pos

                # This will be done only if elements were added to the coverage list
                if coverage:
                    # Write the last entries for the sequence if needed
                    if pos_prec < seq_len-1:
                        for i in range(pos_prec+1, seq_len-1):
                            coverage.append(0)

                    # Create the graph with Utilities.fill_between_graph
                    try:
                        fill_between_graph (
                            X = [i+1 for i in range (len(coverage))],
                            Y = coverage,
                            basename = "{}_{}_{}".format(outpath,ref_name,seq_name),
                            img_type = "svg",
                            title = ("Coverage of reads over {} from {}".format (seq_name, ref_name)),
                            xlabel = 'Position', ylabel = 'Count',
                            xsize = 50, ysize = 10, dpi = 150)
                        covgraph_list.append("{}_{}_{}.svg".format(outpath,ref_name,seq_name))

                    except ImportError as E:
                        print(E)
                        print("Cannot create the required CovGraph file. Skip to the next step")

        return covgraph_list

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class CoverageDecoy(object):
    """
    Decoy class implementing no behaviour
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    def __init__ (self, *args, **kwargs):
        """
        Decoy init method
        """
        pass

    def __repr__(self):
        msg = "COVERAGE DECOY\n"
        return msg

    def __str__(self):
        return "\n<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def make (self, *args, **kwargs):
        """
        Decoy make method
        """
        print ("\tNo file to be generated")
