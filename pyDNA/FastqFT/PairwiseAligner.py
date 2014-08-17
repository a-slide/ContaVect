#~~~~~~~GLOBAL IMPORTS~~~~~~~#
# Third party package import
from Bio import pairwise2

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class PairwiseAligner(object):
    """
    @class  PairwiseAligner
    @brief  Require the third party package Biopython
    """
    # TODO write a wrapper for an external faster SW aligner
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __repr__(self):
        return self.__str__() + self.get_report()

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def __init__ (self, max_alignment=5, match=2, mismatch=-2, open_gap=-2, extend_gap=-2, cutoff=1.5):
        """
        @param max_alignment Maximal number of alignement per read to output
        @param match Gain value in case match
        @param mismatch Penality value in case mismatch
        @param open_gap Penality value in case opening a gap
        @param extend_gap Penality value in case extending a gap
        @param cutoff Raw SW score divided by the lenght of the adapter
        """

        # Store parameters in object variables
        pairwise2.MAX_ALIGNMENTS = max_alignment
        self.match = match
        self.mismatch = mismatch
        self.open_gap = open_gap
        self.extend_gap = extend_gap
        self.cutoff = cutoff

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def find_match (self, adapter, sequence):
        """
        Pairwise alignment using the buildin Biopython function pairwise2
        """

        # Perform a local SW alignment with specific penalties
        match_list = pairwise2.align.localms(
            adapter,
            sequence,
            self.match,
            self.mismatch,
            self.open_gap,
            self.extend_gap)

        #for seqA, seqB, score, begin, end in match_list:
            #if score/len(adapter) > self.cutoff:
                #print("{}\n{}\nScore:{}\tBegin:{}\tEnd:{}".format(seqA, seqB, score, begin, end))

        # Return begin and end position if a match has a score higher than the cutoff value
        return [[begin, end] for seqA, seqB, score, begin, end in match_list if score/len(adapter) > self.cutoff]

    def get_report (self):
        """
        """
        report = "====== PAIRWISE ALIGNER ======\n\n"
        report += "  Parameters of pairwise2\n"
        report += "  Maximal number of alignement : {}\n".format(pairwise2.MAX_ALIGNMENTS)
        report += "  Match bonus : {}\n".format(self.match)
        report += "  Mismatch Penality : {}\n".format(self.mismatch)
        report += "  Gap open Penality : {}\n".format(self.open_gap)
        report += "  Gap extend Penality : {}\n".format(self.extend_gap)
        report += "  Score cutoff : {}\n".format(self.cutoff)

        return report
