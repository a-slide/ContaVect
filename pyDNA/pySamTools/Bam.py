#~~~~~~~GLOBAL IMPORTS~~~~~~~#

# Third party packages import
import pysam

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class BamMaker(object):
    """
    @class BamMaker
    @brief Class creating bam, sam and bam index from a collection of pysam.AlignedRead and
    a pysam pysam.Samfile.header
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__ (self, sort=True, make_bam=True, make_index=True, make_sam=False):
        """
        Create a BamMaker object
        @param sort If true read will be sorted according to sam/bam specifications
        @param make_bam If true a bam file will be generated (should remain true)
        @param make_index If true a bam index will be generated (make_bam needed)
        @param make_sam If true a sam file will be generated in addition to the bam file
        """
        # Creating object variables
        self.sort = sort
        self.make_bam = make_bam
        self.bam = ""
        self.make_index = make_index
        self.bai = ""
        self.make_sam = make_sam
        self.sam = ""

    def __repr__(self):
        msg = "\tBAM MAKER\n"

        if not self.make_bam and not self.make_sam:
            msg += "\t\tNo output requested\n"
            return msg

        msg+= "\t\tSort reads : {}\n".format(str(self.sort))
        msg+= "\t\tOutput requested :"
        if self.make_bam:
            msg+= "\tBam"
        if self.make_index:
            msg+= "\tBam_Index"
        if self.make_sam:
            msg+= "\tSam"
        msg+= "\n"
        return msg

    def __str__(self):
        return "\n<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def get(self, key):
        return self.__dict__[key]

    def set(self, key, value):
        self.__dict__[key] = value

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def make (self, header, read_col, outpath="./out", ref_name = "ref"):
        """
        @param header Multilevel dictionnary containing SAM/BAM header informations. See
        pysam.Samfile for more details
        @param read_col Ordered dictionnary of pysam.AlignedRead grouping reads in lists associated
        to the name of a sequence same as in the sam header. Alternatively a simple flat list of
        reads can be provided but in this case no sorting will be performed.
        @param outpath Basename of the path where to output files
        """

        if not self.make_bam and not self.make_sam:
            print ("\tNo bam/sam/bai files to be generated")
            return

        # Sort the reads if a dictionnary was given and self.sort is true
        if type(read_col) == dict:
            if self.sort: # flat sorted list
                print ("\tSorting reads...")
                read_list = self._sort_read(header, read_col)
                header.to_dict()['HD'] = {'VN': '1.5', 'SO' : 'coordinate'}
            else: # flat unsorted list
                read_list = [item for sublist in list(read_dict.values()) for item in sublist]
        # else just use the list
        else:
            read_list = read_col

        # Ouput bam if requested
        if self.make_bam:
            print ("\tCreate a bam file...")
            self.bam = "{}_{}.bam".format(outpath, ref_name)
            with pysam.Samfile(self.bam, "wb", header=header) as bamfile:
                for read in read_list:
                    bamfile.write(read)
             # Ouput bai if requested
            if self.make_index:
                print ("\tCreate a bam index file...")
                pysam.index(self.bam)
                self.bai = self.bam+".bai"

         # Ouput sam if requested
        if self.make_sam:
            print ("\tCreate a sam file...")
            self.sam = "{}_{}.sam".format(outpath, ref_name)
            with pysam.Samfile(self.sam, "wh", header=header) as samfile:
                for read in read_list:
                    samfile.write(read)

    def _sort_read (self, header, read_dict):
        """
        For coordinate sort, the major sort key is the RNAME field, with order defined
        by the order of @SQ lines in the header. The minor sort key is the POS field.
        """
        # Iterate over the sequence in the bam header to respect the order needed for a sorted bam
        read_list= []
        for seqline in header['SQ']:
            seqname = seqline['SN']

            # If the header entry correspond to to an entry in the seq dict sort the list associated
            # and extend
            if seqname in read_dict:
                read_dict[seqname].sort(key = lambda x: x.pos)
                read_list.extend(read_dict[seqname])

        return read_list

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class BamDecoy(object):
    """
    @class BamDecoy
    @brief Decoy class implementing no behaviour
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    def __init__ (self, *args, **kwargs):
        """
        Decoy init method
        """
        pass

    def __repr__(self):
        msg = "BAM DECOY\n"
        return msg

    def __str__(self):
        return "\n<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def make (self, *args, **kwargs):
        """
        Decoy make method
        """
        print ("\tNo bam/sam/bai file to be generated")
