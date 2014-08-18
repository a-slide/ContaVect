#~~~~~~~GLOBAL IMPORTS~~~~~~~#

# Local Package import
from pyDNA.Utilities import DNA_reverse_comp

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class AdapterTrimmer(object):
    """
    Search matches of a list of adapter in a reference sequence as a Biopython SeqReccord object
    If matches are found return the longer interval of the reference without matches is returned
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__ (self, Aligner, adapters, min_read_len=0.6, min_match_len=0.8, min_match_score=1.4, find_rc = False):
        """
        @param Aligner Wrapper object for pairwise alignement. The aligner needs to accept a
        query and a reference DNA string and return a match object with at least 2 fields
        ref_begin and m.ref_end indicating the starrt and end position along the reference
        @param adapters List of DNA base string corresponding to adapters to be trimmed
        @param min_read_len Fraction of read lenth = minimal size of fragment after trimming
        @param min_match_len Minimal fraction of adapter len that needs to be aligned on the target
        @param min_match_score Minimal score per base for the alignment of adapter and read
        @param find_rc If true will also search for the reverse complementary sequence of the adapter
        """
        #Store object variables
        self.min_read_len = min_read_len
        self.min_match_len = min_match_len
        self.min_match_score = min_match_score
        self.Aligner = Aligner

        # Import a list of adapters and add the reverse complements of adapters to the list
        self.adapter_list = []
        for id, seq in enumerate (adapters, start=1):
            self.adapter_list.append({
                "id": str(id), "seq": seq, "count": 0,
                "min_score": int(self.min_match_score * len(seq)),
                "min_len": int(self.min_match_len * len(seq))})
            if find_rc:
                self.adapter_list.append({
                    "id": str(id)+'#',"seq": DNA_reverse_comp(seq), "count": 0,
                    "min_score": int(self.min_match_score * len(seq)),
                    "min_len": int(self.min_match_len * len(seq))})

        # Initialize generic counters
        self.seq_untrimmed = 0
        self.seq_trimmed = 0
        self.base_trimmed = 0
        self.len_pass = 0
        self.len_fail = 0
        self.run = False

    def __repr__(self):
        msg = "ADAPTER TRIMMER\n"
        msg += "  List of adapters imported for trimming\n"

        for a in self.adapter_list:
            msg += "  id: {}\tSequence: {}\tMin score: {}\tMin len: {}\n".format(
                a['id'], a['seq'], a['min_score'], a['min_len'])

        if self.run:
            msg += "  Sequences untrimmed : {}\n".format(self.seq_untrimmed)
            msg += "  Sequences trimmed : {}\n".format(self.seq_trimmed)
            msg += "  DNA base trimmed : {}\n".format(self.base_trimmed)
            msg += "  Fail len filtering: {}\n".format(self.len_fail)
            msg += "  Pass len filtering : {}\n".format(self.len_pass)
            msg += "  Total pass : {}\n\n".format(self.len_pass+self.seq_untrimmed)
        return msg

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def get(self, key):
        return self.__dict__[key]

    def set(self, key, value):
        self.__dict__[key] = value

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def trimmer (self, record):
        """
        Trim reference sequence by finding matches of adapters with the Aligner object and the
        longuest interval without adapters with _longer_interval
        @param record A BioPython seqRecord object containing the subject reference sequence to be
        trimmed
        @return If no matches were found the original sequence. If matches were found a trimmed
        record if the fraction of lenght remaining after trimming is above min_read_len and
        elsewhere nothing
        """
        self.run = True
        match_list = []
        len_rec = len(record)

        # Set a new reference sequence into the Aligner
        self.Aligner.set_ref(str(record.seq))

        for a in self.adapter_list:
            # Find match of the adapter along the current read
            match = self.Aligner.align(a['seq'], a["min_score"], a["min_len"])

            # if a match was found = increment the counter and append the match to the match list
            if match:
                #print ("Adapter found : {}\tScore : {}\tCigar : {}\tMatchLen : {}".format(
                #a.id, match.score, match.cigar_string, match.query_end-match.query_begin))
                a["count"] += 1
                match_list.append(match)

        # In case no match were found, the sequence doesn't need to be modify
        if not match_list:
            self.seq_untrimmed += 1
            return record

        # Else find the longer interval without adaptor matches
        start, end = self._longer_interval (match_list, len_rec)

        # Update counters
        self.seq_trimmed += 1
        self.base_trimmed += (len_rec-(end-start))

        # Return a slice of the reccord corresponding to the longer interval
        if end-start >= int(self.min_read_len*len_rec):
            self.len_pass +=1
            return record[start:end]
        # Or None if smaller than min_size
        else:
            self.len_fail +=1
            return None

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _longer_interval(self, match_list, len_seq):
        """
        Find the first larger interval that do not overlapp any matches in match list.
        This strategy allow to use an unsorted list of match but will be highly memory consummming
        for large reference.
        """
        # Initialize a list of boolean to False of the same size as the read
        coverage = [False for i in range(len_seq)]

        # Flag positions overlapped by a read by changing the boolean to True
        for m in match_list:
            for i in range (m.ref_begin, m.ref_end):
                coverage[i] = True

        # Read through the list to find the longer inteval between True flags
        start_max = end_max = inter_max = start = inter = 0
        for i in range(len_seq):
            if coverage[i]:
                start = i+1
                inter = 0
            else:
                inter += 1
                if inter > inter_max:
                    inter_max = inter
                    start_max = start
                    end_max = i

        #print ("Longer interval = {} [{}:{}]".format(inter_max, start_max+1, end_max-1))
        return start_max, end_max
