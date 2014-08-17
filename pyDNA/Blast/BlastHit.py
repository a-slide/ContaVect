#~~~~~~~GLOBAL IMPORTS~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class BlastHit(object):
    """
    @class  BlastHit
    @brief  Object oriented class containing informations of one blast hit
    The following instance field are accessible :
    * q_id : Query sequence name
    * s_id : Subject sequence name
    * identity : % of identity in the hit
    * length : length of the hit
    * mis : Number of mismatch in the hit
    * gap : Number of gap in the hit
    * q_orient : Orientation of the query along the hit
    * q_start : Hit start position of the query
    * q_end : Hit end position of the query
    * s_orient : Orientation of the subject along the hit
    * s_start : Hit start position of the subject
    * s_end : Hit end position of the subject
    * evalue : E value of the alignement
    * bscore : Bit score of the alignement
    A class list is used to track all instances generated.
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
    def count_total (self):
        """
        @return Overall number of BlastHit object in Instance list
        """
        return (len(self.Instances))

    @ classmethod
    def stat_per_ref (self):
        """
        @return Number of BlastHit object in Instance list sorted by reference subject sequence
        """
        d = {}
        for hit in self.Instances:
            if hit.s_id in d:
                d[hit.s_id][0] += 1
                d[hit.s_id][1] += hit.length
            else:
                d[hit.s_id] = [1, hit.length]
        return d

    @ classmethod
    def get (self):
        """
        @return The list of all BlastHit object generated
        """
        return self.Instances

    @ classmethod
    def get_ref (self, ref):
        """
        @param ref Name of a reference sequence in the subject database
        @return The list of all BlastHit object generated for this reference
        """
        return [hit for hit in self.Instances if hit.s_id == "ref"]

    @ classmethod
    def reset_list (self):
        """
        Reset the instance tracking list (Usefull after
        """
        self.Instances = []
        self.id_count = 0

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self, q_id, s_id, identity, length, mis, gap, q_start, q_end, s_start, s_end, evalue, bscore):
        """
        Create a BlastHit object which is automatically added to the class tracking instance list
        The object with the following parameters are required for object initialisation
        @param  q_id    Query sequence name
        @param  s_id    Subject sequence name
        @param  identity    % of identity in the hit
        @param  length  length of the hit
        @param  mis Number of mismatch in the hit
        @param  gap Number of gap in the hit
        @param  q_start Hit start position of the query
        @param  q_end   Hit end position of the query
        @param  s_start Hit start position of the subject
        @param  s_end   Hit end position of the subject
        @param  evalue  E value of the alignement
        @param  bscore Bit score of the alignement
        """
        
        self.id = self.next_id()
        self.q_id = q_id
        self.s_id = s_id
        self.identity = float(identity)
        self.length = int(length)
        self.mis = int(mis)
        self.gap = int(gap)
        self.evalue = float(evalue)
        self.bscore = float(bscore)

        # Autoadapt start and end so that start is always smaller than end
        self.q_start = int(q_start) if int(q_start) < int(q_end) else int(q_end)
        self.q_end = int(q_end) if int(q_start) < int(q_end) else int(q_start)
        self.s_start = int(s_start) if int(s_start) < int(s_end) else int(s_end)
        self.s_end = int(s_end) if int(s_start) < int(s_end) else int(s_start)

        # Orientation of the query and subject along the hit. True if positive
        self.q_orient = int(q_start) < int(q_end)
        self.s_orient = int(s_start) < int(s_end)

        # Add the instance to the class instance tracking list
        self.Instances.append(self)

    def __repr__(self):
        msg = "HIT {}".format(self.id)
        msg += "\tQuery\t{}:{}-{}({})\n".format(self.q_id, self.q_start, self.q_end, "+" if self.q_orient else "-")
        msg += "\tSubject\t{}:{}-{}({})\n".format(self.s_id, self.s_start, self.s_end, "+" if self.q_orient else "-")
        msg += "\tLenght : {}\tIdentity : {}%\tEvalue : {}\tBit score : {}\n".format(self.length, self.identity, self.evalue, self.bscore)
        return (msg)

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)
