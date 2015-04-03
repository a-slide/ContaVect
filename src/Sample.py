# -*- coding: utf-8 -*-

"""
@package    Sekator
@brief      Helper class for Sekator to represent Samples
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

# Local Package import
from pyDNA.FileUtils import is_gziped, is_readable_file

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Sample(object):
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##################################################################################################### THE CLASS MIGHT STORE A REPORT???

    #~~~~~~~CLASS FIELDS~~~~~~~#

    SAMPLE_NAMES = []

    #~~~~~~~CLASS METHODS~~~~~~~#

    @ classmethod
    def ADD_TO_SAMPLE_NAMES(self, name):
        self.SAMPLE_NAMES.append(name)

    #~~~~~~~FUNDAMENTAL METHODS~~~~~~~#

    def __init__ (self, name, R1_path, R2_path):

        # test fasta readability
        is_readable_file(R1_path)
        is_readable_file(R2_path)

        # Create self variables
        self.name = name
        self.R1_path = R1_path
        self.R2_path = R2_path

        assert self.name not in self.SAMPLE_NAMES, "Sample name <{}> is duplicated".format(self.name)

        self.ADD_TO_SAMPLE_NAMES(self.name)

    # Fundamental class methods str and repr
    def __str__(self):
        msg = "SAMPLE CLASS\n\tParameters list\n"
        # list all values in object dict in alphabetical order
        keylist = [key for key in self.__dict__.keys()]
        keylist.sort()
        for key in keylist:
            msg+="\t{}\t{}\n".format(key, self.__dict__[key])
        return (msg)

    def __repr__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)
