#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

"""
@package    FastqWriter
@brief      Helper class for Sample
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

from io import open as iopen
from zlib import decompressobj as zlib_decompressobj
from zlib import MAX_WBITS as zlib_MAX_WBITS

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Ungzip(object):
    """
    This class decompress gziped files much more quickly than the gzip package from python standard
    library. It was adapted from the code corner blog m.wolf@code-corner.de
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    def __init__(self, buffer_size=1024*1024*8):
        self.dobj = zlib_decompressobj(16+zlib_MAX_WBITS) #16+zlib.MAX_WBITS -> zlib can decompress gzip
        self.decomp = []
        self.lines = []
        self.buffer_size = buffer_size

    def open(self, filename):
        self.fhwnd = iopen(filename, "rb")
        self.eof = False

    def close(self):
        self.fhwnd.close()
        self.dobj.flush()
        self.decomp = []

    def decompress(self):
        raw = self.fhwnd.read(self.buffer_size)
        if not raw:
            self.eof = True
            self.decomp.insert(0, self.dobj.flush())

        else:
            self.decomp.insert(0, self.dobj.decompress(raw))

    def readline(self):

        out_str = []

        while True:
            if len(self.lines) > 0:
                return self.lines.pop() + "\n"

            elif len(self.decomp) > 0:
                out = self.decomp.pop()
                arr = out.split("\n")

                if len(arr) == 1:
                    out_str.append(arr[0])

                else:
                    self.decomp.append(arr.pop())
                    arr.reverse()
                    out_str.append(arr.pop())
                    self.lines.extend(arr)

                    out_str.append("\n")
                    return "".join(out_str)

            else:
                if self.eof: break
                self.decompress()

        if len(out_str) > 0:
            return "".join(out_str)
