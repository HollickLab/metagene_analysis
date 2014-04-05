#!/usr/bin/python
"""Metagene class for metagene_counts.py

Requires:
    python 2 (https://www.python.org/downloads/)

Joy-El R.B. Talbot Copyright (c) 2014

The MIT License (MIT)

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files 
(the "Software"), to deal in the Software without restriction, 
including without limitation the rights to use, copy, modify, merge, 
publish, distribute, sublicense, and/or sell copies of the Software, 
and to permit persons to whom the Software is furnished to do so, 
subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

import re
import datetime
import subprocess
import math

from metageneMethods import confirm_integer
from MetageneError import MetageneError

class Metagene(object):
    """Create a metagene composed of a padded interval of interest.
    
    Attributes:
        feature_interval -- interval 
            value: non-zero positive integer
        padding -- dictionary of padding values
            keys: 
                'Upstream' -- left-side interval padding 
                    value: positive integer
                'Downstream' -- right-side interval padding 
                    value: positive integer
        length -- paddings + interval 
             value: non-zero positive integer
    
    Methods:
        print_full -- Return metagene positions relative to interval start. 
                      (instead of default summary)
    """
##TODO: add functionality for the start and end of feature_interval to be constant and internal interval to vary
##TODO: add functionality for negative paddings!!

    # restrict attributes for each instance
    __slots__ = ['feature_interval','padding', 'length']
    
    
    def __init__(self, interval=1, padding_upstream=0, padding_downstream=0):
        """Return metagene instance defined by interval and padding sizes.
        
        Keyword arguments:
        interval -- length of interval (default 1)
        padding_upstream -- length of upstream padding (default 0)
        padding_downstream -- length of downstream padding (default 0)
        """
        if confirm_integer(interval, "Interval", minimum=1):
            self.feature_interval = int(interval)
        self.padding = {'Upstream':None, 'Downstream':None}
        if confirm_integer(padding_upstream, "Upstream padding", minimum=0):
            self.padding['Upstream'] = int(padding_upstream)
        if confirm_integer(padding_downstream, "Downstream padding", minimum=0):
            self.padding['Downstream'] = int(padding_downstream)
        self.length = (self.padding['Upstream'] +
                       self.feature_interval + 
                       self.padding['Downstream'])

    # end __init__ function
    
    
    def __str__(self):
        return "Upstream:{} -- Interval:{} -- Downstream:{}\tLength:{}".format(self.padding['Upstream'], self.feature_interval, self.padding['Downstream'], self.length)

    def print_full(self, pretty=False):
        """Return metagene positions relative to interval start as 0.
        
        Keyword arguments:
        pretty -- return human readable version (default False)
        """
        
        output = ""
        # add metagene schematic and position numbers 
        # (relative to feature start as zero)
        if pretty: 
            # ---up---int--down- labeling
            output += "{0:15s}\t\t".format('Metagene')
            for i in range(self.padding['Upstream']):
                output += "---up-"
            for i in range(self.feature_interval):
                output += "--int-"
            for i in range(self.padding['Downstream']):
                output += "-down-"
            output += "\n"
            
            # ---up---int--down-  
            #    -1     0     1   relative position labeling
            output += "{0:15s}:\t".format('Position')
            for i in range(self.padding['Upstream'], 0, -1):
                output += "{0:5d},".format(0-i)
            for i in range(self.feature_interval):
                output += "{0:5d},".format(i)
            for i in range(self.padding['Downstream']):
                output+= "{0:5d},".format(i + self.feature_interval)
            output = output[:-1] + "\n"
        
        else:
            # comma-delimited position output
            # suitable header for metagene_bin.py input files
            output += "{},{}".format('Feature','Orientation:Gap')
            for i in range(self.padding['Upstream'], 0, -1):
                output += ",{}".format(0-i)
            for i in range(self.feature_interval):
                output += ",{}".format(i)
            for i in range(self.padding['Downstream']):
                output+= ",{}".format(i + self.feature_interval)
            output += "\n"
        return output
# end Metagene class
