#!/usr/bin/python
'''Metagene class for metagene_counts.py

Requires:
    python 2 (https://www.python.org/downloads/)

Joy-El R.B. Talbot Copyright (c) 2014

The MIT License (MIT)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
'''

import re, datetime, subprocess, math
from MetageneError import MetageneError

class Metagene(object):
    '''A Metagene is an object representing an interval of interest (length >= 1) 
    with padding on either side (lengths >= 0).  
    
    Attributes:
        feature_interval : interval of interest (int >= 1)
        padding          : dictionary of padding values (int >= 0)
          'Upstream'     : padding upstream (to left) of feature_interval
          'Downstream'   : padding downstream (to right) of feature_interval
        length           : length of entire metagene object
        '''
##TODO: add functionality for the start and end of feature_interval to be constant and internal interval to vary
##TODO: add functionality for negative paddings!!

    # restrict attributes for each instance
    __slots__ = ['feature_interval','padding', 'length']
    
    
    def __init__(self, interval, padding_upstream, padding_downstream):
        '''Initiate lengths from the interval and padding.'''
        
        # assign interval; must be INT >= 1
        try:
            if interval == int(interval) and int(interval) >= 1:
                self.feature_interval = int(interval)
            else:
                raise MetageneError(interval, "Interval must be an interger greater than zero")
        except ValueError as err:
            raise MetageneError(interval, "Interval must be an integer greater than zero")
        
        # assign padding; must be INT >= 0  
        self.padding = {'Upstream':0, 'Downstream':0} # set defaults
        for pad in [(padding_upstream, "Upstream"), (padding_downstream, "Downstream")]:
            try:
                if pad[0] == int(pad[0]) and int(pad[0]) >= 0:
                    self.padding[pad[1]] = int(pad[0])
                else:
                    raise MetageneError(pad[0], "Padding values must be positive integers")
            except ValueError as err:
                raise MetageneError(pad[0], "Padding values must be positive integers")   
             
        self.length = self.feature_interval + self.padding['Upstream'] + self.padding['Downstream']

    # end __init__ function
    
    
    def __str__(self):
        return "Upstream:{} -- Interval:{} -- Downstream:{}\tLength:{}".format(self.padding['Upstream'], self.feature_interval, self.padding['Downstream'], self.length)

    def print_full(self, pretty=False):
        '''Represents the metagene as positions relative to the interval start as 0'''
        
        output = ""
        if pretty: # add metagene schematic and position numbers (relative to feature start as zero)
            output += "{0:15s}\t\t".format('Metagene')
            for i in range(self.padding['Upstream']):
                output += "---up-"
            for i in range(self.feature_interval):
                output += "--int-"
            for i in range(self.padding['Downstream']):
                output += "-down-"
            output += "\n"
        
            output += "{0:15s}:\t".format('Position')
            for i in range(self.padding['Upstream'], 0, -1):
                output += "{0:5d},".format(0-i)
            for i in range(self.feature_interval):
                output += "{0:5d},".format(i)
            for i in range(self.padding['Downstream']):
                output+= "{0:5d},".format(i + self.feature_interval)
            output = output[:-1] + "\n"
        else:
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
