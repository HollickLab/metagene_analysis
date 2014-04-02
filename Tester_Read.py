#!/usr/bin/python
'''Tester_Feature.py to test the Feature.py class.

Requires:
    python 2 (https://www.python.org/downloads/)

Based on Perl code by Karl F. Erhard, Jr Copyright (c) 2011
Extended and modified to Python by Joy-El R.B. Talbot Copyright (c) 2014

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

from Read import Read 
from MetageneError import MetageneError


def test():
    '''Tests of Read class'''
        
    print "\n**** Testing the Read class ****\n"
        
    chromosome_converter = {"1":"chr1", "2":"chr2"}
    print "\t  with chromosome conversions:\t{}\n".format(chromosome_converter)   
        
    try:
        samline1 = "read_1	24	chr1	200	255	3S2M4N3M2X3M	*	0	0	bbbbbbbbbbbbb	bbbbbbbbbbbbb	XA:i:0	MD:Z:40	NH:i:50  NA:i:10"
        print "\twith SAM line:\t{}".format(samline1)
        read1 = Read.create_from_sam(samline1, chromosome_converter)
    except MetageneError as err:
        print "**FAILED**\tCreate Read object from SAM line ?\t\t{}".format(err)
    else:
        print "PASSED\tCreate Read object from SAM line ?"
        print "\t  {}".format(read1[1])

            
    try:
        samline1 = "read_1	0	chr1	200	255	*	*	0	0	bbbbbbbbbbbbb	bbbbbbbbbbbbb	XA:i:0	MD:Z:40	NH:i:50  NA:i:10"
        print "\twith SAM line:\t{}".format(samline1)
        read1 = Read.create_from_sam(samline1, chromosome_converter)
    except MetageneError as err:
        print "**FAILED**\tCreate Read object from SAM line ?\t\t{}".format(err)
    else:
        print "PASSED\tCreate Read object from SAM line ?"
        print "\t  {}".format(read1[1])

                
##TODO finish complete testing of Feature class
    print "\n##TODO finish complete testing of Read class\n"
                
    print "\n**** End of Testing the Read class ****\n"  
# End of test function

if __name__=="__main__":
    test()
