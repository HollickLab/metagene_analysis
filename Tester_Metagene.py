#!/usr/bin/python
'''Tester_Metagene.py to test the Metagene class

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
from Metagene import Metagene
from MetageneError import MetageneError

def test():
    '''Test error handling of Metagene class'''
        
    print "\n**** Testing the Metagene class ****\n"
    # can a metagene be created
    print "\tCommand:\tmetagene = Metagene(11,4,4)"
    metagene = Metagene(11,4,4)
    if str(metagene) == "Upstream:4 -- Interval:11 -- Downstream:4\tLength:19":
        print "PASSED\tCreated a valid metagene ?\t\t{}".format(metagene)
    else:
        print "**FAILED**\tCreated a valid metagene ?"
        print "\tExpected Metagene:\tUpstream:4 -- Interval:11 -- Downstream:4\tLength:19"
        print "\tCreated  Metagene:\t{}".format(metagene)
            
    # create a metagene with different up and downstream padding
    print "\tCommand:\tmetagene = Metagene(10,2,4) "
    metagene = Metagene(10,2,4) 
    if str(metagene) == "Upstream:2 -- Interval:10 -- Downstream:4\tLength:16":
        print "PASSED\tCreated a valid metagene with differing padding ?\t\t{}".format(metagene)
    else:
        print "**FAILED**\tCreated a valid metagene with differing padding ?"
        print "\tExpected Metagene:\tUpstream:2 -- Interval:10 -- Downstream:4\tLength:16"
        print "\tCreated  Metagene:\t{}".format(metagene)
       
    # catch errors from non-desired inputs
    tests = [ ((10.2, 4, 4), "Caught non-integer float interval ?"),
              (("ten", 4, 4), "Caught non-integer string interval ?"),
              ((0, 3, 3), "Caught interval of zero error ?"),
              ((10, -3, 2), "Caught negative padding error ?"),
              ((10, 4.3, 4), "Caught non-integer float padding ?"),
              ((10, 4, "four"), "Caught non-integer string padding ?") ]
        
    for test in tests:
        try: 
            print "\tCommand:\tmetagene = Metagene({})".format(test[0])
            metagene = Metagene(*test[0]) # the * unpacks the list
        except MetageneError as err:
            print "PASSED\t{}\t\t{}".format(test[1], err)
        else:
            print "**FAILED**\t{}".format(test[1])
        
    print "\n**** End of Testing the Metagene class ****\n"
# end of Metagene.test staticmethod 

if __name__ == "__main__":
    test()
