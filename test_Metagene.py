#!/usr/bin/python
'''test_Metagene.py to test the Metagene class

Requires:
    python 2 (https://www.python.org/downloads/)
    nose 1.3 (https://nose.readthedocs.org/en/latest/)

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
import re

from nose.tools import raises

from Metagene import Metagene
from MetageneError import MetageneError

good_input = {}
bad_input = {}

def setup():
    """Create fixtures with expected results"""
    
    good_input['even_padding'] = (9,3,3)
    good_input['uneven_padding'] = (11,2,5)
    good_input['zero_padding'] = (10,0,0)
    
    bad_input['float_interval'] = (10.2, 4, 4)
    bad_input['string_interval'] = ("ten", 4, 4)
    bad_input['zero_interval'] = (0, 3, 3)
    bad_input['negative_padding'] = (10, -3, 2)
    bad_input['float_padding'] = (10, 4, 4.2)
    bad_input['string_padding'] = (10, 4, "four")
    
def test_create_metagene():
    for test in good_input:
        yield (check_create_metagene, test, good_input[test])
              
def check_create_metagene(test, values):
    metagene = Metagene(*values)
    length = sum(values)
    (interval, upstream, downstream) = values
    expected = "Upstream:{} -- Interval:{} -- Downstream:{}\tLength:{}".format(upstream, interval, downstream, length)
    test_description =  "\nTest:    \t{}\n".format(test)
    test_description += "Expected:\t{}\n".format(expected)
    test_description += "Metagene:\t{}\n".format(metagene)
    assert str(metagene) == expected, "{}Error:   \tMetagene does not match expected.".format(test_description)

def test_catch_bad_input():
    for test in bad_input:
        yield (check_catch_bad_input, test, bad_input[test])
           
@raises(MetageneError)
def check_catch_bad_input(test, values):
    print Metagene(*values)
    
def test_print_metagene():
    for check in (check_print_metagene_plain, check_print_metagene_pretty):
        for test in good_input:
            yield (check, test, good_input[test])

def check_print_metagene_plain(test, values):
    expected = (-values[1], values[0] + values[2] - 1)
    metagene = Metagene(*values)
    plain_print = metagene.print_full().strip()
    plain_print_parts = plain_print.split(",")
    new_range = (int(plain_print_parts[2]), int(plain_print_parts[-1]))
    print "\n\tOutput:\n\t{}".format(plain_print)
    test_description =  "\nTest:    \t{}\n".format(test)
    test_description += "Expected:\t{}\n".format(expected)
    test_description += "Range:   \t{}\n".format(new_range)
    test_description += "Output:  \t{}\n".format(plain_print)
    assert new_range == expected, "{}Error:   \tPrinted metagene does not match expected.".format(test_description)
    
def check_print_metagene_pretty(test, values):
    metagene = Metagene(*values)    
    pretty_print = metagene.print_full(pretty=True).strip()
    new_values = (len(re.findall('int', pretty_print)),
                  len(re.findall('up', pretty_print)),
                  len(re.findall('down', pretty_print)))
    print "\n\tOutput:\n\t{}".format("\n\t".join(pretty_print.split("\n")))
    test_description =  "\nTest:    \t{}\n".format(test)
    test_description += "Expected:\t{}\n".format(values)
    test_description += "Range:   \t{}\n".format(new_values)
    test_description += "Output:  \t{}\n".format(pretty_print)
    assert new_values == values, "{}Error:   \tPrinted metagene does not match expected.".format(test_description)



