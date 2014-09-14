#!/usr/bin/python
"""test_Read.py to test the Read class.

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
"""

from nose.tools import raises

from Read import Read
from MetageneError import MetageneError

##TODO: test set_sam_tag method
##TODO: test set_chromosome_sizes

cigar_string = {}
bad_cigar_string = {}
bitwise_flag = {}
bad_bitwise_flag = {}

good_input = {}
bad_input = {}
chromosome_conversion = {"1": "chr1", "2": "chr2"}


def setup():
    """Create fixtures"""
    # define cigar strings; value: ((args for build_positions), expected_result)
    cigar_string['full_match'] = ((1, "10M", "*"), [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    cigar_string['insertion'] = ((1, "5M4I5M", "*"), [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    cigar_string['deletion'] = ((1, "5M4D5M", "*"), [1, 2, 3, 4, 5, 10, 11, 12, 13, 14])
    cigar_string['gapped_match'] = ((1, "5M3N5M", "*"), [1, 2, 3, 4, 5, 9, 10, 11, 12, 13])
    cigar_string['softclipped_match'] = ((4, "3S5M", "*"), [4, 5, 6, 7, 8])
    cigar_string['hardclipped_match'] = ((4, "3H5M3H", "*"), [4, 5, 6, 7, 8])
    cigar_string['padded_match'] = ((1, "3P5M", "*"), [4, 5, 6, 7, 8])
    cigar_string['mismatch'] = ((1, "5=1X3=", "*"), [1, 2, 3, 4, 5, 6, 7, 8, 9])
    cigar_string['no_cigar_match'] = ((1, "*", "aaaaa"), [1, 2, 3, 4, 5])
    bad_cigar_string['unknown_length'] = ((1, "*", "*"), "raise MetageneError")
    bad_cigar_string['illegal_cigar'] = ((1, "5M4B", "*"), "raise MetageneError")
    bad_cigar_string['misordered_cigar'] = ((1, "M5N4M5", "*"), "raise MetageneError")

    # define bitwise flags; value: ((args for parse_sam_bitwise_flag), expected_result(count?, reverse_complemented?))
    bitwise_flag['unmapped'] = ((int("0b000000000100", 2),), (False, False))
    bitwise_flag['unmapped_withflags'] = ((int("0b100111011101", 2),), (False, True))
    bitwise_flag['plus_strand'] = ((int("0b000000000000", 2),), (True, False))
    bitwise_flag['minus_strand'] = ((int("0b000000010000", 2),), (True, True))
    bitwise_flag['multiple_segments'] = ((int("0b000000000001", 2),), (True, False))
    # try various default and user-changed boolean flags
    bitwise_flag['count_secondary_alignment'] = ((int("0b000100000000", 2),), (True, False))
    bitwise_flag['skip_secondary_alignment'] = (
        (int("0b000100000000", 2), False, False, False, True, False, False), (False, False))
    bitwise_flag['skip_failed_quality_control'] = ((int("0b001000000000", 2),), (False, False))
    bitwise_flag['count_failed_quality_control'] = (
        (int("0b001000000000", 2), True, True, False, True, False, False), (True, False))
    bitwise_flag['skip_PCR_optical_duplicate'] = ((int("0b010000000000", 2),), (False, False))
    bitwise_flag['count_PCR_optical_duplicate'] = (
        (int("0b010000000000", 2), True, False, True, True, False, False), (True, False))
    bitwise_flag['count_supplementary_alignment'] = ((int("0b100000000000", 2),), (True, False))
    bitwise_flag['skip_supplementary_alignment'] = (
        (int("0b100000000000", 2), True, False, False, False, False, False), (False, False))
    bitwise_flag['count_only_start_success'] = (
        (int("0b000001000001", 2), True, False, False, True, True, False), (True, False))
    bitwise_flag['count_only_start_fail'] = (
        (int("0b000000000001", 2), True, False, False, True, True, False), (False, False))
    bitwise_flag['count_only_end_success'] = (
        (int("0b000010000001", 2), True, False, False, True, False, True), (True, False))
    bitwise_flag['count_only_end_fail'] = (
        (int("0b000000000001", 2), True, False, False, True, False, True), (False, False))
    bad_bitwise_flag['count_only_both'] = (
        (int("0b000011000001", 2), True, False, False, True, True, True), ("Raise MetageneError",))

    # define good and bad samline inputs
    good_input['no_tags'] = (0, "chr1", 200, "10M", 10, 1, 1, "+")
    good_input['plus_strand_match'] = (0, "chr1", 200, "10M", 10, 2, 4, "+")
    good_input['minus_strand_match'] = (16, "chr1", 200, "10M", 10, 2, 4, "-")
    good_input['no_match'] = (4, "*", 0, "*", 10, 1, 1, ".")

    sample = ["NA:i:4\tNH:i:4",
              "NA:i:4\tNH:i:4",
              "NA:i:4\tNH:i:4",
              "NA:i:4\tNH:i:4",
              "NA:i:4\tNH:i:4",
              "NA:i:4\tNH:i:4",
              "NA:i:4\tNH:i:4",
              "NA:i:4\tNH:i:4",
              "NA:i:4\tNH:i:4",
              "NA:i:4\tNH:i:4"]

    Read.process_set_sam_tag(sample, count_tag=True, tag_regex='NA:i:(\d+)')
    Read.process_set_sam_tag(sample, count_tag=True, tag_regex='NH:i:(\d+)')


def test_build_positions():
    for test in cigar_string:
        yield (check_build_positions, test, cigar_string[test])


def check_build_positions(test, (values, expected)):
    position_array = Read.build_positions(*values)
    test_description = "\nTest:    \t{}\n".format(test)
    test_description += "Expected:\t{}\n".format(expected)
    test_description += "Position:\t{}\n".format(position_array)
    assert position_array == expected, "{}Error:   \tDid not create the expected position array.".format(
        test_description)


def test_catch_bad_cigar_input():
    for test in bad_cigar_string:
        yield (check_catch_bad_cigar_input, test, bad_cigar_string[test])


@raises(MetageneError)
def check_catch_bad_cigar_input(test, (values, expected)):
    print Read.build_positions(*values)


def test_parse_sam_bitwise_flag():
    for test in bitwise_flag:
        yield (check_parse_sam_bitwise_flag, test, bitwise_flag[test])


def check_parse_sam_bitwise_flag(test, (values, expected)):
    bitwise_result = Read.parse_sam_bitwise_flag(*values)
    test_description = "\nTest:    \t{}\n".format(test)
    test_description += "Expected:\t{}\n".format(expected)
    test_description += "Position:\t{}\n".format(bitwise_result)
    assert bitwise_result == expected, "{}Error:   \tDid not parse bitwise flag as expected.".format(test_description)


def test_catch_bad_bitwise_input():
    for test in bad_bitwise_flag:
        yield (check_catch_bad_bitwise_input, test, bad_bitwise_flag[test])


@raises(MetageneError)
def check_catch_bad_bitwise_input(test, (values, expected)):
    print Read.parse_sam_bitwise_flag(*values)


def build_samline(bitcode, chromosome, start, cigar, length, abundance, mappings):
    """Return a SAM format line"""
    string = "a" * length
    return "read\t{}\t{}\t{}\t255\t{}\t*\t0\t0\t{}\t{}\tNH:i:{}\tNA:i:{}".format(
        bitcode,
        chromosome,
        start,
        cigar,
        string,
        string,
        mappings,
        abundance)


def test_create_read():
    for test in good_input:
        yield (check_create_read, test, good_input[test])


def check_create_read(test, values):
    # create expected result
    if int(values[0]) == 4:
        expected = "Non-aligning read"
    else:
        start = int(values[2])
        end = int(values[2]) + int(values[4]) - 1
        if values[7] == "-":
            start = end
            end = int(values[2])
        expected = "Read at {0}:{1}-{2} on {3} strand; counts for {4:2.3f}:".format(
            values[1],  # chromosome
            start,
            end,
            values[7],  # strand
            float(values[5]) / float(values[6]))  # abundance / mappings
    # build input to test
    samline = build_samline(*values[0:-1])  # exclude final value
    (created, read) = Read.create_from_sam(samline, chromosome_conversion.values(), count_method='all')
    output = str(read).split("\t")[0]
    # create description in case test fails
    test_description = "\nTest:    \t{}\n".format(test)
    test_description += "Abundance:\t{}\n".format(Read.has_sam_tag["NA"])
    test_description += "Mappings:\t{}\n".format(Read.has_sam_tag["NH"])
    test_description += "Sam Line:\t{}\n".format(samline)
    test_description += "Expected:\t{}\n".format(expected)
    test_description += "Position:\t{}\n".format(output)
    assert output == expected, "{}Error:   \tDid not create expected read.".format(test_description)


def test_catch_bad_input():
    for test in bad_input:
        yield (check_catch_bad_input, test, bad_input[test])


@raises(MetageneError)
def check_catch_bad_input(test, samline):
    print Read(sam_line)
            
                 


            

