#!/usr/bin/python
'''The second step of metagene_analysis, metagene_binning.py compiles read
abundance over genomic features to create the input for metagene_windows.py.
Please see README for full details and examples.

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

import sys, re, datetime, subprocess, math
import argparse		# to parse the command line arguments
from MetageneError import MetageneError
from Metagene import Metagene

PROGRAM = "metagene_binning.py"
VERSION = "0.1.2"
UPDATED = "140403 JRBT"

def get_arguments():
    '''Collect and parse information from the user's command line arguments.'''
    
    date = datetime.datetime.now().strftime('%y%m%d-%H%M%S')
    
    parser = argparse.ArgumentParser(description=
    '''The second step of metagene_analysis, metagene_binning.py compiles read
abundance over genomic features to create the input for metagene_windows.py.
Please see README for full details and examples.

Requires:
    python 2 (https://www.python.org/downloads/)
    ''')

    parser.add_argument("-v","--version",
                        action = 'version',
                        version = "{} {}\tUpdated {}".format(PROGRAM, VERSION, UPDATED))
    parser.add_argument("-i","--input",
                        help = "Counts file from metagene_counts.py",
                        metavar = 'COUNTS_FILE',
                        required = True,
                        action = "append")
    parser.add_argument("-o", "--output_prefix",
                        help = "Prefix for output files",
                        default = "binned")
    
    parser.add_argument("--window_size",
                        help = "Size of windows for binning; default = 10",
                        type = int,
                        default = 10)                  
    parser.add_argument("--step_size",
                        help = "Step size for binning; default = 10",
                        type = int,
                        default = 10)
    parser.add_argument("--separate_groups",
                        help = "Output orientation and gap types to separate files",
                        action = "store_true")
    
    return parser.parse_args()


def build_output_filenames(infile, prefix, window_size, step_size, separate_groups):
    
    outfiles = {}
    
    # determine orientation:gap combos
    groups = []
    with open(infile) as inf:
        header = inf.readline()
        group = inf.readline().strip().split(",")[1]
        while group not in groups:
            groups.append(group)
            group = inf.readline().strip().split(",")[1]
    
    basefile = "{}.{}bpX{}bp".format(prefix, window_size, step_size)
    if separate_groups:
        for g in groups:
            outfiles[g] = "{}.{}.csv".format(basefile, "_".join(g.split(":")))
    else:
        for g in groups:
            outfiles[g] = "{}.all.csv".format(basefile)
    
    return outfiles

def read_chunk(file_obj,chunk_size):
    '''Read in file by chunk_size chunks returning one line at a time.'''
    # get first chunk
    chunk = file_obj.read(chunk_size)
    # continue looping until a chunk is just EOF (empty line)
    while chunk:
        chunk_list = chunk.split("\n")
        # yield all but last, potentially incomplete line
        for c in chunk_list[:-1]:
            yield c
        # add incomplete line to beginning of next chunk read
        chunk = chunk_list[-1] + file_obj.read(chunk_size)

def metagene_binning():
    '''Main program for metagene_binning.py'''
   
    arguments = get_arguments()

    for infile in arguments.input:
        print "Processing file:\t{}".format(infile)
        
        # returns a dict of file names with keys of orientation:gap_counting
        output_files = build_output_filenames(infile, arguments.output_prefix, arguments.window_size, arguments.step_size, arguments.separate_groups)
        
        # Write header to each file
        for output in output_files.values():
            with open(output, 'w') as outf:
                outf.write("Gene,Orientation,Gapped,Window,Inclusive_Start,Inclusive_End,Abundance\n")

        with open(infile, 'r') as inf:
            header = inf.readline().strip().split(",")
            positions = header[2:] # positions relative to gene start
            
            for counts_line in read_chunk(inf, 1024):
                counts_parts = counts_line.strip().split(",")
                counts = counts_parts[2:]
                length = len(counts)
                (orientation, gap) = counts_parts[1].split(":")
                output = "{},{},{}".format(counts_parts[0], orientation, gap)
                
                window = 0
                exclusive_end = arguments.window_size
                
                while exclusive_end <= length:
                    inclusive_start = exclusive_end - arguments.window_size
                    
                    coverage = 0.0

                    for i in range(inclusive_start, exclusive_end):
                        coverage += float(counts[i])
                    
                    with open(output_files[counts_parts[1]], 'a') as outf:
                        outf.write("{},{},{},{},{}\n".format(output,window,positions[inclusive_start], positions[exclusive_end - 1], coverage))

                    window += 1
                    exclusive_end += arguments.step_size

if __name__ == "__main__":
    metagene_binning()

