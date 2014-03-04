#!/usr/bin/python
'''metagene_windows.py merges alignments with genomic features to create a 
"metagene" style summary of coverage over all of the genomic features. Please 
see README for full details and examples.

Based on Perl code by Karl F. Erhard, Jr Copyright (c) 2011
Modified to Python by Joy-El R.B. Talbot Copyright (c) 2014

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

import sys, re, datetime, subprocess
import argparse		# to parse the command line arguments

PROGRAM = "metagene_windows.py"
VERSION = "0.1.2"
UPDATED = "140303 JRBT"

def get_arguments(program, version, update):
    '''Collect and parse information from the user's command line arguments.'''
    
    date = datetime.datetime.now().strftime('%y%m%d-%H%M%S')
    
    parser = argparse.ArgumentParser(description=
    '''metagene_windows.py merges alignments with genomic features to create a 
"metagene" style summary of coverage over all of the genomic features. Please 
see README for full details and examples.''')

    parser.add_argument("-v","--version",
                        action = 'version',
                        version = "{} {}\tUpdated {}".format(program, version, update))
    parser.add_argument("-a","--alignment",
                        help = "Alignment file, must be an indexed BAM file",
                        metavar = 'BAM',
                        required = True)
    parser.add_argument("-f","--feature",
                        help = "Feature file, must be a GFF file",
                        metavar = 'GFF',
                        required = True)
    parser.add_argument("-o", "--output_prefix",
                        help = "Prefix for output files",
                        default = "{}.metagene.".format(date))
    parser.add_argument("-e", "--end",
                        help = "Flag to analyze the end of features.",
                        action = 'store_true')
    parser.add_argument("-s", "--start",
                        help = "Flag to analyze the start of features.",
                        action = 'store_true')
    parser.add_argument("-i", "--interval",
                        help = "Interval around reference position to query, default = 1000",
                        default = 1000,
                        type = int)
    parser.add_argument("-c","--chromosome_names",
                        help = "Chromosome conversion file, if the names are different between alignment and feature files.",
                        metavar = 'TAB')
    
    arguments = parser.parse_args()
    
    # make sure either start or end but not both are chosen for analysis
    if arguments.start and arguments.end:
        sys.exit("You can only analyze the start (-s) or the end (-e) of the gene at one time.")
    elif not(arguments.start) and not(arguments.end):
        arguments.start = True
        print "WARNING: neither start (-s) or end (-e) of feature was chosen for analysis. Analyzing start by default."
                         
    return arguments

def get_chromosome_sizes(bamfile):
    '''Uses samtools view -H to get the header information and returns a 
    dictionary of the chromosome names and sizes.'''
    
    chromosome_sizes = {}

##Need to add error handling for subprocess!!    
    header = subprocess.check_output(['samtools', 'view', "-H", bamfile]).split("\n")
    
    for line in header:
        if line[0:3] == "@SQ":
            # parse out chromosome information from @SQ lines
            # format example (tab-delimited):
            # @SQ   SN:chromosome_name  LN:chromosome_size
            line_parts = line.split("\t")
            chromosome_sizes[line_parts[1][3:]] = int(line_parts[2][3:])
    
    return chromosome_sizes

def get_chromosome_conversions(tabfile, bam_chromosomes):
    '''Import TAB delimited conversion table for the chromosome names in the 
    feature file (column 0) and in the alignment file (column 1).  Return said
    table as a dictionary with the feature file chromosome names as keys and the 
    alignment file chromosome names as values.'''

##Need to add checking of conversion file!!
    
    conversion_table = {}
    
    if arguments.chromosome_names:
        infile = open(tabfile)
    
        rows = infile.read().strip().split("\n")
        for r in rows:
            if r[0] != "#": # don't process comments
                row_parts = r.split("\t")
                conversion_table[row_parts[0]] = row_parts[1]
    
        infile.close()
    else:
        for c in bam_chromosomes:
            conversion_table[c] = c
    
    return conversion_table
     

if __name__ == "__main__":
    arguments = get_arguments(PROGRAM, VERSION, UPDATED)
    #print "Current arguments: \n{}".format(arguments)
    
    # confirm BAM file and extract chromosome sizes
    chromosomes = get_chromosome_sizes(arguments.alignment)
    #print "Current chromosomes: \n{}".format(chromosomes)
    
    # create chromosome conversion dictionary for GFF to BAM
    chromosome_conversion_table = get_chromosome_conversions(arguments.chromosome_names, chromosomes.keys())
    print "Current conversion table: \n{}".format(chromosome_conversion_table)
    
    
    
   
    

