#!/usr/bin/python
"""The first step of metagene_analysis, metagene_count.py compiles read abundance
over genomic features to create the input for metagene_windows.py. Please 
see README for full details and examples.

Requires:
    python 2 (https://www.python.org/downloads/)
    samtools (http://sourceforge.net/projects/samtools/files/) 

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
"""

import sys
import re
import datetime
import subprocess
import math
import argparse		# to parse the command line arguments
import timeit # for calculating run times...

# import classes
from MetageneError import MetageneError
from Metagene import Metagene
from Feature import Feature
from Read import Read

from metageneMethods import run_pipe
from metageneMethods import read_chunk

PROGRAM = "metagene_count.py"
VERSION = "0.1.2"
UPDATED = "140406 JRBT"

##TODO: better commandline argument parsing!!! And end user display...    
def get_arguments():
    """Collect and parse information from the user's command line arguments."""
    
    date = datetime.datetime.now().strftime('%y%m%d-%H%M%S')
    
    parser = argparse.ArgumentParser(description=
                                     """The first step of metagene_analysis, metagene_count.py compiles read abundance
                                     over genomic features to create the input for metagene_windows.py. Please
                                     see README for full details and examples.

                                     Requires:
                                     python 2 (https://www.python.org/downloads/),
                                     samtools (http://sourceforge.net/projects/samtools/files/)
                                     """)

    parser.add_argument("-v","--version",
                        action = 'version',
                        version = "{} {}\tUpdated {}".format(PROGRAM, VERSION, UPDATED))
    parser.add_argument("-a","--alignment",
                        help = "Alignment file, must be an indexed BAM file",
                        metavar = 'BAM',
                        required = True)
    parser.add_argument("-f","--feature",
                        help = "Feature file, must be a BED or GFF file",
                        metavar = 'BED/GFF',
                        required = True)
    parser.add_argument("-o", "--output_prefix",
                        help = "Prefix for output files",
                        default = "{}.metagene".format(date))
    
    parser.add_argument("--feature_count",
                        help = "Examine only the start, end or all of a feature",
                        choices = ['start','end','all'],
                        default = 'all')                    
    parser.add_argument("--count_method",
                        help = "Count the start, end or all of a read",
                        choices = ['start','end','all'],
                        default = 'all')
    parser.add_argument("--count_partial_reads",
                        help = "Include reads that only partially overlap with region of interest",
                        action = "store_true",
                        default = False)
    
    parser.add_argument("--padding",
                        help = "Padding in nt to add around the feature, default = 1000",
                        default = 1000,
                        type = int,
                        metavar = 'NT')
    parser.add_argument("--interval_size",
                        help = "Normalized size of the feature, default = 1000",
                        default = 1000,
                        type = int,
                        metavar = 'NT')
    parser.add_argument("--interval_variable",
                        help = "Allow variable interval size; will prevent true metagene analysis",
                        action = "store_true")
    parser.add_argument("--ignore_strand",
                        help = "Do not count reads by strand",
                        action = "store_true")
                       
                        
    parser.add_argument("--extract_abundance",
                        help = "Extract abundance information from NA:i:## tag in BAM file",
                        action = 'store_true')
    parser.add_argument("--extract_mappings",
                        help = "Extract number of mappings from NH:i:## tag in BAM file (required for hits normalization)",
                        action = 'store_true')
    parser.add_argument("--uniquely_mapping",
                        help = "Flag if reads are all uniquely mapping",
                        action = 'store_true')

    parser.add_argument("-c","--chromosome_names",
                        help = "Chromosome conversion file (feature_chromosome {tab} alignment_chromosome)",
                        metavar = 'TAB',
                        default = "None")
    
    parser.add_argument("--count_splicing",
                        help = "Count reads as spliced or unspliced",
                        action = 'store_true',
                        default = False)
                        
                        
    parser.add_argument("--include_reads",
                        help = "Include reads with these features, repeat tag up to 4 times. Hint: can ignore if BAM column 2 < 256",
                        choices = ['secondary_alignment', 'failed_quality_control', 'PCR_duplicate', 'supplementary_alignment'],
                        action = 'append',
                        default = [])
                       
    
    arguments = parser.parse_args()
    
    # adjust internal_size if only the start or end of the feature will be counted
    if arguments.feature_count != 'all':
        arguments.interval_size = 1
    
    arguments.count_secondary_alignments = False
    arguments.count_failed_quality_control = False
    arguments.count_PCR_optical_duplicate = False
    arguments.count_supplementary_alignment = False
    
    if len(arguments.include_reads) >= 1:
        if 'secondary_alignment' in arguments.include_reads:
            arguments.count_secondary_alignments = True
        
        if 'failed_quality_control' in arguments.include_reads:
            arguments.count_failed_quality_control = True
        
        if 'PCR_duplicate' in arguments.include_reads:
            arguments.count_PCR_optical_duplicate = True
        
        if 'supplementary_alignment' in arguments.include_reads:
            arguments.count_supplementary_alignment = True
        
    return arguments


           
def metagene_count():
    '''Chain of command for metagene_count analysis.'''
    arguments = get_arguments()
    # confirm BAM file and extract chromosome sizes
    Read.set_chromosome_sizes(arguments.alignment)
##TODO: create a list of chromosomes to analyze and/or exclude
    # create chromosome conversion dictionary for feature (GFF/BED) to alignment (BAM)
    Feature.set_chromosome_conversion(arguments.chromosome_names, Read.chromosome_sizes.keys()) 
    
    # define has_abundance and has_mappings tags for Read class
    Read.set_sam_tag(arguments.extract_abundance, arguments.alignment, "NA:i:(\d+)")
    Read.set_sam_tag(arguments.extract_mappings, arguments.alignment, "NH:i:(\d+)")
          
    # define the metagene array shape (left padding, start, internal, end, right padding)
    # metagene = padding ---- internal region ---- padding 
    try:
        metagene = Metagene(arguments.interval_size, arguments.padding, arguments.padding)
        print "Metagene definition:\t{}".format(metagene)
    except MetageneError as err:
        print err
        raise MetageneError("Unable to create the metagene template")
    
    try:
        Feature.set_format(arguments.feature) # assign file format for the feature file
        print "Reading feature file as {} format".format(Feature.format)
    except MetageneError as err:
        print err
        raise MetageneError("Unable to create the feature object")
    
    # print out the header line...
    if not arguments.interval_variable:
        with open("{}.metagene_counts.csv".format(arguments.output_prefix), 'w') as output_file:
            output_file.write("# Metagene:\t{}\n".format(metagene)) # define for plotting later
            output_file.write(metagene.print_full())
     
    # for each feature
    with open(arguments.feature, 'r') as feature_file:
        for feature_line in read_chunk(feature_file, 1024):
            if feature_line[0] != "#": # skip comment lines
                # change creation with feature_method
                feature = Feature.create(arguments.feature_count, metagene, feature_line, arguments.count_splicing, arguments.ignore_strand)
                
                # pull out sam file lines; it is important to use Feature.get_samtools_region(chromosome_lengths) rather
                # than Feature.get_chromosome_region() because only the first ensures that the interval does not
                # extend beyond the length of the chromosome which makes samtools view return no reads
                (runPipe_worked, sam_sample) = run_pipe(['samtools view {} {}'.format(
                        arguments.alignment,
                        feature.get_samtools_region())])
                if runPipe_worked:
                    for samline in sam_sample:
                        if len(samline) > 0:
                            # create Read feature
                            (created_read, read) = Read.create_from_sam(samline, 
                                                                        Feature.chromosome_conversion.values(),
                                                                        arguments.count_method, 
                                                                        arguments.uniquely_mapping,
                                                                        arguments.ignore_strand,
                                                                        arguments.count_secondary_alignments,
                                                                        arguments.count_failed_quality_control,
                                                                        arguments.count_PCR_optical_duplicate,
                                                                        arguments.count_supplementary_alignment)

                            # count read (if it exists)
                            if created_read:
                                feature.count_read(read, arguments.count_method, arguments.count_splicing, arguments.count_partial_reads, arguments.ignore_strand)

                    # output the resulting metagene
                    with open("{}.metagene_counts.csv".format(arguments.output_prefix), 'a') as output_file:
                        output_file.write("{}\n".format(feature.print_metagene(interval_override=arguments.interval_variable)))
                    
                else:
                    raise MetageneError("Could not pull chromosomal region {} for feature {} from BAM file {}.".format(feature.get_chromosome_region(), feature.name, arguments.alignment))

    
if __name__ == "__main__":
    
    start_time = timeit.default_timer()
    try:
        metagene_count()    
    except MetageneError as err:
        print "\n{}\n".format(err)
        print "Aborting metagene analysis..."
        sys.exit()
    end_time = timeit.default_timer()
    print "Run time:\t{}".format(end_time - start_time)
    

