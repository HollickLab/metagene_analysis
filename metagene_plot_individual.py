#!/usr/bin/python
'''The third step of metagene_analysis, metagene_plot.py uses R to create the
metagene plot (as a PDF) and its associated statistics.
Please see README for full details and examples.

Requires:
    python 2 (https://www.python.org/downloads/)
    R (http://cran.us.r-project.org/)

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

import sys
import subprocess
import re
import datetime
import os
import argparse		# to parse the command line arguments

PROGRAM = "metagene_plot_individual.py"
VERSION = "0.1.0"
UPDATED = "140406 JRBT"

def get_arguments():
    '''Collect and parse information from the user's command line arguments.'''
    
    date = datetime.datetime.now().strftime('%y%m%d-%H%M%S')
    
    parser = argparse.ArgumentParser(description=
    '''The third step of metagene_analysis, metagene_plot.py uses R to create the
metagene plot (as a PDF) and its associated statistics.
Please see README for full details and examples.

Requires:
    python 2 (https://www.python.org/downloads/)
    R (http://cran.us.r-project.org/)
    ''')

    parser.add_argument("-v","--version",
                        action = 'version',
                        version = "{} {}\tUpdated {}".format(PROGRAM, VERSION, UPDATED))
    parser.add_argument("-a","--fileset_a",
                        help = "Bin files from metagene_bin.py",
                        metavar = 'BIN_FILE',
                        required = True,
                        action = 'append')
    parser.add_argument("-b","--fileset_b",
                        help = "Bin file from metagene_bin.py",
                        metavar = 'BIN_FILE',
                        required = True,
                        action = 'append')                   
    parser.add_argument("-o", "--output_prefix",
                        help = "Prefix for output files",
                        required = True)
    
    parser.add_argument("--normalization",
                        help = "Normalization value for file_a; eg total reads would result in reads per million",
                        type = int)
    
    parser.add_argument("--individual_splicing",
                        help = "Plot gapped vs ungapped reads for a single gene",
                        action ="store_true")

    arguments = parser.parse_args()
       
    return arguments
    
if __name__ == "__main__":
    arguments = get_arguments()
    
    try:
        parsed_file_a = re.findall('.(\d+)bpX(\d+)bp.([a-zA-Z]+)_([a-zA-Z]+).csv\Z', arguments.fileset_a[0])[0]
        parsed_file_b = re.findall('.(\d+)bpX(\d+)bp.([a-zA-Z]+)_([a-zA-Z]+).csv\Z', arguments.fileset_b[0])[0]
    except IndexError as err:
        raise MetageneError(err, "You must specify two files for each group -a and -b")
    
    # extract metagene information from first file
    with open(arguments.fileset_a[0]) as inf:
        metagene = re.split('[\s-]+', inf.readline().strip())
        metagene_parts = {}
        for part in metagene:
            search = re.search('([A-Za-z]+):(\d+)', part)
            if search is not None:
                metagene_parts[search.group(1)] = int(search.group(2))
        upstream_start = -metagene_parts['Upstream']
        upstream_end = -1
        interval_start = 0
        interval_end = metagene_parts['Interval'] - 1
        downstream_start = metagene_parts['Interval']
        downstream_end = metagene_parts['Interval'] + metagene_parts['Downstream'] - 1
        
    path_to_script = os.path.dirname(os.path.realpath(__file__))

    if arguments.individual_splicing:
        path_to_script += "/plot_splicing_individual.R"
        subprocess.call(['Rscript', 
                     path_to_script, 
                     str(arguments.fileset_a[0]),     # gapped.file
                     str(arguments.fileset_b[0]),     # ungapped.file
                     str(arguments.normalization), # normalization
                     str(downstream_start),        # gene.length
                     str(arguments.output_prefix)])# output.prefix
                     
   
    print "Finished plotting"
