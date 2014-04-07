#!/usr/bin/python
'''The third step of metagene_analysis, metagene_plot_only.py uses R to create 
the metagene plot (as a PDF).
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

PROGRAM = "metagene_plot_only.py"
VERSION = "0.1.0"
UPDATED = "140407 JRBT"

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
    parser.add_argument("-d","--data_set",
                        help = "comma-dilimited values of file.1.sense(or ungapped),file.1.antisene(or gapped),normalization.factor,color(for plotting),name(for plot legend)",
                        metavar = 'DATA_SET',
                        required = True,
                        action = 'append')                
    parser.add_argument("-o", "--output_prefix",
                        help = "Prefix for output files",
                        required = True)
    parser.add_argument("--feature_counted",
                        help = "Name of feature examined, eg TSS, Start, End, Gene, Intron",
                        required = True)
    arguments = parser.parse_args()
       
    return arguments
    
if __name__ == "__main__":
    arguments = get_arguments()
    
    total_sets = len(arguments.data_set)
    data_sets = []
    for data in arguments.data_set:
        for part in data.split(","):
            data_sets.append(str(part))
    
    
    # identify window.size, top and bottom labels from input filenames 
    try:
        (window_size, top1, top2) = re.findall('.(\d+)bpX\d+bp.([a-zA-Z]+)_([a-zA-Z]+).csv\Z', data_sets[0])[0]
        (bottom1, bottom2) = re.findall('.\d+bpX\d+bp.([a-zA-Z]+)_([a-zA-Z]+).csv\Z', data_sets[1])[0]
    except IndexError as err:
        raise MetageneError(err, "You must specify two files in each data_set -d option")
    
    if top1 == bottom1:
        top = top2
        bottom = bottom2
    else:
        top = top1
        bottom = bottom1
    
    # extract metagene information from first file
    with open(data_sets[0]) as inf:
        metagene = re.split('[\s-]+', inf.readline().strip())
        metagene_parts = {}
        for part in metagene:
            search = re.search('([A-Za-z]+):(\d+)', part)
            if search is not None:
                metagene_parts[search.group(1)] = int(search.group(2))
        total_start = -metagene_parts['Upstream']
        interval_start = 0
        interval_end = metagene_parts['Interval'] - 1
        total_end = metagene_parts['Interval'] + metagene_parts['Downstream'] - 1
        
    path_to_script = os.path.dirname(os.path.realpath(__file__))

    path_to_script += "/plot_only.R"
    call = ['Rscript', "--vanilla", "--verbose",
                     path_to_script, 
                     str(arguments.output_prefix), 
                     str(arguments.feature_counted),
                     str(window_size),
                     str(interval_start),
                     str(interval_end),
                     str(total_start),
                     str(total_end),
                     str(top),
                     str(bottom),
                     str(total_sets)]
    for data in data_sets:
        call.append(str(data))
    
    
    subprocess.call(call)
                         
   
    print "Finished plotting"
