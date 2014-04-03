#!/usr/bin/python
'''The first step of metagene_analysis, metagene_count.py compiles read abundance
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
'''

import sys, re, datetime, subprocess, math
import argparse		# to parse the command line arguments

# import classes
from MetageneError import MetageneError
from Metagene import Metagene
from Feature import Feature
from Read import Read

PROGRAM = "metagene_count.py"
VERSION = "0.1.2"
UPDATED = "140402 JRBT"

    
def get_arguments():
    '''Collect and parse information from the user's command line arguments.'''
    
    date = datetime.datetime.now().strftime('%y%m%d-%H%M%S')
    
    parser = argparse.ArgumentParser(description=
    '''The first step of metagene_analysis, metagene_count.py compiles read abundance
over genomic features to create the input for metagene_windows.py. Please 
see README for full details and examples.

Requires:
    python 2 (https://www.python.org/downloads/), 
    samtools (http://sourceforge.net/projects/samtools/files/)
    ''')

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
                        
    parser.add_argument("--extract_abundance",
                        help = "Extract abundance information from NA:i:## tag in BAM file",
                        action = 'store_true')
    parser.add_argument("--extract_mappings",
                        help = "Extract number of mappings from NH:i:## tag in BAM file (required for hits normalization)",
                        action = 'store_true')

    parser.add_argument("-c","--chromosome_names",
                        help = "Chromosome conversion file (feature_chromosome {tab} alignment_chromosome)",
                        metavar = 'TAB',
                        default = "None")
    
    parser.add_argument("--count_splicing",
                        help = "Count reads as spliced or unspliced",
                        action = 'store_true')
                        
    parser.add_argument("--include_reads",
                        help = "Include reads with these features, repeat tag up to 4 times. Hint: can ignore if BAM column 2 < 256",
                        choices = ['secondary_alignment', 'failed_quality_control', 'PCR_duplicate', 'supplementary_alignment'],
                        action = 'append')
    
    arguments = parser.parse_args()
    
    # adjust internal_size if only the start or end of the feature will be counted
    if arguments.feature_count != 'all':
        arguments.interval_size = 1
    
    arguments.count_secondary_alignments = False
    arguments.count_failed_quality_control = False
    arguments.count_PCR_optical_duplicate = False
    arguments.count_supplementary_alignment = False
    
    if 'secondary_alignment' in arguments.include_reads:
        arguments.count_secondary_alignments = True
        
    if 'failed_quality_control' in arguments.include_reads:
        arguments.count_failed_quality_control = True
        
    if 'PCR_duplicate' in arguments.include_reads:
        arguments.count_PCR_optical_duplicate = True
        
    if 'supplementary_alignment' in arguments.include_reads:
        arguments.count_supplementary_alignment = True
        
    return arguments

def get_chromosome_sizes(bamfile):
    '''Uses samtools view -H to get the header information and returns a 
    dictionary of the chromosome names and sizes.'''
    
    chromosome_sizes = {}

##TODO: Add error handling for subprocess!!    
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
    
    if tabfile != "None":
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

def runPipe(cmds):
    '''runPipe function is from danizgod's post at stackoverflow exchange: 
    http://stackoverflow.com/questions/9655841/python-subprocess-how-to-use-pipes-thrice
    
    Usage: runPipe(['ls -1','head -n 2', 'head -n 1'])'''
    
    try: 
        p = subprocess.Popen(cmds[0].split(' '), stdin = None, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        prev = p
        for cmd in cmds[1:]:
            p = subprocess.Popen(cmd.split(' '), stdin = prev.stdout, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            prev = p
        stdout, stderr = p.communicate()
        p.wait()
        returncode = p.returncode
    except Exception, e:
        stderr = str(e)
        returncode = -1
    if returncode == 0:
        return (True, stdout.strip().split('\n'))
    else:
        return (False, stderr)

def has_sam_tag(bamfile, tag):
    '''Checks for a particular tag in the sam lines.'''
    
    (runPipe_worked, sam_sample) = runPipe(['samtools view {}'.format(bamfile), 'head -n 10'])
    if runPipe_worked:
        num_tags = 0
        for sam_line in sam_sample:
            if re.search(tag,sam_line) != None:
                num_tags += 1
        if num_tags == 10:
            return True
        else:
            return False
    else:
        raise MetageneError(sam_sample, "Checking the bam file failed with error: {}".format(sam_sample))
    
def determine_format(feature_file):
    '''Distinguish between BED and GFF file types.
    BED: chromosome, start, end, name, score, strand
    GFF: chromosome, label, label, start, end, ?, strand, ?, name/misc (often ; delimited)
    
    Distinguishing points: columns 1 & 2 ints and 5 (if it exists) is +, -, or . ==> BED
                           columns 3 & 4 ints and 6 is +, -, or . and 7 or more columns ==> GFF'''
    
    counts = {'BED':0, 'BED_SHORT':0, 'GFF':0, 'UNKNOWN':0}
    header = 0
    total = 0
    
    try:
        with open(feature_file, 'r') as infile:
            for line in infile.readlines(50): # read first 50 bytes
                total += 1
                if line[0] == "#":
                    header += 1
                elif re.search('\A\S+\t\d+\t\d+\t\S+\t\S+\t[+.-]\s+',line) != None:
                    counts['BED'] += 1
                elif re.search('\A\S+\t\S+\t\S+\t\d+\t\d+\t\S+\t[+.-]\t\S+\t\S+\s+', line) != None:
                    counts['GFF'] += 1
                elif re.search('\A\S+\t\d+\t\d+', line) != None:
                    counts['BED_SHORT'] += 1
                else:
                    counts['UNKNOWN'] += 1
    except IOError as err:
        infile.close()
        raise MetageneError(err, "Could not open the feature file.")
    
    # require that at least 80% of the sampled lines are classified the same to auto-determine
    values = list(counts.values())
    keys = list(counts.keys())
    max_key = keys[values.index(max(values))]
    if max(values) >= 0.8 * (total - header) and max_key != "UNKNOWN":
        return max_key
    else:
        raise MetageneError(feature_file, "Could not determine the format of the feature file.")
    
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
  
           
def metagene_count():
    '''Chain of command for metagene_count analysis.'''

    arguments = get_arguments()
    #print "Current arguments: \n{}".format(arguments)
    
    # confirm BAM file and extract chromosome sizes
    chromosomes = get_chromosome_sizes(arguments.alignment)
    #print "Current chromosomes: \n{}".format(chromosomes)
    
    # if extract_abundance or extract_mappings confirm that appropriate tags exist
    if arguments.extract_abundance and not has_sam_tag(arguments.alignment, "NA:i:\d+"):
        raise MetageneError(arguments.extract_abundance, "Your alignment file does not have the required NA:i:## abundance tags\nAdjust the alignment file or remove the --extract_abundance tag to treat each read as an abundance of 1")
    
    if arguments.extract_mappings and not has_sam_tag(arguments.alignment, "NH:i:\d+"):
        raise MetageneError(arguments.extract_mappings, "Your alignment file does not have the required NH:i:## mappings tags\nAdjust the alignment file or remove the --extract_mappings tag to treat each read as unique and remove hits-normalization")
       
    # create chromosome conversion dictionary for feature (GFF/BED) to alignment (BAM)
    chromosome_conversion_table = get_chromosome_conversions(arguments.chromosome_names, chromosomes.keys())   
    #print "Current conversion table: \n{}".format(chromosome_conversion_table)
    
    # define the metagene array shape (left padding, start, internal, end, right padding)
    # metagene = padding ---- internal region ---- padding 
    try:
        metagene = Metagene(arguments.interval_size, arguments.padding, arguments.padding)
        print "Metagene definition:\t{}".format(metagene)
    except MetageneError as err:
        print err
        raise MetageneError(err, "Unable to create the metagene template")
    
    try:
        feature_format = determine_format(arguments.feature) # distinguish between BED, BED_SHORT, and GFF formats
        print "Reading feature file as {} format".format(feature_format)
    except MetageneError as err:
        print err
        raise MetageneError(err, "Unable to create the feature object")

              
    # for each feature
    
    with open(arguments.feature, 'r') as feature_file:
        for feature_line in read_chunk(feature_file, 1024):
            if feature_line[0] != "#": # skip comment lines
                # change creation with feature_method
                feature = Feature.create(feature_format, arguments.feature_count, metagene, feature_line, chromosome_conversion_table)
                
                ##print feature
                
                # pull out sam file lines; it is important to use Feature.get_samtools_region(chromosome_lengths) rather
                # than Feature.get_chromosome_region() because only the first ensures that the interval does not
                # extend beyond the length of the chromosome which makes samtools view return no reads
                (runPipe_worked, sam_sample) = runPipe(['samtools view {} {}'.format(arguments.alignment,feature.get_samtools_region(chromosomes))])
                if runPipe_worked:
                    for samline in sam_sample:
                        if len(samline) > 0:
                            # create Read feature
                            (created_read, read) = Read.create_from_sam(samline, 
                                                                        chromosome_conversion_table, 
                                                                        arguments.count_method, 
                                                                        arguments.extract_abundance, 
                                                                        not(arguments.extract_mappings),
                                                                        arguments.count_secondary_alignments,
                                                                        arguments.count_failed_quality_control,
                                                                        arguments.count_PCR_optical_duplicate,
                                                                        arguments.count_supplementary_alignment)
                            ##print read
                            # count read (if it exists)
                            if created_read:
                                feature.count_read(read, arguments.count_method)
                                ##print Feature.__str__(feature,counts_only=True)
                    # output the resulting metagene
                    with open("{}.metagene_counts.csv".format(arguments.output_prefix), 'a') as output_file:
                        ##print "Finished counting reads..."
                        ##print Feature.__str__(feature,counts_only=True)
                        ##print feature.print_metagene(pretty=True)
                        output_file.write("{}\n".format(feature.print_metagene()))
                    
                else:
                    ##print sam_sample
                    raise MetageneError(sam_sample, "Could not pull chromosomal region {} for feature {} from BAM file {}.".format(feature.get_chromosome_region(), feature.name, arguments.alignment))
      
   

    
if __name__ == "__main__":
    try:
        metagene_count()    
    except MetageneError as err:
        print "\n{}\n".format(err)
        print "Aborting metagene analysis..."
        sys.exit()
    

