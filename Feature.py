#!/usr/bin/python
'''Feature class for metagene_counts.py.
Child of Metagene class.

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
import re
import math

from Metagene import Metagene
from MetageneError import MetageneError
from Read import Read

from metageneMethods import confirm_integer
from metageneMethods import runPipe

class Feature(Metagene):
    '''A Feature is a Metagene object representing an interval of interest  
    with padding on either side, where the positions are defined by chromosomal
    (nucleotide) positions.  
    
    Features can each be different and are used to 
    count aligning reads; then each feature is standardized into a Metagene 
    shape for comparison between features and summarization into a metagene profile.  
    
    Attributes:
        (Inherited)
        feature_interval : interval of interest (int >= 1)
        padding          : dictionary of padding values (int >= 0)
          'Upstream'     : padding upstream (to left) of feature_interval
          'Downstream'   : padding downstream (to right) of feature_interval
        length           : length of entire feature object
        
        (Feature-specific)
        name             : name of feature; can also be used to store additional
                           information for parsing by user; avoid white-space and
                           commas
        chromosome       : chromosome; matches alignment (read) file chromosome name
        strand           : strand of feature: + = Watson; - = Crick; . = Unknown
        metagene_length  : length of the metagene feature interval on which we will map 
                           the final feature (non-padding) counts
        position_array   : array of chromosome positions covered by the feature
                           (including padding); orientated so index 0 is always
                           the start of the feature and last index is the end of 
                           the feature; for Crick-strand feature:
                           position_array[0] > position_array[len(position_array)]
        counts_array     : dictionary of counting objects, value array length is
                           the same as the length of the position_array;
                           gapped/ungapped/allreads can be added to strand information
          'sense'        : counts for sense reads (both read and feature strand is the same)
          'antisense'    : counts for antisense reads (opposite of sense)
          'unstranded'   : counts for when strand is not present or being ignored
          
          'gapped'       : counts for reads with gaps relative to the feature
          'ungapped'     : counts for reads without gaps relative to the feature
          'allreads'     : counts for reads with or without gaps
          
    
    Methods:
        get_chromosome_region()
        get_samtools_region(chromosome_lengths_dictionary)
        print_metagene(pretty=False)   
        adjust_to_metagene(feature_array, verbose=False)
        count_read(read_object, count_method, count_gaps)
    
    Class Methods:
        create(file_format, count_method, metagene_object, feature_line, chromosome_conversion_table)
        create_from_bed(count_method, metagene_object, feature_line, chromosome_conversion_table, short=False)
        create_from_gff(count_method, metagene_object, feature_line, chromosome_conversion_table)

        '''

    __slots__ = ['name', 'chromosome', 'strand', 'metagene_length','counts_array', 'position_array']
    # inherits feature_interval, padding, and length from Metagene
    
    format = "Unknown" # format of feature file (current options handled are BED, SHORT_BED, and GFF)
    previously_warned_start_bigger_than_end = False
    chromosome_conversion = {}
    
    def __init__(self, count_method, metagene_object, name, chromosome, start, end, strand, gap_counting=False):
        '''Not normally called directly; use Feature.create(file_format, count_method, 
        metagene_object, feature_line, chromosome_conversion_table) to call indirectly.
        
        Define a new feature with an interval (represents feature length), 
        up and downstream padding (defined by metagene_object), and genomic 
        (1-based) start and end positions.
        
        Once defined here, the start and end represent the true start and end of
        the feature.  Therefore, if a - strand (Crick strand) feature the start
        will be larger than the end.''' 
        chromosome = Feature.chromosome_conversion[chromosome] # convert to BAM-like chromosome designation
        if (confirm_integer(start, "Start", minimum=1, maximum=Read.chromosome_sizes[chromosome]) and 
                confirm_integer(end, "End", minimum=1, maximum=Read.chromosome_sizes[chromosome])):
           start = int(start)
           end = int(end)
        
        # Define feature-specific metagene where feature_interval respresents 
        # the length of the feature NOT the length of the final metagene interval
        if count_method == 'all':
            interval = (end - start + 1) # length of feature
        else:
            interval = 1 # length of the start (or end) of feature
            
        Metagene.__init__(self,interval, metagene_object.padding['Upstream'], metagene_object.padding['Downstream'])
        self.name = name
        self.chromosome = chromosome  
        self.strand = strand
        self.metagene_length = metagene_object.feature_interval
        
        # define counts_array dictionary 
        # key: orientation:gap_counts string 
        #      where orientation = {'unstranded', 'sense', 'antisense'}
        #            gap_counts  = {'ungapped', 'gapped, 'allreads'}
        #      'ungapped' + 'gapped' = 'allreads'
        #      'sense' + 'antisense' = 'unstranded'
        # 
        # values: arrays of self.length initialized to 0
        if self.strand != "+" and self.strand != "-":
            self.strand = "."
            orientation = ['unstranded']
        else:
            orientation = ['sense', 'antisense']
        if gap_counting:
            gap_counts = ['ungapped', 'gapped']
        else:
            gap_counts = ['allreads']
            
        self.counts_array = {}
        for o in orientation:
            for g in gap_counts:
                self.counts_array["{}:{}".format(o,g)] = []
                for p in range(self.length):
                    #self.counts_array["{}:{}".format(o,g)].append(decimal.Decimal(0.0))
                    self.counts_array["{}:{}".format(o,g)].append(0)
                    
        # define position_array
        # values  : chromosomal 1-based nucleotide positions in 5' to 3' 
        #           orientation WITH RESPECT TO THE FEATURE
        # Example :
        #       + strand:   [10,11,12,13,14,15]
        #       - strand:   [15,14,13,12,11,10] 
        # so position_array[0] is always the start of the feature (with upstream padding) 
        #    position_array[-1] is always the end of the feature (with downstream padding)
        self.position_array = [] 
        if self.strand == "-": 
            # chromosome start = feature end
            # chromosome end   = feature start
            if count_method == 'start':
                start = end 
            elif count_method == 'end':
                end = start
            region_start = start - self.padding['Downstream'] # start is really end
            region_end = end + self.padding['Upstream'] # end is really start
            positions = range(region_start, region_end + 1) # inclusive list
            positions.reverse()
        else:
            if count_method == 'start':
                end = start # set both start and end to the start value
            elif count_method == 'end':
                start = end # set both start and end to the end value
            region_start = start - self.padding['Upstream'] 
            region_end = end + self.padding['Downstream'] 
            positions = range(region_start, region_end + 1) # inclusive list
        
        self.position_array = positions            
    # end Feature.__init__ function
            
    
    def __str__(self, counts_only=False):
        '''Returns pretty graphic of feature information and current counts.'''
        output = ""

        # skip if counts_only option is enabled
        if not(counts_only):
            output += "{} at {} on {} strand\n".format(self.name, self.get_chromosome_region(), self.strand)
        
        # create up(stream), int(erval) and down(stream) labels for each position
        output += "\t\t\t\t"
        for i in range(self.padding['Upstream']):
            output += "---up-"
        for i in range(self.feature_interval):
            output += "--int-"
        for i in range(self.padding['Downstream']):
            output += "-down-"
        output += "\n"
        
        # print out position information
        output += "{0:20s}:\t".format('Position')
        for i in self.position_array:
            output += "{0:5d},".format(i)
        output = output[:-1] + "\n"
        
        # print out counts information           
        for orientation in sorted(self.counts_array.keys(), reverse=True):
            output += "{0:20s}:\t".format(orientation)
            for i in self.counts_array[orientation]:
                output += "{0:>5s},".format("{0:3.2f}".format(i))
            output = output[:-1] + "\n"
        return output
    # end of Feature.__str__ function            
                    
    
    def get_chromosome_region(self):
        '''Return position interval for samtools view (chromosome: start-end (1-based))'''
        if self.strand == "-":
            # order from smaller to larger chromosome position
            return ("{}:{}-{}".format(self.chromosome, self.position_array[-1], self.position_array[0]))
        else:
            return ("{}:{}-{}".format(self.chromosome, self.position_array[0], self.position_array[-1]))
    
    
    def get_samtools_region(self, chromosome_lengths_dictionary=Read.chromosome_sizes):
        '''Return a position interval valid for use in samtools view program.
        
        Same as get_chromosome_region(), but adjust start to 0 and 
        end to chromosome length if they exceed chromosome boundaries.'''
        
        if self.strand == "-":
            start = self.position_array[-1]
            end = self.position_array[0]
        else:
            start = self.position_array[0]
            end = self.position_array[-1]
        
        if start < 1:
            start = 1  
        try: 
            if end > chromosome_lengths_dictionary[self.chromosome]:
                end = chromosome_lengths_dictionary[self.chromosome]
        except:
            raise MetageneError("Could not find chromosome {} in length dictionary".format(self.chromosome))
        
        return ("{}:{}-{}".format(self.chromosome, start, end))
    # end of Feature.get_samtools_region function    
    
    
    def print_metagene(self, pretty=False, header=False):
        '''Converts counts_array data to finalized metagene profiles for printing
        
        Standard printing is in comma-delimited lines for input into metagene_analysis.py
        Pretty printing (pretty=True) gives a human readable, if potentially super long, version'''
        
        final_metagenes = {}
        
        if header:
            metagene = Metagene(self.metagene_length, self.padding['Upstream'], self.padding['Downstream'])
            output = metagene.print_full(pretty)
        else:
            output = ""
        
        # process each subset grouping    
        for subset in sorted(self.counts_array, reverse=True):
            # break counts_array into sections -> upstream padding, interval_feature, and downstream padding
            upstream_counts = self.counts_array[subset][0:self.padding['Upstream']]
            interval_counts = self.counts_array[subset][self.padding['Upstream'] : self.padding['Upstream'] + self.feature_interval]
            downstream_counts = self.counts_array[subset][self.padding['Upstream'] + self.feature_interval : len(self.counts_array[subset])]
        
            # compress (or expand) interval_counts to match the size of the internal metagene
            metagene_interval_counts = self.adjust_to_metagene(interval_counts)
            
            if pretty:
                output += "{0:15s}:\t".format(subset)
                for i in upstream_counts:
                    output += "{0:>5s},".format("{0:3.2f}".format(i)) # keep 2 decimal places in the outputted float
                for i in metagene_interval_counts:
                    output += "{0:>5s},".format("{0:3.2f}".format(i))
                for i in downstream_counts:
                    output += "{0:>5s},".format("{0:3.2f}".format(i))  
                output = output[:-1] + "\n"
            else:    
                # build output
                output += "{},{}".format(self.name, subset)
                for p in upstream_counts:
                    output += ",{0:0.3f}".format(p) # keep 3 decimal places in the outputted float
                for p in metagene_interval_counts:
                    output += ",{0:0.3f}".format(p)
                for p in downstream_counts:
                    output += ",{0:0.3f}".format(p)
                output += "\n"
        
        return output.strip() # remove trailing "\n"
    # end of print_metagene function     
         

    def adjust_to_metagene(self, feature_array, verbose=False):
        '''Expand or collapse the counts data from interval_array into a metagene
        array via the given shrink factor. -- basically a smoothing function :-)'''

        metagene_array = []
        
        # initialize metagene_count and remaining_metagene
        metagene_count = 0.0 #decimal.Decimal(0.0)
        shrink_factor = len(feature_array) / float(self.metagene_length) #decimal.Decimal(len(feature_array)) / decimal.Decimal(self.metagene_length)
        remaining_metagene_bin = shrink_factor
        
        loop_index = 0 # loop index for verbose option
        
        for bin in feature_array: # ensure all data is moved to metagene by looping  through entire interval_array
            # Ideally add in units of 1 (1 bin to 1 metagene_array  position) 
            # unless not possible then start dealing with fractional bins
            
            remaining_feature_bin = 1.0 #decimal.Decimal(1.0) # reset remaining feature for new bin
            
            if verbose:
                print "\n  Loop {}:".format(loop_index)
                print "    Feature Count :\t{}".format(bin)
                print "    Metagene Count:\t{}".format(metagene_count)
                print "    Metagene Bin  :\t{}".format(remaining_metagene_bin)
                print "    Feature Bin   :\t{}".format(remaining_feature_bin)
                loop_index += 1
                while_index = 0 # while loop index
                      
            while remaining_feature_bin > 0: #decimal.Decimal(0.0):
                # keeping adding from this bin until its empty
                
                if verbose:
                    while_index += 1
                    print "    While loop {}:".format(while_index)
                    print "      Feature Count :\t{}".format(bin)
                    print "      Metagene Count:\t{}".format(metagene_count)
                    print "      Metagene Bin  :\t{}".format(remaining_metagene_bin)
                    print "      Feature Bin   :\t{}".format(remaining_feature_bin)
                                
                if remaining_feature_bin <= remaining_metagene_bin:
                    if verbose:
                        print "      Add Remaining Feature Bin:\t{}".format(bin * remaining_feature_bin)    
                    
                    # add entire remaining feature to metagene
                    metagene_count += (bin * float(remaining_feature_bin)) #(decimal.Decimal(bin) * remaining_feature_bin)
                    # adjust bin counters 
                    remaining_metagene_bin -= remaining_feature_bin 
                    remaining_feature_bin = 0 #decimal.Decimal(0.0)
                else:
                    if verbose:
                        print "      Add Remaining Metagene Bin:\t{}".format(bin * remaining_metagene_bin)
                    
                    # add entire remaining_metagene_bin amount of feature to metagene
                    metagene_count += (bin * float(remaining_metagene_bin)) #(decimal.Decimal(bin) * remaining_metagene_bin)
                    # adjust bin counters
                    remaining_feature_bin -= remaining_metagene_bin
                    remaining_metagene_bin = 0 #decimal.Decimal(0.0)
                
                if verbose:
                    print "Remaining_metagene_bin:\t{}".format(remaining_metagene_bin)                    
                # check to see if new metagene bin is ready to be added to the metagene_array
                if remaining_metagene_bin == 0: #decimal.Decimal(0.0):
                    if verbose:
                        print "      Add Count to Metagene Array:\t{}".format(metagene_count)
                    metagene_array.append(metagene_count)
                    metagene_count = 0.0 #decimal.Decimal(0.0)
                    remaining_metagene_bin = shrink_factor
            # end of while loop through current feature bin
        # end of for loop through feature array
        
        if len(metagene_array) < self.metagene_length:
            # print out final metagene that was missed
            if verbose:
                print "      Add Count to Metagene Array:\t{}".format(metagene_count)
            metagene_array.append(metagene_count)
        
        if verbose:                                    
            print "\n  Final Metagene:\t{}".format(metagene_array)
        return metagene_array
    # end of adjust_to_metagene
    
       
    def count_read(self, read_object, count_method, count_gaps=False, count_partial_reads=False):
        '''Add a read object to the sense or antisense counts_array. Requires strand
        options of "+", "-" or "."
        
        Only stranded reads (+ or -) can be counted on stranded features.
        
        Unstranded Features will count in the + direction, but ignore read strand. 
        
        count_gaps=False default to not separate gapped and ungapped reads into different tallies
        count_partial_reads=False default to ignore reads that only partially overlap with the feature'''
        
        # determine orientation (and if countable)
        if self.strand == ".":
            subset = 'unstranded'
        elif read_object.strand != ".":
            if self.strand == read_object.strand:
                subset = 'sense'
            else:
                subset = 'antisense'
        else: 
            raise MetageneError("Can not count unstranded reads on stranded features.")
        
        # determine gap status
        if count_gaps:
            if abs(read_object.position_array[0] - read_object.position_array[-1]) + 1 > len(read_object.position_array):
                # if calculated length > actual length then gapped
                subset += ":gapped"
            else:
                subset += ":ungapped"            
        else:
            subset += ":allreads"
        
        # do we care if the read fully fits?
        if not count_partial_reads: # yes
            # does the read extend beyond the window?
            if read_object.position_array[0] not in self.position_array or read_object.position_array[-1] not in self.position_array: # yes
                # don't count anything then
                return False
        
        # can count if they are on the same chromosome
        if self.chromosome == read_object.chromosome:
            # get positions from read to potentially count
            positions_to_count = []
        
            if count_method == 'start':
                positions_to_count.append(read_object.position_array[0])
            elif count_method == 'end':
                positions_to_count.append(read_object.position_array[-1])
            elif count_method == 'all':
                positions_to_count = read_object.position_array
            else:
                raise MetageneError("Unrecognizable counting method.  Valid options are 'start', 'end', and 'all'")
        
            for p in positions_to_count:
                # make sure it overlaps with the Feature
                if p in self.position_array:
                    self.counts_array[subset][self.position_array.index(p)] += (read_object.abundance / float(read_object.mappings)) #(decimal.Decimal(read_object.abundance) / decimal.Decimal(read_object.mappings))        
    # end of count_read function        
    
    
    #******** creating Feature objects from diffent feature file formats (eg BED and GFF) ********#
    @classmethod
    def set_format(cls,feature_file):
        '''Determines and assigns file format for the feature file.  
        
        Currently, distinguishes between BED and GFF file types.
        BED:        chromosome, start, end, name, score, strand    ==> BED6 or BED12 format (ignores last 6 BED12 columns)
        BED_SHORT:  chromosome, start, end, name(name is optional) ==> BED3 or BED4 format
        GFF:        chromosome, label, label, start, end, ?, strand, ?, name/misc (often ; delimited)
    
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
            raise MetageneError("Could not open the feature file.")
    
        # require that at least 80% of the sampled lines are classified the same to auto-determine
        values = list(counts.values())
        keys = list(counts.keys())
        max_key = keys[values.index(max(values))]
        if max(values) >= 0.8 * (total - header) and max_key != "UNKNOWN":
            cls.format = max_key
        else:
            raise MetageneError("Could not determine the format of the feature file.")
        
        return True
    # end of set_format classmethod        
        
    @classmethod
    def create(cls, count_method, metagene, feature_line):
        '''Calls individual creation methods based on format.'''
        if Feature.format == "GFF":
            return Feature.create_from_gff(count_method, metagene, feature_line)
        elif Feature.format == "BED": # BED6 or BED12
            return Feature.create_from_bed(count_method, metagene, feature_line)
        elif Feature.format == "BED_SHORT":
            return Feature.create_from_bed(count_method, metagene, feature_line, short=True)
        else:
            raise MetageneError("Could not determine the format of features in the feature file")
            
    @classmethod
    def create_from_bed(cls, count_method, metagene_object, bed_line, short=False):
        '''Return a Feature object created from the BED line'''
        
        bed_parts = bed_line.strip().split("\t")
        
        if int(bed_parts[1]) > int(bed_parts[2]):
            if not short and bed_parts[5] == "-":
                if not cls.previously_warned_start_bigger_than_end:
                    print "WARNING: Minus strand start values are bigger than end values.\nConverting to appropriate BED format, assuming that the first (column 2) value is 0-based."
                    cls.previously_warned_start_bigger_than_end = True
                start = int(bed_parts[2]) # 1-based already
                end = int(bed_parts[1]) + 1 # convert to 1-based
            else:
                raise MetageneError("Start value is larger than end value.\nBED format requires the start to be less than the end value")
        else:
            start = int(bed_parts[1]) + 1 # convert to 1-based
            end = int(bed_parts[2])

        if short:
            if len(bed_parts) < 4:
                name = "Unknown_name"
            else:
                name = bed_parts[3]
                
            return (Feature(count_method, 
                        metagene_object, 
                        name,  # name
                        bed_parts[0], # chromosome name
                        start, # start 1-based
                        end, # end 1-based
                        ".")) # strand unknown
        else:
            return (Feature(count_method, 
                        metagene_object, 
                        bed_parts[3],  # name
                        bed_parts[0], # chromosome name
                        start, # start 1-based
                        end, # end 1-based
                        bed_parts[5])) # strand
    # end of create_from_bed classmethod


    @classmethod
    def create_from_gff(cls, count_method, metagene_object, gff_line):
        '''Return a Feature object created from the GFF line'''
        
        gff_parts = gff_line.strip().split("\t")
        
        # ensure there are no commas in the name line
        name = ";".join(gff_parts[8].split(","))
        
        # check if start < end
        if int(gff_parts[3]) > int(gff_parts[4]):
            if gff_parts[6] == "-":
                if not cls.previously_warned_start_bigger_than_end:
                    print "WARNING: Minus strand start values are bigger than end values.\nConverting to appropriate GFF format."
                    cls.previously_warned_start_bigger_than_end = True
                start = int(gff_parts[4])
                end = int(gff_parts[3])
            else:
                raise MetageneError("Start value is larger than end value.\nGFF format requires the start to be less than the end value")
        else:
            start = int(gff_parts[3])
            end = int(gff_parts[4])

        return (Feature(count_method, 
                        metagene_object, 
                        name,  # name (potentially messy)
                        gff_parts[0], # chromosome name
                        start, # start 1-based
                        end, # end 1-based
                        gff_parts[6])) # strand
    # end of create_from_gff classmethod

    @classmethod
    def set_chromosome_conversion(cls, tabfile, bam_chromosomes):
        '''Import TAB delimited conversion table for the chromosome names in the 
        feature file (column 0) and in the alignment file (column 1).  Return said
        table as a dictionary with the feature file chromosome names as keys and the 
        alignment file chromosome names as values.'''
        if tabfile != "None":
            try:
                with open(tabfile) as infile:
                    return cls.process_set_chromosome_conversion(infile.read().strip().split("\n"))
            except IOError as err:
                infile.close()
                raise MetageneError("Could not open the conversion table file.")
        else:
            return cls.process_set_chromosomes_conversion(bam_chromosomes, use_bam_chromosomes=True)
    
    @classmethod
    def process_set_chromosome_conversion(cls, tabfile_lines, use_bam_chromosomes=False):
        if use_bam_chromosomes:
            for row in tabfile_lines:
                cls.chromosome_conversion[row] = row
        else:
            for row in tabfile_lines:
                if row[0] != "#": # don't process comments
                    row_parts = row.split("\t")
                    cls.chromosome_conversion[row_parts[0]] = row_parts[1]
        return True
        
# end Feature class    
