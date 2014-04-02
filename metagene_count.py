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

PROGRAM = "metagene_count.py"
VERSION = "0.1.2"
UPDATED = "140402 JRBT"

class Error(Exception):
    pass

##TODO: add more detail to MetageneErrors (maybe subclasses of error types...)
class MetageneError(Error):
    '''For errors specific to metagene analysis.'''
    def __init__(self, obj, message):
        self.obj = obj
        self.message = message
    
    
    def __str__(self):
        return "MetageneError: {}".format(self.message)


class Metagene(object):
    '''A Metagene is an object representing an interval of interest (length >= 1) 
    with padding on either side (lengths >= 0).  
    
    Attributes:
        feature_interval : interval of interest (int >= 1)
        padding          : dictionary of padding values (int >= 0)
          'Upstream'     : padding upstream (to left) of feature_interval
          'Downstream'   : padding downstream (to right) of feature_interval
        length           : length of entire metagene object
    
    Static Method:
        confirm_int(value, name)
        test()
        '''

##TODO: add functionality for negative paddings!!


    # restrict attributes for each instance
    __slots__ = ['feature_interval','padding', 'length']
    
    
    def __init__(self, interval, padding_upstream, padding_downstream):
        '''Initiate lengths from the interval and padding.'''
        
        # assign interval; must be INT >= 1
        if Metagene.confirm_int(interval, "interval") and int(interval) >= 1:
            self.feature_interval = int(interval)
        else:
            raise MetageneError(interval, "Minimum interval length is 1.")
        
        # assign padding; must be INT >= 0  
        self.padding = {'Upstream':0, 'Downstream':0} # set defaults
        for pad in [(padding_upstream, "Upstream"), (padding_downstream, "Downstream")]:
            if Metagene.confirm_int(pad[0], "{} Padding".format(pad[1])) and int(pad[0]) >= 0:
                self.padding[pad[1]] = int(pad[0])
            else:
                raise MetageneError(pad[0], "Padding values must be positive")   
             
        self.length = self.feature_interval + self.padding['Upstream'] + self.padding['Downstream']

    # end __init__ function
    
    
    def __str__(self):
        return "(Upstream--Interval--Downstream) = {}--{}--{} ({})".format(self.padding['Upstream'], self.feature_interval, self.padding['Downstream'], self.length)
        
    
    @staticmethod
    def confirm_int(value, name):
        try: 
            if float(value) % 1 == 0:
                return True
            else:
                raise MetageneError(value, "{} ({}) must be an integer".format(name, value))         
        except ValueError:
            raise MetageneError(value, "{} ({}) must be an integer".format(name, value))
    # end of Metagene.confirm_int staticmethod
            
    
    @staticmethod
    def test():
        '''Test error handling of Metagene class'''
        
        print "\n**** Testing the Metagene class ****\n"
        # can a metagene be created
        print "\tCommand:\tmetagene = Metagene(11,4,4)"
        metagene = Metagene(11,4,4)
        if str(metagene) == "(Upstream--Interval--Downstream) = 4--11--4 (19)":
            print "PASSED\tCreated a valid metagene ?\t\t{}".format(metagene)
        else:
            print "**FAILED**\tCreated a valid metagene ?"
            print "\tExpected Metagene:\t(Upstream--Interval--Downstream) = 4--11--4 (19)"
            print "\tCreated  Metagene:\t{}".format(metagene)
            
        # create a metagene with different up and downstream padding
        print "\tCommand:\tmetagene = Metagene(10,2,4) "
        metagene = Metagene(10,2,4) 
        if str(metagene) == "(Upstream--Interval--Downstream) = 2--10--4 (16)":
            print "PASSED\tCreated a valid metagene with differing padding ?\t\t{}".format(metagene)
        else:
            print "**FAILED**\tCreated a valid metagene with differing padding ?"
            print "\tExpected Metagene:\t(Upstream--Interval--Downstream) = 2--10--4 (16)"
            print "\tCreated  Metagene:\t{}".format(metagene)
        
        # catch errors from non-desired inputs
        tests = [ ((10.2, 4, 4), "Caught non-integer float interval ?"),
                  (("ten", 4, 4), "Caught non-integer string interval ?"),
                  ((0, 3, 3), "Caught interval of zero error ?"),
##TODO: Make it possible to have negative padding
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

# end Metagene class

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
##TODO: clean up division system for counts_array dict
        counts_array     : dictionary of counting objects, value array length is
                           the same as the length of the position_array;
                           gapped/ungapped/allreads can be added to strand information
          'sense'        : counts for sense reads (both read and feature strand is the same)
          'antisense'    : counts for antisense reads (opposite of sense)
          'unstranded'   : counts for when strand is not present or being ignored
          
          'gapped'       : counts for reads with gaps relative to the feature
          'ungapped'     : counts for reads without gaps relative to the feature
          'allreads'     : counts for reads with or without gaps
          
    
    Functions:
        get_chromosome_region()
        get_samtools_region(chromosome_lengths_dictionary)
        print_metagene(pretty=False)   
        adjust_to_metagene(feature_array, verbose=False)
        count_read(read_object, count_method, count_gaps)
    
    Class methods:
        create(file_format, count_method, metagene_object, feature_line, chromosome_conversion_table)
        create_from_bed(count_method, metagene_object, feature_line, chromosome_conversion_table, short=False)
        create_from_gff(count_method, metagene_object, feature_line, chromosome_conversion_table)
    
    Static methods:
        test()
    
        '''
    
    __slots__ = ['name', 'chromosome', 'strand', 'metagene_length','counts_array', 'position_array']
    # inherits feature_interval, padding, and length from Metagene
    
    
    def __init__(self, count_method, metagene_object, name, chromosome, start, end, strand, gap_counting=False):
        '''Not normally called directly; use Feature.create(file_format, count_method, 
        metagene_object, feature_line, chromosome_conversion_table) to call indirectly.
        
        Define a new feature with an interval (represents feature length), 
        up and downstream padding (defined by metagene_object), and genomic 
        (1-based) start and end positions.
        
        Once defined here, the start and end represent the true start and end of
        the feature.  Therefore, if a - strand (Crick strand) feature the start
        will be larger than the end.''' 
        
               
        if Metagene.confirm_int(start, "Start") and Metagene.confirm_int(end, "End"):
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
                
            for p in positions: 
                self.position_array.append(p)            
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
    
    
    def get_samtools_region(self,chromosome_lengths_dictionary):
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
            raise MetageneError(chromosome_lengths_dictionary, 
                    "Could not find chromosome {} in length dictionary".format(self.chromosome))
        
        return ("{}:{}-{}".format(self.chromosome, start, end))
    # end of Feature.get_samtools_region function    
    
    
    def print_metagene(self, pretty=False):
        '''Converts counts_array data to finalized metagene profiles for printing
        
        Standard printing is in comma-delimited lines for input into metagene_analysis.py
        Pretty printing (pretty=True) gives a human readable, if potentially super long, version'''
        
        final_metagenes = {}
        
        output = ""
        if pretty: # add metagene schematic and position numbers (relative to feature start as zero)
            output += "{0:15s}\t\t".format(self.name)
            for i in range(self.padding['Upstream']):
                output += "---up-"
            for i in range(self.metagene_length):
                output += "--int-"
            for i in range(self.padding['Downstream']):
                output += "-down-"
            output += "\n"
        
            output += "{0:15s}:\t".format('Position')
            for i in range(self.padding['Upstream'], 0, -1):
                output += "{0:5d},".format(0-i)
            for i in range(self.metagene_length):
                output += "{0:5d},".format(i)
            for i in range(self.padding['Downstream']):
                output+= "{0:5d},".format(i + self.metagene_length)
            output = output[:-1] + "\n"
        
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
##TODO: convert to using Decimal library for float math!
        
        metagene_array = []
        
        # initialize metagene_count and remaining_metagene
        metagene_count = 0.0
        shrink_factor = len(feature_array) / float(self.metagene_length)
        remaining_metagene_bin = shrink_factor
        
        loop_index = 0 # loop index for verbose option
        
        for bin in feature_array: # ensure all data is moved to metagene by looping  through entire interval_array
            # Ideally add in units of 1 (1 bin to 1 metagene_array  position) 
            # unless not possible then start dealing with fractional bins
            
            remaining_feature_bin = 1.0 # reset remaining feature for new bin
            
            if verbose:
                print "\n  Loop {}:".format(loop_index)
                print "    Feature Count :\t{}".format(bin)
                print "    Metagene Count:\t{}".format(metagene_count)
                print "    Metagene Bin  :\t{}".format(remaining_metagene_bin)
                print "    Feature Bin   :\t{}".format(remaining_feature_bin)
                loop_index += 1
                while_index = 0 # while loop index
                      
            while remaining_feature_bin > 0:
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
                    metagene_count += (bin * remaining_feature_bin)
                    # adjust bin counters 
                    remaining_metagene_bin -= remaining_feature_bin 
                    remaining_feature_bin = 0
                else:
                    if verbose:
                        print "      Add Remaining Metagene Bin:\t{}".format(bin * remaining_metagene_bin)
                    
                    # add entire remaining_metagene_bin amount of feature to metagene
                    metagene_count += (bin * remaining_metagene_bin)
                    # adjust bin counters
                    remaining_feature_bin -= remaining_metagene_bin
                    remaining_metagene_bin = 0
                
                if verbose:
                    print "Remaining_metagene_bin:\t{}".format(remaining_metagene_bin)                    
                # check to see if new metagene bin is ready to be added to the metagene_array
                if remaining_metagene_bin == 0:
                    if verbose:
                        print "      Add Count to Metagene Array:\t{}".format(metagene_count)
                    metagene_array.append(metagene_count)
                    metagene_count = 0.0
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
    
       
    def count_read(self, read_object, count_method, count_gaps=False):
        '''Add a read object to the sense or antisense counts_array. Requires strand
        options of "+", "-" or "."
        
        Only stranded reads (+ or -) can be counted on stranded features.
        
        Unstranded Features will count in the + direction, but ignore read strand. '''

##TODO - maybe?: option to require full fit of read within the feature 
        
        # determine orientation (and if countable)
        if self.strand == ".":
            subset = 'unstranded'
        elif read_object.strand != ".":
            if self.strand == read_object.strand:
                subset = 'sense'
            else:
                subset = 'antisense'
        else: 
            raise MetageneError(read_object.strand, "Can not count unstranded reads on stranded features.")
        
        # determine gap status
        if count_gaps:
            if abs(read_object.position_array[0] - read_object.position_array[-1]) + 1 > len(read_object.position_array):
                # if calculated length > actual length then gapped
                subset += ":gapped"
            else:
                subset += ":ungapped"            
        else:
            subset += ":allreads"
        
        # confirm that they are on the same chromosome
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
                raise MetageneError(count_method, "Unrecognizable counting method.  Valid options are 'start', 'end', and 'all'")
        
            for p in positions_to_count:
                # make sure it overlaps with the Feature
                if p in self.position_array:
                    self.counts_array[subset][self.position_array.index(p)] += (float(read_object.abundance) / float(read_object.mappings))        
    # end of count_read function        
    
    
    #******** creating Feature objects from diffent feature file formats (eg BED and GFF) ********#
    @classmethod
    def create(cls, format, count_method, metagene, feature_line, chromosome_conversion_table):
        '''Calls individual creation methods based on format.'''
        if format == "GFF":
            return Feature.create_from_gff(count_method, metagene, feature_line, chromosome_conversion_table)
        elif format == "BED":
            return Feature.create_from_bed(count_method, metagene, feature_line, chromosome_conversion_table)
        elif format == "BED_SHORT":
            return Feature.create_from_bed(count_method, metagene, feature_line, chromosome_conversion_table, short=True)
            
            
    @classmethod
    def create_from_bed(cls, count_method, metagene_object, bed_line, chromosome_conversion, short=False):
        '''Return a Feature object created from the BED line'''
        
        bed_parts = bed_line.strip().split("\t")
        
        if short:
            if len(bed_parts) < 4:
                name = "Unknown_name"
            else:
                name = bed_parts[3]
                
            return (Feature(count_method, 
                        metagene_object, 
                        name,  # name
                        chromosome_conversion[bed_parts[0]], # alignment style chromosome name
                        int(bed_parts[1]) + 1, # start 1-based
                        int(bed_parts[2]), # end 1-based
                        ".")) # strand unknown
        else:
            return (Feature(count_method, 
                        metagene_object, 
                        bed_parts[3],  # name
                        chromosome_conversion[bed_parts[0]], # alignment style chromosome name
                        int(bed_parts[1]) + 1, # start 1-based
                        int(bed_parts[2]), # end 1-based
                        bed_parts[5])) # strand
    # end of create_from_bed classmethod


    @classmethod
    def create_from_gff(cls, count_method, metagene_object, gff_line, chromosome_conversion):
        '''Return a Feature object created from the GFF line'''
        
        gff_parts = gff_line.strip().split("\t")
        
        # ensure there are no commas in the name line
        name = ";".join(gff_parts[8].split(","))
        return (Feature(count_method, 
                        metagene_object, 
                        name,  # name (potentially messy)
                        chromosome_conversion[gff_parts[0]], # alignment style chromosome name
                        int(gff_parts[3]), # start 1-based
                        int(gff_parts[4]), # end 1-based
                        gff_parts[6])) # strand
    # end of create_from_gff classmethod
    

    @staticmethod
    def test():
        '''Tests of Feature class'''
        
        print "\n**** Testing the Feature class ****\n"  
        print "**** Testing Feature creation ****\n"
        
        correct_features = {'bed':{}, 'gff':{} }
        correct_features['bed']['all'] = "[17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42]"
        correct_features['bed']['start'] = "[17, 18, 19, 20, 21, 22, 23]"
        correct_features['bed']['end'] = "[36, 37, 38, 39, 40, 41, 42]"
        correct_features['gff']['all'] = "[43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8]"
        correct_features['gff']['start'] = "[43, 42, 41, 40, 39, 38, 37]"
        correct_features['gff']['end'] = "[14, 13, 12, 11, 10, 9, 8]"
         
        chromosome_converter = {"1":"chr1", "2":"chr2"}
                                
        for method in ['all','start','end']:
            print "\nTesting feature_count option: ****{}****".format(method)
            
            if method == 'all':
                metagene = Metagene(10,4,2)
                print "\t  with Metagene:\t{}".format(metagene)
                print "\t  with chromosome conversions:\t{}".format(chromosome_converter)
            else:
                metagene = Metagene(1,4,2)
                print "\t  with Metagene:\t{}".format(metagene)
                print "\t  with chromosome conversions:\t{}".format(chromosome_converter)
        
        
            # create feature from BED line
            try:
                bedline = "{}\t{}\t{}\t{}\t{}\t{}\n".format(1,20,40,"first",44,"+")
                print "\t  with BED line:\t{}".format(bedline.strip())
                feature1 = Feature.create_from_bed(method, metagene, bedline, chromosome_converter)
                if str(feature1.position_array) != correct_features['bed'][method]: 
                    print "**FAILED**\t  Create Feature from BED line ?"
                    print "\t  Desired positions:\t{}".format(correct_features['bed'][method])
                    print "\t  Created positions:\t{}".format(feature1.position_array)          
            except MetageneError as err:
                print "**FAILED**\t  Create Feature from BED line ?"
            else:
                print "PASSED\t  Create Feature from BED line ?\t\t{}".format(feature1.get_chromosome_region())
            
            # create feature from GFF line
            try:
                gffline = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(2,"test","gene",10,39,".","-",".","second")
                print "\t  with GFF line:\t{}".format(gffline.strip())
                feature2 = Feature.create_from_gff(method,metagene, gffline, chromosome_converter)
                if str(feature2.position_array) != correct_features['gff'][method]: 
                    print "**FAILED**\t  Create Feature from GFF line ?\t**FAIL**"
                    print "\t  Desired positions:\t{}".format(correct_features['gff'][method])
                    print "\t  Created positions:\t{}".format(feature2.position_array)
            except MetageneError as err:
                print "**FAILED**\t  Create Feature from GFF line ?"
            else:
                print "PASSED\t  Create Feature from GFF line ?\t\t{}".format(feature2.get_chromosome_region())

##TODO finish complete testing of Feature class
        print "\n##TODO finish testing of Feature class creation\n"
        
        print "\n**** Testing counting and maniputlation ****\n"
        
        expected = { 'all':{}, 'start':{}, 'end':{} }    
        #  Positions in metagene:                           17    18     19   20  21-22,23-24,25-26,27-28,29-30,31-32,33-34,35-36,37-38,39-40,  41,   42
        expected['all'] = {   'all':"first,sense:allreads,3.000,3.000,0.000,0.000,0.000,0.000,0.000,0.000,2.000,4.000,4.000,0.000,0.000,2.000,2.000,0.000\nfirst,antisense:allreads,0.000,0.000,0.000,0.000,0.000,1.000,1.000,1.000,1.000,1.000,0.000,0.000,0.000,0.000,0.000,1.000",
                            'start':"first,sense:allreads,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,2.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000\nfirst,antisense:allreads,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.500,0.000,0.000,0.000,0.000,0.000,0.000",
                              'end':"first,sense:allreads,0.000,3.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,2.000,0.000\nfirst,antisense:allreads,0.000,0.000,0.000,0.000,0.000,0.500,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,1.000" }
        #  Positions in metagene:                           17    18    19    20   [21]   22    23
        expected['start'] = { 'all':"first,sense:allreads,3.000,3.000,0.000,0.000,0.000,0.000,0.000\nfirst,antisense:allreads,0.000,0.000,0.000,0.000,0.000,0.000,0.500",
                            'start':"first,sense:allreads,0.000,0.000,0.000,0.000,0.000,0.000,0.000\nfirst,antisense:allreads,0.000,0.000,0.000,0.000,0.000,0.000,0.000",
                              'end':"first,sense:allreads,0.000,3.000,0.000,0.000,0.000,0.000,0.000\nfirst,antisense:allreads,0.000,0.000,0.000,0.000,0.000,0.000,0.500" }
        #  Positions in metagene:                           36    37    38    39   [40]   41    42
        expected['end'] = {   'all':"first,sense:allreads,0.000,0.000,0.000,0.000,2.000,2.000,0.000\nfirst,antisense:allreads,0.000,0.000,0.000,0.000,0.000,0.000,1.000",
                            'start':"first,sense:allreads,0.000,0.000,0.000,0.000,0.000,0.000,0.000\nfirst,antisense:allreads,0.000,0.000,0.000,0.000,0.000,0.000,0.000",
                              'end':"first,sense:allreads,0.000,0.000,0.000,0.000,0.000,2.000,0.000\nfirst,antisense:allreads,0.000,0.000,0.000,0.000,0.000,0.000,1.000" }                     
                               
        metagene = { 'all':Metagene(10,4,2),
                     'start':Metagene(1,4,2),
                     'end':Metagene(1,4,2)}  
                                   
        for method in ['all','start','end']:
            if method == 'all':
                print "\t  with Metagene:\t{}".format(metagene[method])
                print "\t  with chromosome conversions:\t{}".format(chromosome_converter)
            else:
                print "\t  with Metagene:\t{}".format(metagene[method])
                print "\t  with chromosome conversions:\t{}".format(chromosome_converter)
        
            print "\nTesting feature_count option: ****{}****".format(method)
            feature_line = "{}\t{}\t{}\t{}\t{}\t{}\n".format(1,20,40,"first",44,"+")
            feature1 = Feature.create_from_bed(method, metagene[method], feature_line, chromosome_converter)
            print "\tFeature:\t{}".format(feature1.position_array)
        
            reads = []
            reads.append(Read("chr1", "+", 3, 1, [10,11,12,13,14,15,16,17,18]))
            reads.append(Read("chr1", "-", 1, 2, [23,24,25,26,27,28,29,30,31,32]))
            reads.append(Read("chr1", "+", 4, 2, [30,31,32,33,34,40,41]))
            reads.append(Read("chr1", "-", 1, 1, [42,43,44,45,46,47,48,49,50]))
            
            reads.append(Read("chr1", "+", 10, 1, [51,52,53,54,55]))
            reads.append(Read("chr2", "+", 10, 1, [18,19,20,21,22,23,24,25]))
        
            # starting count
            for count_method in ['all','start','end']:
                print "\nTesting count_method option: ****{}****".format(count_method)

                output = "{}\n".format(feature1)
                                
                for r in reads:
                    output += "{}\n".format(r)
                    feature1.count_read(r, count_method) 
                    output += "{}\n".format(feature1)
                
                output += feature1.print_metagene(pretty=True)
                if str(feature1.print_metagene()).strip() == str(expected[method][count_method]).strip():
                    print "PASSED\tCreated correct metagene with feature method {} and count method {} ?".format(method, count_method)
                else:
                    print "**FAILED**\tCreated correct metagene with feature method {} and count method {} ?".format(method, count_method)
                    print "\tExpected:\n{}".format(expected[method][count_method])
                    print "\tActual  :\n{}".format(feature1.print_metagene())
                    print "\tSummary of run:\n{}".format(output)
                feature1 = Feature.create_from_bed(method, metagene[method], feature_line, chromosome_converter) # zero out counter for next round
        
        try:
            unstranded_read = Read("chr1", ".", 10, 1, [18,19,20,21,22,23,24,25])
            feature1.count_read(unstranded_read, 'all')
        except MetageneError as err:
            print "PASSED\tCaught unstranded read on stranded count ?\t\t".format(err)
        else:
            print "**FAILED**\tCaught unstranded read on stranded count ?"
        
        try:
            feature_line = "{}\t{}\t{}\t{}\t{}\t{}\n".format(1,20,40,"first",44,".")
            feature1 = Feature.create_from_bed(method, metagene[method], feature_line, chromosome_converter)
            unstranded_read = Read("chr1", ".", 10, 1, [18,19,20,21,22,23,24,25])
            feature1.count_read(unstranded_read, 'all')
        except MetageneError as err:
            print "**FAILED**\tAllowed unstranded read on unstranded count ?\t\t".format(err)
        else:
            print "PASSED\tAllowed unstranded read on unstranded count ?"   
        
        
        print "\n**** Testing adjust_to_metagene ****\n"
        
        chromosome_converter = {"1":"chr1", "2":"chr2"}
        
        # ((metagene_tupple),(feature_tupple),expected_result_string, message_string)
        tests = [((8,2,2),(16,8,24,4),'8.000,8.000,4.000,4.000,12.000,12.000,2.000,2.000',"Expand to metagene ?"),
                 ((4,2,2),(6,8,6,2,4,4,2,4,24,8),'17.000,9.000,8.000,34.000',"Contract to metagene ?"),
                 ((4,2,2),(2.5,4,(10.0/3),10,11,7.3,4),'5.500,9.333,17.825,9.475', "Contract with messy floats ?"),
                 ((3,2,2),(2.5,4,(10.0/3),10,11,7.3,4),'7.611,19.555,14.967', "Contract with other messy floats ?")]
                 
        
        for t in tests:
            metagene = Metagene(*t[0])
            print "\t{}".format(metagene)
            feature_line = "{}\t{}\t{}\n".format(1,0,len(t[1]))
            feature = Feature.create_from_bed('all', metagene, feature_line, chromosome_converter, short=True)
            adjusted_feature = ""
            for f in feature.adjust_to_metagene(t[1]):
                adjusted_feature += "{0:0.3f},".format(f)
            if adjusted_feature[:-1] == t[2]:
                print "PASSED\t{}".format(t[3])
            else:
                print "**FAILED**\t{}".format(t[3])
                print "\tExpected:\t{}".format(t[2])
                print "\tActual  :\t{}".format(adjusted_feature[:-1])
                print "\tOriginal:\t{}".format(feature.adjust_to_metagene(t[1]))


        print "\n**** End of Testing the Feature class ****\n"
    # end of Feature.test method
    
# end Feature class    


class Read():
    '''A Read is an object representing a sequencing read as positions along 
    the chromosome to which it aligns.  
    
    Attributes:
        chromosome     : chromosome 
        strand         : strand of read '+', '-', or '.'
        position_array : array of chromosome (1-based) positions covered by the
                         read; ordered from 5'-most read position (start; index = 0) 
                         to 3'-most read position (end; index = -1); gaps in 
                         read alignment are represented by gaps position_array
                         values, but not gaps in the array itself
        abundance      : number of times the read is present; 1 or extract from NA:i:## tag
        mappings       : number of potential alignment positions; 1 or extract from NH:i:## tag
        has_mappings   : boolean - if True then mappings value came from NH:i:## tag (or user)
                                   if False then mappings value set to 1 by script
    
    
        '''
    
    __slots__ = ['chromosome','strand','position_array','abundance','mappings','has_mappings']
    
    def __init__(self, chromosome, strand, abundance, mappings, positions):
        self.position_array = []
            
        self.strand = strand
        if self.strand != "+" and self.strand != "-":
            self.strand = "."
        elif self.strand == "-":
            positions.reverse()
            
        for p in positions:
            self.position_array.append(p)
        
        self.chromosome = chromosome
        
        if Read.confirm_int(abundance, "Abundance") and int(abundance) >= 0:
            self.abundance = int(abundance)
        else:
            raise MetageneError(abundance, "Abundance must be greater than or equal to 0")
            
        if mappings == "Unknown":
            self.has_mappings = False
            self.mappings = 1
        elif Read.confirm_int(mappings, "Mappings") and int(mappings) > 0:
            self.has_mappings = True
            self.mappings = int(mappings)
        else:
            raise MetageneError(mappings, "Mappings must be greater than or equal to 0")
    # End of __init__
    
    
    def __str__(self):
        return "Read at {0}:{1}-{2} on {3} strand; counts for {4:1.3f}:\t\t{5}".format(self.chromosome, self.position_array[0], self.position_array[-1], self.strand, float(self.abundance)/self.mappings, str(self.position_array))
    
    
    @classmethod
    def create_from_sam(cls, sam_line, chromosome_conversion, extract_abundance=False, unique=False):
        '''Create a Read object from a bamfile line, requires that the chromosome 
        is in the chromosome_conversion dictionary
        
        extract_abundance = True --> extract abundance value from NA:i:## tag
                          = False --> assign value of 1 (uncollapsed reads were used)
                          
        unique = True --> force mappings to 1; 
                          only unique alignments represented (or desire to avoid hits-normalization)
               = False --> extract number of alignments from NH:i:## tag'''
        
        sam_parts = sam_line.split("\t")

        # assign chromosome
        if sam_parts[2] not in chromosome_conversion.values():
            raise MetageneError(sam_parts[2], "Read chromosome {} is not in the analysis set".format(sam_parts[2]))
        else:
            chromosome = sam_parts[2]
        
        # parse bitwise flag
        (multiple_flag, unmapped_flag, reversed_flag) = Read.parse_sam_bitwise_flag(int(sam_parts[1]))
        
## TODO Add multiple-mapping functionality 
        if multiple_flag:
            raise MetageneError(sam_line, "Can not parse sam lines with mapping in multiple segments")
           
        
        if not unmapped_flag: # only process mapping reads
            # assign mappings
            if unique:
                mappings = 1
            else: # try to extract mappings from NH:i:## tag
                try:
                    mappings = int(re.search('NH:i:(\d+)', sam_line).group(1))
                except AttributeError:
                    mappings = "Unknown"
        
            # assign abundance either from NA:i:## tag or as 1 (default)
            if extract_abundance:
                try: 
                    abundance = int(re.search('NA:i:(\d+)', sam_line).group(1))
                except AttributeError:
                    raise MetageneError(sam_line, "Could not extract the abundance tag")
            else:
                abundance = 1
           
            # assign strand and positions
            if reversed_flag: # Crick or Minus strand
                strand = "-" 
            else: # Watson or Plus strand
                strand = "+"  
        
            # create genomic positions for read (start, cigar_string, sequence)
            positions = Read.build_positions(int(sam_parts[3]), sam_parts[5], sam_parts[9])
        
            return (True, Read(chromosome, strand, abundance, mappings, positions))
        else:
            return (False, "Non-aligning read")
    # end of create_from_sam
    
   
    @classmethod
    def parse_sam_bitwise_flag(cls, decimal_flag):
        '''Pulls bitwise flags for multiple mapping (bit 0x1), unmapped (bit 0x4), 
        and reversed sequences (bit 0x10) according to the samtools manual.
        
        binary string: .... 0000 0000
        multi-mapping flag          ^ string[-1]
        unmapped flag             ^   string[-3]
        reversed flag          ^      string[-5] (ignoring space)           
        '''
        
## TODO Part of add multiple-mapping functionality -- extend bitwise flag parsing to the multiple mapping flags
        binary_flag = bin(decimal_flag)[2:].zfill(8) # removes "0b" prefix and fills from left out to 8 positions
        
        if int(binary_flag[-1]) == 1:
            multiple = True
        else:
            multiple = False
        
        if int(binary_flag[-3]) == 1:
            unmapped = True
        else:
            unmapped = False
        
        if int(binary_flag[-5]) == 1:
            reverse = True
        else:
            reverse = False
        
        return (multiple,unmapped,reverse)


    @classmethod
    def build_positions(cls, start, cigar, seq):    
        '''Parse through a cigar string to return the genomic positions that are
        covered by the read.  Starts at the left-most 1-based start location'''
        
        array = []
        position = start
        
        # sometime the cigar value is "*", in which case assume a perfect match
        if cigar == "*":
            for i in range(len(seq)):
                array.append(position)
                position += 1
            return array
        
        # cigar_codes adapted using information 
        # from samtools specs (samtools.github.io/hts-specs/SAMv1.pdf)
        cigar_codes = { 'M':True, # alignment match (can be either sequence match or mismatch)
                        'I':False, # insertion to the reference
                        'D':False, # deletion from the reference
                        'N':False, # skipped region from the reference
                        'S':False, # soft clipping (clipped sequences present in SEQ)
                        'H':False, # hard clipping (clipped sequences NOT present in SEQ)
                        'P':False, # padding (silent deletion from padded reference)
                        '=':True, # sequence match
                        'X':True } # sequence mismatch
                        
        advance_position = { 'M':True, # alignment match (can be either sequence match or mismatch)
                             'I':False, # insertion to the reference
                             'D':True, # deletion from the reference
                             'N':True, # skipped region from the reference
                             'S':False, # soft clipping (clipped sequences present in SEQ) <-- start position at first aligned base (Heng et al 2009 Bioinformatics)
                             'H':False, # hard clipping (clipped sequences NOT present in SEQ)
                             'P':True, # padding (silent deletion from padded reference) <-- sort of like a skipped region
                             '=':True, # sequence match
                             'X':True } # sequence mismatch
        
        nucleotides = re.findall('(\d+)',cigar)
        codes = re.split('\d+',cigar)[1:]

        # loop through nucleotide values
        for i in range(len(nucleotides)):
            # iterate nt times adding 1 to start each time
            for j in range(int(nucleotides[i])):
                if cigar_codes[codes[i]]:
                    array.append(position)
                if advance_position[codes[i]]:
                    position += 1
        return array  
          
    # end of build_positions
    
    
    @staticmethod
    def confirm_int(value, name):
        try: 
            if float(value) % 1 == 0:
                return True
            else:
                raise MetageneError(value, "{} ({}) must be an integer".format(name, value))         
        except ValueError:
            raise MetageneError(value, "{} ({}) must be an integer".format(name, value))
    # end of confirm_int     
    
    
    @staticmethod
    def test():
        '''Tests of Read class'''
        
        print "\n**** Testing the Read class ****\n"
        
        chromosome_converter = {"1":"chr1", "2":"chr2"}
        print "\t  with chromosome conversions:\t{}\n".format(chromosome_converter)   
        
        try:
            samline1 = "read_1	24	chr1	200	255	3S2M4N3M2X3M	*	0	0	bbbbbbbbbbbbb	bbbbbbbbbbbbb	XA:i:0	MD:Z:40	NH:i:50  NA:i:10"
            print "\twith SAM line:\t{}".format(samline1)
            read1 = Read.create_from_sam(samline1, chromosome_converter)
        except MetageneError as err:
            print "**FAILED**\tCreate Read object from SAM line ?\t\t{}".format(err)
        else:
            print "PASSED\tCreate Read object from SAM line ?"
            print "\t  {}".format(read1[1])

            
        try:
            samline1 = "read_1	0	chr1	200	255	*	*	0	0	bbbbbbbbbbbbb	bbbbbbbbbbbbb	XA:i:0	MD:Z:40	NH:i:50  NA:i:10"
            print "\twith SAM line:\t{}".format(samline1)
            read1 = Read.create_from_sam(samline1, chromosome_converter)
        except MetageneError as err:
            print "**FAILED**\tCreate Read object from SAM line ?\t\t{}".format(err)
        else:
            print "PASSED\tCreate Read object from SAM line ?"
            print "\t  {}".format(read1[1])

                
##TODO finish complete testing of Feature class
        print "\n##TODO finish complete testing of Read class\n"
                
        print "\n**** End of Testing the Read class ****\n"  
         
# end of Read class
    
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
    
    arguments = parser.parse_args()
    
    # adjust internal_size if only the start or end of the feature will be counted
    if arguments.feature_count != 'all':
        arguments.interval_size = 1
                             
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
                            (created_read, read) = Read.create_from_sam(samline, chromosome_conversion_table, arguments.extract_abundance, not(arguments.extract_mappings))
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
    # testing classes and methods
    Metagene.test()
    Feature.test()
    Read.test()
    
    # actual run...
    try:
        metagene_count()    
    except MetageneError as err:
        print "\n{}\n".format(err)
        print "Aborting metagene analysis..."
        sys.exit()
    

