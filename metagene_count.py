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

import sys, re, datetime, subprocess
import argparse		# to parse the command line arguments

PROGRAM = "metagene_count.py"
VERSION = "0.1.2"
UPDATED = "140327 JRBT"

class Error(Exception):
    pass

##TODO add more detail to MetageneErrors (maybe subclasses of error types...)
class MetageneError(Error):
    '''For illegal/illogical metagene attributes.'''
    def __init__(self, obj, message):
        self.obj = obj
        self.message = message
    
    def __str__(self):
        return "MetageneError: {}".format(self.message)


class Metagene(object):
    '''Feature defined by an interval of interest and padding around said interval (will be 0-based). Regardless of padding side (left or right), negative padding values will shorten the interval, positive padding values will extend the interval. Minimum metagene length is 1.'''

    ##TODO add functionality for negative paddings!!
    # restrict attributes for each instance
    __slots__ = ['feature_interval','padding', 'length']
    
    def __init__(self, interval, padding_upstream, padding_downstream):
        '''Initiate lengths from the interval and padding.'''
        
        # assign interval; must be INT > 0
        if Metagene.confirm_int(interval, "interval") and int(interval) > 0:
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
           
        # confirm resulting metagene has at least a length of 1  
        self.length = self.feature_interval + self.padding['Upstream'] + self.padding['Downstream']
        if self.length < 1:
            raise MetageneError(length, 
                                "Invalid final length to metagene (interval + \
                                 upstream padding + downstream padding = {})".format(self.length))

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
    # end of confirm_int
            
    @staticmethod
    def test_metagene():
        '''Tests of Metagene class'''
        
        print "\n**** Testing the Metagene class ****\n"
        # create a metagene with different up and downstream padding
        first_metagene = Metagene(10,2,4) # uuiiiiiiiiiidddd = metagene
        if (first_metagene.length == 16):
            print "Calculate correct length ?\tTRUE"
        else:
            print "Calculate correct length ?\t**** FAILED ****"
            print "first_metagene.length is {}".format(first_metagene.length)
        
        try: 
            metagene = Metagene(10.2, 4, 4)
        except MetageneError as err:
            print "Caught non-integer interval ?\tTRUE\t\t{}".format(err)
        else:
            print "Caught non-integer interval ?\t**** FAILED ****"
                    
        try: 
            metagene = Metagene(0,3,3)
        except MetageneError as err:
            print "Caught interval of zero error ?\tTRUE\t\t{}".format(err)
        else:
            print "Caught interval of zero error ?\t**** FAILED ****"
        
        try:
            ##TODO Make it possible to have negative padding!
            metagene = Metagene(10,-3,-2)
        except MetageneError as err:
            print "Caught negative padding error ?\tTRUE\t\t{}".format(err)
        else:
            print "Caught negative padding error ?\t**** FAILED ****"
        
        try: 
            metagene = Metagene(10, 4.3, 4)
        except MetageneError as err:
            print "Caught non-integer padding ?\tTRUE\t\t{}".format(err)
        else:
            print "Caught non-integer padding ?\t**** FAILED ****"
        
        try: 
            metagene = Metagene(10, "four", 4)
        except MetageneError as err:
            print "Caught non-integer padding ?\tTRUE\t\t{}".format(err)
        else:
            print "Caught non-integer padding ?\t**** FAILED ****"
                
        
        print "\n**** End of Testing the Metagene class ****\n"

# end Metagene class

class Feature(Metagene):
    '''Defines genomic features in array format analogous to metagene with padding around the interval of interest (genomic feature).'''
    
    __slots__ = ['name', 'chromosome', 'strand', 'shrink_factor','counts_array', 'position_array']
    # inherits feature_interval, padding, and length from Metagene
    
    def __init__(self, count_method, metagene_object, name, chromosome, start, end, strand, gap_counting=False):
        '''Define a new feature with an interval (represents feature length), up and downstream padding (defined by metagene_object), and genomic (1-based) start and end positions.
        
        Once defined here, the start and end represent the true start and end of the feature.  Therefore, if a - strand (Crick strand) feature the start will be larger than the end.''' 
        
        
        if self.confirm_int(start, "Start") and self.confirm_int(end, "End"):
            start = int(start)
            end = int(end)
            # Define feature-specific metagene where feature_interval respresents 
            # the length of the feature NOT the length of the final metagene interval
            if count_method == 'all':
                Metagene.__init__(self,(end - start + 1), metagene_object.padding['Upstream'], metagene_object.padding['Downstream'])
            else:
                Metagene.__init__(self,1, metagene_object.padding['Upstream'], metagene_object.padding['Downstream'])
            self.shrink_factor = self.feature_interval / float(metagene_object.feature_interval)
                
            self.strand = strand
            if self.strand != "+" and self.strand != "-":
                self.strand = "."
                self.counts_array = { 'unstranded':[] }
            elif gap_counting:
                self.counts_array = { 'ungapped':[], 'gapped':[] }
            else:
                self.counts_array = { 'sense':[], 'antisense':[] }
            
            # initialize counts array with zeros
            self.position_array = []
            # go from left-most genomic position to right-most genomic position adding
            # those values to the position array
            # + strand:   [10,11,12,13,14,15]
            # - strand:   [15,14,13,12,11,10] 
            # so array[0] is always the start and array[-1] is always the end
             
            if self.strand == "-": 
                if count_method == 'start':
                    start = end # set both start and end to the end value (which is really the start)
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
                for orientation in self.counts_array:
                    self.counts_array[orientation].append(0)
                self.position_array.append(p)
            
        self.name = name
        self.chromosome = chromosome           

    # end __init__ function
            
    def __str__(self, counts_only=False):
        output = ""

        if not(counts_only):
            output += "{} at {} on {} strand\n".format(self.name, self.get_chromosome_region(), self.strand)
        
        output += "\t\t\t"
        for i in range(self.padding['Upstream']):
            output += "---up-"
        for i in range(self.feature_interval):
            output += "--int-"
        for i in range(self.padding['Downstream']):
            output += "-down-"
        output += "\n"
        
        output += "{0:15s}:\t".format('Position')
        for i in self.position_array:
            output += "{0:5d},".format(i)
        output = output[:-1] + "\n"
                   
        for orientation in sorted(self.counts_array.keys(), reverse=True):
            output += "{0:15s}:\t".format(orientation)
            for i in self.counts_array[orientation]:
                output += "{0:>5s},".format("{0:3.2f}".format(i))
            output = output[:-1] + "\n"
        return output
                    
                    
    
    def get_chromosome_region(self):
        '''Return position interval for samtools view (chromosome: start-end (1-based))'''
        if self.strand == "-":
            # order for samtools must be smaller to larger position
            return ("{}:{}-{}".format(self.chromosome, self.get_region_end(), self.get_region_start()))
        else:
            return ("{}:{}-{}".format(self.chromosome, self.get_region_start(), self.get_region_end()))
    
    def get_samtools_region(self,chromosome_lengths):
        '''Return a position interval valid for use in samtools view program.'''
        
        if self.strand == "-":
            start = self.get_region_end()
            end = self.get_region_start()
        else:
            start = self.get_region_start()
            end = self.get_region_end()
        
        if start < 1:
            start = 1
            
        if end > chromosome_lengths[self.chromosome]:
            end = chromosome_lengths[self.chromosome]
        
        return ("{}:{}-{}".format(self.chromosome, start, end))
            
    def get_region_start(self):
        '''Upstream most position including padding'''
        return self.position_array[0]
    
    def get_region_end(self):
        '''Downstream most position including padding'''
        return self.position_array[-1]
    
    def print_metagene(self, metagene_length, pretty=False):
        '''Converts counts_array data to finalized metagene profiles for printing'''
        
        final_metagenes = {}
        
        output_line = ""
        if pretty:
            output_line += "{0:15s}\t\t".format(self.name)
            for i in range(self.padding['Upstream']):
                output_line += "---up-"
            for i in range(metagene_length):
                output_line += "--int-"
            for i in range(self.padding['Downstream']):
                output_line += "-down-"
            output_line += "\n"
        
            output_line += "{0:15s}:\t".format('Position')
            for i in range(self.padding['Upstream'], 0, -1):
                output_line += "{0:5d},".format(0-i)
            for i in range(metagene_length):
                output_line += "{0:5d},".format(i)
            for i in range(self.padding['Downstream']):
                output_line += "{0:5d},".format(i + metagene_length)
            output_line = output_line[:-1] + "\n"
            
        for orientation in sorted(self.counts_array, reverse=True):
            # break counts_array into sections -> upstream padding, interval_feature, and downstream padding
            upstream_counts = self.counts_array[orientation][0:self.padding['Upstream']]
            interval_counts = self.counts_array[orientation][self.padding['Upstream'] : self.padding['Upstream'] + self.feature_interval]
            downstream_counts = self.counts_array[orientation][self.padding['Upstream'] + self.feature_interval : len(self.counts_array[orientation])]
        
            # compress (or expand) interval_counts to match the size of the internal metagene
            metagene_interval_counts = self.adjust_to_metagene(interval_counts, self.shrink_factor)
            # ensure metagene_interval is the correct length
            if len(metagene_interval_counts) - metagene_length != 0:
                # if the last (extra) bin is really really small then just remove it as the number is a result of float arithmetic issues
                if metagene_interval_counts[-1] < 0.00001:
                    metagene_interval_counts = metagene_interval_counts[:-1] 
                else:
                    raise MetageneError(metagene_interval_counts, "The metagene interval is not the correct length of {}:\nMetagene interval:\t{}".format(metagene_length, metagene_interval_counts))
            
            if pretty:
                output_line += "{0:15s}:\t".format(orientation)
                for i in upstream_counts:
                    output_line += "{0:>5s},".format("{0:3.2f}".format(i))
                for i in metagene_interval_counts:
                    output_line += "{0:>5s},".format("{0:3.2f}".format(i))
                for i in downstream_counts:
                    output_line += "{0:>5s},".format("{0:3.2f}".format(i))  
                output_line = output_line[:-1] + "\n"
            else:    
                # build output
                output_line += "{},{}".format(self.name, orientation)
                for p in upstream_counts:
                    output_line += ",{0:0.3f}".format(p) # keep 3 decimal places in the outputted float
                for p in metagene_interval_counts:
                    output_line += ",{0:0.3f}".format(p)
                for p in downstream_counts:
                    output_line += ",{0:0.3f}".format(p)
                output_line += "\n"
        
        return output_line.strip() # remove trailing "\n"
         
    def adjust_to_metagene(self, feature_array, shrink_factor):
        '''Expand or collapse the counts data from interval_array into a metagene
        array via the given shrink factor.'''
        ##TODO: convert to using Decimal library for float math!
        
        # uncomment ## lines to have debugging function
        
        metagene_array = []
        
        # initialize metagene_count and remaining_metagene
        metagene_count = 0.0
        remaining_metagene_bin = shrink_factor
        
        ##for loop, bin in enumerate(feature_array): # ensure all data is moved to metagene by looping  through entire interval_array
        for bin in feature_array: # ensure all data is moved to metagene by looping  through entire interval_array
            # Ideally add in units of 1 (1 bin to 1 metagene_array  position) 
            # unless not possible then start dealing with fractional bins
            
            remaining_feature_bin = 1.0 # reset remaining feature for new bin
            
            ##print "\n  Loop {}:".format(loop)
            ##print "    Feature Count :\t{}".format(bin)
            ##print "    Metagene Count:\t{}".format(metagene_count)
            ##print "    Metagene Bin  :\t{}".format(remaining_metagene_bin)
            ##print "    Feature Bin   :\t{}".format(remaining_feature_bin)
            ##i = 0
                      
            while remaining_feature_bin > 0:
                # keeping adding from this bin until its empty
                
                ##i += 1
                ##print "    While loop {}:".format(i)
                ##print "      Feature Count :\t{}".format(bin)
                ##print "      Metagene Count:\t{}".format(metagene_count)
                ##print "      Metagene Bin  :\t{}".format(remaining_metagene_bin)
                ##print "      Feature Bin   :\t{}".format(remaining_feature_bin)
                                
                if remaining_feature_bin <= remaining_metagene_bin:
                    ##print "      Add Remaining Feature Bin:\t{}".format(bin * remaining_feature_bin)    
                    
                    # add entire feature to metagene
                    metagene_count += (bin * remaining_feature_bin)
                    # adjust bin counters
                    remaining_metagene_bin -= remaining_feature_bin 
                    remaining_feature_bin = 0
                else:
                    ##print "      Add Remaining Metagene Bin:\t{}".format(bin * remaining_metagene_bin)
                    
                    # add entire remaining_metagene_bin amount of feature to metagene
                    metagene_count += (bin * remaining_metagene_bin)
                    # adjust bin counters
                    remaining_feature_bin -= remaining_metagene_bin
                    remaining_metagene_bin = 0
                
                ##print "Remaining_metagene_bin:\t{}".format(remaining_metagene_bin)                    
                # check to see if new metagene bin is ready to be added to the metagene_array
                if remaining_metagene_bin == 0:
                    ##print "      Add Count to Metagene Array:\t{}".format(metagene_count)
                    metagene_array.append(metagene_count)
                    metagene_count = 0.0
                    remaining_metagene_bin = shrink_factor
            # end of while loop through current feature bin
        # end of for loop through feature array
        
        if metagene_count != 0.0:
            # print out final metagene that was missed
            ##print "      Add Count to Metagene Array:\t{}".format(metagene_count)
            metagene_array.append(metagene_count)
            metagene_count = 0.0
                                            
        ##print "\n  Final Metagene:\t{}".format(metagene_array)
        return metagene_array
    
    def count_read(self, read_object, count_method):
        '''Add a read object to the sense or antisense counts_array. Requires strand
        options of "+", "-" or "."
        
        Only stranded reads (+ or -) can be counted on stranded features.
        
        Unstranded Features will count in the + direction, but ignore read strand. '''

        ##TODO: option to require full fit of read within the feature 
        ##TODO: option to count gapped vs non-gapped instead of strand (or in addition to strands)
        
        # determine orientation (and if countable)
        if self.strand == ".":
            orientation = 'unstranded'
        elif read_object.strand != ".":
            if self.strand == read_object.strand:
                orientation = 'sense'
            else:
                orientation = 'antisense'
        else: 
            raise MetageneError(read_object.strand, "Can not count unstranded reads on stranded features.")
        
        # confirm that they are on the same chromosome
        if self.chromosome == read_object.chromosome:
            # get positions from read to potentially count
            positions_to_count = []
        
            if count_method == 'start':
                positions_to_count.append(read_object.get_start_position())
            elif count_method == 'end':
                positions_to_count.append(read_object.get_end_position())
            elif count_method == 'all':
                positions_to_count = read_object.position_array
            else:
                raise MetageneError(count_method, "Unrecognizable counting method.  Valid options are 'start', 'end', and 'all'")
        
            for p in positions_to_count:
                # make sure it overlaps with the Feature
                if p in self.position_array:
                    self.counts_array[orientation][self.position_array.index(p)] += (float(read_object.abundance) / float(read_object.mappings))
        
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
                        ".")) # strand
        else:
            return (Feature(count_method, 
                        metagene_object, 
                        bed_parts[3],  # name
                        chromosome_conversion[bed_parts[0]], # alignment style chromosome name
                        int(bed_parts[1]) + 1, # start 1-based
                        int(bed_parts[2]), # end 1-based
                        bed_parts[5])) # strand

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

    @staticmethod
    def test_feature():
        '''Tests of Feature class'''
        
        print "\n**** Testing the Feature class ****\n"  
        
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
                print "  with Metagene:\t{}".format(metagene)
                print "  with chromosome conversions:\t{}".format(chromosome_converter)
            else:
                metagene = Metagene(1,4,2)
                print "  with Metagene:\t{}".format(metagene)
                print "  with chromosome conversions:\t{}".format(chromosome_converter)
        
        
            # create feature from BED line
            try:
                bedline = "{}\t{}\t{}\t{}\t{}\t{}\n".format(1,20,40,"first",44,"+")
                print " with BED line:\t{}".format(bedline.strip())
                feature1 = Feature.create_from_bed(method, metagene, bedline, chromosome_converter)
                if str(feature1.position_array) != correct_features['bed'][method]: 
                    print "  Create Feature from BED line ?\t**** FAILED *****"
                    print "  Desired positions:\t{}".format(correct_features['bed'][method])
                    print "  Created positions:\t{}".format(feature1.position_array)          
            except MetageneError as err:
                print "  Create Feature from BED line ?\t**** FAILED ****"
            else:
                print "  Create Feature from BED line ?\tTRUE\t\t{}".format(feature1.get_chromosome_region())
            
            # create feature from GFF line
            try:
                gffline = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(2,"test","gene",10,39,".","-",".","second")
                print " with GFF line:\t{}".format(gffline.strip())
                feature2 = Feature.create_from_gff(method,metagene, gffline, chromosome_converter)
                if str(feature2.position_array) != correct_features['gff'][method]: 
                    print "  Create Feature from GFF line ?\t**** FAILED *****"
                    print "  Desired positions:\t{}".format(correct_features['gff'][method])
                    print "  Created positions:\t{}".format(feature2.position_array)
            
            except MetageneError as err:
                print "  Create Feature from GFF line ?\t**** FAILED ****"
            else:
                print "  Create Feature from GFF line ?\tTRUE\t\t{}".format(feature2.get_chromosome_region())
       
       

        ##TODO finish complete testing of Feature class
        print "\n##TODO finish complete testing of Feature class\n"
                
        print "\n**** End of Testing the Feature class ****\n"
    
    @staticmethod
    def test_counting_methods():
        '''Specifically test the counting and summary methods (all points after
        Feature creation).'''
        
        chromosome_converter = {"1":"chr1", "2":"chr2"}
        
        for method in ['all','start','end']:
            if method == 'all':
                metagene = Metagene(10,4,2)
                print "  with Metagene:\t{}".format(metagene)
                print "  with chromosome conversions:\t{}".format(chromosome_converter)
            else:
                metagene = Metagene(1,4,2)
                print "  with Metagene:\t{}".format(metagene)
                print "  with chromosome conversions:\t{}".format(chromosome_converter)
        
            print "\nTesting feature_count option: ****{}****".format(method)
            feature_line = "{}\t{}\t{}\t{}\t{}\t{}\n".format(1,20,40,"first",44,"+")
            feature1 = Feature.create_from_bed(method, metagene, feature_line, chromosome_converter)
            print "Feature:\t{}".format(feature1.position_array)
        
            reads = []
            reads.append(Read("chr1", "+", 3, 1, [10,11,12,13,14,15,16,17,18]))
            reads.append(Read("chr1", "-", 1, 2, [23,24,25,26,27,28,29,30,31,32]))
            reads.append(Read("chr1", "+", 4, 2, [30,31,32,33,34,40,41]))
            reads.append(Read("chr1", "-", 1, 1, [42,43,44,45,46,47,48,49,50]))
            reads.append(Read("chr1", "+", 10, 1, [51,52,53,54,55]))
            reads.append(Read("chr2", "+", 10, 1, [18,19,20,21,22,23,24,25]))
        
            # starting count
            #TODO: figure out wierd 11th interval bin in feature_method 'end' + count_method 'all'
            for count_method in ['all','start','end']:
                print "\nTesting count_method option: ****{}****".format(count_method)

                output = "\n{0:15s}:\t".format('Feature')
                for i in feature1.position_array:
                    output += "{0:3d},".format(i)
                print output[:-1] 
                   
                for orientation in ['sense','antisense']:
                    output = "{0:15s}:\t".format(orientation)
                    for i in feature1.counts_array[orientation]:
                        output += "{0:3d},".format(i)
                    print output[:-1]
                
                for r in reads:
                    print "\nRead ({},{},NA={},NH={},count={})  :\t{}".format(r.chromosome,r.strand,r.abundance, r.mappings, float(r.abundance) / r.mappings, r.position_array)
                    feature1.count_read(r, count_method) 
                    output = "{0:15s}:\t".format('Feature')
                    for i in feature1.position_array:
                        output += "{0:5d},".format(i)
                    print output[:-1] 
                    
                    for orientation in ['sense','antisense']:
                        output = "{0:15s}:\t".format(orientation)
                        for i in feature1.counts_array[orientation]:
                            output += " {0:2.2f},".format(i)
                        print output[:-1]
        
                print "\nFinal Metagene:\n{}\n".format(feature1.print_metagene(metagene.feature_interval))
                feature1 = Feature.create_from_bed(method, metagene, feature_line, chromosome_converter) # zero out counter for next round
        
        try:
            unstranded_read = Read("chr1", ".", 10, 1, [18,19,20,21,22,23,24,25])
            feature1.count_read(unstranded_read, 'all')
        except MetageneError as err:
            print "Caught unstranded read on stranded count ?\tTRUE\t\t".format(err)
        else:
            print "Caught unstranded read on stranded count ?\t**** FAILED ****"
        
        try:
            unstranded_read = Read("chr1", ".", 10, 1, [18,19,20,21,22,23,24,25])
            feature1.count_read(unstranded_read, 'all')
        except MetageneError as err:
            print "Allowed unstranded read on unstranded count ?\t**** FAILED ****\t\t".format(err)
        else:
            print "Allowed unstranded read on unstranded count ?\tTRUE"
            
        # try on feature with out a strand
        metagene = Metagene(10,4,2)
        print "  with Metagene:\t{}".format(metagene)
        print "  with chromosome conversions:\t{}".format(chromosome_converter)
        feature_line = "{}\t{}\t{}\t{}\t{}\t{}\n".format(1,20,40,"first",44,"u")
        
        
        for count_method in ['all','start','end']:
            feature1 = Feature.create_from_bed('all',metagene, feature_line, chromosome_converter)
            print "\nCount_method = {}".format(count_method)
            for orientation in feature1.counts_array:
                print "{}:\t{}".format(orientation, feature1.counts_array[orientation])
            
            for r in reads:
                print "\nRead:\t{}".format(r.position_array)
                feature1.count_read(r, count_method) 
                for orientation in feature1.counts_array:
                    print "{}:\t{}".format(orientation, feature1.counts_array[orientation])
        
            print "\nFinal Metagene:\n{}\n".format(feature1.print_metagene(metagene.feature_interval))
      

    # end test_counting_methods function
    
    @staticmethod
    def test_adjust_to_metagene_method():
    
        print "\nTesting adjust_to_metagene:\n"
        print "  with feature array: [16, 8, 24, 4] expanding to metagene array of length 8 (shrink = 0.5)"
        metagene = Metagene(8,2,2)
        chromosome_converter = {"1":"chr1", "2":"chr2"}
        
        feature_line = "{}\t{}\t{}\t{}\t{}\t{}\n".format(1,20,24,"first",44,"+")
        feature1 = Feature.create_from_bed(metagene, feature_line, chromosome_converter)
    
        print "Resulting Metagene:\t{}".format(feature1.adjust_to_metagene([16,8,24,4], 0.5))
    
        print "  with feature array: [6,8,6,2,4,4,2,4,24,8] collapsing to metagene array of length 4 (shrink = 2.5)"
        metagene = Metagene(4,2,2)
        chromosome_converter = {"1":"chr1", "2":"chr2"}
        
        feature_line = "{}\t{}\t{}\t{}\t{}\t{}\n".format(1,20,30,"first",44,"+")
        feature2 = Feature.create_from_bed(metagene, feature_line, chromosome_converter)
    
        print "Resulting Metagene:\t{}".format(feature2.adjust_to_metagene([6,8,6,2,4,4,2,4,24,8], 2.5))

        print "\nthe unpretty float math trials:"
        print "  with feature array: [2.5,4,(10.0/3),10,11,7.3,4] collapsing to metagene array of length 4 (shrink = 7/4 = 1.75)"
        metagene = Metagene(4,2,2)
        chromosome_converter = {"1":"chr1", "2":"chr2"}
        
        feature_line = "{}\t{}\t{}\t{}\t{}\t{}\n".format(1,20,27,"first",44,"+")
        feature2 = Feature.create_from_bed(metagene, feature_line, chromosome_converter)
    
        print "  Resulting Metagene: {}".format(feature2.adjust_to_metagene([2.5,4,(10.0/3),10,11,7.3,4], 1.75))
        print "  Expected Metagene:  [5.500, 9.333, 17.825, 9.475]" # 3 decimal place accuracy
        
        print "\nthe unpretty float math trials:"
        print "  with feature array: [2.5,4,(10.0/3),10,11,7.3,4] collapsing to metagene array of length 3 (shrink = 7/3 ~ 2.333)"
        metagene = Metagene(3,2,2)
        chromosome_converter = {"1":"chr1", "2":"chr2"}
        
        feature_line = "{}\t{}\t{}\t{}\t{}\t{}\n".format(1,20,27,"first",44,"+")
        feature2 = Feature.create_from_bed(metagene, feature_line, chromosome_converter)
    
        print "  Resulting Metagene: {}".format(feature2.adjust_to_metagene([2.5,4,(10.0/3),10,11,7.3,4], (7.0 / 3)))
        print "  Expected Metagene:  [7.611, 19.555, 14.967]" # 3 decimal place accuracy
    # end of test_adjust_to_metagene_method function

                 
# end Feature class    


class Read():
    '''Defines positions of a read'''
    
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
        return "Read at {0}:{1}-{2} on {3} strand; counts for {4:1.3f}:\t\t{5}".format(self.chromosome, self.get_start_position(), self.get_end_position(), self.strand, float(self.abundance)/self.mappings, str(self.position_array))
    
    def get_start_position(self):
        return self.position_array[0]
        
    def get_end_position(self):
        return self.position_array[-1]
        
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
    def create_from_sam(cls, sam_line, chromosome_conversion, extract_abundance=False, unique=False):
        '''Create a Read object from a bamfile line, requires that the chromosome 
        is in the chromosome_conversion dictionary'''
        
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
           
        # assign mappings
        if unmapped_flag: # non-mapping read
            raise MetageneError(sam_line, "Can not create a read for an unmapped sequence")    
        elif unique:
            mappings = 1
        # try to extract mappings from NH:i:## tag
        else:
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
        
        return Read(chromosome, strand, abundance, mappings, positions)
    # end of create_from_sam
    
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
    def test_read():
        '''Tests of Read class'''
        
        print "\n**** Testing the Read class ****\n"
        
        chromosome_converter = {"1":"chr1", "2":"chr2"}
        print "  with chromosome conversions:\t{}\n".format(chromosome_converter)   
        
        try:
            samline1 = "read_1	24	chr1	200	255	3S2M4N3M2X3M	*	0	0	bbbbbbbbbbbbb	bbbbbbbbbbbbb	XA:i:0	MD:Z:40	NH:i:50  NA:i:10"
            print "with SAM line:\t{}".format(samline1)
            read1 = Read.create_from_sam(samline1, chromosome_converter)
        except MetageneError as err:
            print "Create Read object from SAM line ?\t**** FAILED ****"
        else:
            print "Create Read object from SAM line ?\tTRUE"
            print "  {}".format(read1)
            print "  Location:\t\t{}:{}-{}".format(read1.chromosome, read1.get_start_position(), read1.get_end_position() )
            print "  Strand:\t\t{}".format(read1.strand)
            print "  Abundance:\t\t{}".format(read1.abundance)
            print "  Number Mappings:\t{}".format(read1.mappings)
            print
            
        try:
            samline1 = "read_1	0	chr1	200	255	*	*	0	0	bbbbbbbbbbbbb	bbbbbbbbbbbbb	XA:i:0	MD:Z:40	NH:i:50  NA:i:10"
            print "with SAM line:\t{}".format(samline1)
            read1 = Read.create_from_sam(samline1, chromosome_converter)
        except MetageneError as err:
            print "Create Read object from SAM line ?\t**** FAILED ****"
        else:
            print "Create Read object from SAM line ?\tTRUE"
            print "  {}".format(read1)
            print "  Location:\t\t{}:{}-{}".format(read1.chromosome, read1.get_start_position(), read1.get_end_position() )
            print "  Strand:\t\t{}".format(read1.strand)
            print "  Abundance:\t\t{}".format(read1.abundance)
            print "  Number Mappings:\t{}".format(read1.mappings)
            print
                
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
                elif re.search('\A\S+\t\d+\t\d+\t\S+\t\S+\t[+.-]\s+\Z',line) != None:
                    counts['BED'] += 1
                elif re.search('\A\S+\t\S+\t\S+\t\d+\t\d+\t\S+\t[+.-]\t\S+\t\S+\s+\Z', line) != None:
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
                
                print feature
                
                # pull out sam file lines; it is important to use Feature.get_samtools_region(chromosome_lengths) rather
                # than Feature.get_chromosome_region() because only the first ensures that the interval does not
                # extend beyond the length of the chromosome which makes samtools view return no reads
                (runPipe_worked, sam_sample) = runPipe(['samtools view {} {}'.format(arguments.alignment,feature.get_samtools_region(chromosomes))])
                if runPipe_worked:
                    for samline in sam_sample:
                        if len(samline) > 0:
                            # create Read feature
                            read = Read.create_from_sam(samline, chromosome_conversion_table, unique=not(arguments.extract_mappings))
                            print read
                            # count read
                            feature.count_read(read, arguments.count_method)
                            print Feature.__str__(feature,counts_only=True)
                    # output the resulting metagene
                    with open("{}.metagene_counts.csv".format(arguments.output_prefix), 'a') as output_file:
                        print "Finished counting reads..."
                        print Feature.__str__(feature,counts_only=True)
                        print feature.print_metagene(metagene.feature_interval, pretty=True)
                        output_file.write("{}\n".format(feature.print_metagene(metagene.feature_interval)))
                    
                else:
                    print sam_sample
                    raise MetageneError(sam_sample, "Could not pull chromosomal region {} for feature {} from BAM file {}.".format(feature.get_chromosome_region(), feature.name, arguments.alignment))
      
   

    
if __name__ == "__main__":
    # testing classes and methods
    #Metagene.test_metagene()
    #Feature.test_feature()
    #Read.test_read()
    
    #Feature.test_counting_methods()
    #Feature.test_adjust_to_metagene_method()
    
    # actual run...
    try:
        metagene_count()    
    except MetageneError as err:
        print "\n{}\n".format(err)
        print "Aborting metagene analysis..."
        sys.exit()
    

