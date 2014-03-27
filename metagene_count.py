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
            print "Caught non-integer interval ?\tTRUE"
            print "\t{}".format(err)
        else:
            print "Caught non-integer interval ?\t**** FAILED ****"
                    
        try: 
            metagene = Metagene(0,3,3)
        except MetageneError as err:
            print "Caught interval of zero error ?\tTRUE"
            print "\t{}".format(err)
        else:
            print "Caught interval of zero error ?\t**** FAILED ****"
        
        try:
            ##TODO Make it possible to have negative padding!
            metagene = Metagene(10,-3,-2)
        except MetageneError as err:
            print "Caught negative padding error ?\tTRUE"
            print "\t{}".format(err)
        else:
            print "Caught negative padding error ?\t**** FAILED ****"
        
        try: 
            metagene = Metagene(10, 4.3, 4)
        except MetageneError as err:
            print "Caught non-integer padding ?\tTRUE"
            print "\t{}".format(err)
        else:
            print "Caught non-integer padding ?\t**** FAILED ****"
        
        try: 
            metagene = Metagene(10, "four", 4)
        except MetageneError as err:
            print "Caught non-integer padding ?\tTRUE"
            print "\t{}".format(err)
        else:
            print "Caught non-integer padding ?\t**** FAILED ****"
                
        
        print "\n**** End of Testing the Metagene class ****\n"

# end Metagene class

class Feature(Metagene):
    '''Defines genomic features in array format analogous to metagene with padding around the interval of interest (genomic feature).'''
    
    __slots__ = ['name', 'chromosome', 'strand', 'shrink_factor','counts_array', 'position_array']
    # inherits feature_interval, padding, and length from Metagene
    
    def __init__(self, metagene_object, name, chromosome, start, end, strand, gap_counting=False):
        '''Define a new feature with an interval (represents feature length), up and downstream padding (defined by metagene_object), and genomic (1-based) start and end positions.
        
        Once defined here, the start and end represent the true start and end of the feature.  Therefore, if a - strand (Crick strand) feature the start will be larger than the end.''' 
        
        
        if self.confirm_int(start, "Start") and self.confirm_int(end, "End"):
            start = int(start)
            end = int(end)
            # Define feature-specific metagene where feature_interval respresents 
            # the length of the feature NOT the length of the final metagene interval
            Metagene.__init__(self,(end - start + 1), metagene_object.padding['Upstream'], metagene_object.padding['Downstream'])
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
                region_start = start - self.padding['Downstream'] # start is really end
                region_end = end + self.padding['Upstream'] # end is really start
                positions = range(region_start, region_end + 1) # inclusive list
                positions.reverse()
            else:
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
    
    def get_chromosome_region(self):
        '''Return position interval for samtools view (chromosome: start-end (1-based))'''
        if self.strand == "-":
            # order for samtools must be smaller to larger position
            return ("{}:{}-{}".format(self.chromosome, self.get_region_end(), self.get_region_start()))
        else:
            return ("{}:{}-{}".format(self.chromosome, self.get_region_start(), self.get_region_end()))
            
    def get_region_start(self):
        '''Upstream most position including padding'''
        return self.position_array[0]
    
    def get_region_end(self):
        '''Downstream most position including padding'''
        return self.position_array[-1]
    
    def print_metagene(self):
        '''Converts counts_array data to finalized metagene profiles for printing'''
        
        final_metagenes = {}
        
        output_line = ""
        for orientation in self.counts_array:
            # break counts_array into sections -> upstream padding, interval_feature, and downstream padding
            upstream_counts = self.counts_array[orientation][0:self.padding['Upstream']]
            interval_counts = self.counts_array[orientation][self.padding['Upstream'] : self.padding['Upstream'] + self.feature_interval]
            downstream_counts = self.counts_array[orientation][self.padding['Upstream'] + self.feature_interval : len(self.counts_array[orientation])]
        
            # compress (or expand) interval_counts to match the size of the internal metagene
            metagene_interval_counts = self.adjust_to_metagene(interval_counts, self.shrink_factor)
            
            # build output
            output_line += "{},{}".format(self.name, orientation)
            for p in upstream_counts:
                output_line += ",{0:0.3f}".format(p) # keep 3 decimal places in the outputted float
            for p in metagene_interval_counts:
                output_line += ",{0:0.3f}".format(p)
            for pi in downstream_counts:
                output_line += ",{0:0.3f}".format(p)
            output_line += "\n"
        
        return output_line.strip() # remove trailing "\n"
         
    def adjust_to_metagene(self, feature_array, shrink_factor):
        '''Expand or collapse the counts data from interval_array into a metagene
        array via the given shrink factor.'''
        
        metagene_array = []
        
        # initialize metagene_count and remaining_metagene
        metagene_count = 0.0
        remaining_metagene_bin = shrink_factor
        
        for loop, bin in enumerate(feature_array): # ensure all data is moved to metagene by looping  through entire interval_array
            # Ideally add in units of 1 (1 bin to 1 metagene_array  position) 
            # unless not possible then start dealing with fractional bins
            
            remaining_feature_bin = 1.0 # reset remaining feature for new bin
            
            print "\n  Loop {}:".format(loop)
            print "    Feature Count :\t{}".format(bin)
            print "    Metagene Count:\t{}".format(metagene_count)
            print "    Metagene Bin  :\t{}".format(remaining_metagene_bin)
            print "    Feature Bin   :\t{}".format(remaining_feature_bin)
            i = 0
                      
            while remaining_feature_bin > 0:
                # keeping adding from this bin until its empty
                i += 1
                print "    While loop {}:".format(i)
                print "      Feature Count :\t{}".format(bin)
                print "      Metagene Count:\t{}".format(metagene_count)
                print "      Metagene Bin  :\t{}".format(remaining_metagene_bin)
                print "      Feature Bin   :\t{}".format(remaining_feature_bin)
                                
                if remaining_feature_bin <= remaining_metagene_bin:
                    print "      Add Remaining Feature Bin:\t{}".format(bin * remaining_feature_bin)    
                    # add entire feature to metagene
                    metagene_count += (bin * remaining_feature_bin)
                    # adjust bin counters
                    remaining_metagene_bin -= remaining_feature_bin 
                    remaining_feature_bin = 0
                else:
                    print "      Add Remaining Metagene Bin:\t{}".format(bin * remaining_metagene_bin)
                    # add entire remaining_metagene_bin amount of feature to metagene
                    metagene_count += (bin * remaining_metagene_bin)
                    # adjust bin counters
                    remaining_feature_bin -= remaining_metagene_bin
                    remaining_metagene_bin = 0
                
                print "Remaining_metagene_bin:\t{}".format(remaining_metagene_bin)                    
                # check to see if new metagene bin is ready to be added to the metagene_array
                if remaining_metagene_bin == 0:
                    print "      Add Count to Metagene Array:\t{}".format(metagene_count)
                    metagene_array.append(metagene_count)
                    metagene_count = 0.0
                    remaining_metagene_bin = shrink_factor
            # end of while loop through current feature bin
        # end of for loop through feature array
        
        if metagene_count != 0.0:
            # print out final metagene that was missed
            print "      Add Count to Metagene Array:\t{}".format(metagene_count)
            metagene_array.append(metagene_count)
            metagene_count = 0.0
                                            
        print "\n  Final Metagene:\t{}".format(metagene_array)
   
    
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
    def create_from_bed(cls, metagene_object, bed_line, chromosome_conversion):
        '''Return a Feature object created from the BED line'''
        
        bed_parts = bed_line.strip().split("\t")
        
        return (Feature(metagene_object, 
                        bed_parts[3],  # name
                        chromosome_conversion[bed_parts[0]], # alignment style chromosome name
                        int(bed_parts[1]) + 1, # start 1-based
                        int(bed_parts[2]), # end 1-based
                        bed_parts[5])) # strand

    @classmethod
    def create_from_gff(cls, metagene_object, gff_line, chromosome_conversion):
        '''Return a Feature object created from the GFF line'''
        
        gff_parts = gff_line.strip().split("\t")
        
        return (Feature(metagene_object, 
                        gff_parts[8],  # name (potentially messy)
                        chromosome_conversion[gff_parts[0]], # alignment style chromosome name
                        int(gff_parts[3]), # start 1-based
                        int(gff_parts[4]), # end 1-based
                        gff_parts[6])) # strand

    @staticmethod
    def test_feature():
        '''Tests of Feature class'''
        
        print "\n**** Testing the Feature class ****\n"
        
        metagene = Metagene(10,4,2)
        chromosome_converter = {"1":"chr1", "2":"chr2"}
        print "  with Metagene:\t{}".format(metagene)
        print "  with chromosome conversions:\t{}".format(chromosome_converter)     
        
        # create feature from BED line
        try:
            bedline = "{}\t{}\t{}\t{}\t{}\t{}\n".format(1,20,40,"first",44,"+")
            print "\nwith BED line:\t{}".format(bedline)
            feature1 = Feature.create_from_bed(metagene, bedline, chromosome_converter)           
        except MetageneError as err:
            print "  Create Feature from BED line ?\t**** FAILED ****"
        else:
            print "  Create Feature from BED line ?\tTRUE"   
            print "  BEDline feature:\t{}".format(feature1)
            print "  Positions ({}): {}".format(len(feature1.position_array),feature1.position_array)
            print "  Chromosome region:\t{}".format(feature1.get_chromosome_region())
            
        # create feature from GFF line
        try:
            gffline = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(2,"test","gene",10,39,".","-",".","second")
            print "\nwith GFF line:\t{}".format(gffline)
            feature2 = Feature.create_from_gff(metagene, gffline, chromosome_converter)
            
        except MetageneError as err:
            print "  Create Feature from GFF line ?\t**** FAILED ****"
        else:
            print "  Create Feature from GFF line ?\tTRUE"   
            print "  GFFline feature:\t{}".format(feature2)
            print "  Positions ({}): {}".format(len(feature2.position_array),feature2.position_array)
            print "  Chromosome region:\t{}".format(feature2.get_chromosome_region())
       
        

        ##TODO finish complete testing of Feature class
        print "\n##TODO finish complete testing of Feature class\n"
                
        print "\n**** End of Testing the Feature class ****\n"
    
    @staticmethod
    def test_counting_methods():
        '''Specifically test the counting and summary methods (all points after
        Feature creation).'''
        
        metagene = Metagene(10,4,2)
        chromosome_converter = {"1":"chr1", "2":"chr2"}
        
        feature_line = "{}\t{}\t{}\t{}\t{}\t{}\n".format(1,20,40,"first",44,"+")
        feature1 = Feature.create_from_bed(metagene, feature_line, chromosome_converter)
        print "Feature:\t{}".format(feature1.position_array)
        
        reads = []
        reads.append(Read("chr1", "+", 3, 1, [10,11,12,13,14,15,16,17,18]))
        reads.append(Read("chr1", "-", 1, 2, [23,24,25,26,27,28,29,30,31,32]))
        reads.append(Read("chr1", "+", 4, 2, [30,31,32,33,34,40,41]))
        reads.append(Read("chr1", "-", 1, 1, [42,43,44,45,46,47,48,49,50]))
        reads.append(Read("chr1", "+", 10, 1, [51,52,53,54,55]))
        reads.append(Read("chr2", "+", 10, 1, [18,19,20,21,22,23,24,25]))
        
        # starting count
        print "Count_method = 'all'"
        for orientation in feature1.counts_array:
            print "{}:\t{}".format(orientation, feature1.counts_array[orientation])
            
        for r in reads:
            print "\nRead:\t{}".format(r.position_array)
            feature1.count_read(r, 'all') 
            for orientation in feature1.counts_array:
                print "{}:\t{}".format(orientation, feature1.counts_array[orientation])
        
        print "\nCount_method = 'start'"
        # reset feature1 to fresh (normally the feature is lost after counting)
        feature1 = Feature.create_from_bed(metagene, feature_line, chromosome_converter)
        
        for r in reads:
            feature1.count_read(r, 'start') 
        
        for orientation in feature1.counts_array:
            print "{}:\t{}".format(orientation, feature1.counts_array[orientation])
        
        print "\nCount_method = 'end'"
        # reset feature1 to fresh (normally the feature is lost after counting)
        feature1 = Feature.create_from_bed(metagene, feature_line, chromosome_converter)
        
        for r in reads:
            feature1.count_read(r, 'end') 
        
        for orientation in feature1.counts_array:
            print "{}:\t{}".format(orientation, feature1.counts_array[orientation])
        
        try:
            unstranded_read = Read("chr1", ".", 10, 1, [18,19,20,21,22,23,24,25])
            feature1.count_read(unstranded_read, 'all')
        except MetageneError as err:
            print "Caught unstranded read on stranded count ?\tTRUE"
            print err
        else:
            print "Caught unstranded read on stranded count ?\t**** FAILED ****"
        
            
        # try on feature with out a strand
        feature_line = "{}\t{}\t{}\t{}\t{}\t{}\n".format(1,20,40,"first",44,"u")
        feature1 = Feature.create_from_bed(metagene, feature_line, chromosome_converter)
        
        print "\nCount_method = 'all'"
        for orientation in feature1.counts_array:
            print "{}:\t{}".format(orientation, feature1.counts_array[orientation])
            
        for r in reads:
            print "\nRead:\t{}".format(r.position_array)
            feature1.count_read(r, 'all') 
            for orientation in feature1.counts_array:
                print "{}:\t{}".format(orientation, feature1.counts_array[orientation])
        
        print "\nCount_method = 'start'"
        # reset feature1 to fresh (normally the feature is lost after counting)
        feature1 = Feature.create_from_bed(metagene, feature_line, chromosome_converter)
        
        for r in reads:
            feature1.count_read(r, 'start') 
        
        for orientation in feature1.counts_array:
            print "{}:\t{}".format(orientation, feature1.counts_array[orientation])
        
        print "\nCount_method = 'end'"
        # reset feature1 to fresh (normally the feature is lost after counting)
        feature1 = Feature.create_from_bed(metagene, feature_line, chromosome_converter)
        
        for r in reads:
            feature1.count_read(r, 'end') 
        
        for orientation in feature1.counts_array:
            print "{}:\t{}".format(orientation, feature1.counts_array[orientation])
        
        try:
            unstranded_read = Read("chr1", ".", 10, 1, [18,19,20,21,22,23,24,25])
            feature1.count_read(unstranded_read, 'all')
        except MetageneError as err:
            print "Allowed unstranded read on unstranded count ?\t**** FAILED ****"
            print err
        else:
            print "Allowed unstranded read on unstranded count ?\tTRUE"
    # end test_counting_methods function
    
    @staticmethod
    def test_adjust_to_metagene_method():
    
        print "\nTesting adjust_to_metagene:\n"
        print "  with feature array: [16, 8, 24, 4] expanding to metagene array of length 8 (shrink = 0.5)"
        metagene = Metagene(8,2,2)
        chromosome_converter = {"1":"chr1", "2":"chr2"}
        
        feature_line = "{}\t{}\t{}\t{}\t{}\t{}\n".format(1,20,24,"first",44,"+")
        feature1 = Feature.create_from_bed(metagene, feature_line, chromosome_converter)
    
        feature1.adjust_to_metagene([16,8,24,4], 0.5)
    
        print "  with feature array: [6,8,6,2,4,4,2,4,24,8] collapsing to metagene array of length 4 (shrink = 2.5)"
        metagene = Metagene(4,2,2)
        chromosome_converter = {"1":"chr1", "2":"chr2"}
        
        feature_line = "{}\t{}\t{}\t{}\t{}\t{}\n".format(1,20,30,"first",44,"+")
        feature2 = Feature.create_from_bed(metagene, feature_line, chromosome_converter)
    
        feature2.adjust_to_metagene([6,8,6,2,4,4,2,4,24,8], 2.5)

        print "\nthe unpretty float math trials:"
        print "  with feature array: [2.5,4,(10.0/3),10,11,7.3,4] collapsing to metagene array of length 4 (shrink = 7/4 = 1.75)"
        metagene = Metagene(4,2,2)
        chromosome_converter = {"1":"chr1", "2":"chr2"}
        
        feature_line = "{}\t{}\t{}\t{}\t{}\t{}\n".format(1,20,27,"first",44,"+")
        feature2 = Feature.create_from_bed(metagene, feature_line, chromosome_converter)
    
        feature2.adjust_to_metagene([2.5,4,(10.0/3),10,11,7.3,4], 1.75)
        print "  Expected Metagene: [5.500, 9.333, 17.825, 9.475]" # 3 decimal place accuracy
        
        print "\nthe unpretty float math trials:"
        print "  with feature array: [2.5,4,(10.0/3),10,11,7.3,4] collapsing to metagene array of length 3 (shrink = 7/3 ~ 2.333)"
        metagene = Metagene(3,2,2)
        chromosome_converter = {"1":"chr1", "2":"chr2"}
        
        feature_line = "{}\t{}\t{}\t{}\t{}\t{}\n".format(1,20,27,"first",44,"+")
        feature2 = Feature.create_from_bed(metagene, feature_line, chromosome_converter)
    
        feature2.adjust_to_metagene([2.5,4,(10.0/3),10,11,7.3,4], (7.0 / 3))
        print "  Expected Metagene: [7.611, 19.555, 14.967]" # 3 decimal place accuracy
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
        return str(self.position_array)
    
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
    def create_from_sam(cls, sam_line, chromosome_conversion, unique=False):
        '''Create a Read object from a bamfile line, requires that the chromosome 
        is in the chromosome_conversion dictionary'''
        
        sam_parts = sam_line.split("\t")

        # assign chromosome
        if sam_parts[2] not in chromosome_conversion.values():
            raise MetageneError(sam_parts[2], "Read chromosome is not in the analysis set")
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
        try: 
            abundance = int(re.search('NA:i:(\d+)', sam_line).group(1))
        except AttributeError:
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
                        default = "{}.metagene.".format(date))
    
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

    parser.add_argument("-c","--chromosome_names",
                        help = "Chromosome conversion file (feature_chromosome {tab} alignment_chromosome)",
                        metavar = 'TAB',
                        default = "None")
    
    parser.add_argument("--count_splicing",
                        help = "Count reads as spliced or unspliced (ignores strand)",
                        action = 'store_true')
    
    arguments = parser.parse_args()
    
    # adjust internal_size if only the start or end of the feature will be counted
    if arguments.feature_count != 'all':
        arguments.interal_size = 1
    
                         
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


def metagene_count():
    '''Chain of command for metagene_count analysis.'''
    ##TODO: finish the flow of command
    arguments = get_arguments()
    #print "Current arguments: \n{}".format(arguments)
    
    # confirm BAM file and extract chromosome sizes
    chromosomes = get_chromosome_sizes(arguments.alignment)
    #print "Current chromosomes: \n{}".format(chromosomes)
    
    # create chromosome conversion dictionary for feature (GFF/BED) to alignment (BAM)
    if arguments.chromosome_names != "None":
        # create dict of 
        chromosome_conversion_table = get_chromosome_conversions(arguments.chromosome_names, chromosomes.keys())
    else:
        # dummy dict of BAM-defined chromosome names in both cases
        chromosome_conversion_table = {}
        for c in chromosomes:
            chromosome_conversion_table[c] = c
        
    #print "Current conversion table: \n{}".format(chromosome_conversion_table)
    
    # define the metagene array shape (left padding, start, internal, end, right padding)
    # metagene = padding ---- internal region ---- padding 
    try:
        metagene = Metagene(arguments.interval_size, arguments.padding, arguments.padding)
    except MetageneError as err:
        print err
        sys.exit()
    
    # testing metagene object    
    print metagene.padding_upstream, metagene.feature_interval, metagene.padding_downstream, metagene.length
    print metagene.get_interval_start()
    print metagene.get_interval_end()
    
    feature_format = "gff" # get_format(arguments.feature) 
      
    # for each feature
    for feature_line in read_features(arguments.feature):
        feature = Feature.create(feature_format, metagene, feature_line, chromosome_conversion_table)# create an empty metagene
        
        # define genomic positions
        feature_region = define_genomic(feature, arguments.padding)
        
        for read in get_reads(feature_region):
            # read = array of genomic positions to tally back to the feature_metagene
            feature_metagene = add_read(read, feature_metagene)
            
        # output feature_metagene          
   

    
if __name__ == "__main__":
    get_arguments()
    Metagene.test_metagene()
    Feature.test_feature()
    Read.test_read()
    
    Feature.test_counting_methods()
    Feature.test_adjust_to_metagene_method()
    
#    metagene_count()    
   
    

