#!/usr/bin/python
'''The first step of metagene_analysis, metagene_count.py compiles read abundance
over genomic features to create the input for metagene_windows.py. Please 
see README for full details and examples.

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
UPDATED = "140326 JRBT"

class Error(Exception):
    pass

##TODO fix error handling for MetageneErrors: print error and exit!
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
    
    def __init__(self, metagene_object, name, chromosome, start, end, strand):
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
            # initialize counts array with zeros
            self.counts_array = []
            self.position_array = []
            # go from left-most genomic position to right-most genomic position adding
            # those values to the position array
            # + strand:   [10,11,12,13,14,15]
            # - strand:   [15,14,13,12,11,10] 
            # so array[0] is always the start and array[-1] is always the end
            
            if strand == "-": 
                region_start = start - self.padding['Downstream'] # start is really end
                region_end = end + self.padding['Upstream'] # end is really start
                positions = range(region_start, region_end + 1) # inclusive list
                positions.reverse()
            else:
                region_start = start - self.padding['Upstream'] 
                region_end = end + self.padding['Downstream'] 
                positions = range(region_start, region_end + 1) # inclusive list
                
            for p in positions: 
                self.counts_array.append(0)
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
    
    def in_counts_array(self, genomic_position):
        if genomic_position in self.position_array:
            return True
        else:
            return False
            
    def get_counts_array_index(self, genomic_position):
        '''Returns the position in the counts_array represented by 1-based genomic position'''
        if genomic_position in self.position_array:
            return (self.position_array.index(genomic_position))
        else:
            raise MetageneError(genomic_position, "Position is not in the feature metagene counting region")
    
    ##TODO method to shrink interval to metagene size
    ##TODO method to print metagene
    ##TODO method to count in feature
    
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
         
# end Feature class    


class Read():
    '''Defines positions of a read'''
    
    __slots__ = ['chromosome','strand','position_array','abundance','mappings','has_mappings']
    
    def __init__(self, chromosome, start, strand, abundance, mappings, positions):
        if Read.confirm_int(start, "Start"):
            self.position_array = []
            
            self.strand = strand
            if self.strand == "-":
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
    
    @classmethod
    def build_positions(cls, start, cigar, seq):    
        '''Parse through a cigar string to return the genomic positions that are
        covered by the read.  Starts at the left-most 1-based start location'''
        
        array = []
        position = start
        
        # sometime the cigar value is "*", in which case assume a perfect match
        if cigar == "*":
            for i in len(seq):
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
        
        # assign start as left-most position (1-based)
        start = int(sam_parts[3])
        
        # create genomic positions for read
        positions = Read.build_positions(start, sam_parts[5], sam_parts[9])
        
        return Read(chromosome, start, strand, abundance, mappings, positions)
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
            
            

def metagene_count():
    '''Chain of command for metagene_count analysis.'''
    
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
   
    
    
    
def get_arguments():
    '''Collect and parse information from the user's command line arguments.'''
    
    date = datetime.datetime.now().strftime('%y%m%d-%H%M%S')
    
    parser = argparse.ArgumentParser(description=
    '''The first step of metagene_analysis, metagene_count.py compiles read abundance
over genomic features to create the input for metagene_windows.py. Please 
see README for full details and examples.''')

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






#def count_reads_around_features(arguments, chromosomes, chromosome_conversion):
    '''For each feature count the number of alignments (5' most base of each
    alignment only is counted) at each position {INTERVAL} nucleotides around
    the reference point-- start {-s} or end {-e} of the feature.
    
    Brief note about coordinate systems:
    The reference position is 0 with + and - {INTERVAL} positions on either side.
    However, the array storing the data runs from 0 to 2*{INTERVAL} with a length
    of 2*{INTERVAL} + 1.
    
    {INTERVAL} = 5
    Feature Index: -5 -4 -3 -2 -1  0  1  2  3  4  5    length = 11 
      Array Index:  0  1  2  3  4  5  6  7  8  9  10   length = 11
     Genome Index:  10 11 12 13 14 15 16 17 18 19 20   length = 11
      
    Therefore: array_index(reference point) = {INTERVAL}
    and the reference point is the ({INTERVAL} + 1)-ith value in the list
    
    More importantly, if the reference point is really 15 (1-based position) then
    genome_index(start) = genome_index(reference point) - ({INTERVAL} + 1) + 1
                        = 15 - (5+1) + 1 = 10
        (treating reference point as end and INTERVAL + 1 as length)
    genome_index(end)   = genome_index(start) + ({INTERVAL}*2 + 1) - 1
                        = 10 + (5*2 + 1) - 1 = 20
        (where the length = {INTERVAL} * 2 + 1 )'''
    

 

    
if __name__ == "__main__":
    Metagene.test_metagene()
    Feature.test_feature()
    
    samline = "read_1	24	chr1	250	255	10M40N20M	*	0	0	aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa	aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa	XA:i:0	MD:Z:40	NH:i:50  NA:i:10"
    
    print Read.create_from_sam(samline,{"1":"chr1"}).strand
    

#    metagene_count()    
   
    

