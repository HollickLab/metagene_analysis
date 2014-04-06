#!/usr/bin/python
'''Read class for metagene_counts.py.

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
import subprocess

from MetageneError import MetageneError
from metageneMethods import confirm_integer
from metageneMethods import runPipe

##TODO: add support for different alignment times (eg. bigwig or bigbed?)
class Read():
    """Create a read represented by a list of chromosomal positions.
    
    Class Attributes:
        has_sam_tag -- file level information about SAM tags
            value: dictionary; key = two-letter SAM tag (eg. "NA", "NH")
                               value = boolean
        cigar_codes -- description of CIGAR alignment options
            value: dictionary; key = single character CIGAR code
                               value = (boolean, boolean) for (counting, advancing position)
            reference: Descriptions from Sequence Alignment/Map Format Secification, 
                       28 Feb 2014; version 7fd84b0 from https://github.com/samtools/hts-specs 
                               
    Attributes:
        chromosome -- chromosome 
            value: matches a key in the chromosome_conversion dict
        strand -- strand relative to the chromosome
            value: '+' | '-' | '.'
        position_array -- chromosome positions (1-based) covered by the read
            value: array of integers from 5' to 3' relative to the read itself
                   with gaps in the read represented as jumps in the positions
                   start of read (5'-most base) = position_array[0]
                   end of read (3'-most base) = position_array[-1]
        abundance -- read count represented by the alignment line
            value: non-zero positive integer; 1 or extracted from NA:i:## tag
        mappings -- count of potentially alignment positions
            value: non-zero positive integer; 1 or extracted from NH:i:## tag
    
    Methods:
        create_from_sam -- create read object from SAM/BAM line
        parse_sam_bitwise_flag -- return countable and reverse_complement booleans
        build_positions -- build positions array from CIGAR alignment
        set_sam_tag -- add key:value pairs to has_sam_tag class dictionary 
    """
    
    __slots__ = ['chromosome','strand','position_array','abundance','mappings']
    
    # Keeps track of presence/absence of certain SAM file tags 
    # (defined by get_sam_tag classmethod)
    # Key = two-letter code; value = boolean
    has_sam_tag = {}
    
    # How to count[0] and advance[1] chromosome position by CIGAR code
    cigar_codes = { 'M':(True, True),  # alignment match (can be either sequence match or mismatch)
                    'I':(False, False), # insertion to the reference
                    'D':(False, True), # deletion from the reference
                    'N':(False, True), # skipped region from the reference
                    'S':(False, False), # soft clipping (clipped sequences present in SEQ)
                    'H':(False, False), # hard clipping (clipped sequences NOT present in SEQ)
                    'P':(False, True), # padding (silent deletion from padded reference)
                    '=':(True, True),  # sequence match
                    'X':(True, True) } # sequence mismatch  
    
    def __init__(self, chromosome, strand, abundance, mappings, positions):
        """Create read object. Invoke with a constructor rather than directly.
        
        Keyword Arguments:
        chromosome -- reference sequence name, often a chromosome
        strand -- read strand relative to the reference ['+'|'-'|'.']
        abundance -- identical reads represented by the same alignment
        mappings -- alignment positions for this read
        positions -- array of 1-based chromosome positions the read covers """
        self.chromosome = chromosome
        if strand != "+" and strand != "-":
            self.strand = "."
        else:
            self.strand = strand
        
        if self.strand == "-":
            positions.reverse()
        self.position_array = positions
        
        if confirm_integer(abundance, "Abundance", minimum=1):
            self.abundance = int(abundance)
        
        if confirm_integer(mappings, "Alignments", minimum=1):
            self.mappings = int(mappings)
    # End of __init__
    
    def __str__(self):
        return "Read at {0}:{1}-{2} on {3} strand; counts for {4:2.3f}:\t\t{5}".format(self.chromosome, self.position_array[0], self.position_array[-1], self.strand, float(self.abundance)/self.mappings, str(self.position_array))
    
    @classmethod
    def create_from_sam(cls, sam_line, 
                             chromosome_conversion, 
                             count_method, 
                             unique=False,
                             count_secondary_alignments=True,
                             count_failed_quality_control=False,
                             count_PCR_optical_duplicate=False,
                             count_supplementary_alignment=True):
        """Create a Read object from a BAM or SAM line.
        
        Keyword Arguments:
        sam_line -- raw line from SAM file or 'samtools view BAM_file' output
        chromosome_conversion -- dictionary of chromosome names
        count_method -- how to count the read ['all'|'start'|'end']
        unique -- boolean for mapping; if True mappings = 1 (default False)
        count_secondary_alignments -- process secondary alignment reads (default True)
        count_failed_quality_control -- process failed quality control reads (default False)
        count_PCR_optical_duplicate -- process PCR/optical duplicate reads (default False)
        count_supplementary_alignment -- process supplementary alignment reads (default True)
        """
        sam_parts = sam_line.split("\t")      
        if count_method == 'start':
            count_only_start = True
            count_only_end = False
        elif count_method == 'end':
            count_only_start = False
            count_only_end = True
        else:
            count_only_start = False
            count_only_end = False
        # parse bitwise flag
        (countable, reverse_complement) = Read.parse_sam_bitwise_flag(int(sam_parts[1]), 
                                                                      count_secondary_alignments, 
                                                                      count_failed_quality_control,
                                                                      count_PCR_optical_duplicate,
                                                                      count_supplementary_alignment,
                                                                      count_only_start,
                                                                      count_only_end)
        if countable: 
            # assign chromosome
            if sam_parts[2] not in chromosome_conversion.values():
                raise MetageneError("Read chromosome {} is not in the analysis set".format(sam_parts[2]))
            else:
                chromosome = sam_parts[2]
            # assign mappings
            if unique:
                mappings = 1
            # try to extract mappings from NH:i:## tag
            elif 'NH' in cls.has_sam_tag and cls.has_sam_tag['NH']: 
                try:
                    mappings = int(re.search('NH:i:(\d+)', sam_line).group(1))
                except AttributeError:
                    raise MetageneError("Could not determine number of mappings")
            else:
                raise MetageneError("Could not determine number of mappings")
        
            # assign abundance either from NA:i:## tag or as 1 (default)
            if 'NA' in cls.has_sam_tag and cls.has_sam_tag['NA']:
                try: 
                    abundance = int(re.search('NA:i:(\d+)', sam_line).group(1))
                except AttributeError:
                    raise MetageneError("Could not extract the abundance tag")
            else:
                abundance = 1
            # assign strand and positions
            if reverse_complement: # Crick or Minus strand
                strand = "-" 
            else: # Watson or Plus strand
                strand = "+"  
        
            # create genomic positions for read (start, cigar_string, sequence)
            positions = Read.build_positions(int(sam_parts[3]), sam_parts[5], sam_parts[9])
        
            return (countable, Read(chromosome, strand, abundance, mappings, positions))
        else:
            return (countable, "Non-aligning read")
    # end of create_from_sam
    
    @classmethod
    def parse_sam_bitwise_flag(cls, flags, 
                                    count_secondary_alignments=True, 
                                    count_failed_quality_control=False,
                                    count_PCR_optical_duplicate=False,
                                    count_supplementary_alignments=True,
                                    count_only_start=False,
                                    count_only_end=False):
        """Parse bitwise flag and return (countable, reverse_complemented) booleans.
        
        Keyword Arguments:
        flags -- decimal number representing bitwise flag
        count_secondary_alignments -- count reads with bitwise flag 0x100 (default True)
        count_failed_quality_control -- count reads with bitwise flag 0x200 (default False)
        count_PCR_optical_duplicate -- count reads with bitwise flag 0x400 (default False)
        count_supplementary_alignment -- count reads with bitwise flag 0x800 (default True)
        count_only_start -- count only alignment start reads; bitwise flag 0x40 (default False)
        count_only_end -- count only alignment end reads; bitwise flag 0x80 (default False)
       
        Understanding the SAM/BAM bitwise flag:  
        Based on description from Sequence Alignment/Map Format Secification, 
        28 Feb 2014; version 7fd84b0 from https://github.com/samtools/hts-specs
        
        Bit(hex)  binary_representation
        0x1       0000 0000 0001  template in multiple segments; if 0, MUST ignore string[-2,-4,-6,-7,-8]
        0x2       0000 0000 0010  each segment is properly aligned (according to aligner)
        0x4       0000 0000 0100  unmapped flag read MUST ignore string[-2,-5,-9,-12] and prev_string[-6]
        0x8       0000 0000 1000  next segment in template is unmapped MUST ignore string[-6]
        0x10      0000 0001 0000  reverse complement of seq (- strand)
        0x20      0000 0010 0000  next segment is reverse complemented
        0x40      0000 0100 0000  first segment of template
        0x80      0000 1000 0000  last segment of template
         *** User input can flip these default behaviors... ***
        0x100     0001 0000 0000  secondary alignment                   
        0x200     0010 0000 0000  not passing quality controls         
        0x400     0100 0000 0000  PCR or optical duplicate              
        0x800     1000 0000 0000  supplementary alignment               
        """        
        # make sure that only one of count_only_start or count_only_end is true
        if count_only_start and count_only_end:
            raise MetageneError("You can not count only the start and only the end, choose one or neither")
        if (flags & 0x10) == 0:
            reverse_complement = False
        else:
            reverse_complement = True
            
        # Is the read countable?
        # Does it map? 
        if (flags & 0x4) == 0x4: # flag is set and read is unmapped
            return (False,reverse_complement)
        # Is it a secondary alignment and do we care?
        elif (flags & 0x100) == 0x100 and not count_secondary_alignments:
            return (False,reverse_complement)  
        # Did it fail the quality control and do we care?
        elif (flags & 0x200) == 0x200 and not count_failed_quality_control:
            return (False,reverse_complement) 
        # Is it a PCR or optical duplicate and do we care?
        elif (flags & 0x400) == 0x400 and not count_PCR_optical_duplicate:
            return (False,reverse_complement)
        # Is it a supplementary alignment and do we care?
        elif (flags & 0x800) == 0x800 and not count_supplementary_alignments:
            return (False,reverse_complement)
        # Do we care about counting only the start or end? and does it matter (because part of a multi-segment template)?
        elif (count_only_start or count_only_end) and (flags & 0x1) == 0x1:
            # Do we care about the start and does this segment contain the start?
            if count_only_start and (flags & 0x40) == 0x40:
                return (True,reverse_complement)
            # Do we care about the end and does this segment contain the end?
            elif count_only_end and (flags & 0x80) == 0x80:
                return (True,reverse_complement)
            else:
                return (False,reverse_complement)
        else:
            # Made it through everything that could negate counting the read, so count it!     
            return (True, reverse_complement)

    @classmethod
    def build_positions(cls, start, cigar, seq):    
        """Return array of 1-based positions ordered relative to the chromosome.
        
        Keyword Arguments:
            start -- start of read relative to the chromosome: 1-based, left-most position
            cigar -- string representation of alignment (discussed below)
            seq -- sequence of the read
        """
        array = []
        position = int(start)
        # sometime the cigar value is "*", in which case assume a perfect match
        if cigar == "*":
            if seq != "*":
                for i in range(len(seq)):
                    array.append(position)
                    position += 1
                return array
            else:
                raise MetageneError("Unable to determine alignment length")

        # separate CIGAR string into nucleotide counts and CIGAR codes
        nucleotides = re.findall('(\d+)',cigar)
        codes = re.split('\d+',cigar)[1:]
        # loop through nucleotide values
        for i in range(len(nucleotides)):
            # iterate nt times adding 1 to start each time
            for j in range(int(nucleotides[i])):
                try:
                    if cls.cigar_codes[codes[i]][0]:
                        array.append(position)
                except KeyError:
                    raise MetageneError("Incorrect CIGAR string")
                if cls.cigar_codes[codes[i]][1]:
                    position += 1
        return array  
    # end of build_positions

    @classmethod
    def set_sam_tag(cls, count_tag, bamfile_name, tag_regex):
        """Add key:value pair to class variable: has_sam_tag.
        
        Keyword Arguments:
        count_tag -- boolean on whether to count with this tag
        bamfile_name -- file to query for tag
        tag_regex -- regular expression for the tag (eg 'NA:i:(\d+)')
        """
        tag = tag_regex.split(":")[0]
        (runPipe_worked, sam_sample) = Read.runPipe(['samtools view {}'.format(bamfile_name), 'head -n 10'])
        if runPipe_worked:
            num_tags = 0
            for sam_line in sam_sample:
                if re.search(tag_regex,sam_line) != None:
                    num_tags += 1
            if num_tags == 10:
                has_sam_value = True
            else:
                has_sam_value = False
                if count_tag:
                    raise MetageneError("Your alignment file does not have the required {} tag.".format(tag))
            cls.has_sam_tag[tag] = has_sam_value
            return True
        else:
            raise MetageneError("Checking the bam file failed with error: {}".format(sam_sample))  
# end of Read class
