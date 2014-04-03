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

import re, subprocess
from Metagene import Metagene
from Feature import Feature
from MetageneError import MetageneError

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
    
    # Keeps track of presence/absence of certain SAM file tags (defined by get_sam_tag classmethod)
    # Key = two-letter tag code (from SAM format specs)
    # Value = boolean 
    has_sam_tag = {}
    
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
        
        try:
            if abundance == int(abundance) and int(abundance) >= 0:
                self.abundance = int(abundance)
            else:
                raise MetageneError(abundance, "Abundance must be a positive integer")
        except ValueError as err:
            raise MetageneError(abundance, "Abundance must be a positive integer")
        
        try:    
            if mappings == "Unknown":
                self.has_mappings = False
                self.mappings = 1
            elif mappings == int(mappings) and int(mappings) > 0:
                self.has_mappings = True
                self.mappings = int(mappings)
            else:
                raise MetageneError(mappings, "Mapping must be a positive integer")
        except ValueError as err:
            raise MetageneError(mappings, "Mapping must be a positive integer")
    # End of __init__
    
    
    def __str__(self):
        return "Read at {0}:{1}-{2} on {3} strand; counts for {4:1.3f}:\t\t{5}".format(self.chromosome, self.position_array[0], self.position_array[-1], self.strand, float(self.abundance)/self.mappings, str(self.position_array))
    
    
    @classmethod
    def create_from_sam(cls, sam_line, 
                             chromosome_conversion, 
                             count_method, 
                             extract_abundance=False, 
                             unique=False,
                             count_secondary_alignments=True,
                             count_failed_quality_control=False,
                             count_PCR_optical_duplicate=False,
                             count_supplementary_alignment=True):
        '''Create a Read object from a bamfile line, requires that the chromosome 
        is in the chromosome_conversion dictionary
        
        extract_abundance = True --> extract abundance value from NA:i:## tag
                          = False --> assign value of 1 (uncollapsed reads were used)
                          
        unique = True --> force mappings to 1; 
                          only unique alignments represented (or desire to avoid hits-normalization)
               = False --> extract number of alignments from NH:i:## tag
        
        Default options for counting reads with certain SAM flags       
        count_secondary_alignments=True,
        count_failed_quality_control=False,
        count_PCR_optical_duplicate=False,
        count_supplementary_alignment=True'''
        
        sam_parts = sam_line.split("\t")

        # assign chromosome
        if sam_parts[2] not in chromosome_conversion.values():
            raise MetageneError(sam_parts[2], "Read chromosome {} is not in the analysis set".format(sam_parts[2]))
        else:
            chromosome = sam_parts[2]
        
        # parse bitwise flag
        count_only_start = False
        count_only_end = False
        if count_method == 'start':
            count_only_start = True
        elif count_method == 'end':
            count_only_end = True
        
        (countable, reverse_complement) = Read.parse_sam_bitwise_flag(int(sam_parts[1]), 
                                                                      count_secondary_alignments, 
                                                                      count_failed_quality_control,
                                                                      count_PCR_optical_duplicate,
                                                                      count_supplementary_alignment,
                                                                      count_only_start,
                                                                      count_only_end)
        
        if countable: # only process countable reads
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
            if reverse_complement: # Crick or Minus strand
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
    def parse_sam_bitwise_flag(cls, flags, 
                                    count_secondary_alignments=True, 
                                    count_failed_quality_control=False,
                                    count_PCR_optical_duplicate=False,
                                    count_supplementary_alignment=True,
                                    count_only_start=False,
                                    count_only_end=False):
        '''Parses bitwise flag to determine how to handle the read.
        Default values for count flags, pass flag name and opposite boolean to switch: 
            count_secondary_alignments      = True 
            count_failed_quality_control    = False
            count_PCR_optical_duplicate     = False
            count_supplementary_alignment   = True
            count_only_start                = False*
            count_only_end                  = False*
            * only one of these can be True at a time
          
        Based on description from Sequence Alignment/Map Format Secification, 
        28 Feb 2014; version 7fd84b0 from https://github.com/samtools/hts-specs
        
        Bit(hex)    binary string: 0000 0000 0000
        0x1         string[-1]                  ^  template in multiple segments; if 0, MUST ignore string[-2,-4,-6,-7,-8]
        0x2         string[-2]                 ^   each segment is properly aligned (according to aligner)
        0x4         string[-3]                ^    unmapped flag read MUST ignore string[-2,-5,-9,-12] and prev_string[-6]
        0x8         string[-4]               ^     next segment in template is unmapped MUST ignore string[-6]
        
        0x10        string[-5]             ^       reverse complement of seq (- strand)
        0x20        string[-6]            ^        next segment is reverse complemented
        0x40        string[-7]           ^         first segment of template
        0x80        string[-8]          ^          last segment of template
        
        *** User input can flip these default behaviors... ***
        0x100       string[-9]        ^            secondary alignment                   
        0x200       string[-10]      ^             not passing quality controls         
        0x400       string[-11]     ^              PCR or optical duplicate              
        0x800       string[-12]    ^               supplementary alignment               
        
          
        Return tuple of booleans for (all start as True):
        (Countable?, Reverse Complemented?)

        '''
        
        # make sure that only one of count_only_start or count_only_end is true
        if count_only_start and count_only_end:
            raise MetageneError((count_only_start, count_only_end), "You can not count only the start and only the end, choose one or neither")
        
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
        '''Parse through a cigar string to return the genomic positions that are
        covered by the read.  Starts at the left-most 1-based start location
        
        CIGAR codes and explainations from: 
        Descriptions from Sequence Alignment/Map Format Secification, 
        28 Feb 2014; version 7fd84b0 from https://github.com/samtools/hts-specs'''
        
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
    def set_sam_tag(cls, count_tag, bamfile_name, tag_regex):
        '''Sets has_sam_tag value for a particular tag'''
    
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
                    raise MetageneError(count_tag, "Your alignment file does not have the required {} tag.".format(tag))
        
            cls.has_sam_tag[tag] = has_sam_value

        else:
            raise MetageneError(sam_sample, "Checking the bam file failed with error: {}".format(sam_sample))  
        return True
    

    @staticmethod
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
# end of Read class
