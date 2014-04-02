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
# end of Read class
