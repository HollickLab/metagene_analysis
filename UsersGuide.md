metagene_analysis - User's Guide
================================

Step 1: metagene_count.py
=========================

Quick Start
-----------
Required input:
    1. Alignment file 
        - MUST be a sorted and indexed BAM file
        - The index file (.bai) must be in the same directory as the BAM file (.bam)
    2. Feature file 
        - Either BED or GFF format
    3. Chromosome conversion file 
        - TAB file if chromosome names vary between alignment and feature files

Basic command:
    python metagene_count.py --alignment alignement_file.bam --feature feature_file.bed [--chromosome_names conversion_file.tab] [--output_prefix my_counts] 
    
Additional Options
------------------

###**--feature_count [ all | start | end ]**

    Specify how each feature is processed and counted.

    Default: -- all

    Options:
        * all   -- each position from feature's start to end are included 
        * start -- only the feature's start position (1 nt) is included
        * end   -- only the feature's end position (1 nt) is included

    Additional notes:
        1. Start and end positions are relative to the feature with the start 
           being the 5'-most base of the feature and end being the 3'-most base.
           Therefore, a minus (Crick) strand feature's start will be larger than
           its end. 
        2. Output metagene will only scale the feature to the --interval_size if
           --feature_count all was set. For start or end counting only the 
           --interval_size is automatically 1.

--count_method [ all | start | end ]
    Specify how each read (alignment) is processed and counted.
    Default -- all
    Options
        1. all   -- each position from the read's start to end are included
        2. start -- only the read's start position (5'-most base) is included
        3. end   -- only the read's end position (3'-most base) is included
    Additional notes:
        1. Start and end positions are relative to the read with the start 
           being the 5'-most base of the feature and end being the 3'-most base.
           Therefore, a minus strand read's start will be larger than its end. 
        2. Counting is done in reads, not nucleotides. To convert the counts 
           under option --count_method all to reads, the counts are divided by
           the length of the read that maps to the feature (gaps for instance 
           are excluded).

--count_partial_reads
    
