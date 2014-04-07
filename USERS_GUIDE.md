metagene_analysis - User's Guide
================================

Step 1: metagene_count.py
=========================

Quick Start
-----------

####Required input:

    1. Alignment file 
        - MUST be a sorted and indexed BAM file
        - The index file (.bai) must be in the same directory as the BAM file (.bam)
    2. Feature file 
        - Either BED or GFF format
    3. Chromosome conversion file 
        - TAB file if chromosome names vary between alignment and feature files

####Basic command:

    python metagene_count.py --alignment alignement_file.bam --feature feature_file.bed [--chromosome_names conversion_file.tab] [--output_prefix my_counts] 

####Output:

    Comma-delimited file for metagene_bin.py
    
    Nameing: output_prefix.metagene_counts.csv, where output_prefix is either the
    prefix you defined or the default "date.metagene"
    
    
Additional Options
------------------

####**--feature_count [ all | start | end ]**

    Specify how each feature is processed and counted.

    Default: all

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

####**--count_method [ all | start | end ]**

    Specify how each read (alignment) is processed and counted.
    
    Default: all
    
    Options:
        * all   -- each position from the read's start to end are included
        * start -- only the read's start position (5'-most base) is included
        * end   -- only the read's end position (3'-most base) is included
        
    Additional notes:
        1. Start and end positions are relative to the read with the start 
           being the 5'-most base of the feature and end being the 3'-most base.
           Therefore, a minus strand read's start will be larger than its end. 
        2. Counting is done in reads, not nucleotides. To convert the counts 
           under option --count_method all to reads, the counts are divided by
           the length of the read that maps to the feature (gaps for instance 
           are excluded).

####**--count_partial_reads**

    Flag to include in count reads only partially aligning with feature 
    
    Default: False (require that reads completely fall in feature (+/- padding) boundary)

####**--padding INT**

    Define size of padding (in nucleotides) added to each side of the feature
    
    Default: 1000

####**--interval_size INT**

    Define the size of the metagene feature region
    
    Default: 1000 (unless --feature_count = [start | end] then default is 1, see --feature_count for details)
    
####**--interval_variable**

    Flag to prevent shrinking/expansion of feature to metagene interval_size.  
    *Only works for downstream analysis with a feature_file of 1 feature* 
    
    Default: False

####**--ignore_strand**

    Combine + and - strand reads into an unstranded count
    
    Default: False
    
    Additional notes:
        *  Normally, the metagene_count.py script extracts strand information from
           both features and reads (alignments). Following the rules:
             1. If the features are unstranded all reads (regardless of strand 
                information) will be counted as unstranded.
             2. If the features are stranded and reads are unstranded 
                the program will error.
             3. If the features are stranded and the reads are stranded, 
                then stranded counts will be produced.
             4. Whenever strand information is available it will at least be 
                used to determine the start and end positions of the read or 
                feature. Start/end of unstranded reads or features are determined,
                but not necessarily counted (see rules 1 and 2) as though from 
                the + strand.
        *  --ignore_strand allows stranded reads to be counted on stranded features
           without separating the counts into distinct tallies as is normally done
           under rule 3 above.

####**--extract_abundance**

    Flag to pull abundance information from BAM alignment file
    
    Requires: "NA:i:##" tag in the BAM file lines; ## = abundance 
    
    Default: False; and each alignment line in BAM file is counted with abundance = 1
    
####**--extract_mappings**
    
    Flag to pull total number of mappings (to bowtie index) from BAM alignment file;
    will result in hits-normalized (abundance / mappings) counts.
    
    Requires: "NH:i:##" tag in the BAM file lines; ## = number of mappings
      
    Default: False; and each alignment line is assumed to have mappings = 1

####**--uniquely_mapping**

    Flag to assert that each alignment in BAM aligned uniquely to the reference
    
    Default: False

####**--chromosome_names chromosome_conversion_file.tab**

    Add a conversion to account for name differences between feature and alignment files.
    
    Default: None
    
    Format Requirements:
        1. Tab ("\t")-delimited file
        2. First column is the feature file chromosome name
        3. Second column is the alignment file chromosome name
            - must be the same as the name in the BAM header (SN: column in @SQ 
              lines from command: samtools view -H bamfile.bam)

####**--count_splicing**
    
    Flag to separate reads into gapped (potentially spliced) and ungapped (potentially
    unspliced and other ungapped reads like fully intronic, fully exonic, fully
    intergenic even) categories with separate output lines.
    
    Default: False

####**--include_reads [ secondary_alignment | failed_quality_control | PCR_duplicate | supplementary_alignment ]**

    Flags to change how different types of reads are handled
    
    Default: 
        * secondary_alignment:     True
        * failed_quality_control:  False
        * PCR_duplicate:           False
        * supplementary_alignment: True
        
    Additional notes:
        1. These four options may be represented in the bitwise flag of your BAM
           alignment file.
        2. If the second column of your BAM alignment lines is always less than 256
           then these tags have not been accounted for (or are not present) in 
           by your aligner program.
        3. See Samtools documentation for more details on the bitwise flag:
           Sequence Alignment/Map Format Secification, 28 Feb 2014; version 7fd84b0
           from https://github.com/samtools/hts-specs


Step 2: metagene_bin.py
=======================

Quick Start
-----------
Required input:
    1. Output from metagene_count.py

Basic command:
    python metagene_bin.py --input prefix.metagene_count.csv [--output_prefix my_bins] [--separate_groups]

    --separate_groups is technically optional, but required to flow output into the
    metagene_plot*.py commands
    
Output:
    Comma-delimited file for metagene_plot*.py commands.
    
    NOTE: to be parsed by the metagene_plot*.py, the output names *must end* in 
    the pattern: 
        * '.\d+bpX\d+bp.[A-Za-z]+_[A-Za-z]+.csv' 
        * example: '.10bpX10bp.sense_allreads.csv'
    
    Nameing: output_prefix.metagene_counts_prefix.10bpX10bp.orientation_gap.csv
    
Additional Options
------------------

####**--input metagene_count_output_file**
    
    Input file
    
    Additional notes:
        1. Multiple files can be processed at once by repeating the --input file 
           command with each new file name

####**--output_prefix prefix**

    Output prefix to append to resulting file names
    
    Default: binned
    
####**--window_size INT**
    
    Size of the windows (in nucleotides) to bin metagene_counts
    
    Default: 10
    
####**--step_size INT**

    Size of the step (in nucleotides) to slide the windows.
    
    Default: 10
    
    Additional note:
        1. To have non-overlapping windows: set step_size = window_size
        2. Decrease step_size relative to window_size to further smooth the plots

####**--separate_groups**

    Flag to separate different counting groups (eg sense vs antisense or 
    gapped vs ungapped) into different files. REQUIRED for downstream plotting 
    with metagene_plot*.py
    
    Default: False


Step 3: metagene_plot*.py
=========================

Note: These scripts are much more specific in their requirements. If you prefer, 
please feel free to check out their underlying .R script files for tips to guide
you in custom plotting of your binned metagenes with R. 

Quick Start
-----------
Required input:
    1. Output from metagene_bin.py; DO NOT CHANGE file names (see above)!

Output:
    PDF file of a plot
    
    NOTE: The rendering of the plot may sometimes result in legends overlapping
    your data. While I hope to address this in the future, the best work-around 
    for now is to import the PDF plot into either Adobe Illustrator or the free
    Inkscape (http://www.inkscape.org/en/). If you ungroup the image enough times
    the text of the legend will become editable and you can move it to a new
    location. These tools are also a great option to turn your data from 
    metagene_analysis into summary figures and adjust line sizes and fonts for
    publication in your favorite journal. :-)

The rest of the details vary by which metagene_plot*.py is chosen.  

Options:
    1. metagene_plot_only.py  
        - plot many metagenes (of same size) at once
    2. metagene_plot_individual.py  
        - plotting output of metagene_count.py --interval_variable
    3. metagene_plot_with_statistics.py
        - similar to plot_only, but requires 2 sets of sense/antisense data and 
          returns a plot with red horizontal bars to indicate (strand specific)
          significantly different pairwise Welch's t-tests passing the 
          Holm-Bonferoni cutoff (alpha = 0.05).  Also returns tables of the 
          statistics performed.

Step 3a: metagene_plot_only.py
==============================

Quick Start
-----------

Required input:
    1. Output from metagene_bin.py; DO NOT CHANGE file names (see above)!

Basic command:
    python metagene_plot_only.py --output_prefix my_plot --feature_counted [ TSS | Start | End | other_string ] --data_set binned.10bpX10bp.sense_allreads.csv,binned.10bpX10bp.antisense_allreads.csv,total_reads_in_library,color(for plotting),name(for legend)
    
    --data_set entries must be comma-delimited with NO spaces, and can be repeated
    for as many data_sets as you want
    
    order of data_set entry:
        1. binned file for the top strand
        2. binned file for the bottom strand
        3. normalization value (eg. total reads mapped, or total reads)
        4. color to draw plot lines (eg. blue, black, red, green) 
           see http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf for options
        5. name of the data_set (will be used in the legend)
    
Output:
    PDF file of a plot

Additional Options
------------------

####**--feature_counted [ TSS | Start | End | othertext ]**

    This text will appear in the x-axis label of your plot. It is also used to 
    determine how the grey feature box in drawn (see below).

    Drawing the grey feature box:
        * Shows the feature region of the metagene, but if you only counted the
          Start or End (--feature_count [ Start | End ]) of your features in 
          the metagene_count.py step, then you may want to extend the box to
          represent the entire feature region.
        * Options TSS or Start:
            - Draws the box from the Start of the feature to the right end of the 
              metagene (encompasses the downstream padding)
        * Option End:
            - Draws the box from the left end of the metagene to the End of the 
              feature (encompasses the upstream padding)
        * Option anything else:
            - Draws the box from the Start of the feature to the End of the 
              feature (length will vary and may be as little as 1 unit if only
              the Start of End of the feature were originally counted)

Step 3b: metagene_plot_individual.py
==============================

Quick Start
-----------

Required input:
    1. Output from metagene_bin.py; DO NOT CHANGE file names (see above)!
       *MUST* originate from a metagene_count.py --interval_variable run with a
       *single* feature!!

Basic command:
    python metagene_plot_individual.py --output_prefix my_plot --normalization INT --fileset_a binned.10bpX10bp.unstranded_gapped.csv --fileset_b binned.10bpX10bp.unstranded_ungapped.csv --individual_splicing
        
Output:
    PDF file of a plot
    
Step 3c: metagene_plot_with_statistics.py
==============================

Quick Start
-----------

Required input:
    1. Output from metagene_bin.py; DO NOT CHANGE file names (see above)!

Basic command:
    python metagene_plot_with_statistics.py --output_prefix my_plot --fileset_a binned.set1.10bpX10bp.sense_allreads.csv --fileset_a binned.set1.10bpX10bp.antisense_allreads.csv --normalization_a INT --fileset_b binned.set2.10bpX10bp.sense_allreads.csv --fileset_b binned.set2.10bpX10bp.antisense_allreads.csv --normalization_b INT --feature_counted [ TSS | Start | End | other_string ]
    
Output:
    1. PDF file of a plot
    2. CSV files of statistics for sense and antisense strands
    

