metagene_analysis 
=================
Updated: 07 April 2014 by Joy-El R.B. Talbot

Description
===========
Create custom plots that summarize read alignments patterns over user-defined
genomic features. Requires alignment files (indexed BAM format) and feature files
(currently BED or GFF format) which can describe any type of feature expressed
by start and end chromosomal coordinates. Please see the USERS_GUIDE.md for 
more details and a list of metagene counting options.

Requirements
============
1. python 2.7.x (https://www.python.org/downloads/)
2. samtools 0.1.18+ (http://sourceforge.net/projects/samtools/files/)
3. R 2.15.2 (http://cran.us.r-project.org/)

Installation
============
1. Unzip download from github
2. Make sure that samtools and R are in your $PATH
3. Start by running metagene_count.py (then metagene_bin.py, then your flavor of metagene_plot_*.py)

Acknowledgements
================
The initial step in metagene_analysis (metagene_counts.py) is based on 
original Perl code by Karl F. Erhard, Jr Copyright (c) 2011.

License
=======
Joy-El R.B. Talbot Copyright (c) 2014

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
