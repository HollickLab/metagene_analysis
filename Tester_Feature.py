#!/usr/bin/python
'''Tester_Feature.py to test the Feature.py class.

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

from Feature import Feature
from Metagene import Metagene
from MetageneError import MetageneError
from Read import Read 

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

        # create feature from GFF line with start and end swapped
        try:
            gffline = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(2,"test","gene",39,10,".","-",".","second")
            print "\t  with GFF line:\t{}".format(gffline.strip())
            feature2 = Feature.create_from_gff(method,metagene, gffline, chromosome_converter)
            if str(feature2.position_array) != correct_features['gff'][method]: 
                print "**FAILED**\t  Create Feature from GFF line with swapped start and end ?\t**FAIL**"
                print "\t  Desired positions:\t{}".format(correct_features['gff'][method])
                print "\t  Created positions:\t{}".format(feature2.position_array)
        except MetageneError as err:
            print "**FAILED**\t  Create Feature from GFF line with swapped start and end ?"
        else:
            print "PASSED\t  Create Feature from GFF line with swapped start and end ?\t\t{}".format(feature2.get_chromosome_region())
        try:
            gffline = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(2,"test","gene",39,10,".","+",".","second")
            print "\t  with GFF line:\t{}".format(gffline.strip())
            feature2 = Feature.create_from_gff(method,metagene, gffline, chromosome_converter)
            if str(feature2.position_array) != correct_features['gff'][method]: 
                print "**FAILED**\t  Do not create Feature from GFF line with swapped start and end, + strand ?\t**FAIL**"
                print "\t  Desired positions:\t{}".format(correct_features['gff'][method])
                print "\t  Created positions:\t{}".format(feature2.position_array)
        except MetageneError as err:
            print "PASSED\t  Do not create Feature from GFF line with swapped start and end, + strand ?\t\t{}".format(err)
        else:
            print "**FAILED**\t  Do not create Feature from GFF line with swapped start and end, + strand ?\t\t{}".format(feature2.get_chromosome_region())
              
            
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
             ((3,2,2),(2.5,4,(10.0/3),10,11,7.3,4),'7.611,19.556,14.967', "Contract with other messy floats ?")]
                 
        
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


if __name__=="__main__":
    test()    
