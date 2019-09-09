#!/usr/bin/python

from __future__ import absolute_import, division, print_function, unicode_literals
from optparse import OptionParser

import pathlib

#import matplotlib.pyplot as plt
import pandas as pd
#import seaborn as sns
import sys
#!pip install -q tensorflow==2.0.0-beta1
import tensorflow as tf
import re

from tensorflow import keras
from tensorflow.keras import layers



 
#print(pathofconfig);
parser = OptionParser(
"\n\n"+

"\n"+
"Prediction of ctDNA for low coverage data\n"+
"\n"+

"\tcynch [options] [mode] [fof]\n"+
"\n"+
"There are 2 modes:\n"+
"\n"+
"\tcynch train [fof]\n"+
"\tcynch predict [fo]\n"
"\n"+
    "\tThe fof (file of file) as the following format for \"training\"\n"+
    "\t\t[bamfile1] [ctDNA frac 1]\n"+
    "\t\t[bamfile2] [ctDNA frac 2]\n"+
    "\tFor example:\n"+
    "\t\t/data/path/in1.bam 0.59\n"+
    "\t\t/data/path/in2.bam 0.11\n"+

    "\n"+
"\tThe fof has the following format for \"predict\"\n"+
    "\t\t[bamfile1]\n"+
    "\t\t[bamfile2]\n"+
    "\t\t"+
"\tThe program will auto-resume if you have sent commands to the cluster or used parallel\n"+
    
"\n"
);


##parser.add_option("-t", "--theta",       dest="theta",        help="Theta, default is 20",                          default=20,      type="float");
#parser.add_option("-r", "--reference",      dest="ref",          help="Reference genome used for alignment",                 default=None,    type="string");
#parser.add_option("--tmpfolder",            dest="tmpfolder",    help="Temporary folder",                                    default="/tmp/",    type="string");
parser.add_option("-o","--outdir",         dest="resultso",     help="Output directory to store the results",                  default=None,    type="string");

#
#parser.add_option("--samtools",             dest="samtools",     help="Use this version of samtools",                        default="samtools",    type="string");
#parser.add_option("--tabix",                dest="tabix",        help="Use this version of tabix",                           default="tabix",    type="string");
#
#parser.add_option("--map",                  dest="mappability",  help="Use a mappability map in BED format (Recommended)",   default=None,    type="string");
#
#parser.add_option("--hpc",                  dest="hpc",          help="Use high-performance computing (for queueing systems ex: SGE)",          action="store_true");
#parser.add_option("--resume",               dest="resume",       help="Resume by providing the temp directory used",                              type="string");
#parser.add_option("--nice",                 dest="nice",         help="Nice the jobs",                                                          action="store_true");
#parser.add_option("-t"  , "--threads",      dest="threads",      help="Number of threads to use if the local machine is used",  default=1,   type="int");
#
#parser.add_option("--mismap",               dest="mismappingrate", help="Mismapping rate , default is "+str(mismappingrate)+"",            default=mismappingrate, type="float");
#
##parser.add_option(""  , "--branchl",     dest="branchlscale", help="Seq-gen branch scale, default is 0.00045",            default=0.00045, type="float");
##parser.add_option(""  , "--chrlen",      dest="lengthchr",    help="Chromosome length, default is 10kb",              default=10000,   type="int");
##parser.add_option("-c", "--numcont",     dest="numcont",      help="Number of present-day human contaminants, default is 2", default=2, type="int")
#parser.add_option("--winsize",              dest="winsize",      help="Size of genomic windows, default is 1Mbp",             default=1000000, type="int");
#parser.add_option("--gc",                   dest="gcfile",       help="File containing the % of GC per window size (see above) default: "+gcContent+"",                  default=gcContent,    type="string");
#
#parser.add_option("-c", "--conf",           dest="configfile",   help="Configuration for various conditions file to use for species/build default: "+str(pathofconfig)+"", default=pathofconfig, type="string");



#print handle_job("which ls");

(options,args) = parser.parse_args()

if( len(args) < 2 ):
    sys.stderr.write("\nplease use -h to see options.\n");
    sys.exit(1);

if( options.resultso == None):
    sys.stderr.write("\nPlease specify the outdir\n");
    sys.exit(1);

resultso = options.resultso;
if(not resultso.endswith("/")):
    resultso = resultso+"/";




# train

if(args[0] == "train"):

    foffile=args[1];

    foffilesub=re.sub('/','_',foffile);
    foffilesub=re.sub("\ ",'_',foffilesub);
    
    
    logfile = (resultso+foffilesub+"_train.log");
    print(logfile);

    #stage 1: feature extraction
    #stage 2: training+writing model
    
    
    if(not os.path.exists(logfile)):
        #insert size


        #CNV
        fileHandleLC = open ( ""+tfm+"/listcommands_1.txt", 'w' ) ;
        #print commands here
        
        fileHandleLC.close();
        
        
        logfilefp = open(logfile, "w");
        logfile.write("#stage1\n");
        

        
    else:
           
        logfilefp = open(logfile, "r");

        for linelog in logfilefp:
            print(linelog);

        logfilefp.close();
 

        


    
    sys.exit(0);

#predict

if(args[0] == "predict"):

    #print(args[1]);
    #stage 1: feature extraction
    #stage 2: reading model +prediction

    foffile=args[1];

    print(foffile);
    sys.exit(0);



sys.stderr.write("\nUnknown mode.\n");
sys.exit(1);


