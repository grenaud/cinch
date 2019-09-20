#!/usr/bin/python

from __future__ import absolute_import, division, print_function, unicode_literals
from optparse import OptionParser

import pathlib

#import matplotlib.pyplot as plt
import pandas as pd
#import seaborn as sns
import sys
import os, random
import subprocess

#!pip install -q tensorflow==2.0.0-beta1
import tensorflow as tf
import re

from tensorflow import keras
from tensorflow.keras import layers


def which(program):
    sys.stderr.write("Detecting program: "+str(program)+"\n");
    cjob = "type "+str(program);

    sp = subprocess.Popen(["/bin/bash", "-i", "-c", cjob],
                                stdout=subprocess.PIPE, 
                                stderr=subprocess.PIPE)
    

    out, err = sp.communicate();
    errcode  = sp.returncode;
    if(errcode != 0): #has finished but wrong code
        sys.stderr.write("Cannot find program: "+str(program)+" please make sure it is installed\n");
        sys.exit(1);
    
        #print("#"+str(str(out).find("aliased"))+"#");
        #out=str(out);
    if(out.find("aliased") != -1): #was aliased
        return out[out.find(" to ")+5:-2];
    else:
        if(out.find(" is ") != -1):
            return out[out.find(" is ")+4:];
        else:
            sys.stderr.write("Cannot seem to find program: "+str(program)+" please make sure it is installed\n");
            sys.exit(1);
            

    

    return out;

def handle_job(cjob):

    jobcreated=subprocess.Popen(cjob,shell=True,
                                executable="/bin/bash",
                                stdout=subprocess.PIPE, 
                                stderr=subprocess.PIPE); 
    jobcreated.wait()
    #print str(cjob);

    out, err = jobcreated.communicate()
    errcode  = jobcreated.returncode;

    if(err != ""): #has finished but has printed to stderr
        print("Job failed "+cjob+" failed");
        sys.exit(1);
    else:
        print("Job finished succesfully");

    return out;



sys.stderr.write("Detecting program: insize");

pathofexec  = os.path.abspath(sys.argv[0]);
pathofexecarray=pathofexec.split('/')[:-1];

pathofinsize =  ("/".join(pathofexecarray))+"/lib/insertsize/src/insize";


rcmd        = re.sub('\s+','',which("R"));

print(rcmd);

rscmd        = re.sub('\s+','',which("Rscript"));

print(rscmd);
#hmmcopycmd =  rcmd+"  CMD BATCH --vanilla --silent <(echo \"is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]);  is.installed('HMMcopy');\") /dev/stdout";
hmmcopycmd =  rscmd+" -e \"is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]);  is.installed('HMMcopy');\" ";

print(hmmcopycmd);
outputrcmd = handle_job(hmmcopycmd);

if( "[1] FALSE"  in outputrcmd ):
    sys.stderr.write("\nERROR: cannot find package HMMcopy, please install it (see http://bioconductor.org/packages/release/bioc/html/HMMcopy.html)\n");
    sys.exit(1);


if ( not (  "[1] TRUE"  in outputrcmd )):
    sys.stderr.write("\nERROR: not sure if we can find package HMMcopy, contact developers, this is an unknown case\n");
    sys.exit(1);
        

    

#sys.exit(1);


if(not os.path.exists(pathofinsize)):
    sys.stderr.write("\nERROR: The executable file "+pathofinsize+" does the exist, please type make in the main directory\n");
    sys.exit(1);


def which(program):
    sys.stderr.write("Detecting program: "+str(program)+"\n");
    cjob = "type "+str(program);

    sp = subprocess.Popen(["/bin/bash", "-i", "-c", cjob],
                                stdout=subprocess.PIPE, 
                                stderr=subprocess.PIPE)
    

    out, err = sp.communicate();
    errcode  = sp.returncode;
    if(errcode != 0): #has finished but wrong code
        sys.stderr.write("Cannot find program: "+str(program)+" please make sure it is installed\n");
        sys.exit(1);
    
        #print("#"+str(str(out).find("aliased"))+"#");
        #out=str(out);
    if(out.find("aliased") != -1): #was aliased
        return out[out.find(" to ")+5:-2];
    else:
        if(out.find(" is ") != -1):
            return out[out.find(" is ")+4:];
        else:
            sys.stderr.write("Cannot seem to find program: "+str(program)+" please make sure it is installed\n");
            sys.exit(1);
            

    

    return out;



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

    foffilefd = open(foffile, "r");
    bamfiles=[];
    label   =[];

    for linefd in foffilefd:
        fields=linefd.split("\t");
        if(len(fields)!=2):
            sys.stderr.write("\nThe line "+linefd+" does not have 2 tab separated fields\n");
            sys.exit(1);
        bamfiles.append( fields[0] );
        label.append(    fields[1] );#todo check if digit
        
    foffilesub=re.sub('/','_',foffile);
    foffilesub=re.sub("\ ",'_',foffilesub);
    
    
    logfile = (resultso+foffilesub+"_train.log");
    print(logfile);

    #stage 1: feature extraction
    #stage 2: training+writing model
    
    
    if(not os.path.exists(logfile)):#step 1
        #read file of file
        
        #insert size


        #CNV
        fileHandleLC = open ( ""+tfm+"/listcommands_1.txt", 'w' ) ;
        for bami in range(0,len(bamfiles)):
            fileHandleLC.write(pathofinsize+" "+bamfiles[bami]+" |sort -n |uniq -c |gzip > "+options.resultso+"/stage1/"+str(bami)+".isize.gz");
        fileHandleLC.close();
        
        
        logfilefp = open(logfile, "w");
        logfile.write("#-o:"+options.resultso+"\n");
        logfile.write("#fof:"+foffile+"\n");
        logfile.write("#stage1\n");
        
        print("Please run the commands manually either using:");
        print("cat "+tfm+"/listcommands.txt | parallel -j "+str(options.threads));
        print("on the use a batch/queueing system to launch:");
        print("cat "+tfm+"/listcommands.txt | sbatch ...");
        print("");
        print("Once commands are done, rerun with:\n");
        print("cinch.py   -o "+options.resultso+"  train "+foffile);
        print("");

        cmdtolaunch="cat "+tfm+"/listcommands.txt | parallel  -j "+str(options.threads);


        
    else:
        stage=1;
        logfilefp = open(logfile, "r");

        #for linelog in logfilefp:
        linelog = logfilefp.readline();
        if(linelog.startswith("#-o:")):
            options.resultso = linelog[len("#-o:"):len(linelog)];                
        else:
            sys.stderr.write("\nThe line "+linelog+" in "+logfile+" should start with #-o: something went wrong, please delete the log file.\n");
            sys.exit(1);

        linelog = logfilefp.readline();
        if(linelog.startswith("#fof:")):
            if(foffile != linelog[len("#fof:"):len(linelog)]):
                sys.stderr.write("\nThe line "+linelog+" in "+logfile+" should have #fof with the same as the ones provided as arguments: something went wrong, please delete the log file.\n");
                sys.exit(1);
                
        else:
            sys.stderr.write("\nThe line "+linelog+" in "+logfile+" should have #fof: something went wrong, please delete the log file.\n");
            sys.exit(1);

        linelog = logfilefp.readline();
        if(linelog.startswith("#stage1:")):
            stage=2;
        else:
            sys.stderr.write("\nThe line "+linelog+" in "+logfile+" should have with #stage1: something went wrong, please delete the log file.\n");
            sys.exit(1);

        logfilefp.close();

        if(stage==1):
            sys.stderr.write("\nCannot identify stage in "+logfile+" should have with #stage1: something went wrong, please delete the log file.\n");
            sys.exit(1);
                    
        if( stage == 2):
            #parse isize
            for bami in range(0:len(bamfiles)):
                fileisize = options.resultso+"/stage1/"+str(bami)+".isize.gz";
                if(not os.path.exists(fileisize)):
                    sys.stderr.write("\nThe file "+fileisize+" does not exist, please run all commands.\n");
                    sys.exit(1);
                else:
                    fileisizefd = open(fileisize, "r");

                    for lineisfd in fileisizefd:
                        fields=lineisfd.split("\t");
                    #TODO, add ability to read gzipped
                    
            #if(linelog read 


 

        


    
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


