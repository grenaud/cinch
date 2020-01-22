#!/usr/bin/python

from __future__ import absolute_import, division, print_function, unicode_literals
from optparse import OptionParser

import pathlib

#import matplotlib.pyplot as plt

#import seaborn as sns
import sys
import os, random
import subprocess
import gzip
from tqdm import tqdm

import numpy as np
from sklearn.decomposition import NMF, non_negative_factorization
from sklearn.preprocessing import normalize
from sklearn import metrics
from sklearn.metrics import roc_auc_score

from itertools             import combinations 


#!pip install -q tensorflow==2.0.0-beta1
#import pandas as pd
#import tensorflow as tf
import re
#from tensorflow import keras
#from tensorflow.keras import layers

#def normalize_rows(x: numpy.ndarray):
#    return x/numpy.linalg.norm(x, ord=2, axis=1, keepdims=True)


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



        

    

#sys.exit(1);


    return out;


mininsize=1;
maxinsize=700;
wcomponents=2;

#print(pathofconfig);
parser = OptionParser(
"\n\n"+

"\n"+
" ------------------------------   \n"+
    "       __          __       \n"+
"      /  ` | |\ | /  ` |__| \n"+
"      \__, | | \| \__, |  | \n"+
"\n"+
" ------------------------------   \n"+
"\n"
"Prediction of ctDNA for low coverage data\n"+
"\n"+

"\tcynch [options] [mode] [fof]\n"+
"\n"+
"There are 2 modes:\n"+
"\n"+
"\tcynch train [fof]\n"+
    "\tcynch predict [model] [fof]\n"
"\n"+
    "\tThe fof (file of file) as the following format for \"training\"\n"+
    "\t\t[bamfile1] [ctDNA frac 1]\n"+
    "\t\t[bamfile2] [ctDNA frac 2]\n"+
    "\t\t...\n"+

    "\tFor example:\n"+
    "\t\t/data/path/in1.bam 0.59\n"+
    "\t\t/data/path/in2.bam 0.11\n"+

    "\n"+
"\tThe fof has the following format for \"predict\"\n"+
    "\t\t[bamfile1]\n"+
    "\t\t[bamfile2]\n"+
    "\t\t..."+
"\n"+
"\n"+

"\tThe [model] is the *H.dat file created by the training and has the following shape\n"+
    "\t\tcomp1_1 comp1_2 comp1_3 ... \n"+
    "\t\tcomp2_1 comp2_2 comp2_3 ...\n"+
    "\t\t..."+
"\n"+
"\n"+

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
#parser.add_option("--map",                  dest="mappability",  help="Path to the mappability file for the reference genome",   default=None,    type="string");
#parser.add_option("--gc",                   dest="gc",           help="Path to the GC content file  for the reference genome",   default=None,    type="string");
#
#parser.add_option("--hpc",                  dest="hpc",          help="Use high-performance computing (for queueing systems ex: SGE)",          action="store_true");
#parser.add_option("--resume",               dest="resume",       help="Resume by providing the temp directory used",                              type="string");
#parser.add_option("--nice",                 dest="nice",         help="Nice the jobs",                                                          action="store_true");
parser.add_option("-t"  , "--threads",      dest="threads",      help="Number of threads to use if the local machine is used",  default=1,   type="int");
#
#parser.add_option("--mismap",               dest="mismappingrate", help="Mismapping rate , default is "+str(mismappingrate)+"",            default=mismappingrate, type="float");
#
##parser.add_option(""  , "--branchl",      dest="branchlscale", help="Seq-gen branch scale, default is 0.00045",            default=0.00045, type="float");
##parser.add_option(""  , "--chrlen",       dest="lengthchr",    help="Chromosome length, default is 10kb",              default=10000,   type="int");
parser.add_option("--minl",                 dest="minl",         help="Minimum fragment length to use, default is "+str(mininsize),                       default=mininsize,type="int")
parser.add_option("-c",                     dest="numcomp",      help="Number of components, default is 2",                   default=2,        type="int")
parser.add_option("--winsize",              dest="winsize",      help="Size of genomic windows, default is 1Mbp",             default=1000000,  type="int");
parser.add_option(""  , "--chrnum",         dest="numchr",       help="Number of chromosomes",                                default=22,       type="int");
parser.add_option(""  , "--chrprefix",      dest="prechr",       help="Prefix for chromosomes",                               default="chr",    type="string");

#parser.add_option("--gc",                  dest="gcfile",       help="File containing the % of GC per window size (see above) default: "+gcContent+"",                  default=gcContent,    type="string");
#
#parser.add_option("-c", "--conf",           dest="configfile",   help="Configuration for various conditions file to use for species/build default: "+str(pathofconfig)+"", default=pathofconfig, type="string");




#print handle_job("which ls");

(options,args) = parser.parse_args()

if( len(args) < 2 ):
    sys.stderr.write("\nneed a least 2 arguments, please use -h to see options.\n");
    sys.exit(1);


chrarray=[];
for chrn in range(1, options.numchr):
    chrarray.append( options.prechr+str(chrn));

sys.stderr.write("Using chromosomes: "+str( ",".join(chrarray) )+"\n" );

sys.stderr.write("Detecting program: insize");

pathofexec      = os.path.abspath(sys.argv[0]);
pathofexecarray = pathofexec.split('/')[:-1];

pathofinsize =  ("/".join(pathofexecarray))+"/lib/insertsize/src/insize";

if(not os.path.exists(pathofinsize)):
    sys.stderr.write("\nERROR: The executable file "+pathofinsize+" does not exist, please type make in the main directory\n");
    sys.exit(1);

#pathofreadcounter =  ("/".join(pathofexecarray))+"/lib/hmmcopy_utils/bin/readCounter";

#if(not os.path.exists(pathofreadcounter)):
#    sys.stderr.write("\nERROR: The executable file "+pathofreadcounter+" does not exist, please type make in the main directory\n");
#    sys.exit(1);


    
if(options.minl<1 or options.minl>=200 ):
    sys.stderr.write("\nERROR: The minimum length of the fragments:"+str(options.minl)+" cannot be less than 1 or more than 200\n");
    sys.exit(1);

mininsize = options.minl;

if(options.numcomp<=1):
    sys.stderr.write("\nERROR: The number of components  "+str(options.numcomp)+" cannot be less than 1\n");
    sys.exit(1);

if(options.numcomp>=9):
    sys.stderr.write("\nERROR: The number of components  "+str(options.numcomp)+" cannot greater or equal to 9\n");
    sys.exit(1);

wcomponents = options.numcomp;

#rcmd        = re.sub('\s+','',which("R"));
#print(rcmd);

#rscmd        = re.sub('\s+','',which("Rscript"));
#print(rscmd);
#hmmcopycmd =  rcmd+"  CMD BATCH --vanilla --silent <(echo \"is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]);  is.installed('HMMcopy');\") /dev/stdout";
#hmmcopycmd =  rscmd+" -e \"is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]);  is.installed('HMMcopy');\" ";

#print(hmmcopycmd);
#outputrcmd = handle_job(hmmcopycmd);

#if( "[1] FALSE"  in outputrcmd ):
#    sys.stderr.write("\nERROR: cannot find package HMMcopy, please install it (see http://bioconductor.org/packages/release/bioc/html/HMMcopy.html)\n");
#    sys.exit(1);


#if ( not (  "[1] TRUE"  in outputrcmd )):
#    sys.stderr.write("\nERROR: not sure if we can find package HMMcopy, contact developers, this is an unknown case\n");
#    sys.exit(1);

if( options.resultso == None):
    sys.stderr.write("\nPlease specify the outdir using -o outputDirectory/ or --outdir=outputDirectory/\n");
    sys.exit(1);

sys.stderr.write("\n\n\n");
resultso = options.resultso;
if(not resultso.endswith("/")):
    resultso = resultso+"/";

foffile=args[1];    
if(not os.path.exists(foffile)):
    sys.stderr.write("\nERROR: The file of files "+foffile+" does not exist\n");
    sys.exit(1);
pathoffof      = os.path.abspath(foffile);
pathoffofarray = pathoffof.split('/')[:-1];
pathoffofnof   = "/".join(pathoffofarray);


# train
def readFOFlabeled(foffile):
    foffilefd = open(foffile, "r");
    bamfiles  = [];
    label     = [];
    
    for linefd in foffilefd:
        linefd = linefd.strip();
        fields = linefd.split();

        if(len(linefd)==0):
            continue;

        if(len(fields)!=2):
            sys.stderr.write("\nThe line ->"+str(linefd)+"<- does not have 2 columns, please check the expected format using -h\n");
            sys.exit(1);
       
        bamfile = os.path.abspath(fields[0]);
        #print(bamfile);
        if(not os.path.exists(bamfile)):            
            bamfile2 = os.path.abspath(pathoffofnof+"/"+fields[0]);
            if(not os.path.exists(bamfile2)):            
                sys.stderr.write("\nERROR: The BAM file "+bamfile+" or "+bamfile2+" does not exist, make sure you have the correct relative or absolute path\n");
                sys.exit(1);
            else:
                if(not os.path.exists(bamfile2+".bai")):            
                    sys.stderr.write("\nERROR: The BAM file "+bamfile2+" is not indexed, please index the bam files and rerun\n");
                    sys.exit(1);
                    
                bamfiles.append( bamfile2 );
        else:
            if(not os.path.exists(bamfile+".bai")):            
                sys.stderr.write("\nERROR: The BAM file "+bamfile+" is not indexed, please index the bam files and rerun\n");
                sys.exit(1);

            bamfiles.append( bamfile   );
        fracbam = fields[1];
        try:
            fracbam = float(fracbam)
        except ValueError:
            sys.stderr.write("\nThe line ->"+str(len(linefd))+"<- does not have 2 columns where the second column is a floating point\n");
            sys.exit(1);

        if(fracbam>1 or fracbam<0):
            sys.stderr.write("\nThe second column should between 0 and 1, found "+str(fracbam)+" on line "+str(linefd)+"\n");
            sys.exit(1);
        label.append( fracbam );
        
    foffilesub=re.sub('/','_',foffile);
    foffilesub=re.sub("\ ",'_',foffilesub);
    foffilesub=re.sub("\.",'_',foffilesub);

    if(len(bamfiles) == 0):
        sys.stderr.write("\nNo BAM files were found in "+str(foffile)+"\n");
        sys.exit(1);

    
    return bamfiles,label,foffilesub;

def readFOFunlabeled(foffile):

    foffilefd = open(foffile, "r");
    bamfiles  = [];
    
    for linefd in foffilefd:
        linefd = linefd.strip();
        fields = linefd.split();

        if(len(linefd)==0):
            continue;

        if(len(fields)!=1):
            sys.stderr.write("\nThe line ->"+str(linefd)+"<- does not have 1 column, please check the expected format using -h\n");
            sys.exit(1);
       
        bamfile = os.path.abspath(fields[0]);
        #print(bamfile);
        if(not os.path.exists(bamfile)):            
            bamfile2 = os.path.abspath(pathoffofnof+"/"+fields[0]);
            if(not os.path.exists(bamfile2)):            
                sys.stderr.write("\nERROR: The BAM file "+bamfile+" or "+bamfile2+" does not exist, make sure you have the correct relative or absolute path\n");
                sys.exit(1);
            else:
                if(not os.path.exists(bamfile2+".bai")):            
                    sys.stderr.write("\nERROR: The BAM file "+bamfile2+" is not indexed, please index the bam files and rerun\n");
                    sys.exit(1);
                    
                bamfiles.append( bamfile2 );
        else:
            if(not os.path.exists(bamfile+".bai")):            
                sys.stderr.write("\nERROR: The BAM file "+bamfile+" is not indexed, please index the bam files and rerun\n");
                sys.exit(1);

            bamfiles.append( bamfile   );
            
        
    foffilesub=re.sub('/','_',foffile);
    foffilesub=re.sub("\ ",'_',foffilesub);
    foffilesub=re.sub("\.",'_',foffilesub);

    if(len(bamfiles) == 0):
        sys.stderr.write("\nNo BAM files were found in "+str(foffile)+"\n");
        sys.exit(1);


    return bamfiles,foffilesub;


def runStage1(resultso,threads,winsize,foffilesub,bamfiles,logfile):

    sys.stderr.write("\n");
    sys.stderr.write("    #####################################    \n");
    sys.stderr.write("    #   stage 1: feature extraction     #    \n");
    sys.stderr.write("    #####################################    \n");
    sys.stderr.write("\n");

    #CNV

    if not os.path.exists( ""+resultso+"/stage1/"):
        os.mkdir( ""+resultso+"/stage1/", 0755 );
    else:
        sys.stderr.write("\nThe directory "+resultso+"/stage1/"+" already exists\n");            

    fileHandleLC = open ( ""+resultso+"/listcommands_1.txt", 'w' ) ;
    #isize

    #sys.stderr.write("bamfiles"+str(bamfiles));

    for bami in range(0,len(bamfiles)):
        fileHandleLC.write(pathofinsize+"  "+bamfiles[bami]+" |sort --temporary-directory="+resultso+"/stage1/ -n |uniq -c |gzip > "+resultso+"/stage1/"+str(bami)+".isize.gz\n");

    #readcount, todo reenable?
    #if False:
    #    for bami in range(0,len(bamfiles)):
    #        fileHandleLC.write(pathofreadcounter+" -w "+str(winsize)+" -c "+str( ",".join(chrarray) )+" "+bamfiles[bami]+"  > "+resultso+"/stage1/"+str(bami)+".seg\n");
    #fileHandleLC.close();

    #logfile = (resultso+foffilesub+"_train.log");
    logfilefp = open(logfile, "w");
    logfilefp.write("#-o:"+resultso+"\n");
    logfilefp.write("#fof:"+foffile+"\n");
    logfilefp.write("#stage1:\n");
    logfilefp.close();

    print("Please run the commands manually either using:");
    print("  cat "+resultso+"/listcommands_1.txt | parallel -j "+str(threads));
    print("on the use a batch/queueing system to launch:");
    print("  cat "+resultso+"/listcommands_1.txt | sbatch ...");
    print("");
    print("Once commands are done, rerun with:\n");
    print("  cinch.py   -o "+resultso+"  train "+foffile);
    print("");

    cmdtolaunch="cat "+resultso+"/listcommands_1.txt | parallel  -j "+str(threads);
    
def parseIsize(resultso,foffilesub,bamfiles):
    isizedat = (resultso+foffilesub+"_isize.dat");
    isizedatfp = open(isizedat, "w");

    isizedatfp.write("#fileidx");
    for i in range(mininsize,(maxinsize+1)):
        isizedatfp.write( "\t"+str( i  ) );
    isizedatfp.write( "\n" );
    sys.stderr.write("Reading bam files\n");

    #dataAllisize = np.array([])
    dataAllisize = [];

    for bami in tqdm(range(0,len(bamfiles))):
        fileisize = resultso+"/stage1/"+str(bami)+".isize.gz";
        #datatoadd = np.array([]);
        datatoadd = [];

        if(not os.path.exists(fileisize)):
            sys.stderr.write("\nThe file "+fileisize+" does not exist, please run all commands.\n");
            sys.exit(1);
        else:
            isizeCount=[];
            for i in range(0,maxinsize+1):
                isizeCount.append(0);
            #print(fileisize);
            fileisizefd = gzip.open(fileisize, "r");

            for lineisfd in fileisizefd:
                #add filters and store
                fields=lineisfd.split( );

                try:
                    count   = int(fields[0]);
                    isize   = int(fields[1]);

                except ValueError:
                    continue;

                if( (isize<mininsize) or (isize>maxinsize) ):
                    continue;

                isizeCount[ isize ] = count;

            isizedatfp.write(str(bami));
            for i in range(mininsize,(maxinsize+1)):
                isizedatfp.write( "\t"+str( isizeCount[ i ]) );
                datatoadd.append( isizeCount[ i ]  );

            isizedatfp.write( "\n" );

            dataAllisize.append( datatoadd );

    isizedatfp.close();
    sys.stderr.write("\nWriten isize data to "+str(resultso+foffilesub+"_isize.dat")+".\n");
    return dataAllisize;


if(args[0] == "train"):

    #####################################
    #   stage 1: feature extraction     #
    #####################################

    bamfiles,label,foffilesub  =  readFOFlabeled(foffile);
    

    logfile = (resultso+foffilesub+"_train.log");
    #print(logfile);


    #stage 2: training+writing model
    
    
    if(not os.path.exists(logfile)):#step 1
        stage=1;
        #read file of file

        if not os.path.exists( ""+resultso):
            os.mkdir( ""+resultso, 0755 );
        else:
            sys.stderr.write("\nThe directory "+resultso+" already exists\n");            

        #insert size
        runStage1(options.resultso,options.threads,options.winsize,foffilesub,bamfiles,logfile);

        
    else:
        stage=2;
        sys.stderr.write("\n");
        sys.stderr.write("    #####################################    \n");
        sys.stderr.write("    #   stage 2: parsing features       #    \n");
        sys.stderr.write("    #####################################    \n");
        sys.stderr.write("\n");
        sys.stderr.write("Opening logfile:"+logfile+"\n");

        logfilefp = open(logfile, "r");

        #for linelog in logfilefp:
        linelog = logfilefp.readline();
        linelog = linelog.strip();
        
        if(linelog.startswith("#-o:")):
            options.resultso = linelog[len("#-o:"):len(linelog)];                
        else:
            sys.stderr.write("\nThe line "+linelog+" in "+logfile+" should start with #-o: something went wrong, please delete the log file.\n");
            sys.exit(1);

        linelog = logfilefp.readline();
        linelog = linelog.strip();

        if(linelog.startswith("#fof:")):
            if(foffile != linelog[len("#fof:"):len(linelog)]):
                sys.stderr.write("\nThe line "+linelog+" in "+logfile+" should have #fof with the same as the ones provided as arguments: something went wrong, please delete the log file.\n#"+linelog[len("#fof:"):len(linelog)]+"#");
                sys.exit(1);
                
        else:
            sys.stderr.write("\nThe line "+linelog+" in "+logfile+" should have #fof: something went wrong, please delete the log file.\n");
            sys.exit(1);

        linelog = logfilefp.readline();
        linelog = linelog.strip();

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
            #########################
            #      PARSE ISIZE      #
            #########################
            dataAllisize = parseIsize(options.resultso,foffilesub,bamfiles);
            #########################
            #      run HMMcopy      #
            #########################
            #normalize
            #print("dataAllisize");
            #print(dataAllisize.shape);
            #print(dataAllisize);

            #row_sums          = dataAllisize.sum(axis=1);          
            #dataAllisizeNorm  = dataAllisize / row_sums[:, numpy.newaxis]
            
            dataAllisizeNorm = np.array(normalize(dataAllisize,norm='l1'));
            #print(dataAllisize[0]);
            #print(dataAllisizeNorm[0]);

            #nmf = NMF(n_components=2, init='random', random_state=0)
            

            #W = nmf.fit_transform(dataAllisizeNorm);
            #H = nmf.components_;
            sys.stderr.write("\nRunning NMF...");
            W, H, n_iter = non_negative_factorization(dataAllisizeNorm, n_components=wcomponents,init='random', random_state=0,solver="mu",beta_loss=1,max_iter=1000)
            sys.stderr.write("done\n");

            #print("W");
            #print(W);
            
            hfile = (resultso+foffilesub+"_"+str(wcomponents)+"_H.dat");
            hfilefp = open(hfile, "w");
            
            rows = H.shape[0]
            cols = H.shape[1]
            for i in range(0, rows):
                hfilefp.write( str(H[i,0]) );
                for j in range(1, cols):
                    hfilefp.write( "\t"+str(H[i,j]) );
                hfilefp.write( "\n" );
            hfilefp.close();


            sys.stderr.write("Writen H matrix  "+hfile+"\n");

            WNorm = np.array(normalize(W,norm='l1'));
            #print(WNorm);


            wfile = (resultso+foffilesub+"_"+str(wcomponents)+"_W.dat");
            wfilefp = open(wfile, "w");
            
            rows = W.shape[0]
            cols = W.shape[1]
            for i in range(0, rows):
                wfilefp.write( str(W[i,0]) );
                for j in range(1, cols):
                    wfilefp.write( "\t"+str(W[i,j]) );
                wfilefp.write( "\n" );
            wfilefp.close();

            sys.stderr.write("Writen raw W matrix  "+wfile+"\n");

            wfile = (resultso+foffilesub+"_"+str(wcomponents)+"_Wnorm.dat");
            wfilefp = open(wfile, "w");
            
            rows = WNorm.shape[0]
            cols = WNorm.shape[1]
            for i in range(0, rows):
                wfilefp.write( str(WNorm[i,0]) );
                for j in range(1, cols):
                    wfilefp.write( "\t"+str(WNorm[i,j]) );
                wfilefp.write( "\n" );
            wfilefp.close();

            sys.stderr.write("Writen normalized W matrix  "+wfile+"\n");



            #check correlations
            
            # Get all combination of [1, 2, 3] 
            arraycomb=[]
            for j in range(1,wcomponents):
                comb = combinations(range(0,wcomponents),j) 
                # Print the obtained combutations 
                for i in list(comb): 
                    arraycomb.append(i);

            #print(arraycomb);
            correlationArray=[];
            correlationTOP= float("-inf");
            correlationTOPi= -1;
            
            #for i in arraycomb: #for each possible combination
            auc=False;
            if(len(np.unique(label))==2):
                auc=True;
                sys.stderr.write("Labels are binary, using AUC\n");

            for i in range(0,len(arraycomb)): #for each possible combination
                #print("i="+str(i))
                #add components
                #arraySumW=np.zeros(rows) 

                arraySumW=[];
                #sum the components 
                for j in range(0,rows):
                    #print("j="+str(j)+" "+str(WNorm[j]));
                    #print("j="+str(j)+" "+str(WNorm[j,i]));
                    arraySumW.append( sum( WNorm[ j , arraycomb[i] ] ) );
                    
                #print("i="+str(i))
                #print(bamfiles);
                #print(arraySumW)
                #print(label);
                #correlation

                if(auc):
                    corr=roc_auc_score(label,arraySumW);                    
                    corrFactor = corr;
                else:                 
                    corr=np.corrcoef(arraySumW, label);
                    corrFactor=corr[1][0];
  
                    #print(corrFactor);
                correlationArray.append(corrFactor);
                if(corrFactor>correlationTOP):
                    correlationTOP=corrFactor;
                    correlationTOPi= i;

            sys.stderr.write("Correlation factors:\n");
            sys.stderr.write("components\tcorrelation\n");
            wfilec = (resultso+foffilesub+"_"+str(wcomponents)+"_Wcomp.dat");

            for i in range(0,len(arraycomb)): #for each possible combination
                if(correlationTOPi == i):
                    wfilecfp = open(wfilec, "w");
                    wfilecfp.write( str(arraycomb[i]) );
                    wfilecfp.close();

                    sys.stderr.write(str(arraycomb[i])+"\t"+str(correlationArray[i])+"\tbest\n");
                else:
                    sys.stderr.write(str(arraycomb[i])+"\t"+str(correlationArray[i])+"\n");
                    
            sys.stderr.write("Writen best components  "+wfilec+"\n");

            #########################
            #      run HMMcopy      #
            #########################

            

 

        


    
    sys.exit(0);

#predict

if(args[0] == "predict"):

    #print(args[1]);
    #stage 1: feature extraction
    #stage 2: reading model +prediction

    hmatrix=args[1];
    foffile=args[2];


    #print(foffile);
    bamfiles,foffilesub  =  readFOFunlabeled(foffile);
    

    logfile = (resultso+foffilesub+"_predict.log");
    #sys.stderr.write("Trying to detect "+str(logfile)+"\n");



    #stage 2: training+writing model
    
    
    if(not os.path.exists(logfile)):#step 1
        stage=1;
        #read file of file
        
        #insert size
        runStage1(options.resultso,options.threads,options.winsize,foffilesub,bamfiles,logfile);
    else:

        hfile = hmatrix;
        if(not hfile.endswith("_H.dat")):
            sys.stderr.write("\nThe hfile "+str(hfile)+" should end with _H.dat\n");
            sys.exit(1);

        #detecting # of components

        fieldshf=hfile.split("_");
        try:
            wcomponents=int(fieldshf[len(fieldshf)-2] )
        except ValueError:
            sys.stderr.write("\nTried to infer the # of components used during training using "+str(hfile)+", failed. Make sure the file name ends with #_H.dat where # is the number of components\n");
            sys.exit(1);s

        sys.stderr.write("Model used "+str(wcomponents)+" components\n");
        sys.stderr.write("Opening model file "+str(hfile)+"\n");
        hfilefp = open(hfile, "r");
        Hmat=[];
        for hline in hfilefp:
            fields=hline.split( );
            Hmattoadd=[];

            for i in range(0,len(fields)):
                #print(fields[i]);
                Hmattoadd.append( float(fields[i]) );
            
            Hmat.append( Hmattoadd );

        Hnp = np.array( Hmat );

        #reading correlation
        correlationTOP= -1;

        
        wcompfile   = hmatrix[0:-5]+"Wcomp.dat";
        sys.stderr.write("Opening component file "+str(wcompfile)+"\n");
        
        wcompfilefp = open(wcompfile, "r");
        for wcompline in wcompfilefp:
            correlationTOP=wcompline;
        correlationTOP=correlationTOP.replace("(","") 
        correlationTOP=correlationTOP.replace(")","") 
        correlationTOP=correlationTOP.replace(" ","") 



        correlationTOP=np.fromstring(correlationTOP,dtype=int, sep=',');
        sys.stderr.write("Using best W components: "+str(correlationTOP)+"\n");
        

        #sys.exit(1);

        dataAllisize = parseIsize(options.resultso,foffilesub,bamfiles);
        #print(dataAllisize);
        dataAllisizeNorm = np.array(normalize(dataAllisize,norm='l1'));

        #print(Hnp.shape);
        #print(dataAllisizeNorm.shape);


        W, H, n_iter = non_negative_factorization(dataAllisizeNorm, 
                                                  n_components=wcomponents,
                                                  #W=None, 
                                                  H=Hnp,
                                                  update_H=False,
                                                  init='custom', 
                                                  random_state=0,
                                                  solver="mu",
                                                  beta_loss=1,
                                                  max_iter=200);
        #writing out W    
        wfile = (options.resultso+foffilesub+"_"+str(wcomponents)+"_W.out");
        sys.stderr.write("\nWriten W matrix to "+wfile+"\n");

        wfilefp = open(wfile, "w");

        rows = W.shape[0]
        cols = W.shape[1]
        for i in range(0, rows):
            wfilefp.write( str(W[i,0]) );
            for j in range(1, cols):
                wfilefp.write( "\t"+str(W[i,j]) );
            wfilefp.write( "\n" );
        wfilefp.close();

        Wnorm = np.array(normalize(W,norm='l1'));

        #writing out W norm
        wfile = (options.resultso+foffilesub+"_"+str(wcomponents)+"_Wnorm.out");
        sys.stderr.write("\nWriten normalized W matrix to "+wfile+"\n");

        wfilefp = open(wfile, "w");

        rows = Wnorm.shape[0]
        cols = Wnorm.shape[1]
        for i in range(0, rows):
            wfilefp.write( str(Wnorm[i,0]) );
            for j in range(1, cols):
                wfilefp.write( "\t"+str(Wnorm[i,j]) );
            wfilefp.write( "\n" );
        wfilefp.close();

        #writing out subsampled W norm
        wfile = (options.resultso+foffilesub+"_"+str(wcomponents)+"_Wnormcomp.out");
        sys.stderr.write("\nWriten correlated score matrix to "+wfile+"\n");

        wfilefp = open(wfile, "w");

        rows = Wnorm.shape[0]
        cols = Wnorm.shape[1]
        for i in range(0, rows):
            wfilefp.write( str(sum(Wnorm[i,correlationTOP])) );
            wfilefp.write( "\n" );
        wfilefp.close();



    #end else not stage 1
    
        
        #print(W);
    #print(hmatrix);
    #print(foffile);

    sys.exit(0);



sys.stderr.write("\nUnknown mode.\n");
sys.exit(1);


