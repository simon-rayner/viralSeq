from pypesteps import abstractStep

'''
Created on Dec 18, 2020

@author:     simon rayner
@contact:    simon.rayner@medisin.uio.no
'''
import os
import csv
import subprocess
import shlex
import glob

import pandas as pd
import plotnine as p9


import logging

logger = logging.getLogger(__name__)
INDENT = 6


class StepBAMGCReadCorr(abstractStep.AbstractStep):
    '''
    classdocs
    This calculates the correlation between read coverage and gc content
    for a set of BAM files.
    
    The analysis assumes the read coverage and GC calculation has already
    been performed for the specified BAM files and are in the folders
    specified at runtime (using the -r/--readfolder and -g/gcfolder parameters)
    
    The analysis is done in two steps:
    
        1. SAMTools is used to calculate the read coverage/nt, 
        2. This data to calculate a sliding window average (if requested)
    
    SAMTools is used because it's faster
    
    The user must also specify:
        a sliding window and step size.
        both values need to be specified
        (-w/--window_size & (-s/--step_size)
    
    '''
    CLASSID             = "StepBAMGCReadCorr"
    OUTPUTFOLDER        = "gcreadcorrelation"
    GCFILESHORT         = "-g"
    GCFILELONG          = "--gc_coverage_file"
    READFILESHORT       = "-r"
    READFILELONG        = "--read_coverage_file"
    WINSIZESHORT        = "-w"
    WINSIZELONG         = "--window_size"
    STEPSIZESHORT       = "-s"
    STEPSIZELONG        = "--step_size"
    BAMFILEFOLDERSHORT  = "-b"
    BAMFILEFOLDERLONG   = "--bam_file_folder"
    
    PLOTHEIGHT          = 8
    PLOTWIDTH           = 10
    PLOTUNITS           = 'in'
    PLOTDPI             = 1000
    XVAR                = "pos"
    YVAR                = "readcoverage"
    

    def __init__(self, windowSize=0, stepSize=0, gccoveragefile= "", readcoveragefile="", bamFileFolder=""):
        '''
        Constructor
        '''
        self.windowSize = windowSize
        self.stepSize = stepSize
        self.gcCoverageFile = gccoveragefile
        self.readCoverageFile = readcoveragefile
        self.bamFileFolder = bamFileFolder

        
    def checkInputData(self):
        '''
        build file paths and check all input resources exist
        filepath can be absolute or relative (to Project Root)
        '''
        
        # Check input folders exist
        if( os.path.exists(os.path.join(self.projectRoot, self.inFolder)) is False):
            logging.error("input folder for BAM files <" + os.path.join(self.projectRoot, self.inFolder) + "> not found")
            raise Exception("input folder for BAM files <" + os.path.join(self.projectRoot, self.inFolder) + "> not found")
       
                
        if os.path.exists(self.gcCoverageFile) == False:
            raise RuntimeError ("gc coverage file <" + self.gcCoverageFile + "> not found")
            logging.error("gc coverage file <" + self.gcCoverageFile + "> not found")
        else:
            logging.info(INDENT*'-' + "found gc coverage file <" + self.gcCoverageFile + ">")
        
        if os.path.exists(self.readCoverageFile) == False:
            raise RuntimeError ("read coverage file <" + self.readCoverageFile + "> not found")
            logging.error("read coverage file <" + self.readCoverageFile + "> not found")
        else:
            logging.info(INDENT*'-' + "found read coverage file <" + self.readCoverageFile + ">")
        
        
        ## check window size and step size are positive integers 
        if(self.windowSize > 0 and self.stepSize > 0):
            logging.info(INDENT*'-' + "-Read Coverage will be calculated using a sliding window of <" +\
                          str(self.windowSize) + ">nt and a step size of <" + str(self.stepSize) + ">")
        elif(self.windowSize <= 0 and self.stepSize <= 0):
            logging.info("No sliding window selected")            
        else:
            logging.error("both window and step size must be specified and > 0: (found window size <" \
                            + self.windowSize + " and step size <" + self.stepSize + ">)")
            raise Exception("both window and step size must be specified and > 0: (found window size <" \
                            + self.windowSize + " and step size <" + self.stepSize + ">)")
            
            
        # do the window and step size match for the GC and Read Coverage files?
        if len(self.inputFiles) == 1 & (not self.inputFiles[0]):
            self.inputFiles = glob.glob(os.path.join(self.projectRoot, self.inFolder) + os.path.sep + "*gen__trim_paired__sorted.bam")
        
        
        
    def execute(self):
        '''
        contains the main operations for the step
            1. load gc coverage file
            2. for each BAM file, load read coverage file 
            3. for each position get GC coverage and Read Coverage
        '''
        logger.info(INDENT*'-' + "executing step")
        bamFileFolder = os.path.join(self.projectRoot, self.inFolder)
        resultFolder = os.path.join(self.projectRoot, self.outFolder)
        
        logging.info(INDENT*'-' + "--results will be written to output folder <" + resultFolder + ">")
        if not os.path.exists(resultFolder):
            logging.info(INDENT*'-' + "----folder doesn't exist, creating")
            os.makedirs(resultFolder)        
        
        dfGCcoverage = pd.read_csv(self.gcCoverageFile, skiprows=1, delimiter="\t")
        dfGCcoverage.columns=["ID", "start", "stop", "nothing1", "GCpercent", "nothing2"]
        dfReadCoverage = pd.read_csv(self.readCoverageFile)
        dfGCRC= pd.merge(left=dfReadCoverage, right=dfGCcoverage, how='left', left_on='pos', right_on='start')
        
            # plot GC coverage
        logging.info(INDENT*'-' + "--plotting")
        gcPlotFile = os.path.join(resultFolder, self.projectID + "_w" + str(self.windowSize) + "s" + str(self.stepSize) + self.md5string + ".png")
        gcPlotTitle = self.CLASSID + "_" + self.projectID + "_w" + str(self.windowSize) + "s" + str(self.stepSize)
        
        p = (p9.ggplot(data=dfGCRC,
                   mapping=p9.aes(x='normcoverage',
                                  y="GCpercent", colour="datasource"))
            + p9.geom_point( alpha=0.25, size=0.25) + p9.labs(title=gcPlotTitle) 
        )
        p.save(filename = gcPlotFile, height=self.PLOTHEIGHT, width=self.PLOTWIDTH, dpi=self.PLOTDPI)   
                         
            
            
        logging.info(INDENT*'-' + "finishing")



        
            
        

    def shortDescription(self):
        print('calculate GC coverage for fasta file with an optional sliding window')


    def longDescription(self):
        print('calculate GC coverage for fasta file.')
        print('by default, the GC coverage is calculated at each site')
        print('but a window size and step interval can also be specified to ')
        print('calculate an moving average.')
        print('  window size: -w / --window_size')
        print('    step size: -s / -- step_size')
        print('')      
        print('The output is in BED format. If an output file is not specified, ')
        print('the output file the same as the input file with a bed extension.')
        print('')


    def parseJSON(self, stepJSON):
        '''
        require:  `parameters`, `inFolder`, `inFiles`
        optional: `outFolder`
        '''
        logging.info(INDENT*'-' + "parsing JSON")
        
        logging.debug(INDENT*'-' + "--parsing [" + abstractStep.AbstractStep.PARAMID + "] string")
        self.paramString = stepJSON.get(abstractStep.AbstractStep.PARAMID)
        if(self.paramString is not None):
            self.paramString = stepJSON[abstractStep.AbstractStep.PARAMID]
            logging.debug(INDENT*'-' + "--found [" + abstractStep.AbstractStep.PARAMID + "] string + <" + self.paramString + ">")
        else:
            raise Exception("missing [" + abstractStep.AbstractStep.PARAMID + "] string in JSON line")
        
        logging.debug(INDENT*'-' + "--parsing [" + abstractStep.AbstractStep.INFILESID + "] string")
        self.inputFiles = stepJSON.get(abstractStep.AbstractStep.INFILESID)
        if(self.inputFiles is not None):
            self.inputFiles = stepJSON[abstractStep.AbstractStep.INFILESID]       
            logging.debug(INDENT*'-' + "--found [" + abstractStep.AbstractStep.INFILESID + "] string + <" + "-".join(stepJSON['inFiles']) + ">")
        else:
            raise Exception("missing [" + abstractStep.AbstractStep.INFILESID + "] string in JSON")
        
        logging.debug(INDENT*'-' + "--parsing [" + abstractStep.AbstractStep.INFOLDERID + "] string")
        self.inFolder = stepJSON.get(abstractStep.AbstractStep.INFOLDERID)
        if(self.inFolder is not None):
            self.inFolder = stepJSON[abstractStep.AbstractStep.INFOLDERID]
            logging.debug(INDENT*'-' + "--found [" + abstractStep.AbstractStep.INFOLDERID + "] string + <" + self.inFolder + ">")
        else:
            raise Exception("missing [" + abstractStep.AbstractStep.INFOLDERID + "] string in JSON")
        
        logging.debug(INDENT*'-' + "--parsing [" + abstractStep.AbstractStep.OUTFOLDERID + "] string")
        self.outFolder = stepJSON.get(abstractStep.AbstractStep.OUTFOLDERID)
        if(self.outFolder is not None):
            self.outFolder = stepJSON[abstractStep.AbstractStep.OUTFOLDERID]
            logging.debug(INDENT*'-' + "--found [" + abstractStep.AbstractStep.OUTFOLDERID + "] string + <" + self.outFolder + ">")
        else:
            self.outFolder = self.OUTPUTFOLDER
            logging.debug(INDENT*'-' + "--did not find [" + abstractStep.AbstractStep.OUTFOLDERID + "] string, setting to <" + self.outFolder + ">")
            
            
        logging.info(INDENT*'-' + "finished parsing JSON")

                

    def parseParameterString(self):
        '''
        require: 
            'read coverage folder'
            'gc coverage folder'
            
            'window_size` and `step_size` parameters in the string
             both need to be present (i.e., you can't specify a `window_size` without
             specifying the `step_size` )
            
        '''
        logging.info(INDENT*'-' + "parsing parameters strings")
        params = self.paramString.split(",")
        for param in params:
            
            if self.STEPSIZESHORT in param or self.STEPSIZELONG in param:
                if self.STEPSIZELONG in param:
                    self.stepSize = int(param.split(self.STEPSIZELONG)[1].strip())
                else:
                    self.stepSize = int(param.split(self.STEPSIZESHORT)[1].strip())
                logging.info(INDENT*'-' + "step size set to <" + str(self.stepSize) + ">")
                
            elif self.GCFILESHORT in param or self.GCFILELONG in param:
                if self.GCFILELONG in param:
                    self.gcCoverageFile = param.split(self.GCFILELONG)[1].strip()
                else:
                    self.gcCoverageFile = param.split(self.GCFILESHORT)[1].strip()
                logging.info(INDENT*'-' + "gc coverage file set to <" + self.gcCoverageFile + ">")
            
            elif self.READFILESHORT in param or self.READFILELONG in param:
                if self.READFILELONG in param:
                    self.readCoverageFile = param.split(self.READFILELONG)[1].strip()
                else:
                    self.readCoverageFile = param.split(self.READFILESHORT)[1].strip()
                logging.info(INDENT*'-' + "read coverage folder set to <" + self.readCoverageFile + ">")
            elif self.WINSIZESHORT in param or self.WINSIZELONG in param:
                if self.WINSIZELONG in param:
                    self.windowSize = int(param.split(self.WINSIZELONG)[1].strip())
                else:
                    self.windowSize = int(param.split(self.WINSIZESHORT)[1].strip())
                logging.info(INDENT*'-' + "window size set to <" + str(self.windowSize) + ">")
            
            elif self.STEPSIZESHORT in param or self.STEPSIZELONG in param:
                if self.STEPSIZELONG in param:
                    self.stepSize = int(param.split(self.STEPSIZELONG)[1].strip())
                else:
                    self.stepSize = int(param.split(self.STEPSIZE)[1].strip())
                logging.info(INDENT*'-' + "step size set to <" + str(self.stepSize) + ">")
                
            if self.BAMFILEFOLDERSHORT in param or self.BAMFILEFOLDERLONG in param:
                if self.BAMFILEFOLDERLONG in param:
                    self.bamFileFolder = param.split(self.BAMFILEFOLDERLONG)[1].strip()
                else:
                    self.bamFileFolder = param.split(self.BAMFILEFOLDERLONG)[1].strip()
                logging.info(INDENT*'-' + "bamFileFolder set to <" + str(self.bamFileFolder) + ">")

            
            
        if self.gcCoverageFile == "":
            logging.error("you need to specify a folder containing the results of the GC coverage analysis")
            raise Exception("you need to specify a folder containing the results of the GC coverage analysis")
        if self.readCoverageFile == "":
            logging.error("you need to specify a folder containing the results of the read coverage analysis")
            raise Exception("you need to specify a folder containing the results of the read coverage analysis")
        if self.bamFileFolder == "":
            logging.error("you need to specify a folder containing the BAM files")
            raise Exception("you need to specify a folder containing the BAM files")
        
           

        

        
        