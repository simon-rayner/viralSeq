from pypesteps import abstractStep

'''
Created on Dec 18, 2020

@author:     simon rayner
@contact:    simon.rayner@medisin.uio.no
'''
import os
from Bio import SeqIO
from Bio.SeqUtils import GC
import csv
import pandas as pd
import plotnine as p9


import logging

logger = logging.getLogger(__name__)
INDENT = 6


class StepGCReadCoverage(abstractStep.AbstractStep):
    '''
    classdocs
    This calculates the GC coverage across the specified fasta file
    a sliding window and step size must be specified to calculate average value
    (using -w/--window_size & -s/--step_size)
    
    output is written to an output file in BED format
    
    To do: add parameters to set x axis plot range in GC plot
    '''
    CLASSID             = "StepGCReadCoverage"
    OUTPUTFOLDER        = "gccoverage"
    WINSIZESHORT        = "-w"
    WINSIZELONG         = "--window_size"
    STEPSIZESHORT       = "-s"
    STEPSIZELONG        = "--step_size"
    
    PLOTHEIGHT          = 5
    PLOTWIDTH           = 10
    PLOTUNITS           = 'in'
    PLOTDPI             = 1000
    XVAR                = "pos"
    YVAR                = "gcpercent"
    

    def __init__(self, windowSize=0, stepSize=0):
        '''
        Constructor
        '''

        self.windowSize = windowSize
        self.stepSize = stepSize

        
    def checkInputData(self):
        '''
        build file paths and check all input resources exist
        '''
        
        if( os.path.exists(os.path.join(self.projectRoot, self.inFolder)) is False):
            logging.error("input folder <" + os.path.join(self.projectRoot, self.inFolder) + "> not found")
            raise Exception("input folder <" + os.path.join(self.projectRoot, self.inFolder) + "> not found")
        
        inputFolder = os.path.join(self.projectRoot, self.inFolder)
        
        for inputFile in self.inputFiles:    
            inFile = os.path.join(os.path.join(inputFolder,inputFile))
       
            if os.path.exists(inFile) == False:
                logging.error("input file <" + inFile + "> not found")
                raise RuntimeError ("input file <" + inFile + "> not found")
            else:
                logging.info(INDENT*'-' + "found input file <" + inFile + ">")
        
        ## check window size and step size are positive integers (this must be specified for GC coverage calculation)
        if(self.windowSize > 0 and self.stepSize > 0):
            logging.info("Read Coverage will be calculated using a sliding window of <" +\
                          str(self.windowSize) + "> nt and a step size of <" + str(self.stepSize) + "> nt")
        elif(self.windowSize <= 0 and self.stepSize <= 0):
            logging.error("No sliding window selected")            
            raise Exception("No sliding window selected") 
        else:
            logging.error("both window and step size must be specified and > 0: (found window size <" \
                            + str(self.windowSize) + " and step size <" + str(self.stepSize) + ">)")
            raise Exception("both window and step size must be specified and > 0: (found window size <" \
                            + str(self.windowSize) + " and step size <" + str(self.stepSize) + ">)")       
        
        
        
    def execute(self):
        '''
        contains the main operations for the step
        '''
        logger.info(INDENT*'-' + "executing step")

        inputFolder = os.path.join(self.projectRoot, self.inFolder)
        resultFolder = os.path.join(self.projectRoot, self.outFolder)    
            
        for inputFile in self.inputFiles:
        # for each fasta file
            #   read file
            inFile = os.path.join(os.path.join(inputFolder,inputFile))
            logger.info(INDENT*'-' + "processing file <" + inFile + ">")
            for record in SeqIO.parse(inFile, "fasta"):
                genomeSeq = record.seq
                genomeID  = record.id
                genomeLen = len(genomeSeq)
                logger.info(INDENT*'-' + "-- found record ID <" + genomeID + ">")
                
                # calc GC content
                posInGenome = 0
                gcData = []
                logging.debug(INDENT*'-' + "-- calculating sliding window")
                while posInGenome < genomeLen - self.windowSize:
                    gcData.append({self.XVAR: (posInGenome + self.windowSize/2), self.YVAR: (GC(genomeSeq[posInGenome:posInGenome+self.windowSize]))})
                    posInGenome += self.stepSize
                dfGCdata = pd.DataFrame.from_dict(gcData)
                logging.debug(INDENT*'-' + "-- done")

                # create output filename
                outputFolder = os.path.join(os.path.dirname(inputFile), self.OUTPUTFOLDER)
                logging.info(INDENT*'-' + "--GC results will be written to output folder <" + resultFolder + ">")
                if not os.path.exists(resultFolder):
                    logging.info(INDENT*'-' + "----folder doesn't exist, creating")
                    os.makedirs(resultFolder)
                
                inBaseName = os.path.splitext(os.path.basename(inputFile))[0]
                gcResultsFile = os.path.join(resultFolder, inBaseName + "__w" + str(self.windowSize) + "_s" + str(self.stepSize) + "__" + self.md5string + ".bed")
                logging.info(INDENT*'-' + "--GC output BED file is <" + gcResultsFile + ">")
                logging.info(INDENT*'-' + "--writing")
                
                with open(gcResultsFile, 'wt+') as bedfile:
                    bedwriter = csv.writer(bedfile, delimiter='\t')
                    bedwriter.writerow(["track name=GC percentage description = sliding window " \
                                       + str(self.windowSize) + "nt/step size " + str(self.stepSize) + "nt"])
                    for entry in gcData:                        
                        bedwriter.writerow([genomeID, int(entry['pos']), int(entry['pos']), ".", str(entry['gcpercent']), "."])
                logging.info(INDENT*'-' + "--done")
                
                # plot GC coverage
                logging.info(INDENT*'-' + "--plotting")
                gcPlotFile = os.path.join(resultFolder, inBaseName + "__w" + str(self.windowSize) + "_s" + str(self.stepSize) + "__" + self.md5string + ".png")
                logging.info(INDENT*'-' + "--plot file is <" + gcPlotFile + ">")
                gcPlotTitle = inBaseName + "_w" + str(self.windowSize) + "s" + str(self.stepSize)
                
                p = (p9.ggplot(data=dfGCdata,
                           mapping=p9.aes(x=self.XVAR,
                                          y=self.YVAR, colour=self.YVAR))
                    + p9.geom_point( alpha=0.25, size=0.25) + p9.labs(title=gcPlotTitle) 
                    + p9.scale_x_continuous(name=self.XVAR, limits=[0, 30000] ) + p9.ylab(self.YVAR)
                )
                p.save(filename = gcPlotFile, height=self.PLOTHEIGHT, width=self.PLOTWIDTH, dpi=self.PLOTDPI)   
                             
                
                logging.info(INDENT*'-' + "finishing")



        
            
        

    def shortDescription(self):
        print('calculate GC coverage for fasta file with a sliding window')


    def longDescription(self):
        print('calculate GC coverage for fasta file across a sliding window.')
        print('')
        print('  window size: -w / --window_size')
        print('    step size: -s / -- step_size')
        print('')      
        print('The output is in BED format. If an output file is not specified, ')
        print('the output file the same as the input file with a bed extension.')
        print('')


    def parseJSON(self, stepJSON):

        # check the 
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
        look for 
        1. `window_size` and `step_size` parameters in the string
            both need to be present (i.e., you can't specify a `window_size` without
            specifying the `step_size` )
        2.  `output_file`, otherwise set the output filename to the input with .BED extension
        '''
        logging.info(INDENT*'-' + "parsing parameters strings")
        params = self.paramString.split(",")
        for param in params:
            if self.WINSIZESHORT in param or self.WINSIZELONG in param:
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
            
        
           

        
    def gcPrint(self):
        pass
        
        
        