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

import numpy as np
import pandas as pd
import plotnine as p9
from sklearn import preprocessing
from Bio import SeqIO

import logging

logger = logging.getLogger(__name__)
INDENT = 6


class StepBAMReadCoverage(abstractStep.AbstractStep):
    '''
    classdocs
    This calculates the read coverage for a set of BAM files
    The analysis is done in two steps:
    
        1. SAMTools is used to calculate the read coverage/nt, 
        2. This data to calculate a sliding window average (if requested)
    
    SAMTools is used because it's faster
    
    required parameters
        the FASTA file that was used for the alignment 
        (This is used to get the genome length)
    
    optional parameters
    The user may specify:
        a sliding window and step size can also be specified to calculate average value.
        both values need to be specified
        (-w/--window_size & (-s/--step_size)
    
        the location of the samtools software package (some Python installations have trouble locating installs)
        (-s/--software_location)
    
    To do: generate integrated read coverage plot
    '''
    CLASSID             = "StepBAMReadCoverage"
    OUTPUTFOLDER        = "bamreadcoverage"
    REFFASTASHORT       = "-r"
    REFFASTALONG        = "--ref_fasta"
    WINSIZESHORT        = "-w"
    WINSIZELONG         = "--window_size"
    STEPSIZESHORT       = "-s"
    STEPSIZELONG        = "--step_size"
    SOFTWARELOCSHORT    = "-p"
    SOFTWARELOCLONG     = "--path_to_software"
    
    PLOTHEIGHT          = 3
    PLOTWIDTH           = 10
    PLOTUNITS           = 'in'
    PLOTDPI             = 1000
    XVAR                = "pos"
    YVAR                = "readcoverage"
    SAMCOL0             = "id"
    NORMCOL             = "normcoverage"
    STEPNORM            = "steppednorm"
    

    def __init__(self, windowSize=0, stepSize=0, softwarePath= "samtools", refFastA=""):
        '''
        Constructor
        '''
        self.windowSize = windowSize
        self.stepSize = stepSize
        self.softwarePath = softwarePath
        self.refFastA = refFastA

        
    def checkInputData(self):
        '''
        build file paths and check all input resources exist
        filepath can be absolute or relative (to Project Root)
        '''

        if( os.path.exists(os.path.join(self.projectRoot, self.inFolder)) is False):
            raise Exception("input folder <" + os.path.join(self.projectRoot, self.inFolder) + "> not found")
        
        bamFileFolder = os.path.join(self.projectRoot, self.inFolder)
        
        for inputFile in self.inputFiles:            
            if os.path.exists(os.path.join(bamFileFolder, inputFile)) == False:
                logging.error("input file <" + os.path.exists(os.path.join(bamFileFolder, inputFile)) + "> not found")
                raise RuntimeError ("input file <" + os.path.join(bamFileFolder, inputFile) + "> not found")
            else:
                logging.info(INDENT*'-' + "found input file <" + os.path.join(bamFileFolder, inputFile) + ">")
        
        ## check reference FastA file
        if os.path.exists(self.refFastA) == False:
            logging.error("reference FastA file <" + self.refFastA + "> not found")
            raise RuntimeError ("reference FastA file <" + self.refFastA + "> not found")
        else:
            logging.info(INDENT*'-' + "found reference FastA file <" + self.refFastA + ">")
                        
        ## check window size and step size are positive integers (optional for Read coverage calculation)
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
            
            
        
        
        
        
    def execute(self):
        '''
        contains the main operations for the step
        '''
        logger.info(INDENT*'-' + "executing step")

        # load the genome to get the number of nucleotides
        recordCount = 0
        for record in SeqIO.parse(self.refFastA, "fasta"):
            genomeSeq = record.seq
            genomeID  = record.id
            genomeLen = len(genomeSeq)
            recordCount += 1
            
        if recordCount > 1:
            logging.warn(INDENT*'-' + "----fasta file contains more than one record. Using the last loaded record (#" + recordCount + ")")
            logging.warn(INDENT*'-' + "----Using the last loaded record <" + genomeID + "> which is <" + genomeLen + "> nt" )
                    
        # create dataframe with genomeLen rows
        a = np.arange(self.windowSize/2, 1000, self.stepSize)
        dfAllCSV = pd.DataFrame(np.arange(self.windowSize/2, 1000, self.stepSize), columns = ['nt']).astype(int) # for CSV file
        dfAllPlot = pd.DataFrame(columns=[self.XVAR, self.YVAR, self.NORMCOL]) # for plotting
        
        bamFileFolder = os.path.join(self.projectRoot, self.inFolder)
        resultFolder = os.path.join(self.projectRoot, self.outFolder)
        
        logging.info(INDENT*'-' + "--results will be written to output folder <" + resultFolder + ">")
        if not os.path.exists(resultFolder):
            logging.info(INDENT*'-' + "----folder doesn't exist, creating")
            os.makedirs(resultFolder)        

        
        offset = 1  # the y distance for no SNV in a single sample
        dOffset = 4 # the y distance between successive samples on the plot
        delta = 2   # the y distance for SNV in a single sample
        
        for inputFile in self.inputFiles:
            basename = os.path.splitext(os.path.basename(inputFile))[0]            
        # for each BAM file
            #   1. generate coverage/nt using SAMTools
            #   2. generate sliding window coverage if requested
            
            
            
            # 1. use SAMTools to get read coverage at each base
            bamFile = os.path.join(bamFileFolder, inputFile)
            ntCovFile = os.path.join(resultFolder, os.path.splitext(os.path.basename(bamFile))[0] + "_ntcov.tsv")
            stderrFile = os.path.basename(bamFile) + "_ntcov.stderr"

            command = self.softwarePath + ' depth -a ' + bamFile 
            logging.debug(INDENT*'-' + "--SAMTools command is <"+ command + ">")
            with open(ntCovFile,"wb") as fout, open(stderrFile,"wb") as err:
                process = subprocess.Popen(shlex.split(command), stdout=subprocess.PIPE)

                while True:
                    output = process.stdout.readline()
                    if process.poll()==0:
                        break
                    if output:
                        fout.write(output)
            
                logger.info(INDENT*'-' + "--process finished with return code <" + str(process.poll()) + ">")
            
            
            # 2. if window parameters have been set, calculate sliding window coverage
            if(self.windowSize > 0):
                dfThisBAMCoverage = pd.read_csv(ntCovFile, sep='\t', header=None)
                #dfThisBAMCoverage.columns = [self.SAMCOL0, self.XVAR, self.YVAR]
                meanCovWin = []
                
                posInGenome = int(self.windowSize/2)
                while (posInGenome < len(dfThisBAMCoverage)- self.windowSize/2):
                    meanCovWin.append({self.XVAR: posInGenome, self.YVAR: (dfThisBAMCoverage.iloc[posInGenome:posInGenome+self.windowSize+1][2].mean())})
                    posInGenome += self.stepSize
                dfThisBAMWin = pd.DataFrame.from_dict(meanCovWin)
                
                # calculate mean of column + normalise the reads 
                min_max_scaler = preprocessing.MinMaxScaler()
                #  x= dfThisBAMWin[self.YVAR]
                
                dfThisBAMWin[self.NORMCOL] = min_max_scaler.fit_transform(dfThisBAMWin[[self.YVAR]].astype(float)) # + offset
                dfThisBAMWin['datasource'] = basename

                dfAllCSV = pd.merge(left=dfAllCSV, right=dfThisBAMWin.loc[:, [self.XVAR, self.YVAR, self.NORMCOL]], how='left', left_on='nt', right_on=self.XVAR)
                dfAllCSV = dfAllCSV.rename(
                    columns={self.YVAR:self.YVAR + "_" + basename, self.NORMCOL:self.NORMCOL + "_" + basename, self.XVAR: self.XVAR + "_" + basename,})
                #dfThisBAMWin = dfThisBAMWin.rename(columns={"snvplot":"snvplot_"+ basename, "frqMean":"frqMean_"+ basename, "Pos":"Pos_"+ basename,})

                # x = df[['score']].values.astype(float)
                # min_max_scaler = preprocessing.MinMaxScaler()
                # x_scaled = min_max_scaler.fit_transform(x)
                # df_normalized = pd.DataFrame(x_scaled)
                dfThisBAMWin[self.STEPNORM] = dfThisBAMWin[self.NORMCOL] + offset
                dfAllPlot = dfAllPlot.append(dfThisBAMWin)  
                                              
                # create output filename
                inBaseName = os.path.splitext(os.path.basename(bamFile))[0]
                outputFolder = os.path.join(self.projectRoot, self.outFolder)
                logging.info(INDENT*'-' + "--read coverage results will be written to output folder <" + outputFolder + ">")
                if not os.path.exists(outputFolder):
                    logging.info(INDENT*'-' + "----folder doesn't exist, creating")
                    os.makedirs(outputFolder)
                    
                ntCovFileWinAv = os.path.splitext(os.path.basename(bamFile))[0] + "__w" + str(self.windowSize) + "_s" + str(self.stepSize) + "__" + self.md5string + ".bed"
                ntCovFileWinAv = os.path.join(outputFolder, ntCovFileWinAv)
                logging.info(INDENT*'-' + "--read coverage output BED file is <" + ntCovFileWinAv + ">")
                logging.info(INDENT*'-' + "--writing")
                
                with open(ntCovFileWinAv, 'wt+') as bedfile:
                    bedwriter = csv.writer(bedfile, delimiter='\t')
                    bedwriter.writerow(["track name=read coverage description = sliding window " \
                                       + str(self.windowSize) + "nt/step size " + str(self.stepSize) + "nt"])
                    for index, row in dfThisBAMWin.iterrows():                       
                        bedwriter.writerow([inBaseName, str(row[self.XVAR]), str(row[self.XVAR]), ".", str(row[self.YVAR]), "."])
                logging.info(INDENT*'-' + "--done")
            
            # plot read coverage for this BAM file
            logging.info(INDENT*'-' + "--plotting")
            gcPlotFile = os.path.join(outputFolder, inBaseName + "_w" + str(self.windowSize) + "s" + str(self.stepSize) + self.md5string + ".png")
            gcPlotTitle = self.CLASSID + "_" + inBaseName + "_w" + str(self.windowSize) + "s" + str(self.stepSize)
            
            p = (p9.ggplot(data=dfThisBAMWin,
                       mapping=p9.aes(x=self.XVAR,
                                      y=self.YVAR, colour=self.YVAR))
                + p9.geom_point( alpha=0.1, size=0.25) + p9.labs(title=gcPlotTitle) 
            )
            p.save(filename = gcPlotFile, height=self.PLOTHEIGHT, width=self.PLOTWIDTH, dpi=self.PLOTDPI)   
                         
            offset += dOffset
            
        logging.info(INDENT*'-' + "finishing")

        # write out single file containing normalised read coverage for all files
        allDataAsCSV = os.path.join(resultFolder, self.projectID + "__normreads__" 
                                     + "__w" + str(self.windowSize) + "_s" + str(self.stepSize) + "__"+ self.md5string + ".csv")
        logging.info(INDENT*'-' + "--saving combined data to <" + allDataAsCSV +">")
        dfAllCSV.to_csv(allDataAsCSV)
        plotDataAsCSV = os.path.join(resultFolder, self.projectID + "__normreads__" 
                                     + "__w" + str(self.windowSize) + "_s" + str(self.stepSize) + "__plot__"+ self.md5string + ".csv")
        logging.info(INDENT*'-' + "--saving plot data to <" + allDataAsCSV +">")
        dfAllPlot.to_csv(plotDataAsCSV)
                
        # plot read coverage for all BAM files
        logging.info(INDENT*'-' + "--plotting combined SNV data")
        covPlotFile = os.path.join(resultFolder, self.projectID + "__normreads__" 
                                     + "__w" + str(self.windowSize) + "_s" + str(self.stepSize) + "__"+ self.md5string + ".png")
        logging.info(INDENT*'-' + "--plot file is to <" + covPlotFile +">")

        dfAllPlot[self.XVAR] = pd.to_numeric(dfAllPlot[self.XVAR])
        p = (p9.ggplot(data=dfAllPlot, mapping=p9.aes(x=self.XVAR, y=self.STEPNORM, color='datasource', size = self.XVAR)) \
             + p9.geom_point( alpha=0.1) + p9.scale_size(range = [0, 1]) \
             + p9.labs(title=self.projectID)
        + p9.scale_x_continuous(name=self.XVAR) + p9.ylab(self.YVAR)) 
        #+ p9.scale_x_continuous(name=self.XVAR, breaks=np.arange(0, 30000, 5000), limits=[0, 30000] ) + p9.ylab(self.YVAR)) 

        p.save(filename = covPlotFile, height=self.PLOTHEIGHT, width=self.PLOTWIDTH,  dpi=self.PLOTDPI)   

        
            
        

    def shortDescription(self):
        print('calculate read coverage for fasta file with an optional sliding window')


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
           `window_size` and `step_size` parameters in the string
            both need to be present (i.e., you can't specify a `window_size` without
            specifying the `step_size` 
            
        optional:
        `software_location`
        '''
        logging.info(INDENT*'-' + "parsing parameters strings")
        params = self.paramString.split(",")
        for param in params:
            param = param.strip()
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
                    self.stepSize = int(param.split(self.STEPSIZESHORT)[1].strip())
                logging.info(INDENT*'-' + "step size set to <" + str(self.stepSize) + ">")
                
            elif self.SOFTWARELOCSHORT in param or self.SOFTWARELOCLONG in param:
                if self.SOFTWARELOCLONG in param:
                    self.softwarePath = param.split(self.SOFTWARELOCLONG)[1].strip()
                else:
                    self.softwarePath = param.split(self.SOFTWARELOCSHORT)[1].strip()
                logging.info(INDENT*'-' + "samtools software location set to <" + str(self.softwarePath) + ">")
            elif self.REFFASTASHORT in param or self.REFFASTALONG in param:
                if self.REFFASTALONG in param:
                    self.refFastA = param.split(self.REFFASTALONG)[1].strip()
                else:
                    self.refFastA = param.split(self.REFFASTALONG)[1].strip()
                logging.info(INDENT*'-' + "reference FastA file set to <" + str(self.refFastA) + ">")
            
            
        
           

        

        
        