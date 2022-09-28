from pypesteps import abstractStep

'''
Created on Dec 18, 2020

@author:     simon rayner
@contact:    simon.rayner@medisin.uio.no
'''
import os

from Bio import SeqIO

import pandas as pd
import numpy as np
import plotnine as p9
import glob


import logging

logger = logging.getLogger(__name__)
INDENT = 6


class StepSNVProcessShorahResults(abstractStep.AbstractStep):
    '''
    classdocs
    This step processes the estimated SNV data from a set of BAM files that were generated 
    by the `shorah` software package.
    Shorah is well suited to SNV calling in viruses where you can have variable rates of evolution 
    and read coverage due sequencing efficiency across the alignment.
    
    An analysis generates the following file structure
    
     30210555 Dec 18 09:49 A1_S1.hq.txt.sorted.cor.fas

       -   122880 Dec 17 23:35 corrected
       -   102400 Dec 17 23:25 debug
       -   126976 Dec 17 23:39 freq
       -   118784 Dec 17 23:39 raw_reads
       -   106496 Dec 17 23:35 sampling
       -        0 Dec 18 10:17 snv
       -   118784 Dec 17 23:35 support
       
            28197 Dec 17 09:05 coverage.txt
             6516 Dec 18 09:49 proposed.dat
           932421 Dec 18 09:49 raw_snv.txt
        316609399 Dec 17 09:05 reads.fas
           352124 Dec 18 10:17 shorah.log

    In this step, we only consider the SNV calling summarised in `snv/SNVs_0.010000_final.csv`
    We use the CSV rather than the VCF as it is simpler to parse and the information is identical
    
    To Do:
    Add parameters to allow user to specify nt start and stop
    Add parameters to allow user to specify plot settings (size and dpi)
    Add code to handle tick intervals
     
    '''
    CLASSID             = "StepSNVProcessShorahResults"
    OUTPUTFOLDER        = "processedsnvs"
    REFFASTASHORT       = "-r"
    REFFASTALONG        = "--ref_fasta"
    BAMFILEFOLDERSHORT  = "-b"
    BAMFILEFOLDERLONG   = "--bam_file_folder"

    SNVFILEEND          = "snv/SNVs_0.010000_final.csv"
    PLOTHEIGHT          = 5
    PLOTWIDTH           = 20
    PLOTUNITS           = 'in'
    PLOTDPI             = 300
    XVAR                = "pos"
    YVAR                = "SNVs"
    

    def __init__(self, refFastA="", bamFileFolder=""):
        '''
        Constructor
        '''
        self.refFastA = refFastA
        self.bamFileFolder = bamFileFolder
        pass

        
    def checkInputData(self):
        '''
        build file paths and check all input resources exist
        filepath can be absolute or relative (to Project Root)
        '''
        if( os.path.exists(os.path.join(self.projectRoot, self.inFolder)) is False):
            raise Exception("input folder <" + os.path.join(self.projectRoot, self.inFolder) + "> not found")
        

        
        if len(self.inputFiles) == 1 & (not self.inputFiles[0]):
            self.inputFiles = glob.glob(os.path.join(self.projectRoot, self.bamFileFolder) + os.path.sep + "*gen__trim_paired__sorted.bam")
        for inputFile in self.inputFiles:            
            basename = os.path.splitext(os.path.basename(inputFile))[0]
            resultFolder = os.path.join(self.projectRoot, self.inFolder)
            snv_vcf_file = os.path.join(resultFolder, basename, self.SNVFILEEND)        
                                         
            if os.path.exists(snv_vcf_file) is False:
                #raise RuntimeError ("input file <" + snv_vcf_file + "> not found")
                logging.warn("input file <" + snv_vcf_file + "> not found")
            else:
                logging.info(INDENT*'-' + "found SNV VCF file <" + snv_vcf_file + ">")
        
        ## check reference FastA file
        if os.path.exists(self.refFastA) == False:
            raise RuntimeError ("reference FastA file <" + self.refFastA + "> not found")
            logging.error("reference FastA file <" + self.refFastA + "> not found")
        else:
            logging.info(INDENT*'-' + "found reference FastA file <" + self.refFastA + ">")
            
        
        
        
        
    def execute(self):
        '''
        contains the main operations for the step
        '''
        logger.info(INDENT*'-' + "executing step")
        sourceFolder = os.path.join(self.projectRoot, self.inFolder)
        resultFolder = os.path.join(self.projectRoot, self.outFolder)
        
        logging.info(INDENT*'-' + "--results will be written to output folder <" + resultFolder + ">")
        if not os.path.exists(resultFolder):
            logging.info(INDENT*'-' + "----folder doesn't exist, creating")
            os.makedirs(resultFolder)        
        

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
        dfSNV = pd.DataFrame(list(range(1,genomeLen+1)), columns = ['nt']).astype(int)
        bamFileCount = 1
        
        # The following is for plot cosmetics. 
        offset = 1  # the y distance for no SNV in a single sample
        dOffset = 1 # the y distance between successive samples on the plot
        delta = 2   # the y distance for SNV in a single sample
        
        dfAll = pd.DataFrame(columns=["Pos","frqMean","snvplot","datasource"])
        for inputFile in self.inputFiles:
            
            # for each VCF file
            #    1. Skip first 18 header lines
            #    2. grab POS (#2) and Frq1 (#4), Frq2 (#5) & Frq3 (#6) columns
            #    3. add command to execute shorah with absolute filepaths
            
            # I'm not sure what we gain from working with the VCF files as 
            # it requires more work to parse and the CSV contains the same information.
            # The only advantage is that the VCF format is fixed, we know what each column
            # contains.
            # However, for now, work with CSV as we are just trying to figure out the data
            # 
            basename = os.path.splitext(os.path.basename(inputFile))[0]

            snv_vcf_file = os.path.join(sourceFolder, basename, self.SNVFILEEND) 
            if os.path.exists(snv_vcf_file) is False:
                #raise RuntimeError ("input file <" + snv_vcf_file + "> not found")
                logging.warn("input file <" + snv_vcf_file + "> not found")
                continue
            
            #dfTempVCF = pd.read_csv(snv_vcf_file, delimiter="\t", skiprows=17)
            dfTempCSV = pd.read_csv(snv_vcf_file)
            dfTempCSV['Frq1'] = pd.to_numeric(dfTempCSV['Frq1'], errors='coerce')
            dfTempCSV['Frq2'] = pd.to_numeric(dfTempCSV['Frq2'], errors='coerce')
            dfTempCSV['Frq3'] = pd.to_numeric(dfTempCSV['Frq3'], errors='coerce')
            dfTempCSV['frqMean'] = np.nanmean(dfTempCSV.loc[:, 'Frq1':'Frq3'], axis=1)
    
            dfSNV = pd.merge(left=dfSNV, right=dfTempCSV.loc[:, ['Pos','frqMean']], how='left', left_on='nt', right_on='Pos')
            dfSNV['snvplot'] = np.where(dfSNV['Pos'].isnull(), offset, offset+delta)
            dfSNV = dfSNV.rename(columns={"snvplot":"snvplot_"+ basename, "frqMean":"frqMean_"+ basename, "Pos":"Pos_"+ basename,})
            
            
            dfTemp = dfTempCSV.loc[:, ['Pos','frqMean']]
            dfTemp['snvplot'] = offset+delta
            dfTemp['datasource'] = basename
            dfAll = dfAll.append(dfTemp)

            offset += dOffset
            
            bamFileCount+=1

        # write the unified dataframe as CSV
        plotFileAsCSV = os.path.join(resultFolder, self.projectID + "__SNVs__"+ self.md5string + ".csv")
        logging.info(INDENT*'-' + "--saving combined SNV data to <" + plotFileAsCSV +">")
        dfAll.to_csv(plotFileAsCSV)
        
        
        
        # plot the SNV data
        logging.info(INDENT*'-' + "--plotting combined SNV data")
        snvPlotFile = os.path.join(resultFolder, self.projectID + "__SNVs__"+ self.md5string + ".png")
        logging.info(INDENT*'-' + "--plot file is to <" + snvPlotFile +">")
        
        xAxisStart = 0
        xAxisEnd = genomeLen
        xAxisNoOfTicks = 10
        xAxisInterval = int((xAxisEnd - xAxisStart)/xAxisNoOfTicks)
        dfAll['Pos'] = pd.to_numeric(dfAll['Pos'])
        dfAll['frqMean'] = dfAll['frqMean']*100.0

        p = (p9.ggplot(data=dfAll, mapping=p9.aes(x='Pos', y='snvplot', color='datasource', size = 'frqMean')) \
             + p9.geom_point( alpha=0.1) + p9.scale_size(range = [0, 10]) \
             + p9.labs(title=self.projectID)
        + p9.scale_x_continuous(name=self.XVAR) + p9.ylab(self.YVAR)) 
        #+ p9.scale_x_continuous(name=self.XVAR, breaks=np.arange(0, 30000, 5000), limits=[0, 30000] ) + p9.ylab(self.YVAR)) 

        p.save(filename = snvPlotFile, height=self.PLOTHEIGHT, width=self.PLOTWIDTH,  dpi=self.PLOTDPI)   
        
        logger.info(INDENT*'-' + "done")
        

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
        look for 
        1. `window_size` and `step_size` parameters in the string
            both need to be present (i.e., you can't specify a `window_size` without
            specifying the `step_size` 
        2.  `output_file`, otherwise set the output filename to the input with .BED extension
        '''
        logging.info(INDENT*'-' + "parsing parameters strings")
        params = self.paramString.split(",")
        for param in params:
            if self.REFFASTASHORT in param or self.REFFASTALONG in param:
                if self.REFFASTALONG in param:
                    self.refFastA = param.split(self.REFFASTALONG)[1].strip()
                else:
                    self.refFastA = param.split(self.REFFASTALONG)[1].strip()
                logging.info(INDENT*'-' + "reference FastA file set to <" + str(self.refFastA) + ">")
            
            if self.BAMFILEFOLDERSHORT in param or self.BAMFILEFOLDERLONG in param:
                if self.BAMFILEFOLDERLONG in param:
                    self.bamFileFolder = param.split(self.BAMFILEFOLDERLONG)[1].strip()
                else:
                    self.bamFileFolder = param.split(self.BAMFILEFOLDERLONG)[1].strip()
                logging.info(INDENT*'-' + "bamFileFolder set to <" + str(self.bamFileFolder) + ">")
                
        if self.bamFileFolder == "":
            logging.error("you need to specify a folder containing the BAM files")
            raise Exception("you need to specify a folder containing the BAM files")

        
           

        

        
        