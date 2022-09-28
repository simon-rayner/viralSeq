'''
Created on Jan 5, 2021

@author: simonray
'''

from pypesteps import abstractStep

import os
import logging

from Bio import SeqIO
import subprocess
import shlex
import glob


logger = logging.getLogger(__name__)
INDENT = 6


class StepSliceAndSampleBAM(abstractStep.AbstractStep):
    '''
    classdocs
    This slices and samples the specified list of BAM files
    the BAM file can be: 
    sliced by specifying a start and stop position
    (using -b/--begin & -e/--end parameters)
    
    sampled
    
    output is written to an output file in BED format
    
    To do: add parameters to set x axis plot range in GC plot
    '''
    
    CLASSID             = "StepSliceAndSampleBAM"
    OUTPUTFOLDER        = "subbedbams"
    REFFASTASHORT       = "-r"
    REFFASTALONG        = "--ref_fasta"
    SLICEBEGINSHORT     = "-g"
    SLICEBEGINLONG      = "--begin"
    SLICEENDSHORT       = "-e"
    SLICEENDLONG        = "--end"
    SAMPLEMINSHORT      = "-n"
    SAMPLEMINLONG       = "--sample_min"
    SAMPLEMAXSHORT      = "-x"
    SAMPLEMAXLONG       = "--sample_max"
    SAMPLESTEPSHORT     = "-p"
    SAMPLESTEPLONG      = "--sample_step"
    SAMPLETYPESHORT     = "-t"
    SAMPLETYPELONG      = "--sample_type"
    SAMPLEPERCENT       = "bypercent"
    SAMPLEREADS         = "byreads"
    SOFTWARELOCSHORT    = "-p"
    SOFTWARELOCLONG     = "--path_to_software"   
    BAMFILEFOLDERSHORT  = "-b"
    BAMFILEFOLDERLONG   = "--bam_file_folder"     
    
    PLOTHEIGHT          = 5
    PLOTWIDTH           = 10
    PLOTUNITS           = 'in'
    PLOTDPI             = 1000
    XVAR                = "pos"
    YVAR                = "gcpercent"
    



    def __init__(self, begin = 0, end=0, minS=0, maxS=0, stepS=0, typeS = SAMPLEPERCENT, softwarePath= "samtools", refFastA="", bamFileFolder=""):
        '''
        Constructor
        '''
        self.softwarePath = softwarePath
        self.refFastA = refFastA
        self.beginSlice = begin
        self.endSlice = end
        self.sampleMin = minS
        self.sampleMax = maxS 
        self.sampleStep = stepS
        self.sampleType = typeS
        
        self.sliceBAM = False
        self.sample = False
        
        self.bamFileFolder = bamFileFolder
        
        

    def shortDescription(self):
        print('sliceBAM / sample a list of BAM files')


    def longDescription(self):
        print('sliceBAM / sample a list of BAM files')
        print('')
        print('  sliceBAM begin position: -b / --begin')
        print('    sliceBAM end position: -e / --end')
        print('          sampling min: -n / --sampling_min')
        print('          sampling max: -x / --sampling_max')
        print('         sampling step: -t / --sampling_step')
        print('         sampling type: -p / --sampling_specs <percent|total>')
        print('')      
        print('Where sampling type specifies whether the sampling is ')
        print('in terms of total reads or percentage of reads')
        print('')
        print('For example')
        print('  ---sampling_min 10, -sampling_max 100, --sampling_step 10, --sampling_type percent ')
        print('  specifies sampling from 10% to 100% in steps of 10%')
        print('')
        print('  ---sampling_min 1000, -sampling_max 10000, --sampling_step 1000, --sampling_type total ')
        print('  specifies sampling from 1000 to 10000 read in steps of 1000')
        print('')
        print('')
        
        

    def execute(self):
        '''
        contains the main operations for the step
        '''
        logger.info(INDENT*'-' + "executing step")

        # load the genome to get the number of nucleotides
        recordCount = 0
        genomeID = ""
        for record in SeqIO.parse(self.refFastA, "fasta"):
            genomeSeq = record.seq
            genomeID  = record.id
            genomeLen = len(genomeSeq)
            recordCount += 1
            
        if recordCount > 1:
            logging.warn(INDENT*'-' + "----fasta file contains more than one record. Using the last loaded record (#" + recordCount + ")")
            logging.warn(INDENT*'-' + "----Using the last loaded record <" + genomeID + "> which is <" + genomeLen + "> nt" )

        bamFileFolder = os.path.join(self.projectRoot, self.inFolder)
        resultFolder = os.path.join(self.projectRoot, self.outFolder)
        
        logging.info(INDENT*'-' + "--results will be written to output folder <" + resultFolder + ">")
        if not os.path.exists(resultFolder):
            logging.info(INDENT*'-' + "----folder doesn't exist, creating")
            os.makedirs(resultFolder)        
            
        cmds = []
        
        for inputFile in self.inputFiles:
            basename = os.path.splitext(os.path.basename(inputFile))[0]   
            bamFile = os.path.join(bamFileFolder, inputFile)
            
            sampleSize = self.sampleMin
            while(sampleSize < self.sampleMax):
                   
                # for each BAM file:
    
                #   1. sample the BAM file
                #   2. index the output
                #   3. slice the BAM file
                #   4. index the output 
                #   5. clean up temporary files
                
                #   1. sample BAM file at sampleSize % and pipe the output for sorting
                #      samtools view -s 0.15 -b test.bam|samtools sort -o test__sp_p15_so.bam
                sampledBasename = basename + "__sp_" + str(sampleSize) + "_so"
                sampledBamFile = os.path.join(resultFolder, sampledBasename + ".bam")
                cmd1 = self.softwarePath + ' view -s ' + str(float(sampleSize)/100.0) + " -b " + bamFile\
                 + " | " + self.softwarePath + ' sort -o ' + sampledBamFile
                logging.debug(INDENT*'-' + "--SAMTools sample/sort command is <"+ cmd1 + ">")

                
                #   2. index the sorted file
                #      samtools index test__sp_p15_so.bam
                cmd2 = self.softwarePath + ' index ' + sampledBamFile 
                logging.debug(INDENT*'-' + "--SAMTools index command is <"+ cmd2 + ">")

                
                #   3. sample sliced the file, pipe the output for sorting.
                #      samtools view -hb test__sample_p15_so.bam "NC_045512.2:2-100" > test__sp_p15_so__sl_2-100.bam 
                slicedBasename = sampledBasename + "__sl_" + str(sampleSize) + "_sorted"
                slicedBamFile = os.path.join(resultFolder, slicedBasename + ".bam")
                sliceString = ' "' + genomeID + ":" + str(self.beginSlice) + "-" + str(self.endSlice)+ '" '
                cmd3 = self.softwarePath + ' view -hb ' + sampledBamFile + sliceString + " > " + slicedBamFile
                logging.debug(INDENT*'-' + "--SAMTools slice command is <"+ cmd3 + ">")

                 
                #   4. index the output.
                #      samtools index test__sp_p15_so__sl_2-100.bam
                cmd4 = self.softwarePath + ' index ' + slicedBamFile 
                logging.debug(INDENT*'-' + "--SAMTools command is <"+ cmd4 + ">")
                    
                cmds = cmds + [cmd1, cmd2, cmd3, cmd4]
                sampleSize += self.sampleStep
                
                
        shellFile = os.path.join(self.projectRoot, self.outFolder, self.projectID + "_" + "BAM_sampling" + ".sh")
        if len(cmds) > 0:
            with open(shellFile, 'w') as shfile:
                for cmd in cmds:
                    print(cmd, file=shfile)
        
        logger.info(INDENT*'-' + "done")
                

    def checkInputData(self):
        '''
        build file paths and check all input resources exist
        '''
        
        ## check reference FastA file
        if os.path.exists(self.refFastA) == False:
            raise RuntimeError (INDENT*"-" + "reference FastA file <" + self.refFastA + "> not found")
            logging.error(INDENT*"-" + "reference FastA file <" + self.refFastA + "> not found")
        else:
            logging.info(INDENT*'-' + "found reference FastA file <" + self.refFastA + ">")

        ## check begin and end values for slicing
        # it isn't necessary to specify both. If begin is missing it us set to 1, 
        # if end is missing set to the reference seq length
        # (SAM uses 1 base location)       
        if(self.beginSlice > 0):
            logging.info(INDENT*"-" + "Slice begins at <" + str(self.beginSlice) + "> nt")
        else:
            logging.info(INDENT*"-" + "No sliceBAM begin position specified")            

        if(self.endSlice > 0):
            logging.info(INDENT*"-" + "Slice ends at <" + str(self.endSlice) + "> nt")
        else:
            logging.info(INDENT*"-" + "No sliceBAM end position specified")            
            
        if self.beginSlice > 0 or self.endSlice > 0:
            self.sliceBAM = True
            logging.info(INDENT*"-" + "BAM files will be sliced")
        
        # check sampleType, min, max and step values
        
        # sampleType defaults to bypercent
        if(self.sampleType==self.SAMPLEPERCENT):
            logger.info(INDENT*"-" + "sampling type set to " + self.SAMPLEPERCENT)
        elif(self.sampleType==self.SAMPLEREADS):
            logger.info(INDENT*"-" + "sampling type set to " + self.SAMPLEREADS)
        else:
            logger.error(INDENT*"-" + "unrecognised sampling type. Options are < "\
                         + self.SAMPLEPERCENT + "|" + self.SAMPLEREADS + ">")
            raise Exception(INDENT*"-" + "unrecognised sampling type. Options are < "\
                         + self.SAMPLEPERCENT + "|" + self.SAMPLEREADS + ">")
            
        # max, min, step
        # perform some basic checks on value to ensure they make sense
        # less than zero
        if(self.sampleMin < 0):
            logger.error(INDENT*"-" + "sampling min < 0 (" + str(self.sampleMin + ")"))            
            raise Exception(INDENT*"-" + "sampling min < 0 (" + str(self.sampleMin + ")"))             
        
        if(self.sampleMax < 0):
            logger.error(INDENT*"-" + "sampling max < 0 (" + str(self.sampleMax + ")"))            
            raise Exception(INDENT*"-" + "sampling max < 0 (" + str(self.sampleMax + ")"))             
        
        if(self.sampleStep < 0):
            logger.error(INDENT*"-" + "sampling step < 0 (" + str(self.sampleStep + ")"))            
            raise Exception(INDENT*"-" + "sampling step < 0 (" + str(self.sampleStep + ")"))             
                
        # max < min, etc
        if(self.sampleMax < self.sampleMin):
            logger.error(INDENT*"-" + "sampling max < sampling min ("\
                         + str(self.sampleMax) +  "<" + str(self.sampleMin) + ")")            
            raise Exception(INDENT*"-" + "sampling max < sampling min ("\
                         + str(self.sampleMax) +  "<" + str(self.sampleMin) + ")")            
        if(self.sampleStep > (self.sampleMax - self.sampleMin)):
            logger.warn(INDENT*"-" + "sampling step > (sampling max - sampling max)"\
                         + self.sampleStep + "> (" + str(self.sampleMax) +  " - " + str(self.sampleMin) + ")")

        # if the user selected sampling reads by percent we can check values are less than 100
        if(self.sampleType==self.SAMPLEPERCENT):
            if(self.sampleMin > 100):
                logger.error(INDENT*"-<" + self.SAMPLEPERCENT + "> was selected but sampling min > 100 (" \
                                + str(self.sampleMin) + ")")  
                raise Exception(INDENT*"-<" + self.SAMPLEPERCENT + "> was selected but sampling min > 100 (" \
                                + str(self.sampleMin) + ")")  
                      
            if(self.sampleMax > 100):
                logger.error(INDENT*"-<" + self.SAMPLEPERCENT + "> was selected but sampling max > 100 (" \
                                + str(self.sampleMax) + ")")  
                raise Exception(INDENT*"-<" + self.SAMPLEPERCENT + "> was selected but sampling max > 100 (" \
                                + str(self.sampleMax) + ")")  
                
            if(self.sampleStep > 100):
                logger.error(INDENT*"-<" + self.SAMPLEPERCENT + "> was selected but sampling step > 100 (" \
                                + str(self.sampleStep) + ")")  
                raise Exception(INDENT*"-<" + self.SAMPLEPERCENT + "> was selected but sampling step > 100 (" \
                                + str(self.sampleStep) + ")")  

        if(self.sampleMax != 0): # probably should have a more thorough check for this.
            self.sample = True
            logging.info(INDENT*"-" + "BAM files will be sampled")
            
            
        # finally, check the specified BAM files exist
        if( os.path.exists(os.path.join(self.projectRoot, self.inFolder)) is False):
            logging.error("input folder <" + os.path.join(self.projectRoot, self.inFolder) + "> not found")
            raise Exception("input folder <" + os.path.join(self.projectRoot, self.inFolder) + "> not found")
        
        inputFolder = os.path.join(self.projectRoot, self.inFolder)
        if len(self.inputFiles) == 1 & (not self.inputFiles[0]):
            self.inputFiles = glob.glob(os.path.join(self.projectRoot, self.bamFileFolder) + os.path.sep + "*gen__trim_paired__sorted.bam")
        
        for inputFile in self.inputFiles:    
            inFile = os.path.join(os.path.join(inputFolder,inputFile))
       
            if os.path.exists(inFile) == False:
                logging.error(INDENT*"-" + "--input file <" + inFile + "> not found")
                raise RuntimeError ("input file <" + inFile + "> not found")
            else:
                logging.info(INDENT*'-' + "--found input file <" + inFile + ">")
        
        

    def parseJSON(self, stepJSON):

        '''
        check all the required parameters are present
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

            OUTPUTFOLDER        = "subbedbams"
            REFFASTASHORT       = "-r"
            REFFASTALONG        = "--ref_fasta"
            SLICEBEGINSHORT     = "-b"
            SLICEBEGINLONG      = "--begin"
            SLICEENDSHORT       = "-e"
            SLICEENDLONG        = "--end"
            SAMPLEMINSHORT      = "-n"
            SAMPLEMINLONG       = "--sample_min"
            SAMPLEMAXSHORT      = "-x"
            SAMPLEMAXLONG       = "--sample_max"
            SAMPLESTEPSHORT     = "-p"
            SAMPLESTEPLONG      = "--sample_step"
            SAMPLETYPESHORT     = "-t"
            SAMPLETYPELONG      = "--sample_type"
    
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
                     
            elif self.SOFTWARELOCSHORT in param or self.SOFTWARELOCLONG in param:
                if self.SOFTWARELOCLONG in param:
                    self.softwarePath = param.split(self.SOFTWARELOCLONG)[1].strip()
                else:
                    self.softwarePath = param.split(self.SOFTWARELOCSHORT)[1].strip()
                logging.info(INDENT*'-' + "samtools software location set to <" + str(self.softwarePath) + ">")

            elif self.SLICEBEGINSHORT in param or self.SLICEBEGINLONG in param:
                if self.SLICEBEGINLONG in param:
                    self.beginSlice = int(param.split(self.SLICEBEGINLONG)[1].strip())
                else:
                    self.beginSlice = int(param.split(self.SLICEBEGINSHORT)[1].strip())
                logging.info(INDENT*'-' + "sliceBAM begin set to <" + str(self.beginSlice) + ">")
            
            elif self.SLICEENDSHORT in param or self.SLICEENDLONG in param:
                if self.SLICEENDLONG in param:
                    self.endSlice = int(param.split(self.SLICEENDLONG)[1].strip())
                else:
                    self.endSlice = int(param.split(self.SLICEENDSHORT)[1].strip())
                logging.info(INDENT*'-' + "sliceBAM end set to <" + str(self.endSlice) + ">")
            
            elif self.SAMPLEMINSHORT in param or self.SAMPLEMINLONG in param:
                if self.SAMPLEMINLONG in param:
                    self.sampleMin = int(param.split(self.SAMPLEMINLONG)[1].strip())
                else:
                    self.sampleMin = int(param.split(self.SAMPLEMINSHORT)[1].strip())
                logging.info(INDENT*'-' + "sample min set to <" + str(self.sampleMin) + ">")
            
            elif self.SAMPLEMAXSHORT in param or self.SAMPLEMAXLONG in param:
                if self.SAMPLEMAXLONG in param:
                    self.sampleMax = int(param.split(self.SAMPLEMAXLONG)[1].strip())
                else:
                    self.sampleMax = int(param.split(self.SAMPLEMAXSHORT)[1].strip())
                logging.info(INDENT*'-' + "sample max set to <" + str(self.sampleMax) + ">")
            
            elif self.SAMPLESTEPSHORT in param or self.SAMPLESTEPLONG in param:
                if self.SAMPLESTEPLONG in param:
                    self.sampleStep = int(param.split(self.SAMPLESTEPLONG)[1].strip())
                else:
                    self.sampleStep = int(param.split(self.SAMPLESTEPSHORT)[1].strip())
                logging.info(INDENT*'-' + "sample step set to <" + str(self.sampleStep) + ">")
            
            elif self.SAMPLETYPESHORT in param or self.SAMPLETYPELONG in param:
                if self.SAMPLETYPELONG in param:
                    self.sampleType = param.split(self.SAMPLETYPELONG)[1].strip()
                else:
                    self.sampleType = param.split(self.SAMPLETYPESHORT)[1].strip()
                logging.info(INDENT*'-' + "sample type set to <" + self.sampleType + ">")

            elif self.BAMFILEFOLDERSHORT in param or self.BAMFILEFOLDERLONG in param:
                if self.BAMFILEFOLDERLONG in param:
                    self.bamFileFolder = param.split(self.BAMFILEFOLDERLONG)[1].strip()
                else:
                    self.bamFileFolder = param.split(self.BAMFILEFOLDERLONG)[1].strip()
                logging.info(INDENT*'-' + "bamFileFolder set to <" + str(self.bamFileFolder) + ">")            
            
        if self.bamFileFolder == "":
            logging.error("you need to specify a folder containing the BAM files")
            raise Exception("you need to specify a folder containing the BAM files")        
           
        
        