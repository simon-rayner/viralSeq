from pypesteps import abstractStep

'''
Created on Dec 18, 2020

@author:     simon rayner
@contact:    simon.rayner@medisin.uio.no
'''
import os

import logging

logger = logging.getLogger(__name__)
INDENT = 6


class StepSNVgenerateShorahCmds(abstractStep.AbstractStep):
    '''
    classdocs
    This generates the shell commands needed to run the `shorah` software package
    for estimating SNVs from a BAM file. Shorah is well suited to SNV calling in viruses 
    where you can have variable rates of evolution and read coverage due sequencing efficiency
    As it can take several hours to execute a single script, script execution doesn't take place
    inside the step.
    
    The user needs to specify the reference genome (in fasta format) that was used to align the reads 
    (-r/--ref_fasta)
    
    optional parameters
    The user may specify:
        to split the scripts into subscripts (to execute the SNV calling in parallel) 
        (g-/no_of_groups)
        
        the location of the shorah software (some Python installations have trouble locating installs)
        (-s/--software_location)
    
     
    '''
    CLASSID             = "StepSNVgenerateShorahCmds"
    OUTPUTFOLDER        = "bamsnvanalysis"
    REFFASTASHORT       = "-r"
    REFFASTALONG        = "--ref_fasta"
    NOOFGRPSSHORT       = "-g"
    NOOFGRPSLONG        = "--no_of_groups"
    SOFTWARELOCSHORT    = "-p"
    SOFTWARELOCLONG     = "--path_to_software"

    
    # the following constants have no meaning in this step
    PLOTHEIGHT      = 3
    PLOTWIDTH       = 10
    PLOTUNITS       = 'in'
    PLOTDPI         = 1000
    XVAR            = "pos"
    YVAR            = "readcoverage"
    

    def __init__(self, refFastA="", noOfGroups=1, softwarePath="shorah"):
        '''
        Constructor
        '''
        self.refFastA = refFastA
        self.noOfGroups = noOfGroups
        self.softwarePath = softwarePath

        
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
                raise RuntimeError ("input file <" + os.path.join(bamFileFolder, inputFile) + "> not found")
                logging.error("input file <" + os.path.exists(os.path.join(bamFileFolder, inputFile)) + "> not found")
            else:
                logging.info(INDENT*'-' + "found input file <" + os.path.join(bamFileFolder, inputFile) + ">")
        
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
        bamFileFolder = os.path.join(self.projectRoot, self.inFolder)
        resultFolder = os.path.join(self.projectRoot, self.outFolder)
        
        logging.info(INDENT*'-' + "--results will be written to output folder <" + resultFolder + ">")
        if not os.path.exists(resultFolder):
            logging.info(INDENT*'-' + "----folder doesn't exist, creating")
            os.makedirs(resultFolder)        
        
        cmds = []
        bamFileCount = 1
        groupCount = 0
        groupSize = int(len(self.inputFiles)/self.noOfGroups)
        
        for inputFile in self.inputFiles:
            
        # for each BAM file
            #    1. get basename and create subfolder with this name within the output folder
            #    2. add command to cd into this subfolder
            #    3. add command to execute shorah with absolute filepaths
            basename = os.path.splitext(os.path.basename(inputFile))[0]
            runFolder = os.path.join(resultFolder, basename)
            cmd1 = "mkdir " + runFolder
            logger.debug(INDENT*'-' + "-- cmd 1 is : " + cmd1)
            cmd2 = "cd " + runFolder
            logger.debug(INDENT*'-' + "-- cmd 2 is : " + cmd2)
            cmd3 = self.softwarePath + " shotgun -b " + os.path.join(self.projectRoot, self.inFolder, inputFile) + " -f " + self.refFastA
            logger.debug(INDENT*'-' + "-- cmd 3 is : " + cmd3)
            
            cmds = cmds + [cmd1, cmd2, cmd3]
            
            bamFileCount+=1
            if(bamFileCount > groupSize and groupCount < self.noOfGroups):
                # write command
                shellFile = os.path.join(self.projectRoot, self.outFolder, self.projectID + "_" + str(groupCount) + ".sh")
                with open(shellFile, 'w') as shfile:
                    for cmd in cmds:
                        print(cmd, file=shfile)
                
                bamFileCount = 1
                groupCount += 1
                cmds = []
            
        shellFile = os.path.join(self.projectRoot, self.outFolder, self.projectID + "_" + str(groupCount) + ".sh")
        if len(cmds) > 0:
            with open(shellFile, 'w') as shfile:
                for cmd in cmds:
                    print(cmd, file=shfile)
        
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
            
            elif self.NOOFGRPSLONG in param or self.NOOFGRPSSHORT in param:
                if self.NOOFGRPSLONG in param:
                    self.noOfGroups = int(param.split(self.NOOFGRPSLONG)[1].strip())
                else:
                    self.noOfGroups = int(param.split(self.NOOFGRPSSHORT)[1].strip())
                logging.info(INDENT*'-' + "no of groups set to <" + str(self.noOfGroups) + ">")
                
            elif self.SOFTWARELOCSHORT in param or self.SOFTWARELOCLONG in param:
                if self.SOFTWARELOCLONG in param:
                    self.softwarePath = param.split(self.SOFTWARELOCLONG)[1].strip()
                else:
                    self.softwarePath = param.split(self.SOFTWARELOCSHORT)[1].strip()
                logging.info(INDENT*'-' + "shorah software path set to <" + self.softwarePath + ">")
            
        
           

        

        
        