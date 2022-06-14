from pypesteps import abstractStep

'''
Created on Dec 18, 2020

@author:     simon rayner
@contact:    simon.rayner@medisin.uio.no
'''

import logging

logger = logging.getLogger(__name__)
INDENT = 6


class StepExit(abstractStep.AbstractStep):
    '''
    classdocs
    this is an empty step. It is used to terminate the pipe
    '''
    CLASSID             = "StepExit"
    OUTPUTFOLDER        = "exitpipe"
    

    def __init__(self):
        '''
        Constructor
        '''
        pass

        
    def checkInputData(self):
        '''
        nothing to check in this step
        '''
        pass
        
        
    def execute(self):
        '''
        nothing to execute in this step
        '''
        logger.info(INDENT*'-' + "executing EXIT step")
        logger.info(INDENT*'-' + "nothing to do")

        

    def shortDescription(self):
        print('instructs the pipe to terminate')


    def longDescription(self):
        print('instructs the pipe to terminate')
        print('This step can be inserted in the JSON file to test a single step.')
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
        
        
        