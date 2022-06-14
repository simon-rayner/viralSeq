'''
Created on Dec 19, 2020

@author: simonray
'''

import json
import logging

from fairpype import fairpypeConstants
from pypesteps.stepFactory import StepFactory
from pypesteps.stepExit import StepExit


logger = logging.getLogger(__name__)
INDENT = 4


class VirusProject(object):
    '''
    Stores the information about a virus project
    '''
    global vStepFactory

    def __init__(self):
        '''
        Constructor
        '''
        self.projectFile = ""
        self.projectID = ""
        self.pipelineName = ""

        self.projectSteps = []
        self.stepsToExecute = []
        
        self.md5string = " "
        

        

        
        
    def loadProjectFile(self):
        '''
        load project parameters from JSON file
        '''
        logging.debug(INDENT*"-" + "load ProjectFile <" + self.projectFile + ">")
        with open(self.projectFile) as f:
            projectJSON = json.load(f)
            
        self.projectID = projectJSON[fairpypeConstants.FairPypeConstants.PROJECT_ID]
        self.projectRoot = projectJSON[fairpypeConstants.FairPypeConstants.PROJECT_ROOT]
        self.pipelineName = projectJSON[fairpypeConstants.FairPypeConstants.PIPELINE_NAME]
        self.projectSteps = projectJSON[fairpypeConstants.FairPypeConstants.PROJECT_STEPS]
        logging.info(INDENT*"-" + "ProjectID is    <" + self.projectID + ">")
        logging.info(INDENT*"-" + "ProjectRoot is  <" + self.projectRoot + ">")
        logging.info(INDENT*"-" + "PipelineName is <" + self.pipelineName + ">")
        
        
    def buildStepSet(self):
        '''
        process the steps to translate them to executable format
        '''
        logging.debug(INDENT*"-" + "create StepFactory")
        vStepFactory = StepFactory()
        
        # first register the available steps
        vStepFactory.register_step(StepExit.CLASSID, StepExit)
        for step in self.projectSteps:
            thisStep = vStepFactory.get_step(step['step'])
            if thisStep.CLASSID == StepExit.CLASSID:
                logger.info("Found step [" + StepExit.CLASSID + "] - Exiting")
                return
            
            logging.info(INDENT*"-" + "--found step [" + thisStep.CLASSID + "]")
            thisStep.projectRoot = self.projectRoot
            thisStep.projectID = self.projectID
            thisStep.md5string = self.md5string
            thisStep.jsonString = step
            thisStep.parseJSON(step)
            thisStep.parseParameterString()
            thisStep.checkInputData()
            thisStep.execute()
            self.stepsToExecute.append(thisStep)
            logging.debug(self.projectSteps)
            
        
    def registerSteps(self):
        
        pass
        
        