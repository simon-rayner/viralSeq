'''
Created on Dec 20, 2020

@author: simonray
'''
import logging
logger = logging.getLogger(__name__)
from pypesteps.stepGCReadCoverage import StepGCReadCoverage
from pypesteps.stepSliceAndSampleBAM import StepSliceAndSampleBAM
from pypesteps.stepBAMReadCoverage import StepBAMReadCoverage
from pypesteps.stepSNVGenerateShorahCmds import StepSNVgenerateShorahCmds
from pypesteps.stepSNVProcessShorahResults import StepSNVProcessShorahResults
from pypesteps.stepBAMGCReadCorr import StepBAMGCReadCorr
from pypesteps.stepExit import StepExit

class StepFactory:
    '''
    creates an instance of a Step, derived from the AbstractStep class
    '''


    def __init__(self):
        '''
        Constructor
        '''
        self._creators = {}
        
        
    def get_step(self, classID):
        if classID == StepSliceAndSampleBAM.CLASSID:
            return StepSliceAndSampleBAM()
        if classID == StepGCReadCoverage.CLASSID:
            return StepGCReadCoverage()
        if classID == StepBAMReadCoverage.CLASSID:
            return StepBAMReadCoverage()
        if classID == StepSNVgenerateShorahCmds.CLASSID:
            return StepSNVgenerateShorahCmds()
        if classID == StepSNVProcessShorahResults.CLASSID:
            return StepSNVProcessShorahResults()
        if classID == StepBAMGCReadCorr.CLASSID:
            return StepBAMGCReadCorr()
        if classID == StepExit.CLASSID:
            return StepExit()
        
       
        logging.error("didn't recognise the specified step ID <" + classID + ">")
        raise ValueError(format)



    def register_step(self, stepType, creator):
        '''
        at this point, I'm not even clear if i need this
        '''
        self._creators[stepType] = creator

    def listAvailableSteps(self):
        '''
        The user needs to be able to get a list of all available Steps
        Need to determine the best way to do this.
        '''

#factory = SerializerFactory()