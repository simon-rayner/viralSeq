'''
Created on Dec 18, 2020

@author: simonray
'''
from abc import ABC, abstractmethod

class AbstractStep(ABC):
    '''
    this is an Abstract class inherited by all other step classes.
    '''
    PARAMID         = "parameters"
    INFOLDERID      = "inFolder"
    INFILESID       = "inFiles"
    OUTFOLDERID     = "outFolder"

    def __init__(self, params):
        '''
        Constructor
        '''
        self.projectFile = ""
        self.projectID = ""
        self.projectRoot = ""
        self.jsonString = ""
        self.inFolder = ""
        self.inputFiles = []
        self.outFolder = ""
        self.outputFiles = []
        self.paramString = ""        
        self.md5string = ""
        
    @abstractmethod
    def shortDescription(self):
        pass
    
    @abstractmethod
    def longDescription(self):
        pass
    
    @abstractmethod
    def parseParameterString(self):
        '''
        parse out the JSON string specified for this step
        For example:
        
            {"step": "stepname",  "parameters": "", "inFolder": "", "inFiles": ["", ""], "outFolder": ""}, 
    
        Where 
          "step":  is the step name
          "parameters": is a string of parameters for the step
          "inFolder":  is where the input files for the step are located
          "inFiles" :  is a list of input to be processed by the step
          "outFolder": is an optional parameter that specifies where the output files will be saved
                       If this isn't specified, the default folder name specified for this step will be used.
        '''
        pass
      
    
    @abstractmethod
    def parseJSON(self):
        pass
    
    @abstractmethod
    def checkInputData(self):
        pass    
    
    @abstractmethod
    def execute(self):
        pass    