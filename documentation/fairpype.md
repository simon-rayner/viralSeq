# Viral NGS Analysis

## Python pipeline

### Introduction
This should be written to make it easy to add new analysis and perform automated analysis of datasets.

Thus, it similar to the `java` pipeline but written in `Python` to make it easier to work on

It uses the same approach of specifying analysis steps in `JSON` format, and each analysis implemented as a separate step specified in a `Python` module

### Logging
We use the `logger` package to write output. There should be no print statements anywhere in the code.

### Exceptions
logging.error + raise Exception 

### MD5 String
```
    md5String = hashlib.md5(b"CBGAMGOUS").hexdigest()
```

### Filenames
you can specify just a filename, in which case the step will use the default foldername, or you can specify a folder as part of the filename (e.g. `folder/filename`) to 


### Analysis Steps
the steps inherit the abstract class `AbstractStep` to ensure a basic structure.

So far, this requires that the following methods are overridden

```
    @abstractmethod
    def shortDescription(self):
        pass
    
    
    @abstractmethod
    def longDescription(self):
        pass
    
    @abstractmethod
    def checkInputData(self):
        pass    

    @abstractmethod
    def execute(self):
        pass    
      
``` 
The first two attempt to provide the user with basic information about the function of the step, and what data is input/output

`checkInputData` performs some basic check of input data (such as making sure the input files exist) so that we can throw a more informative exception if something goes wrong

`execute` is where the step carries out its primary function. The idea behind this is that it makes it easier to understand the code if it needs to be modified.

The steps are placed in the `fairpype` package

### Naming convention
The step module has to begin with the string `step`. This is to allow the code to find all available steps. This allows us to do things such as presenting the user with a list of all the steps and a description of function.

### Use of CONSTANTS in a class
e.g. `self.SNVFILEEND`

## Creating a step
`class StepFactory:`

need to add a case to `StepFactory:get_step()`

```
    def get_step(self, classID):
        if classID == StepGenomeGCCoverage.CLASSID:
            return StepGenomeGCCoverage()
        if classID == StepBAMReadCoverage.CLASSID:
            return StepBAMReadCoverage()
       
        logging.error("didn't recognise the specified step ID <" + classID + ">")
        raise ValueError(format)
```

For example, when creating a new step called `StepSNVgenerateShorahCmds` the following needs to be added.

First of all, we need to tell the `StepFactory` class that we have a new step, and where to find it. So at, the top of the file we add the following line

```
from pypesteps.stepSNVGenerateShorahCmds import StepSNVgenerateShorahCmds
```
then, in the `get_step` method, we have to add a check for this new class

```p
    def get_step(self, classID):
        if classID == StepGenomeGCCoverage.CLASSID:
            return StepGenomeGCCoverage()
        if classID == StepBAMReadCoverage.CLASSID:
            return StepBAMReadCoverage()
  ==>   if classID == StepSNVgenerateShorahCmds.CLASSID:
  ==>       return StepSNVgenerateShorahCmds()            
       
        logging.error("didn't recognise the specified step ID <" + classID + ">")
        raise ValueError(format)
```

## calling a step in the pipeline

The main processing is carried out in the `VirusProject.buildStepSet()`

The four methods that need to be implemented are `parseJSON`,`parseParameterString`,`checkInputData` & `execute` 

```
    for step in self.projectSteps:
        thisStep = vStepFactory.get_step(step['step'])
        logging.info(INDENT*"-" + "--found step [" + thisStep.CLASSID + "]")
        thisStep.projectRoot = self.projectRoot
        thisStep.projectID = self.projectID
        thisStep.md5string = self.md5string
        thisStep.jsonString = step
  ==>   thisStep.parseJSON(step)
  ==>   thisStep.parseParameterString()
  ==>   thisStep.checkInputData()
  ==>   thisStep.execute()
        self.stepsToExecute.append(thisStep)
        logging.debug(self.projectSteps)
```

### parseJSON()
this shouldn't require any modifications, unless you want to add custom parameters to the JSON, which is probably a bad idea

(I should probably try to derive a class from the `AbstractStep` that implements this method to save the effort of implementing it for every new step)

### parseParameterString()
This is where you parse out the parameters specified in the `parameters:` field in the `JSON` entry. This is how you customise the step to allow the user to specify input data that is specific to the step. 

For example, in the `StepSNVgenerateShorahCmds` step, I am generating a shell script file containing the commands needed to run the `shorah` software package for SNV calling for a list of BAM files. To do this, I need the user to specify the location of the reference FastA file (using the `-r/--ref_fasta` argument) that was used in the alignment. 

Additionally, the user can specify two optional parameters. 

`-g/--no_of_groups` to tell the program to split the shell script into several smaller scripts (because it can take several hours to perform SNV calling on a single file, so it will be faster to run in parallel on multiple processors)

`-s/--software_path` to tell the program the absolute path to the `shorah` executable (on some `Python` installs there seems to be a problem locating the executable, even when located in a publicly accessible location such as `/usr/local/bin`)

`parseParameterString()` needs to check the specified values of these parameters. For example, for the `-r/--ref_fasta` parameter, does the FastA file exist at the specified location?

For other steps, you may be asking the user to specify numerical parameters. For example, when calculating a sliding window to calculate read coverage, we need a `window_size` and `step_size`. In this case, it is necessary to check that *both* these values are specified and they are both > 0. See `pypesteps.stepBAMReadCoverage.checkInputData()` for an example