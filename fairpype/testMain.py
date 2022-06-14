'''
Created on Dec 18, 2020

@author: simonray
'''
import sys
import logging
from pypesteps import stepGenomeGCCoverage


def main(argv):
    
    # Initialize logging
    handler = logging.StreamHandler(sys.stdout)
    logging.basicConfig(level=logging.DEBUG)
    
    fileh = logging.FileHandler('logfile.txt', 'a')
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fileh.setFormatter(formatter)
    
    log = logging.getLogger()  # root logger
    for hdlr in log.handlers[:]:  # remove all old handlers
        log.removeHandler(hdlr)
    log.addHandler(fileh)      # set the new handler  
    log.addHandler(handler)  
    logging.info('This will get logged')     
    
    testFA = "/Users/simonray/Dropbox/DResearch/Simmi/NGS_analysis/SARS-CoV-2_NC_045512.fa"
    
    stepGC = stepGenomeGCCoverage.StepBAMReadCoverage(windowSize=100, inputFile = testFA)
    stepGC.checkInputData()
    stepGC.execute()
    stepGC.gcPrint()
    print("done")


if __name__ == '__main__':
    main(sys.argv[1:])
