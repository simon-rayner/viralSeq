#!/usr/local/bin/python2.7
# encoding: utf-8
'''
fairpype.virusPipe -- 

fairpype.virusPipe is a simple framework for stepwise analysis of virus data

It defines classes_and_methods

@author:     Simon Rayner

@copyright:  2020 Oslo University Hospital. All rights reserved.

@license:    license

@contact:    simon.rayner@medisin.uio.no
@deffield    updated: Updated
'''

import sys
import os
from datetime import datetime
import hashlib
import logging

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

from fairpype import virusProject
from anaconda_navigator.utils.py3compat import getcwd



__all__ = []
__version__ = 0.1
__date__ = '2020-12-19'
__updated__ = '2020-12-19'

DEBUG = 1
TESTRUN = 0
PROFILE = 0

class CLIError(Exception):
    '''Generic exception to raise and log different fatal errors.'''
    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: %s" % msg
    def __str__(self):
        return self.msg
    def __unicode__(self):
        return self.msg
    
    
virusProject = virusProject.VirusProject()

def main(argv=None): # IGNORE:C0111
    
    if argv is None:
        argv = sys.argv

    md5String = hashlib.md5(b"CBGAMGOUS").hexdigest()
    parseArgs(argv)
    initLogger(md5String)
    virusProject.projectFile = projectFile
    virusProject.md5string = md5String
    logging.info("project file is <" + projectFile + ">")
    virusProject.loadProjectFile()
    virusProject.buildStepSet()
    

    
    
def initLogger(md5string):
    
    ''' setup log file based on project name'''
    projectBaseName = os.path.splitext(os.path.basename(projectFile))[0]
    now = datetime.now()
    dt_string = now.strftime("%Y%m%d_%H%M%S")
    logFolder = os.path.join(getcwd(), "logfiles")
    if not os.path.exists(logFolder):
        print("--log folder <" + logFolder + "> doesn't exist, creating")
        os.makedirs(logFolder)   
    logfileName = os.path.join(logFolder, projectBaseName + "__" + dt_string + "__" + md5string +".log")
    handler = logging.StreamHandler(sys.stdout)
    logging.basicConfig(level=logging.DEBUG)
    
    fileh = logging.FileHandler(logfileName, 'a')
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fileh.setFormatter(formatter)
    
    log = logging.getLogger()  # root logger
    for hdlr in log.handlers[:]:  # remove all old handlers
        log.removeHandler(hdlr)
    log.addHandler(fileh)      # set the new handler 
    log.addHandler(handler) 
    logging.info("+" + "*"*78 + "+")   
    logging.info("project log file is <" + logfileName + ">")         
    logging.info("+" + "*"*78 + "+")   
    

def parseArgs(argv):
    
    '''parse out Command line options.'''

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''%s
i
  Created by Simon Rayner on %s.
  Copyright 2020 Oslo University Hospital. All rights reserved.

  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
''' % (program_shortdesc, str(__date__))

    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument("-p", "--projectfile", dest="projectfile", action="store", help="project file in JSON format [default: %(default)s]")
        
        #parser.add_argument("-v", "--verbose", dest="verbose", action="count", help="set verbosity level [default: %(default)s]")
        #parser.add_argument("-i", "--include", dest="include", help="only include paths matching this regex pattern. Note: exclude is given preference over include. [default: %(default)s]", metavar="RE" )
        #parser.add_argument("-e", "--exclude", dest="exclude", help="exclude paths matching this regex pattern. [default: %(default)s]", metavar="RE" )
        #parser.add_argument('-V', '--version', action='version', version=program_version_message)
        #parser.add_argument(dest="paths", help="paths to folder(s) with source file(s) [default: %(default)s]", metavar="path", nargs='+')

        # Process arguments
        args = parser.parse_args()

        global projectFile 
        projectFile = args.projectfile
        #paths = args.paths
        #verbose = args.verbose
        #recurse = args.recurse
        #inpat = args.include
        #expat = args.exclude

        #if verbose > 0:
        #    print("Verbose mode on")
        print(projectFile)

    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception as e:
        print(e)
        if DEBUG or TESTRUN:
            raise(e)
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2

if __name__ == "__main__":
    if DEBUG:
        pass
        #sys.argv.append("-h")
        #sys.argv.append("-v")

    if TESTRUN:
        import doctest
        doctest.testmod()
    if PROFILE:
        import cProfile
        import pstats
        profile_filename = 'fairpype.virusPipe_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())