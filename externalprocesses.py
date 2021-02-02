# -*- coding: utf-8 -*-
#
# Library for running matlab commands and external shell commands in blocking and non-blocking ways.
#
# (c) Inserm U836 2012-2014 - Manik Bhattacharjee
#                 2016-2021 - Francois Tadel
# License GNU GPL v3

import subprocess, os
from string import letters
from tempfile import gettempdir
from os.path import expanduser
from brainvisa.data import neuroHierarchy
from random import choice

# Set the environment for external commands without BrainVisa interference (remove brainvisa-specific paths)
myEnv = os.environ.copy()
for var in myEnv.keys():
    if var.find('BRAINVISA_UNENV_') >= 0:
        realvar = var.split('BRAINVISA_UNENV_')[1]
        print "On remplace %s par %s"%(realvar, var)
        myEnv[realvar] = myEnv[var]
        del myEnv[var]
		
# Windows computer
if os.name == 'nt':
	isWindows = True
	isWindowsWSL = False
# Linux: Check whether the program is executing on a Windows/WSL or native Linux
else:
    isWindows = False
    with open('/proc/version', 'r') as procfile:
		# isWindowsWSL = (procfile.read().find("Microsoft") != -1)
        isWindowsWSL = False

	
def getTmpDir():
    """ Get temporary folder """
    # In Windows/WSL: Use $HOME/tmp instead of /tmp
    if isWindowsWSL:
        tmpDir = expanduser("~/tmp");
        if not os.path.exists(tmpDir):
            os.makedirs(tmpDir)
    else:
        tmpDir = gettempdir()   # Expands to /tmp
    return tmpDir


def getTmpFilePath(extension='txt', isFull=False):
    """
      Generate a random file name in a temp directory

      file extension (txt by default) can be given as argument
      Usage : filepath = getTmpFilePath('jpg')
              dict{'dir','filename','fullpath'} = getTmpFilePath('jpg', True)
    """
    tmpDir = getTmpDir()
    tmpFile = ''.join(choice(letters) for i in xrange(15))
    tmpPath = os.path.join(tmpDir, tmpFile + '.' + extension)
    if isFull:
        return {'dir':tmpDir, 'filename':tmpFile, 'fullpath':tmpPath}
    else:
        return tmpPath


def matlabRun(cmd):
    """
      Runs the provided matlab code and returns the stdoutput
      The call is blocking until command completion.
    """
    # Save matlab script
    matlabCall = saveMatlabCall(cmd)
    # Run code
    [result, errMsg] = subprocess.Popen(matlabCall['code'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, env = myEnv).communicate()
    # Delete temp file
    os.remove(matlabCall['fullpath'])
    # Print error message
    if errMsg:
        print "\n===================================\n" + \
              "Matlab execution returned an error:\n" + errMsg + \
              "\n===================================\n"
    # Return stderr
    return errMsg.splitlines()


def saveMatlabCall(cmd):
    """ Format Matlab call in a cross-platform setup (replace paths if needed) """
    cmd = str(cmd)
    # For WSL/Windows, replace home and temporary folders
    if isWindowsWSL:
        cmd = cmd.replace(expanduser("~")+"/", "L:\\")
    # Get temp file
    tmp = getTmpFilePath('m', True)
    # Write Matlab script in temp folder
    f = open(tmp['fullpath'], 'w')
    f.write("disp('===============================================================');\n")
    f.write("disp('Executing Matlab script: "+tmp['fullpath']+"');\n")
    f.write("disp(' ');\n")
    f.write("fprintf(1, '" + cmd.replace("\\", "\\\\").replace("'", "''").replace('\n','\\n').replace('%','%%') + "\\n');\n")
    f.write("disp(['===============================================================' 10]);\n\n")
    f.write(cmd + "\n\n")
    f.close()
    matlabExe = ['matlab', '-nodesktop', '-r']
    print("Calling script '"+tmp['fullpath']+ "' and with Matlab command : "+repr(matlabExe))
    scriptCall = ["cd '%s';%s"%(formatExternalPath(tmp['dir']),tmp['filename']),]
    return {'code':matlabExe+scriptCall, 'file':f, 'fullpath':tmp['fullpath']}


def runCmd(cmd):
    """ Executes a command and returns the output as an array of strings (the output lines)"""
    return subprocess.Popen(cmd, stdout=subprocess.PIPE, env = myEnv).communicate()[0].splitlines()
    
    
def formatExternalPath(fullpath):
    if isWindowsWSL:
        return str(fullpath).replace(expanduser("~"), "L:").replace("/", "\\")
    else:
        return str(fullpath)


def createItemDirs(item):
    """ Create the directories containing the provided WriteDiskItem and insert them in the BrainVisa database """
    # Copied from brainvisa.in ExecutionContext._processExecution()
    dirname = os.path.dirname( item.fullPath() )
    dir=dirname
    dirs = []
    while not os.path.exists( dir ):
        dirs.append(dir)
        dir=os.path.dirname(dir)
    if dirs:
        try:
            os.makedirs( dirname )
        except OSError, e:
            if not e.errno == errno.EEXIST:
                # filter out 'File exists' exception, if the same dir has
                # been created concurrently by another instance of BrainVisa
                # or another thread
                raise
        for d in dirs:
            dirItem=neuroHierarchy.databases.createDiskItemFromFileName(d, None)


def removeFromDB(file, db=None):
    """
    If the file is a directory, recursive call to remove all its content before removing the directory.
    Corresponding diskitem is removed from the database if it exists.
    Taken from brainvisa-4.3.0/python/brainvisa/data/qt4gui/hierarchyBrowser.py
    """
    # Delete database entry
    if db is None:
        try:
            db=neuroHierarchy.databases.database(neuroHierarchy.databases.getDiskItemFromFileName(file).get("_database"))
        except:
            pass
    if db:
        diskItem=db.getDiskItemFromFileName(file, None)
        if diskItem:
            db.removeDiskItem(diskItem)
    # Delete folder/file
    if os.path.isdir(file):
        for f in os.listdir(file):
            removeFromDB(os.path.join(file, f), db)
        os.rmdir(file)
    elif os.path.exists(file):
        os.remove(file)
    # Delete .minf
    if os.path.exists(file + '.minf'):
        os.remove(file + '.minf')

    