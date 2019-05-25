# -*- coding: utf-8 -*-
#
#
# Library for running matlab commands and external shell commands in blocking and non-blocking ways.
#
# (c) Inserm U836 2012-2014 - Manik Bhattacharjee
# License GNU GPL v3
#

# TODO : redo using QRunnable instead of QThread with QThreadPool.globalInstancer().start(myRunnable)

import subprocess, traceback, os, types, tempfile, time, random, string, sys, os
from soma.qt_gui.qt_backend import QtGui, QtCore
from os.path import expanduser
from brainvisa.data import neuroHierarchy


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
# Linux: Check whether the program is executing on a Windows/WSL or native Linux system
else:
    isWindows = False
    with open('/proc/version', 'r') as procfile:
		# isWindowsWSL = (procfile.read().find("Microsoft") != -1)
        isWindowsWSL = False

# Base call for matlab
# matlabCall = ['matlab', '-nosplash', '-nodisplay', '-r']
matlabCall = ['matlab', '-nodesktop', '-r']

	
def getTmpDir():
    """ Get temporary folder """
    # In Windows/WSL: Use $HOME/tmp instead of /tmp
    if isWindowsWSL:
        tmpDir = expanduser("~/tmp");
        if not os.path.exists(tmpDir):
            os.makedirs(tmpDir)
    else:
        tmpDir = tempfile.gettempdir()   # Expands to /tmp
    return tmpDir

def getTmpFilePath(extension='txt', isFull=False):
    """
      Generate a random file name in a temp directory

      file extension (txt by default) can be given as argument
      Usage : filepath = getTmpFilePath('jpg')
              dict{'dir','filename','fullpath'} = getTmpFilePath('jpg', True)
    """
    tmpDir = getTmpDir()
    tmpFile = ''.join(random.choice(string.letters) for i in xrange(15))
    tmpPath = os.path.join(tmpDir, tmpFile + '.' + extension)
    if isFull:
        return {'dir':tmpDir, 'filename':tmpFile, 'fullpath':tmpPath}
    else:
        return tmpPath

class Executor(QtCore.QThread):
    """ This class executes a shell command in a thread.
    The Executor object must NOT be destroyed until the execution is complete.
    Parameters include :
      :param commandList (for example ['cat', '/tmp/file1', '/tmp/file2']
      :param parent a parent for the QThread (None by default)
      :param objectsToKeep : objets to keep in the Executor object until the end of execution
      :param exitFunc : a callback function that will be called at the end of the command execution

     As this is a QThread object, the object is created, then the non-blocking start() function must be called to start the thread.
     Do NOT call directly the run() function
    """
    def __init__(self, commandList, parent=None, objectsToKeep=None, exitFunc=None):
        QtCore.QThread.__init__(self,parent)
        self.commandList = commandList
        self.objectsToKeep = objectsToKeep
        self.exitFunc = exitFunc

    def run(self):
        """
           Reimplementation of the run function of the QThread.
           This SHOULD NOT BE CALLED DIRECTLY as it would run in the current thread.
           Run self.start() to start the execution in a separate thread
        """
        #pdb.set_trace()
        #self.emit(QtCore.SIGNAL("Started( QString )"),'Command started')
        print "######## Calling process in a QThread ####"+' '.join(self.commandList)+"\n\n\n"
        try:
            lines = subprocess.Popen(self.commandList, stdout=subprocess.PIPE, env = myEnv).communicate()[0].splitlines()
        except:
            print "Erreur lors de l'exécution de la commande : ", sys.exc_info()[0]
            print("Execution failed: \n" + traceback.format_exc());
            lines = ""
        if self.exitFunc:
            self.exitFunc()
        print "QThread finished"
        print "*****************LINES"+repr(lines)+"\n\n\n"
        #self.emit(QtCore.SIGNAL("Finished( QStringList )"),lines)


class PythonExecutor(QtCore.QThread):
    """ This class executes a python function with no arguments (can be created with 'lambda(x=myXvalue,y=myYvalue) : myFunc(x,y)') in a thread.
    The object must NOT be destroyed until the execution is complete.
    After creation, call the start() function to start the execution in a separate thread.
    """
    def __init__(self, func, toKeep = None, parent=None):
        """:param func : The function that will be run in a separate thread
           :param parent : a parent object for the QThread (default is None)
           :param toKeep : object or dictionnary or list of objects to keep a reference to.
        """
        QtCore.QThread.__init__(self,parent)
        self.func = func
        self.out = None
        self.toKeep = toKeep

    def output(self):
        """ Returns the output value of the function when execution is terminated"""
        return self.out

    def kept(self):
        return self.toKeep

    def run(self):
        """
           Reimplementation of the run function of the QThread.
           This SHOULD NOT BE CALLED DIRECTLY as it would run in the current thread.
           Run self.start() to start the execution in a separate thread
        """
        #self.emit(QtCore.SIGNAL("Started( QString )"),'Command started')
        self.out = self.func()
        #print "*****************OUTPUT"+repr(out)+"\n\n\n"

# def matlabIsPresent():
#     """
#         Check if the 'matlab' command is available in the path, by running a very small matlab command
#         Returns True or False
#     """
#     try:
#      	result = subprocess.Popen(matlabCall+['quit;',], stdout=subprocess.PIPE, env = myEnv).communicate()[0].splitlines()
#     except OSError:
#      	print "matlab is not present : OSError"
#       	return False
#     except:
#      	print "matlab is not present or fails !"
#       	return False
#    	return True

def matlabRun(cmd):
    """
      Runs the provided matlab code and returns the stdoutput
      The call is blocking until command completion.
      Use matlabRunNB to get non-blocking command execution
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

def matlabRunNB(cmd, callback = None):
    """ Returns a QThread object that can launch the matlab command by calling its start() function.
    This allows to run a matlab command without blocking the rest of the program.
     A "Finished( QStringList )" Qt signal is emitted when matlab exits.
     :param callback is a callback function that will be called at the end of the execution
    """
    # Save matlab script
    matlabCall = saveMatlabCall(cmd)
    # Create callback function : destroy temporary file and call the provided callback function
    if callback is None:
        cb = lambda:os.remove(matlabCall['fullpath'])
    else:
        def cb():
            os.remove(matlabCall['fullpath'])
            callback()
    return Executor(matlabCall['code'], objectsToKeep=matlabCall['file'], exitFunc = cb)


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
    print("Calling script '"+tmp['fullpath']+ "' and with Matlab command : "+repr(matlabCall))
    scriptCall = ["cd '%s';%s"%(formatExternalPath(tmp['dir']),tmp['filename']),]
    return {'code':matlabCall+scriptCall, 'file':f, 'fullpath':tmp['fullpath']}

def runCmd(cmd):
    """ Executes a command and returns the output as an array of strings (the output lines)"""
    return subprocess.Popen(cmd, stdout=subprocess.PIPE, env = myEnv).communicate()[0].splitlines()
    
def formatExternalPath(fullpath):
    if isWindowsWSL:
        return str(fullpath).replace(expanduser("~"), "L:").replace("/", "\\")
    else:
        return str(fullpath)
    

def removeFromDB(file, db=None):
    """
    If the file is a directory, recursive call to remove all its content before removing the directory.
    Corresponding diskitem is removed from the database if it exists.
    Taken from brainvisa-4.3.0/python/brainvisa/data/qt4gui/hierarchyBrowser.py
    """
    if db is None:
        try:
            db=neuroHierarchy.databases.database(neuroHierarchy.databases.getDiskItemFromFileName(file).get("_database"))
        except:
            pass
    
    if os.path.isdir(file):
        for f in os.listdir(file):
            removeFromDB(os.path.join(file, f), db)
        os.rmdir(file)
    else:
        os.remove(file)
    if db:
        diskItem=db.getDiskItemFromFileName(file, None)
        if diskItem:
            db.removeDiskItem(diskItem)
                
                
    