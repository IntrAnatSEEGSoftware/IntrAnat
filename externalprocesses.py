# -*- coding: utf-8 -*-
#
#
# Library for running matlab commands and external shell commands in blocking and non-blocking ways.
#
# (c) Inserm U836 2012-2014 - Manik Bhattacharjee
# License GNU GPL v3
#

# TODO : redo using QRunnable instead of QThread with QThreadPool.globalInstancer().start(myRunnable)

import subprocess, os, types, tempfile, time, random, string,sys
from PyQt4 import QtCore
import pdb


# Set the environment for external commands without BrainVisa interference (remove brainvisa-specific paths)
myEnv = os.environ.copy()
for var in myEnv.keys():
	if var.find('BRAINVISA_UNENV_') >= 0:
		realvar = var.split('BRAINVISA_UNENV_')[1]
		print "On remplace %s par %s"%(realvar, var)
		myEnv[realvar] = myEnv[var]
		del myEnv[var]

# Base call for matlab
matlabCall = ['matlab', '-nosplash', '-nodisplay','-r']     #essayer sans le -nojvm


def getTmpFilePath(extension='txt'):
	"""
	  Generate a random file name in a temp directory

	  file extension (txt by default) can be given as argument
	  Usage : filepath = getTmpFilePath('jpg')
	"""
	tmpdir = tempfile.gettempdir()
	tmpCmd = ''.join(random.choice(string.letters) for i in xrange(15))
	tmpfile = tmpCmd + '.' + extension
	#print "TMPdir : %s, filename : %s"%(repr(tmpdir), repr(tmpfile))
	fullpath = os.path.join(tmpdir, tmpfile)
	return fullpath

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
			#pdb.set_trace()
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

def matlabIsPresent():
  """
    Check if the 'matlab' command is available in the path, by running a very small matlab command
    Returns True or False
  """
  try:
    result = subprocess.Popen(matlabCall+['quit;',], stdout=subprocess.PIPE, env = myEnv).communicate()[0].splitlines()
  except OSError:
    print "matlab is not present : OSError"
    return False
  except:
    print "matlab is not present or fails !"
    return False
  return True

def matlabRun(cmd):
	"""
	  Runs the provided matlab code and returns the stdoutput
	  The call is blocking until command completion.
	  Use matlabRunNB to get non-blocking command execution
	"""
	if not isinstance(cmd, types.ListType):
		cmd = [cmd,]
		# Because the command may be too long for a command line, let's put it in a temp file and run the file
	tmpdir = tempfile.gettempdir()
	tmpCmd = ''.join(random.choice(string.letters) for i in xrange(15))
	tmpfile = tmpCmd + '.m'
	#print "TMPdir : %s, filename : %s"%(repr(tmpdir), repr(tmpfile))
	fullpath = os.path.join(tmpdir, tmpfile)
	# os.path.join(tempfile.tempdir,''.join(random.choice(string.letters) for i in xrange(15))+'.m')
	f = open(fullpath, 'w')
	f.write(' '.join(cmd) + "\n\n")
	f.close()
	print "Calling matlab in "+tmpdir+" with file "+tmpfile+ " and matlab command : "+repr(matlabCall)
	cmd = ["cd '%s';%s"%(tmpdir,tmpCmd),]
	result = subprocess.Popen(matlabCall+cmd, stdout=subprocess.PIPE, env = myEnv).communicate()[0].splitlines()
	#pdb.set_trace()
	os.remove(fullpath)
	return result

def matlabRunNB(cmd, callback = None):
	""" Returns a QThread object that can launch the matlab command by calling its start() function.
	This allows to run a matlab command without blocking the rest of the program.
	 A "Finished( QStringList )" Qt signal is emitted when matlab exits.
	 :param callback is a callback function that will be called at the end of the execution
	"""
	cmd = [cmd,]
	# Because the command may be too long for a command line, let's put it in a temp file and run the file
	tmpdir = tempfile.gettempdir()
	tmpCmd = ''.join(random.choice(string.letters) for i in xrange(15))
	tmpfile = tmpCmd + '.m'
	#print "TMPdir : %s, filename : %s"%(repr(tmpdir), repr(tmpfile))
	fullpath = os.path.join(tmpdir, tmpfile)
	# os.path.join(tempfile.tempdir,''.join(random.choice(string.letters) for i in xrange(15))+'.m')
	f = open(fullpath, 'w')
	f.write(' '.join(cmd) + "\n\n")
	f.close()
	print "Calling matlab in "+tmpdir+" with file "+tmpfile+ " and matlab command : "+repr(matlabCall)
	cmd = ["cd '%s';%s"%(tmpdir,tmpCmd),]
	# Create callback function : destroy temporary file and call the provided callback function
	if callback is None:
		cb = lambda:os.remove(fullpath)
	else:
		def cb():
			os.remove(fullpath)
			callback()
	#pdb.set_trace()
	return Executor(matlabCall+cmd, objectsToKeep=f, exitFunc = cb)

def runCmd(cmd):
	""" Executes a command and returns the output as an array of strings (the output lines)"""
	return subprocess.Popen(cmd, stdout=subprocess.PIPE, env = myEnv).communicate()[0].splitlines()
