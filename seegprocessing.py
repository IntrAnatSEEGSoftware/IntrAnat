# -*- coding: utf-8 -*-
#
# Functions to process SEEG files, and dictionary of available processing methods for each manip type
#
# (c) Inserm U836 2012-2014 - Manik Bhattacharjee
#
# License GNU GPL v3
#
# Define all functions to process SEEG data like this:
#----------------------------------------------
#  def extractBetaFunction(seegfilepath):
#     try:
#       doSomething(seegfilepath)
#       return True
#     else:
#       return False
#  extractBeta = ('Extract Beta', 'Extracts the beta band from the seeg signal and saves it in the DB', extractBetaFunction)
#---------------------------------------------------
# Then you can add the processing method to 'processors' and add the name of each 'manip' it can be applied to


from externalprocesses import *

# ---------------------   Specific Functions ----------------------------
def localizeFunction(seegfilepath):
	try:
	  matlabCallNB("localize('%s')"%seegfilepath)
	  return True
	except:
		return False

localize = ('Localizer processing', 'Launches JP Lachaux\'s localize matlab process', localizeFunction)

# ---------------------   Linking functions to manips ----------------------------
processors = { localize:['VISU', 'MVEB', 'MLAH'], }
# ----------------------  Dependencies ------------------------------
# if "advancedAnalysis" needs "basicAnalysis", add advancedAnalysis:[basicAnalysis,] to dependencies
dependencies = {localize:[],}
# ------------------------ Generic functions --------------------------
def getProcessingMethods(currentManip = None):
	""" Returns all processing methods availables for the manip name provided. If no manip is provided, returns a dictionnary with all processors"""
	allManips = set([m for manips in processors.itervalues() for m in manips])
	allProcessorsByManip = dict([(m,[p for p in processors if m in processors[p]]) for m in allManips])
	if currentManip is not None:
		return allProcessorsByManip[currentManip]
	else:
		return allProcessorsByManip

def getManipNameFromDirName(dirname):
	return dirname.split('_')[0]
	