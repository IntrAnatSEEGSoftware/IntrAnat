# -*- coding: utf-8 -*-

# Group of functions to read, organize, and convert DICOM files
#
# (c) Inserm U836 2012-2014 - Manik Bhattacharjee
#
# License GNU GPL v3
#

import re, shutil
from PyQt4 import QtCore
from externalprocesses import *


#  Matlab code : convert DICOM files to a NIFTI file with SPM
spm_convert_dicom = "spm_get_defaults; cd('%s'); hdr = spm_dicom_headers(%s);name = '%s';a=spm_dicom_convert(hdr, 'all', 'flat', 'nii');if length(a.files)<1 || isempty(a.files{1}); quit; end;for i=1:length(a.files) if i==1 movefile(a.files{i}, [name '.nii']);else movefile(a.files{i}, [name '_' num2str(i)  '.nii']);end;end;quit; "
	

# Construction de la liste des codes:
# - dcmdump d'un fichier DICOM
# - suppression des lignes sans intéret
# -recherche et remplace avec kwrite mode regexp
#  (\([\da-f]+,[\da-f]+\)).*#[^a-zA-Z]+(\w+) remplacé par '\1':'\2',\\
	
# The DICOM codes that will be extracted from each image.
keepCodes = {\
#  '(0008,0012)':'InstanceCreationDate',\
#	'(0008,0013)':'InstanceCreationTime',\
#	'(0008,0014)':'InstanceCreatorUID',\
#	'(0008,0016)':'SOPClassUID',\
#	'(0008,0018)':'SOPInstanceUID',\
	'(0008,0020)':'StudyDate',\
	'(0008,0021)':'SeriesDate',\
	'(0008,0022)':'AcquisitionDate',\
#	'(0008,0023)':'ContentDate',\
	'(0008,0030)':'StudyTime',\
	'(0008,0031)':'SeriesTime',\
	'(0008,0032)':'AcquisitionTime',\
#	'(0008,0033)':'ContentTime',\
	'(0008,0050)':'AccessionNumber',\
#	'(0008,0060)':'Modality',\
#	'(0008,0070)':'Manufacturer',\
#	'(0008,0080)':'InstitutionName',\
#	'(0008,1010)':'StationName',\
	'(0008,1030)':'StudyDescription',\
	'(0008,103e)':'SeriesDescription',\
#	'(0008,1040)':'InstitutionalDepartmentName',\
	'(0010,0010)':'PatientsName',\
#	'(0010,0020)':'PatientID',\
	'(0010,0030)':'PatientsBirthDate',\
	'(0010,0040)':'PatientsSex',\
#	'(0018,0015)':'BodyPartExamined',\
#	'(0018,0020)':'ScanningSequence',\
	'(0018,0023)':'MRAcquisitionType',\
	'(0018,0050)':'SliceThickness',\
#	'(0018,0080)':'RepetitionTime',\
#	'(0018,0081)':'EchoTime',\
#	'(0018,0083)':'NumberOfAverages',\
#	'(0018,0087)':'MagneticFieldStrength',\
	'(0018,0088)':'SpacingBetweenSlices',\
	'(0018,1030)':'ProtocolName',\
#	'(0018,1250)':'ReceiveCoilName',\
#	'(0018,1310)':'AcquisitionMatrix',\
	'(0018,9073)':'AcquisitionDuration',\
	'(0020,000d)':'StudyInstanceUID',\
	'(0020,000e)':'SeriesInstanceUID',\
	'(0020,0010)':'StudyID',\
#	'(0020,0011)':'SeriesNumber',\
#	'(0020,0012)':'AcquisitionNumber',\
#	'(0020,0013)':'InstanceNumber',\
#	'(0020,0032)':'ImagePositionPatient',\
#	'(0020,0037)':'ImageOrientationPatient',\
#	'(0020,0052)':'FrameOfReferenceUID',\
#	'(0028,0010)':'Rows',\
#	'(0028,0011)':'Columns',\
#	'(0028,0030)':'PixelSpacing',\
#	'(0028,0100)':'BitsAllocated',\
#	'(0028,0101)':'BitsStored',\
#	'(0032,1033)':'RequestingService',\
#	'(0032,1060)':'RequestedProcedureDescription',\
#	'(0040,0241)':'PerformedStationAETitle',\
#	'(0040,0244)':'PerformedProcedureStepStartDate',\
#	'(0040,0245)':'PerformedProcedureStepStartTime',\
#	'(0040,0250)':'PerformedProcedureStepEndDate',\
#	'(0040,0251)':'PerformedProcedureStepEndTime',\
#	'(0040,0253)':'PerformedProcedureStepID',\
#	'(0040,0254)':'PerformedProcedureStepDescription',\
	'(2001,100a)':'SliceNumberMR',\
	'(2001,100b)':'SliceOrientation'}

# Regular expression to recognize a dicom tag from a text line
dcmCodeRE = re.compile(r'^\([\da-f]+,[\da-f]+\)')
# Regular expression to get the content of a DICOM tag
dcmValueRE = re.compile(r'\([\da-f]+,[\da-f]+\)\s+\w+\s+\[?([^\]#]*)\]?')


def getDcmValue(line):
	"""Gets the value of a DICOM tag from a text line
     If there is no such value, this will throw an exception
	"""
	m = dcmValueRE.search(line)
	if m is None:
		print 'No value in DCM line **'+line
	return m.group(1)


def organizeDicomFiles(files, thread = False):
	""" Creates a dictionary of all files provided with DICOM tag data extracted by dcmdump (of DCMTK)
	    :param thread : If true, a Qt signal 'dicomFilesOutput(PyQt_PyObject)' will be emitted on completion with the result dictionary
	    Type of result : dcm={'DCM0001':{'StudyDate':'20140101','SeriesDate':'20140101',...}, 'DCM0002':{'StudyDate':'20140101','SeriesDate':'20140101',...},... }
	"""
	dcm = {}
	goodkeys = []
	for k in keepCodes.keys():
		goodkeys.extend(['+P',k[1:-1]])
		
	for f in files:
		print 'Processing '+f
		dcm[f] = {}
		# Get dcmdump info for the current file f
		lines = runCmd(['dcmdump','-M']+goodkeys+[f, ])
		
		# Python 2.7 -> subprocess.check_output(['dcmdump', f])
		for line in lines:
			k = line[:11] # Get the key from the line
			if k in keepCodes:
				dcm[f][keepCodes[k]] = getDcmValue(line)
	if thread:
	  print "emitting signal"
	  QtCore.QThread.currentThread ().emit(QtCore.SIGNAL('dicomFilesOutput(PyQt_PyObject)'), dcm)
	return dcm

def patientNameDecode(p):
	""" Converts DICOM-style patient name LastName^FirstName^SecondName to LastName_FirstName_SecondName"""
	return p.replace('^', '_').rstrip('_').replace(' ','-')
	
def datetimeDecode(d,t):
  """ Decodes DICOM datetime format 20140101, 120531 to human-readable 'Date: 2014-01-01, heure: 12:05:31'"""
  if len(d) < 8 or len(t)<6:
    print "DICOM Datetime is invalid !"
    return ''
  return str('Date: '+d[:4] +'-'+d[4:6]+'-'+ d[6:] +', heure: '+ t[:2] + ':'+t[2:4]+':'+t[4:])

def dateDecode(d):
  """ Decodes DICOM date format 20140101, to human-readable '2014-01-01' """
  if len(d)<8:
    print "DICOM Datetime is invalid !"
    return ''
  return d[:4] +'-'+d[4:6]+'-'+ d[6:]
  
def getSeries(dcmFiles):
	""" From organizeDicomFiles(f) output, returns a dictionnary of series with the most important DICOM tag values and list of files"""
	series={}
	
	for d,dcm in dcmFiles.iteritems():
		dc = lambda x:dcm[x] if x in dcm else ''
		if dcm['SeriesInstanceUID'] in series:
			series[dcm['SeriesInstanceUID']]['files'].append(d)
		else:
			
			series[dcm['SeriesInstanceUID']] = {'files':[d,],\
				'SeriesDescription':dc('SeriesDescription'),\
				'SeriesDate':dc('SeriesDate'),\
				'SeriesTime':dc('SeriesTime'),\
				'StudyDescription':dc('StudyDescription'),\
				'PatientsName':str(dc('PatientsName')),\
				'PatientsBirthDate':dc('PatientsBirthDate'),\
				#'SliceThickness':dcm['SliceThickness'],\
				#'PixelSpacing':dcm['PixelSpacing'],\
				'StudyInstanceUID':str(dc('StudyInstanceUID'))
				}
	return series

def serie2filename(s):
  """ Generate the filename for the provided serie dictionary"""
  seriesName = s['SeriesDate']+'_'+s['SeriesTime']+'_'+s['SeriesDescription']
  return seriesName.replace(' ', '_')
  
def getStudies(dcmFiles):
	""" Returns a dictionnary of studies available in organizeDicomFiles(f) output"""
	studies={}
	for dcm in dcmFiles.itervalues():
		if dcm['StudyInstanceUID'] not in studies:
			studies[dcm['StudyInstanceUID']]={'name':datetimeDecode(dcm['StudyDate'],dcm['StudyTime']), 'PatientsName':dcm['PatientsName']}
	return studies

def getPatients(series):
	""" Returns a list of patients available in organizeDicomFiles(f) output"""
	patients = []
	for s in series.itervalues():
		if str(s['PatientsName']) not in patients:
			patients.append(str(s['PatientsName']))
	return patients
	
def copySeriesTo(series, path, move = False, nifti = False):
	""" Copy the series listed as a dictionary to a directory (path) and organize them in a patientName/Serie/ hierarchy
	    :param series : the series to copy/move
	    :param path the base directory that will contain the patientName/Serie/ hierarchy
	    :param move If true, files will be move. If false (default), they will be copied.
	    :param nifti : if nifti is true, the nifti files listed in the series parameter will be copied/moved. If not, the dicom files will be copied/moved.
	"""
	for s in series.itervalues():
		try:
			patientPath = os.path.join(path, patientNameDecode(s['PatientsName']))
			os.mkdir(patientPath)
		except OSError, e:
			pass
			#if e.errno != os.errno.EEXIST: # DEPEND DE LA VERSION PYTHON
			#	raise
				
		seriesName = serie2filename(s)
		
		try:
			seriesPath = os.path.join(patientPath, seriesName)
			os.mkdir(seriesPath)
		except OSError, e:
			#if e.errno != os.errno.EEXIST:
				print 'ERREUR : impossible de creer le dossier '+seriesPath+': '+e.strerror
				# On  oublie cette série ou on écrase ?
				#continue
		try:
			if nifti is False:
				files = s['files']
			else:
				if 'nifti' not in s:
					print "Cannot move unavailable Nifti files : use convertSeriesToNifti first !"
					continue
				files = s['nifti']
			
			if move == False:
				for f in files:
					shutil.copyfile(f, seriesPath)
			else:
				for f in files:
					shutil.move(f, seriesPath)
					
		except OSError, e:
			print 'ERREUR : impossible de copier/déplacer '+f+' vers '+seriesPath+': '+e.strerror
		except IOError, e:
		  print 'Erreur IO : impossible de copier/déplacer '+f+' vers '+seriesPath+': '+e.strerror 



def convertSeriesToNifti(series, path):
	"""Convert the series of dicom files to nifti files. Files are named by the serie name and in patientName directories.
	The conversion is done with matlab/SPM and the returned QThread objects must be kept until completion of the conversion."""
	threads = []
	for i, s in series.iteritems():
		try:
			patientPath = os.path.join(path, patientNameDecode(s['PatientsName']))
			os.mkdir(patientPath)
		except OSError, e:
			pass
			#if e.errno != os.errno.EEXIST: # DEPEND DE LA VERSION DE PYTHON
			#	raise

		seriesName = s['SeriesDate']+'_'+s['SeriesTime']+'_'+s['SeriesDescription']
		seriesName = seriesName.replace(' ', '_')
		niftiPath = os.path.join(patientPath, seriesName) # without the '.nii', it will be added in matlab
		series[i]['niftiPath'] = niftiPath + '.nii'
		print "Creating nifti at "+niftiPath
		# spm needs a list of filenames where ALL FILENAMES ARE THE SAME LENGTH (stupid matlab string behaviour grrr)
		# So let's make a stupid array with same-length filenames (add spaces at the end)
		matfiles = s['files']
		maxcars = max([len(s) for s in matfiles])
		matfiles = [s.ljust(maxcars) for s in matfiles]
		mtlb = matlabRunNB(spm_convert_dicom%(patientPath,  repr(matfiles).replace(',', ';'), niftiPath))
		mtlb.start()
		threads.append(mtlb)
		
	return threads
