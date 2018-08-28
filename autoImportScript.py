# -*- coding: utf-8 -*-
#
# (c) INSERM U836 - Manik Bhattacharjee
# Script to import old databases of patients into new BrainVisa hierarchy

import os, json,sys
from soma.qt_gui.qt_backend import QtGui, QtCore

#Init PyQt4
app=QtGui.QApplication(sys.argv)

# Dossier
#root = '/home/lnppeeg/data/LOCA_DATABASE'
root = '/media/My Passport/epiBaseDo/LyonNuit/LOCA_DATABASE'
rootEpi = '/media/My Passport/epiBaseDo/EpiG'

# readPrefs
pr={}
try:
  prfile = open(os.path.join(os.path.expanduser('~'), '.autoImportBrainVisa'),'r')
  pr = json.load(prfile)
  prfile.close()
except:
  pass

###################### LOCA_DATABASE -> un dossier avec les IRM et schémas d'implantation (MRI_2010), un dossier avec la SEEG (2010).
# Dans chaque dossier, un dossier par sujet
dirs = os.listdir(root)
subjects = []
if 'oldsubjects' in pr:
	subjects = pr['oldsubjects']
else:
	for d in dirs:
		if os.path.isdir(os.path.join(root,d)):
			subjects.extend([f for f in os.listdir(os.path.join(root, d)) if os.path.isdir(os.path.join(root,d,f))])

	# Unique subjects, sorted
	subjects = sorted(list(set(subjects)))

	pr['oldsubjects'] = subjects


# Build a correspondance oldnames -> newNames
newNames={}
if "newNames" in pr:
	newNames = pr["newNames"]

widgets={}

def updateNewNames():
	for s in subjects:
		newNames[s] = str(widgets[s].text())
	pr["newNames"] = newNames

dialog = QtGui.QDialog()
lay = QtGui.QVBoxLayout(dialog)
#dialog.setLayout(
button=QtGui.QPushButton("Save", dialog)

for s in subjects:
	lay2 = QtGui.QHBoxLayout()
	lay2.addWidget(QtGui.QLabel(s))
	widgets[s] = QtGui.QLineEdit()
	if s in newNames:
		widgets[s].setText(newNames[s])
	lay2.addWidget(widgets[s])
	lay.addLayout(lay2)

lay.addWidget(button)

dialog.show()
#fin de l'application lorsque toutes les fenêtres sont fermées
app.connect(app,QtCore.SIGNAL("lastWindowClosed()"),app,QtCore.SLOT("quit()"))
#fin de l'application lorsque l'utilisateur clique sur le bouton
app.connect(button, QtCore.SIGNAL("clicked()"),updateNewNames)
#boucle principale de traitement des évènements
app.exec_()

output = None
output2 = None

if 'oldnifti' in pr: # not in pr to use the cache !
	import subprocess
	output = subprocess.Popen(['find',root, '-iname', '*.img', '-o', '-iname', '*.pts','-o', '-iname', '*.nii', '-o', '-iname', 'electrode.txt', '-o', '-iname', 'electrodes.txt','-o', '-iname', 'DICOM', '-o', '-iname', '*.jpg', '-o', '-iname', '*.jpeg', '-o', '-iname', '*.pdf', '-o', '-iname', '*.trc'], 
														stdout=subprocess.PIPE).communicate()[0].split('\n')
	#output.extend(subprocess.Popen(['find',rootEpi, '-iname', '*.img', '-o', '-iname', '*.pts','-o', '-iname', '*.nii', '-o', '-iname', 'electrode.txt', '-o', '-iname', 'electrodes.txt','-o', '-iname', 'DICOM'], 
	#													stdout=subprocess.PIPE).communicate()[0].split('\n'))
	
	## Trier les .nii et .img
	images = [ o for o in output if o[-4:].lower()=='.img'] + [ o for o in output if o[-4:].lower()=='.nii']
	pr['oldnifti'] = images
	# Trier les dossiers DICOM et les étiqueter
	pr['dicom'] = [ o for o in output if o[-5:].upper()=='DICOM']
	## Trier les PTS et les electrode(s).txt
	pr['oldpts'] = [ o for o in output if o[-4:].lower()=='.pts']
	
	pr['oldelectrodestxt'] = [ o for o in output if o[-13:].lower()=='electrode.txt'] + [ o for o in output if o[-14:].lower()=='electrodes.txt']
	pr['impl'] = [ o for o in output if os.path.splitext(o)[1].lower()=='.jpg'] + [ o for o in output if os.path.splitext(o)[1].lower()=='.jpeg']
	pr['trc'] = [ o for o in output if os.path.splitext(o)[1].lower()=='.trc']

# What do we have for each subject ?
subjData=dict([(s,{'dicom':[], 'nifti':[], 'pts':[], 'electrodestxt':[]}) for s in subjects])

for s in subjects:
	subjData[s]['nifti']=[{'path':i,'type':None, 'acq':None} for i in pr['oldnifti'] if i.split(root)[1].split('/')[2] == s and i.split(root)[1].split('/')[3].lower().find('seeg') == -1 and os.path.split(i)[1].find('y_inverse') == -1 ]
	# if os.path.split(i)[1].startswith('r') -> image registered, 'w' pour les images normalizée, s pour smooth
	subjData[s]['dicom']=[{'path':i,'type':None, 'acq':None} for i in pr['dicom'] if i.split(root)[1].split('/')[2] == s]
	subjData[s]['pts']=[{'path':i,'type':None, 'acq':None} for i in pr['oldpts'] if i.split(root)[1].split('/')[2] == s]
	subjData[s]['electrodestxt']=[{'path':i,'type':None, 'acq':None} for i in pr['oldelectrodestxt'] if i.split(root)[1].split('/')[2] == s]
	subjData[s]['impl']=[{'path':i,'type':None, 'acq':None} for i in pr['impl'] if i.split(root)[1].split('/')[2] == s]
	subjData[s]['trc']=[{'path':i,'type':None, 'acq':None} for i in pr['trc'] if i.split(root)[1].split('/')[2] == s]
	subjData[s]['comment'] = "Has %d nifti, %d dicom dirs, %d pts, %d electrodes.txt and %d implantation scans"%(len(subjData[s]['nifti']), len(subjData[s]['dicom']), len(subjData[s]['pts']), len(subjData[s]['electrodestxt']), len(subjData[s]['impl']))
	# If not defined, try to extract name of subject from TRC
	if len(newNames[s]) < 10 and len(subjData[s]['trc']) > 0:
		output = subprocess.Popen(['getTRCpatientData', subjData[s]['trc'][0]['path']], stdout=subprocess.PIPE).communicate()[0].split('\n')
		print "Auto name : %s  --> oldname %s --> newName %s"%(s, newNames[s], str(output[2])+'+++')
		newNames[s] = str(output[2])+'+++'

# On veut -> trouver le nouveau nom, trouver une T1pre, une T1post, une T2pre/post si dispo, une CT..., tous les TRC et les localizers déjà faits et insérer tout ça dans la db.


##### Save preferences

prfile = open(os.path.join(os.path.expanduser('~'), '.autoImportBrainVisa'),'w')
json.dump(pr, prfile)        
prfile.close()
#################### EpiG ############################