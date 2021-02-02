# -*- coding: utf-8 -*-
# Author: Francois Tadel, 2021
# License GNU GPL v3

import sys, os, time, datetime
import ImageImport
from soma import aims
from brainvisa import axon
# from brainvisa.configuration import neuroConfig
# neuroConfig.gui = True
from brainvisa.data import neuroHierarchy
# import brainvisa.registration as registration
# from brainvisa.processes import *
from soma.qt_gui.qt_backend import QtGui, QtCore, uic
from brainvisa.data.readdiskitem import ReadDiskItem


# ===== IMPORT FREESURFER =====
def importFreesurfer(w, subject, fsDir):
    errMsg = []
    # Load patient
    try:
        # Find T1pre for subject
        diT1 = ReadDiskItem('Raw T1 MRI', 'aims readable volume formats', requiredAttributes={'subject':subject, 'modality':'t1mri', 'normalized':'no'})
        rdiT1 = list(diT1._findValues({}, None, False))
        diT1pre = [x for x in rdiT1 if ('T1pre' in x.attributes()['acquisition'])]
        if (len(diT1pre) > 1):
            return [False, ['Multiple T1pre']]
        elif (len(diT1pre) == 0):
            return [False, ['No T1pre']]
        diT1pre = diT1pre[0]
    
        # Load patient
        w.brainvisaPatientAttributes = diT1pre.attributes()
        w.currentProtocol = diT1pre.attributes()['center']
        w.t1pre2ScannerBasedTransform = None
        w.electrodes = []
        
        # Load T1 in anatomist
        if w.diskItems and w.dispObj and w.diskItems['T1pre'] and w.dispObj['T1pre']:
            del w.dispObj['T1pre']
            del w.diskItems['T1pre']
        w.diskItems['T1pre'] = diT1pre
        w.dispObj['T1pre'] = w.a.loadObject(diT1pre)
        # Delete existing transformations, otherwise we can't match the T1pre and FreeSurferT1 scanner-based exactly
        w.deleteNormalizedTransf(w.dispObj['T1pre'])
        
        # Brain center
        if diT1pre.get('brainCenter'):
            w.t1preCenter = diT1pre.get('brainCenter')
        elif diT1pre.get('volume_dimension') and diT1pre.get('voxel_size'):
            volSize = diT1pre.get('volume_dimension')
            voxSize = diT1pre.get('voxel_size')
            w.t1preCenter = [volSize[0]*voxSize[0]/2, volSize[1]*voxSize[1]/2, volSize[2]*voxSize[2]/2];
        else:
            w.t1preCenter = [128, 128, 128]
        
        # Find FreeSurfer atlases
        diFs = ReadDiskItem('FreesurferAtlas', 'aims readable volume formats', requiredAttributes={'subject':subject, 'modality':'freesurfer_atlas'})
        rdiFs = list(diFs._findValues({}, None, False))
        # Destrieux
        diDes = [x for x in rdiFs if ('FreesurferAtlaspre' in x.attributes()['acquisition'])] 
        if diDes:
            w.vol_destrieux = aims.read(str(diDes[0]))
        else:
            w.vol_destrieux = None
        # DKT
        diDKT = [x for x in rdiFs if ('DKT' in x.attributes()['acquisition'])]
        if diDKT:
            w.vol_dkt = aims.read(str(diDKT[0]))
        else:
            w.vol_destrieux = None
        # HCP-MMP1
        diHCP = [x for x in rdiFs if ('HCP-MMP1' in x.attributes()['acquisition'])]
        if diHCP:
            w.vol_hcp = aims.read(str(diHCP[0]))
        else:
            w.vol_destrieux = None
    
        # Load electrodes and contacts
        w.loadElectrodes(subject, w.currentProtocol, isGui=False)
        if not w.electrodes:
            return [False, ['No implantation file']]
    except:
        errMsg += ["Could not load patient"]
        return [False, errMsg]
      
    # Compute parcels
    try:
        elecfile, errMsgCompute = w.computeParcels(diT1pre)
        errMsg += errMsgCompute
        if not elecfile:
            return [False, errMsg]
    except:
        errMsg += ["Could not compute parcels"]
        return [False, errMsg]

    return [True, errMsg]
    
    
# ===== MAIN =====
def main(input_dir, logFilename):
    # Start BrainVISA
    app = QtGui.QApplication(sys.argv)
    axon.initializeProcesses()
    w = ImageImport.ImageImport(app=app, isGui=False)

    # Find FreeSurfer database    
    FsSubjDir = None
    for db in neuroHierarchy.databases.iterDatabases():
        if db.directory.lower().find("freesurfer") != -1:
            FsSubjDir = db.directory
            break
    if not FsSubjDir:
        print "ERROR: No local FreeSurfer database found."
        return
    
    # Find subjects in BrainVISA database
    diSub = ReadDiskItem( 'Subject', 'Directory',requiredAttributes={'_ontology':'brainvisa-3.2.0'}) #, requiredAttributes={'center':'Epilepsy'} )
    rdiSub = list(diSub._findValues( {}, None, False ))
    bvSubjects = dict()
    bvCenter = dict()
    for s in rdiSub:
        sub = s.attributes()['subject'].split(' ')[0]
        bvSubjects[sub] = s.attributes()['subject']
        bvCenter[sub] = s.attributes()['center']
    if not bvSubjects:
        print "ERROR: No subjects to update in the BrainVISA database."
        return
    
    # Find subjects in input dir
    fsSubjects = dict()
    for sub in os.listdir(input_dir):
        if os.path.exists(os.path.join(input_dir, sub, 'freesurfer', 'mri', 'aparc.DKTatlas+aseg.mgz')):
            fsSubjects[sub.split(' ')[0]] = sub
    if not fsSubjects:
        print "ERROR: No FreeSurfer segmentations found in input folder."
        return
    
    # Subjects to import
    subList = sorted(list(set(fsSubjects.keys()) & set(fsSubjects.keys())))
    if not subList:
        print "ERROR: No subjects found both in input folder and BrainVISA database."
        return
    
    # Open log file
    log = open(logFilename, 'wb')
    # Save list of subjects to update
    log.write("Number of patients : %d\n" % len(subList))
    log.write("Folder to import   : " + input_dir + "\n")
    log.write("FreeSurfer database: " + FsSubjDir + "\n\n")
    maxLen = max([len(s) for s in subList])

    # Loop on subjects
    for iSub in range(len(subList)):
        # Write patient name to log
        tstamp = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
        strPatient = ("[%s] Patient #%4d: %" + str(maxLen) + "s...    ") % (tstamp, iSub+1, subList[iSub])
        print "\n" + strPatient
        log.write(strPatient)
        log.flush()
        # Current protocol: based on selected subject
        w.currentProtocol = bvCenter[subList[iSub]]
        # Import FreeSurfer folder
        importDir = os.path.join(input_dir, fsSubjects[subList[iSub]], 'freesurfer')
        isOk, errMsg = w.importFSoutput(bvSubjects[subList[iSub]], w.currentProtocol, importDir, FsSubjDir, False, False)
        # Handling error/success
        if isOk and not errMsg:
            log.write("OK\n")
        elif isOk:
            log.write("OK - WARNING: " + " / ".join(errMsg) + "\n")
        else: 
            log.write("ERROR: " + " / ".join(errMsg) + "\n")
        log.flush()
    # Close log file
    log.close()
    # Quit Qt application
    app.quit()
    del app


# ===== COMMAND LINE =====
if __name__ == "__main__":
    defLog = "log_importfs_" + datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H.%M.%S') + ".txt"
    # Test input parameters
    if (len(sys.argv) < 2) or not os.path.exists(sys.argv[1]):
        print("USAGE: batch_importfs.py input_dir [logfile=$HOME/" + defLog + "]")
        sys.exit(2)
    inputDir = sys.argv[1]
    # Log file
    if (len(sys.argv) == 3):
        logFilename = sys.argv[2]
    else:
        logFilename = os.path.join(os.path.expanduser("~"), defLog)
    # Call processing function
    main(inputDir, logFilename)
    # Close application
    os._exit(0)

