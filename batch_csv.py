# -*- coding: utf-8 -*-
# Author: Francois Tadel, 2018-2021
# License GNU GPL v3

import sys, os, time, datetime
import locateElectrodes
from soma import aims
from brainvisa import axon
from soma.qt_gui.qt_backend import QtGui, QtCore, uic
from brainvisa.data.readdiskitem import ReadDiskItem


# ===== SAVE CSV =====
def generateCsv(w, subject, isCsv, isBids):
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
            w.t1preCenter = [volSize[0]*voxSize[0]/2, volSize[1]*voxSize[1]/2, volSize[2]*voxSize[2]/2]
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
      
    # Generate CSV
    if isCsv:
        try:
            csvfile, errMsgCsv = w.saveCSV(diT1pre)
            errMsg += errMsgCsv
            if not csvfile:
                return [False, errMsg]
        except:
            errMsg += ["Could not save CSV"]
            return [False, errMsg]

    # Generate TSV
    if isBids and w.bidspath:
        try:
            tsvfiles, errMsgCsv = w.saveBidsTsv(diT1pre)
            errMsg += errMsgCsv
            if not tsvfiles:
                return [False, errMsg]
        except:
            errMsg += ["Could not save BIDS TSV"]
            return [False, errMsg]

    return [True, errMsg]
    
    
# ===== MAIN =====
def main(isCsv, isBids, start_index, logFilename):
    # Open log file
    log = open(logFilename, 'wb')
    # Start BrainVISA
    app = QtGui.QApplication(sys.argv)
    axon.initializeProcesses()
    w = locateElectrodes.LocateElectrodes(app=app, loadAll=True, isGui=False)
    # Find available patients in BV database
    rdi = ReadDiskItem( 'Subject', 'Directory',requiredAttributes={'_ontology':'brainvisa-3.2.0'}) #, requiredAttributes={'center':'Epilepsy'} )
    w.allSubjects = list( rdi._findValues( {}, None, False ))
    w.currentProtocol = 'Epilepsy'
    w.subjects = [s.attributes()['subject'] for s in w.allSubjects if 'center' in s.attributes() and s.attributes()['center'] == w.currentProtocol]
    w.subjects = sorted(w.subjects)
    log.write("Number of patients: %d\n\n" % len(w.subjects))
    maxLen = max([len(s) for s in w.subjects])
    # Loop on subjects
    for iSubj in range(start_index - 1, len(w.subjects)):
        # Write patient name to log
        tstamp = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
        strPatient = ("[%s] Patient #%4d: %" + str(maxLen) + "s...    ") % (tstamp, iSubj+1, w.subjects[iSubj])
        print "\n" + strPatient
        log.write(strPatient)
        log.flush()
        # Create CSV
        isOk, errMsg = generateCsv(w, w.subjects[iSubj], isCsv, isBids)
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
    defLog = "log_csv_" + datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H.%M.%S') + ".txt"
    # Help
    if (len(sys.argv) < 2) or ((sys.argv[1] != "--all") and (sys.argv[1] != "--csv-only") and (sys.argv[1] != "--bids-only")):
        print("USAGE: batch_csv.py --all       [start_index=1] [logfile=$HOME/" + defLog + "]")
        print("       batch_csv.py --csv-only  [start_index=1] [logfile=$HOME/" + defLog + "]")
        print("       batch_csv.py --bids-only [start_index=1] [logfile=$HOME/" + defLog + "]")
        sys.exit(2)
    # Get command
    if (sys.argv[1] == "--all"):
        isCsv = True
        isBids = True
    elif (sys.argv[1] == "--csv-only"):
        isCsv = True
        isBids = False
    elif (sys.argv[1] == "--bids-only"):
        isCsv = False
        isBids = True
    # Start index
    if (len(sys.argv) >= 3):
        start_index = int(sys.argv[2])
    else:
        start_index = 1
    # Log file
    if (len(sys.argv) >= 4):
        logFilename = sys.argv[3]
    else:
        logFilename = os.path.join(os.path.expanduser("~"), defLog)
    # Call processing function
    main(isCsv, isBids, start_index, logFilename)
    # Close application
    os._exit(0)
    
    