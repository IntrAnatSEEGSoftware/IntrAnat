# -*- coding: utf-8 -*-
# Author: Francois Tadel
# License GNU GPL v3

import sys, os, json, time, datetime
from collections import  OrderedDict

from locateElectrodes import LocateElectrodes
# BrainVISA/anatomist imports
from soma import aims
from brainvisa import axon
from brainvisa.configuration import neuroConfig
from IPython.config.application import catch_config_error
# from IPython.config.application import catch_config_error
# from operator import pos
neuroConfig.gui = True
from brainvisa.data import neuroHierarchy
import brainvisa.registration as registration
from brainvisa.processes import *
from soma.qt_gui.qt_backend import QtGui, QtCore, uic
from brainvisa.data.readdiskitem import ReadDiskItem
        

def main(computeMni, start_index, logFilename):
    # Open log file
    log = open(logFilename, 'wb')
    
    # Start BrainVISA
    app = QtGui.QApplication(sys.argv)
    axon.initializeProcesses()
    w = LocateElectrodes(app=app, loadAll=True, isGui=False)

    # Find available patients in BV database
    rdi = ReadDiskItem( 'Subject', 'Directory',requiredAttributes={'_ontology':'brainvisa-3.2.0'}) #, requiredAttributes={'center':'Epilepsy'} )
    w.allSubjects = list( rdi._findValues( {}, None, False ))
    w.currentProtocol = 'Epilepsy'
    w.subjects = [s.attributes()['subject'] for s in w.allSubjects if 'center' in s.attributes() and s.attributes()['center'] == w.currentProtocol]
    w.subjects = sorted(w.subjects)

    # Export options
    if computeMni:
        selOptions = [True, True, False, False, True, True, False, False]
    else:
        selOptions = [False, True, False, False, True, True, False, False]
    
    # Loop on subjects
    for iSubj in range(start_index - 1, len(w.subjects)):
        # Write patient name to log
        tstamp = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
        log.write("[%s] Patient %4d/%4d: %18s...    " % (tstamp, iSubj+1, len(w.subjects), w.subjects[iSubj]))
        log.flush()
        
        # Load patient  (capture console output)
        try:
            w.loadPatientWorker(w.subjects[iSubj], thread=None, isGui=False)
        except: 
             log.write("ERROR: Could not load patient.\n")
        # Skip subjects where nothing was loaded
        if not w.dispObj:
            log.write("ERROR: No volume loaded.\n")
            continue
        # Export process  (capture console output)
        try:
            res = w.exportAllWorker(selOptions)
            # Log error message
            if res and res[1]:    # errMsg
                log.write("ERROR: " + " / ".join(res[1]) + "\n")
                continue
        except:
            log.write("ERROR: Could not export patient.\n")
            continue
        # It worked with no errors
        log.write("OK\n")
        log.flush()
        # Unload patient
        w.changePatient()
    # Close log file
    log.close()



# Calling from command line
if __name__ == "__main__":
    # Test input parameters
    if (len(sys.argv) < 2) or ((sys.argv[1] != "--mni-skip") and (sys.argv[1] != "--mni-recompute")):
        print("USAGE: generate_csv.py --mni-skip [start_index] [logfile]")
        print("       generate_csv.py --mni-recompute [start_index] [logfile]")
        sys.exit(2)
    elif (sys.argv[1] == "--mni-skip"):
        computeMni = False
    else:
        computeMni = True
        
    # Start index
    if (len(sys.argv) == 3):
        start_index = int(sys.argv[2])
    else:
        start_index = 1
        
    # Log file
    if (len(sys.argv) == 4):
        logFilename = sys.argv[3]
    else:
        logFilename = os.path.join(os.path.expanduser("~"), "log_csv_export.txt")
        
    main(computeMni, start_index, logFilename)
    
    