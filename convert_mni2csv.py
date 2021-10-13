# -*- coding: utf-8 -*-
# Author: Francois Tadel
# License GNU GPL v3

import sys, os, time, datetime, json

from soma import aims
from brainvisa import axon
from soma.qt_gui.qt_backend import QtGui, QtCore, uic
from brainvisa.data import neuroHierarchy
from brainvisa.data.readdiskitem import ReadDiskItem
from brainvisa.data.writediskitem import WriteDiskItem
from collections import  OrderedDict
import locateElectrodes
from externalprocesses import createItemDirs


def main(mniFile, subject):
    
    # ===== READ MNI COORDINATES =====
    # Read input file
    f = open(mniFile, 'r')
    lines = f.readlines()
    f.close()
    # Create dictionnary
    plot_dict_MNI = OrderedDict()
    for iLine in range(len(lines)):
        # Skip header line
        if (iLine == 0):
            continue
        # Get contact name and positions
        contactInfo = lines[iLine].split()
        # Separate electrode name / contact index (find last letter)
        origName = list(contactInfo[0])
        iLastLetter = None
        for i in reversed(list(range(len(origName)))):
            if not origName[i].isdigit():
                iLastLetter = i
                break
        if iLastLetter is None:
            continue
        # Upper case for all the letters except from "p" that stand for ' (prime)
        for i in range(iLastLetter+1):
            if (i == 0) or (origName[i] != 'p'):
                origName[i] = origName[i].upper()
            elif (i > 0) and (origName[i] == 'p'):
                origName[i] = "'"
        # Format CSV contact name
        cleanName = ''.join(origName[:iLastLetter+1]) + "%02d" % int(''.join(origName[iLastLetter+1:]))
        # Add contact to the list
        plot_dict_MNI[cleanName] = [float(contactInfo[1]), float(contactInfo[2]), float(contactInfo[3])]
    # Sort contact names
    plot_dict_MNI = OrderedDict(sorted(plot_dict_MNI.items()))

    # ===== START LOCATEELECTRODES =====
    # Start BrainVISA
    app = QtGui.QApplication(sys.argv)
    axon.initializeProcesses()
    w = locateElectrodes.LocateElectrodes(app=app, loadAll=True, isGui=False)
    
    # ===== LOAD/CREATE SUBJECT =====
    # Find patient in BV database
    w.currentProtocol = 'Epilepsy'
    rdi = ReadDiskItem( 'Subject', 'Directory', requiredAttributes={'_ontology':'brainvisa-3.2.0', 'subject':subject, 'center':w.currentProtocol})
    diSubj = rdi.findValue({})
    # If subject does not exist: create it
    if not diSubj:
        wdi = WriteDiskItem('Subject', 'Directory', requiredAttributes={'_ontology':'brainvisa-3.2.0', 'subject':subject, 'center':w.currentProtocol})
        diSubj = wdi.findValue({})
        # Create directory that do not exist yet
        createItemDirs(diSubj)
        if not os.path.exists(diSubj.fileName()):
            os.mkdir(diSubj.fileName())
        neuroHierarchy.databases.insertDiskItem(diSubj, update=True)
    
    # Load patient
    w.brainvisaPatientAttributes = diSubj.attributes()
    w.t1pre2ScannerBasedTransform = None
    # w.diskItems['T1pre'] = diSubj
    
    # Compute parcels
    elecfile = w.computeParcels(diSubj, plot_dict_MNI)
    # Generate CSV
    csvfile = w.saveCSV(diSubj, plot_dict_MNI)

    # Quit Qt application
    app.quit()
    del app


# Calling from command line
if __name__ == "__main__":
    # Test input parameters
    if (len(sys.argv) < 3) or not os.path.exists(sys.argv[1]):
        print("USAGE: convert_mni2mcs.py input_mni.txt subject_id")
        sys.exit(2)
    # Run conversion
    main(sys.argv[1], sys.argv[2])
    # Close application
    os._exit(0)
    
    
    