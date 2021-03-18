# -*- coding: utf-8 -*-
# Author: Francois Tadel
# License GNU GPL v3

import sys, os, json
from collections import  OrderedDict

from locateElectrodes import LocateElectrodes
from brainvisa import axon
from soma.qt_gui.qt_backend import QtGui, QtCore, uic


def main(inFile, outFolder):
    
    # Read input file
    f = open(inFile, 'r')
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
        for i in reversed(range(len(origName))):
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
    
    # Get input base name
    fPath, fName = os.path.split(inFile)
    baseName, fExt = os.path.splitext(fName)
    # Create output filenames
    fileEleclabel = os.path.join(outFolder, baseName + ".eleclabel")
    fileCsv = os.path.join(outFolder, baseName + ".csv")   
    
    # Start BrainVISA
    app = QtGui.QApplication(sys.argv)
    axon.initializeProcesses()
    w = LocateElectrodes()

    # Call export functions from locateElectrodes
    w.exportParcels2(True, True, plot_dict_MNI, None, fileEleclabel)
    w.exportCSVdictionaries(False, fileEleclabel, None, fileCsv, plot_dict_MNI)


# Calling from command line
if __name__ == "__main__":
    # Test input parameters
    if (len(sys.argv) < 3) or not os.path.exists(sys.argv[1]) or not os.path.exists(sys.argv[2]):
        print("USAGE: convert_mni2mcs.py input_mni.txt output_dir")
        sys.exit(2)
    main(sys.argv[1], sys.argv[2])
    
    