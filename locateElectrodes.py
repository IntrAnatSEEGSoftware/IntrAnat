#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Localisation graphique des electrodes
# (c) Inserm U836 2012-2021 - Manik Bhattacharjee, Pierre Deman, Francois Tadel 
# License GNU GPL v3
 
# Standard Python imports
import sys, os, pickle, csv, numpy, re, string, time, subprocess, json, io
from numpy import *
from scipy import ndimage
from collections import OrderedDict
from PIL import Image
from collections import Counter

# BrainVISA/anatomist imports
from soma import aims
from brainvisa import axon
from brainvisa.configuration import neuroConfig
neuroConfig.gui = True
from brainvisa.data import neuroHierarchy
import brainvisa.registration as registration
from brainvisa.processes import *
from soma.qt_gui.qt_backend import QtGui, QtCore, uic, QtWidgets
from brainvisa.data.readdiskitem import ReadDiskItem
from brainvisa.data.writediskitem import WriteDiskItem
from brainvisa import anatomist
        
# IntrAnat local imports
from externalprocesses import *
from referentialconverter import ReferentialConverter
from readLabels import readLabels
from readFreesurferLabelFile import readFreesurferLabelFile
from TimerMessageBox import *
from generate_contact_colors import *
from electrode import ElectrodeModel
from bipoleSEEGColors import bipoleSEEGColors
from DialogCheckbox import DialogCheckbox
from progressbar import ProgressDialog


# =============================================================================
# ===== SPM CALLS =============================================================
# =============================================================================
spm12_inverse_y = """try
    addpath(genpath(%s));
    spm('defaults', 'FMRI');
    spm_jobman('initcfg');
    clear matlabbatch;
    matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def = {%s};
    matlabbatch{1}.spm.util.defs.comp{1}.inv.space = {%s};
    matlabbatch{1}.spm.util.defs.out{1}.savedef.ofname = %s;
    matlabbatch{1}.spm.util.defs.out{1}.savedef.savedir.saveusr = {%s};
    spm_jobman('run',matlabbatch);
catch
    disp 'AN ERROR OCCURED'; 
end
quit;"""

# Read deformation field y_<subject>_inverse.nii, apply the vectorf field to scanner-based coordinates of electrodes
# FT 24-Feb-2019: NEW OPTIMIZED VERSION WITHOUT LOOPS
spm12_normalizePoints = """try
    fprintf('Reading inverse MNI transformation... ');
    addpath(genpath(%s));
    P='%s';
    P1=spm_vol([P ',1,1']);
    P2=spm_vol([P ',1,2']);
    P3=spm_vol([P ',1,3']);
    fprintf(' X '); drawnow('update');
    V1=spm_read_vols(P1);
    fprintf(' Y '); drawnow('update');
    V2=spm_read_vols(P2);
    fprintf(' Z\\n'); drawnow('update');
    V3=spm_read_vols(P3);

    %% Apply transformation to electrodes
    PosElectrode = dlmread('%s')';
    %% Get the corresponding coordinates in the transformation volume
    Tinv = inv(P1(1).mat);
    ijk = Tinv * [PosElectrode; ones(1,size(PosElectrode,2))];
    %% Interpolate values for contacts
    wPosElectrode = [...
        interp3(V1, ijk(2,:), ijk(1,:), ijk(3,:), 'linear')', ...
        interp3(V2, ijk(2,:), ijk(1,:), ijk(3,:), 'linear')', ...
        interp3(V3, ijk(2,:), ijk(1,:), ijk(3,:), 'linear')'];
    %% Detect if any point outside of the MNI transformation volume
    iOut = find(any(isnan(wPosElectrode),2));
    %% Fix them with spline interpolation (much slower)
    for i = 1:length(iOut)
        wPosElectrode(iOut(i),:) = [...
            interp3(V1, ijk(2,iOut(i)), ijk(1,iOut(i)), ijk(3,iOut(i)), 'spline')', ...
            interp3(V2, ijk(2,iOut(i)), ijk(1,iOut(i)), ijk(3,iOut(i)), 'spline')', ...
            interp3(V3, ijk(2,iOut(i)), ijk(1,iOut(i)), ijk(3,iOut(i)), 'spline')'];
    end
    %% Write output file
    dlmwrite('%s',wPosElectrode,'precision',18);
catch
    disp 'AN ERROR OCCURED'; 
end
quit;"""


#################### READ ATLAS LABELS #########################################
labels = dict()
(labels['Destrieux'], freesurfer_colormap) = readFreesurferLabelFile('labels/freesurfer_labels.txt', 14175)
labels['DKT'] = labels['Destrieux']
(labels['HCP-MMP1'], hcp_colormap) = readFreesurferLabelFile('labels/hcp_l0r200_labels.txt', 14175)
(labels['VEP'], vep_colormap) = readFreesurferLabelFile('labels/VepFreeSurferColorLut.txt', 72077)
(labels['Lausanne2008-33'], lausanne33_colormap) = readFreesurferLabelFile('labels/lausanne33_labels.txt', 83)
(labels['Lausanne2008-60'], lausanne60_colormap) = readFreesurferLabelFile('labels/lausanne60_labels.txt', 129)
(labels['Lausanne2008-125'], lausanne125_colormap) = readFreesurferLabelFile('labels/lausanne125_labels.txt', 234)
(labels['Lausanne2008-250'], lausanne250_colormap) = readFreesurferLabelFile('labels/lausanne250_labels.txt', 463)
(labels['Lausanne2008-500'], lausanne500_colormap) = readFreesurferLabelFile('labels/lausanne500_labels.txt', 1015)
labels['MarsAtlas'] = readLabels('labels/marsatlas_labels.txt')

# Functions to sort the contacts A1,A2...,A10 and not A1, A10, A2..
def atoi(text):
    return int(text) if text.isdigit() else text
def natural_keys(text):
    """alist.sort(key=natural_keys) sorts in human order"""
    return [ atoi(c) for c in re.split('(\d+)', text) ]


##################### Electrode functions ##################################
def moveElectrode(target, entry, referential, newRef, a, meshes):
    # a is anatomist object
    if entry is None or target is None:
        return

    if newRef is None:
        newRef = a.createReferential()
    transl = target
    i = array([[1,0,0]])
    z = array([target,]) - array([entry,])
    if z[0][1] < 0.001 and z[0][2] < 0.001:# si i est colinéaire avec z, changer i car on n'obtiendra pas un vecteur perpendiculaire aux deux
        i = array([[0,1,0]])
    if linalg.norm(z) == 0:
        return
    z = -z / linalg.norm(z)
    y = cross (i,z)
    y = -y/linalg.norm(y)
    x = cross(y,z)
    m = [x[0][0], y[0][0], z[0][0], x[0][1],y[0][1], z[0][1], x[0][2],  y[0][2], z[0][2]]
    try:
        transf = a.createTransformation(transl+m, origin = newRef, destination = referential)
    except:
        return
    a.assignReferential(newRef, meshes)
    return (newRef, transf)


def createElectrode(target, entry, referential, ana=None, windows=None, model = None, dispMode=None, dispParams=None, isTransparent=False):
    elecModel = ElectrodeModel(ana)
    elecModel.open(model,dispMode, dispParams)
    #elecModel.setDisplayReferential(newRef)
    meshes = elecModel.getAnatomistObjects()
    (newRef, transf) = moveElectrode(target, entry, referential, None, ana, meshes)
    # Make electrodes transparent
    if isTransparent:
        for m in meshes:
            mat = m.getInternalRep().GetMaterial().genericDescription()['diffuse']
            m.setMaterial(diffuse=[mat[0], mat[1], mat[2], 0])
    # Add objects to all windows
    if windows is not None:
        ana.addObjects(meshes, windows)
    return (newRef, transf, elecModel)


def createBipole(target, entry, referential, ana=None, windows=None, model = None, dispMode=None, dispParams=None):
    elecModel = ElectrodeModel(ana)
    elecModel.open(str(model),dispMode, dispParams)
    #elecModel.setDisplayReferential(newRef)
  
    meshes = elecModel.getAnatomistObjects()
    (newRef, transf) = moveElectrode(target, entry, referential, None, ana, meshes)

    if windows is not None:
        ana.addObjects(meshes, windows)
    return (newRef, transf, elecModel)


# Récupération des plots dans l'ordre plot1 -> plot n
def getPlots(elecModel):
    cyls = elecModel.getCylinders()
    return dict((name,cyls[name]) for name in cyls if cyls[name]['type'] == 'Plot')

def getPlotsNames(elecModel):
    cyls = elecModel.getCylinders()
    plots = [n for n in cyls if cyls[n]['type'] == 'Plot']
    return sorted(plots, key=natural_keys)

# Récupération des coordonnées des centres des plots dans le référentiel électrode
def getPlotsCenters(elecModel):
    plots = getPlots(elecModel)
    return dict((p, plots[p]['center']) for p in plots)


# ==========================================================================
# ===== MAIN WINDOW ========================================================
# ==========================================================================

class LocateElectrodes(QtWidgets.QDialog):

    def __init__(self, app=None, loadAll=True, isGui=True):        
        # UI init
        if loadAll == True:
            QtWidgets.QDialog.__init__(self)
            self.ui = uic.loadUi("locateElectrodes.ui", self)
            self.setWindowTitle('Localisation des electrodes - NOT FOR MEDICAL USE')
            self.app = app
        self.nameEdit.setText('A')
    
        # Load the list of protocols, patients and electrode models from BrainVisa
        if loadAll == True:
            self.modalities = ['Raw T1 MRI', 'T2 MRI', 'CT', 'PET', 'Electrode Implantation Coronal Image', 'Electrode Implantation Sagittal Image','fMRI-epile', 'Statistic-Data','FLAIR', 'resection', 'FreesurferAtlas', 'FGATIR']
            # Electrode models
            self.elecModelList = []
            self.elecModelListByName = []
            self.loadFromBrainVisa()

        # Init of variables
        self.dispObj = {} # All displayable objects "t1mri-pre", "t2"...
        self.objtokeep = {} #all object we must keep alive for anatomist but not in other variables
        self.diskItems = {} # For each dispObj, a list of dictionnaries {'refName', 'refObj', 'refId', 'transf'}
        # Coordinates displayed using referential : 'Natif' par defaut
        self.coordsDisplayRef = 'Natif'
        self.referentialCombo.clear()
        self.referentialCombo.addItems(['Natif',])
        self.dispMode = 'real'
        self.dispParams = None
        self.t1pre2ScannerBasedId = None
        self.electrodes = []# {Les objets electrodes avec les coordonnées, les meshes
        self.bipoles = [] #{Les objects bipoles}
        self.electrodeTemplateStubs = [] # Un objet electrode par template disponible dans la base de données (chargé à la demande par getElectrodeTemplates)
        self.contacts = [] # {name:n, number:2, electrode: e, mesh:m}
        self.currentWindowRef = None # Referential used by windows (because pyanatomist AWindow.getReferential is not implemented yet)
        self.threads = [] # List of running threads
        self.t1pre2ScannerBasedTransform = None #Transfo from T1pre native to scanner-based referential (Anatomist Transformation object)
        self.t1preCenter = []
        self.brainvisaPatientAttributes = None # Attributes of a BrainVisa ReadDiskItem MRI of the loaded patient
        self.isModified = False
        self.lastClickedPos = None
        self.lastClickedRef = None
        self.windowContent = {}  # list of objects to display in window for each scenario (MNI, pre, post, etc)
        self.updateComboboxes()
        
        # Start anatomist
        if (loadAll == True):
            self.a = anatomist.Anatomist('-b') #Batch mode (hide Anatomist window)
        # Create Anatomist windows
        if (loadAll == True) and isGui:
            #self.a.config()['setAutomaticReferential'] = 1
            #self.a.config()['commonScannerBasedReferential'] = 1
            self.a.onCursorNotifier.add(self.clickHandler)
            
            # Create 4 windows
            self.wins = []
            for wcont in [self.windowContainer1, self.windowContainer2, self.windowContainer3, self.windowContainer4]:
                # Create anatomist window
                w = self.a.createWindow('Sagittal', options={'hidden': 1})
                #w.setParent(wcont)#, no_decoration=True )
                self.wins.append(w)
                # Add to main window
                wlayout = QtGui.QHBoxLayout(wcont)
                wlayout.setSpacing(0)
                #wlayout.setMargin(0)
                wlayout.addWidget(w.getInternalRep())
                # Configure window further
                w.getInternalRep().menuBar().setVisible(False)
                #w.getInternalRep().layout().setMargin(0)
                #w.getInternalRep().centralWidget().layout().setMargin(0)
                w.getInternalRep().findChild(QtGui.QToolBar,'controls').setVisible(False)
                w.getInternalRep().statusBar().setSizeGripEnabled(False)
                
            # Set correctwindow type (not done at the creation time, otherwise, they get different sizes)
            self.wins[0].getInternalRep().muteAxial()
            self.wins[1].getInternalRep().muteSagittal()
            self.wins[2].getInternalRep().muteCoronal()
            self.wins[3].getInternalRep().mute3D()
            # Hide control window
            self.a.getControlWindow().setVisible(False)
            # By default: use only two views
            self.setWindowNumber(2)

            # ===== SET CALLBACKS =====
            # Patient selection
            self.loadPatientButton.clicked.connect(self.loadPatient)
            self.changePatientButton.clicked.connect(self.changePatient)
            self.patientList.itemDoubleClicked.connect(lambda x:self.loadPatient())
            self.protocolCombo.currentIndexChanged[int].connect(self.updateBrainvisaProtocol)
            self.filterSiteCombo.currentIndexChanged[int].connect(self.filterSubjects)
            self.filterYearCombo.currentIndexChanged[int].connect(self.filterSubjects)
            # Image manipulation
            self.makefusionButton.clicked.connect(self.makeFusion)
            # self.connect(self.generateResectionArray.clicked.connect(self.generateResection)
            # self.connect(self.validateROIresection.clicked.connect(self.ROIResectiontoNiftiResection)
            self.deleteMarsAtlasfiles.clicked.connect(self.deleteMarsAtlasFiles)
            # Electrodes
            self.addElectrodeButton.clicked.connect(self.addElectrode)
            self.removeElectrodeButton.clicked.connect(self.removeElectrode)
            self.nameEdit.editingFinished.connect(self.editElectrodeName)
            self.typeComboBox.currentIndexChanged[str].connect(self.updateElectrodeModel)
            self.targetButton.clicked.connect(self.updateTarget)
            self.entryButton.clicked.connect(self.updateEntry)
            self.electrodeList.currentRowChanged.connect(self.electrodeSelect)
            self.electrodeList.itemDoubleClicked.connect(self.electrodeGo)
            self.contactList.itemClicked.connect(self.contactSelect)
            self.contactList.itemDoubleClicked.connect(self.contactGo)
            self.electrodeSaveButton.clicked.connect(self.saveElectrodes)
            self.electrodeLoadButton.clicked.connect(self.loadElectrodes)
            self.normalizeExportButton.clicked.connect(self.exportAll)
            # Display options
            self.colorConfigButton.clicked.connect(self.configureColors)
            self.dispModeCombo.currentIndexChanged[int].connect(self.updateDispMode)
            self.Clipping_checkbox.clicked.connect(self.clippingUpdate)
            self.electrodeRefCheck.stateChanged.connect(self.updateElectrodeView)
            # self.connect(self.electrodeRefRotationSlider, QtCore.SIGNAL('valueChanged(int)'), self.updateElectrodeViewRotation)
            self.electrodeRefRotationSlider.hide()
            # Show anatomist
            self.showAnatomist.setIcon(QtGui.QIcon('logoAnatomist.png'))
            self.showAnatomist.setIconSize(QtGui.QSize(24, 24))
            self.showAnatomist.clicked.connect(self.toggleAnatomistWindow)
            # Anatomist windows
            self.windowCombo1.currentIndexChanged[str].connect(lambda s: self.updateWindow(0, s, True))
            self.windowCombo2.currentIndexChanged[str].connect(lambda s: self.updateWindow(1 ,s, True))
            self.windowCombo3.currentIndexChanged[str].connect(lambda s: self.updateWindow(2, s, True))
            self.windowCombo4.currentIndexChanged[str].connect(lambda s: self.updateWindow(3 ,s, True))
            self.referentialCombo.currentIndexChanged[str].connect(self.updateCoordsDisplay)
            self.buttonWin1.clicked.connect(lambda : self.setWindowNumber(1))
            self.buttonWin2.clicked.connect(lambda : self.setWindowNumber(2))
            self.buttonWin4.clicked.connect(lambda : self.setWindowNumber(4))
            

            # List of controls to enable when a subject is loaded
            self.widgetsLoaded = [self.loadPatientButton, self.patientList, self.protocolCombo, self.filterSiteCombo, self.filterYearCombo]
            self.widgetsUnloaded = [self.changePatientButton, self.groupManip, self.groupDisplay, self.referentialCombo, self.referentialCombo,\
                                    self.referentialCombo, self.addElectrodeButton, self.removeElectrodeButton, self.nameEdit, \
                                    self.typeComboBox, self.targetButton, self.entryButton, self.electrodeList, self.contactList, \
                                    self.electrodeSaveButton, self.electrodeLoadButton, \
                                    self.windowCombo1, self.windowCombo2, self.windowCombo3, self.windowCombo4, \
                                    self.windowContainer1, self.windowContainer2, self.windowContainer3, self.windowContainer4]
            # Update enabled/disabled controls
            for w in self.widgetsLoaded:
                w.setEnabled(True)
            for w in self.widgetsUnloaded:
                w.setEnabled(False)
                
            # Display warning: Not for medical use
            self.warningMEDIC()
        # Initialize empty GUI if not set yet (for command line calls)
        if not hasattr(self, 'wins'):
             self.wins = []
             self.widgetsLoaded = []
             self.widgetsUnloaded = []  
            
        # ===== LOAD PREFERENCES =====
        if (loadAll == True):
            self.spmpath = None
            self.bidspath = None
            prefpath_imageimport = os.path.join(os.path.expanduser('~'), '.imageimport')
            try:
                if (os.path.exists(prefpath_imageimport)):
                    filein = open(prefpath_imageimport, 'rU')
                    prefs_imageimport = pickle.load(filein)
                    if 'spm' in prefs_imageimport.keys():
                        self.spmpath = prefs_imageimport['spm']
                    if 'bids' in prefs_imageimport.keys():
                        self.bidspath = prefs_imageimport['bids']
                    filein.close()
            except:
                pass

        # Get BrainVISA context
        self.brainvisaContext = defaultContext()
        # Get Transformation Manager
        self.transfoManager = registration.getTransformationManager()
        # Get ReferentialConverter (for Talairach, AC-PC...)
        self.refConv = ReferentialConverter()


    # ==========================================================================
    # ===== WINDOW EVENTS ======================================================
    # ==========================================================================
    def closeEvent(self, event):
        """Function called when there is a request to close the app window."""
        # Ask for confirmation
        if self.isModified:
            reply = QtGui.QMessageBox.question(self, 'Message', "Quit the software without saving?",
                    QtGui.QMessageBox.Yes | QtGui.QMessageBox.No, QtGui.QMessageBox.No)
            isClose = (reply == QtGui.QMessageBox.Yes)
        else:
            isClose = True
        # If user confirmed
        if isClose:
            self.quit()
            if event:
                event.accept()
        elif event:
            event.ignore()
        

    def quit(self):
        """Quits application."""
        # Closes BrainVISA
        axon.cleanup()
        # Quits Qt app
        self.app.quit()

    
    def keyPressEvent( self, event ) :
        if (event.key()==QtCore.Qt.Key_Return):
            if hasattr(self, 'app') and (self.app.focusWidget() == self.nameEdit):
                self.editElectrodeName()
            event.accept()


    # ==========================================================================
    # ===== OTHER FUNCTIONS ====================================================
    # ==========================================================================

    def warningMEDIC(self):
        shortwarning = TimerMessageBox(5,self)
        shortwarning.exec_()


    def toggleAnatomistWindow(self):
        self.a.getControlWindow().setVisible(self.showAnatomist.isChecked())

    def loadFromBrainVisa(self):
        # Find available patients in BV database
        rdi = ReadDiskItem( 'Subject', 'Directory',requiredAttributes={'_ontology':'brainvisa-3.2.0'}  ) #, requiredAttributes={'center':'Epilepsy'} )
        subjects = list( rdi._findValues( {}, None, False ) )
        protocols = list(set([s.attributes()['center'] for s in subjects if 'center' in s.attributes()]))
        # Fill the combo
        self.protocolCombo.clear()
        self.protocolCombo.addItems(sorted(protocols))
        self.allSubjects = subjects
        self.updateBrainvisaProtocol()


    def updateBrainvisaProtocol(self, idx=None):
        """Updates the UI when the selected protocol changes"""
        self.currentProtocol = str(self.protocolCombo.currentText())
        self.subjects = [s.attributes()['subject'] for s in self.allSubjects if 'center' in s.attributes() and s.attributes()['center'] == self.currentProtocol]
        self.patientList.clear()
        self.patientList.addItems(sorted(self.subjects))
        
        # Update the filters
        years = []
        sites = []
        for s in self.subjects:
            subsplit = s.split('_')
            if len(subsplit) > 1:
                # CHUGA subject id format: "GRE_YYYYY_XXXz"
                if len(subsplit[1]) == 4:
                    years += [subsplit[1]]
                    sites += [subsplit[0]]
                # ftract subject id format: "0016GRE_DDMMYYYY"
                elif len(subsplit[1]) == 8:
                    years += [subsplit[1][4:]]
                    sites += [subsplit[0][4:7]]
        years = ['*',] + sorted(list(set(years)))
        sites = ['*',] + sorted(list(set(sites)))
        
        self.filterSiteCombo.clear()
        self.filterSiteCombo.addItems(sites)
        self.filterYearCombo.clear()
        self.filterYearCombo.addItems(years)
        # Loading electrode models
        self.loadElectrodeModels()


    # Loading electrode models
    def loadElectrodeModels(self):
        # rdiEM = ReadDiskItem('Electrode Model', 'Electrode Model format', requiredAttributes={'center':self.currentProtocol})
        rdiEM = ReadDiskItem('Electrode Model', 'Electrode Model format')
        self.elecModelList = list (rdiEM._findValues( {}, None, False ) )
        elecNames = [e.attributes()['model_name'] for e in self.elecModelList]
        self.elecModelListByName = dict((e.attributes()['model_name'], e) for e in self.elecModelList)
        self.typeComboBox.clear()
        listItems = sorted(elecNames)
        listItems += ['[Reload electrode models]']
        self.typeComboBox.addItems(listItems)


    def filterSubjects(self, value=None):
        """Filtering subject list"""
        subs = self.subjects
        if str(self.filterSiteCombo.currentText()) != '*':
            subs = [s for s in subs if ((s.split('_')[0] == str(self.filterSiteCombo.currentText())) or \
                                        (len(s.split('_')[0]) == 7) and (s.split('_')[0][4:7] == str(self.filterSiteCombo.currentText())))]
        if str(self.filterYearCombo.currentText()) != '*':
            subs = [s for s in subs if len(s.split('_')) > 1 and str(self.filterYearCombo.currentText()) in s.split('_')[1]]
        self.patientList.clear()
        self.patientList.addItems(sorted(subs))


    def getT1preMniTransform(self):
        """Returns the path of the MNI transformation and compute it if necessary"""
        # If not T1pre: error
        if 'T1pre' not in self.diskItems:
            print "No T1pre loaded : cannot get MNI transform from it"
            return None
        # Look for an existing y inverse file
        rdi_inv_read = ReadDiskItem('SPM normalization inverse deformation field', 'NIFTI-1 image')
        di_inv_read = rdi_inv_read.findValue(self.diskItems['T1pre'])
        if di_inv_read:
            return di_inv_read.fileName()
        # Look for a y file
        rdi_y = ReadDiskItem('SPM normalization deformation field','NIFTI-1 image')
        di_y = rdi_y.findValue(self.diskItems['T1pre'])
        if not di_y:
            print "SPM MNI deformation field not found: Recompute MNI normalization from ImageImport."
            return None
        # Compute y inverse
        else:
            # Get y_inverse file path
            wdi_inverse = WriteDiskItem('SPM normalization inverse deformation field','NIFTI-1 image')
            dir_yinv_split = str(di_y.fileName()).split('/')
            name_yinverse = dir_yinv_split.pop()[2:]
            dir_yinverse = "/".join(dir_yinv_split)
            di_inv_write = wdi_inverse.findValue(di_y)
            # Compute inverse deformation with SPM12
            errMsg = matlabRun(spm12_inverse_y%("'"+self.spmpath+"'","'"+str(di_y.fileName())+"'","'"+self.dispObj['T1pre'].fileName()+"'","'"+name_yinverse.replace('.nii','_inverse.nii')+"'","'"+dir_yinverse+"'"))
            if errMsg:
                print errMsg
                return None
            # Check that output file exists
            if os.path.exists(di_inv_write.fileName()):
                neuroHierarchy.databases.insertDiskItem(di_inv_write, update=True)
                return di_inv_write.fileName()
            else:
                print "Error: Could not compute SPM MNI inverse deformation field"
                return None


    def changePatient(self):
        # Update enabled/disabled controls
        for w in self.widgetsLoaded:
            w.setEnabled(True)
        for w in self.widgetsUnloaded:
            w.setEnabled(False)
        # Delete all the graphical objects
        self.a.removeObjects(self.a.getObjects(), self.wins)
        # Remove unused referentials
        referentials = self.a.getReferentials()
        for element in referentials:
            if element.getInfos() and (element.getInfos().get('name') not in ('Talairach-MNI template-SPM', 'Talairach-AC/PC-Anatomist')):
                self.a.deleteElements(element)
        # Re-initialize the entire window
        self.electrodeList.clear()
        self.contactList.clear()
        self.currentElectrodes = []
        self.currentContacts = []
        self.__init__(loadAll=False, isGui=True)
  
  
    def loadPatient(self, patient=None):
        # If multiple subjects selected: Select only one
        if (len(self.patientList.selectedItems()) > 1):
            self.patientList.setCurrentItem(self.patientList.selectedItems()[0])
        # Get current patient
        if not patient:
            patient = str(self.patientList.selectedItems()[0].text())
        # Block callbacks
        self.windowCombo1.blockSignals(True)
        self.windowCombo2.blockSignals(True)
        self.windowCombo3.blockSignals(True)
        self.windowCombo4.blockSignals(True)
        # Call loading function
        errMsg = ProgressDialog.call(lambda thr:self.loadPatientWorker(patient, thr), True, self, "Processing...", "Load patient: " + patient)
        #errMsg = self.loadPatientWorker(patient)
        # Update enabled/disabled controls
        for w in self.widgetsLoaded:
            w.setEnabled(False)
        for w in self.widgetsUnloaded:
            w.setEnabled(True)
        # Restore callbacks
        self.windowCombo1.blockSignals(False)
        self.windowCombo2.blockSignals(False)
        self.windowCombo3.blockSignals(False)
        self.windowCombo4.blockSignals(False)
        # Display all
        self.updateAllWindows(True)
        # Display error message
        if errMsg:
            errMsg += ["", "The various images may not be aligned correctly.", 
                       "Close LocateElectrodes, start ImageImport, delete the extra images (the ones added last), then start again LocateEletrodes and try again. ", 
                       "If the images and electrodes are still not aligned correctly, you may need to delete the patient and start over."]
            QtGui.QMessageBox.warning(self, u'Database error', "\n".join(errMsg))
    

    def loadPatientWorker(self, patient, thread=None, isGui=True):
        errMsg = []

        self.t1pre2ScannerBasedTransform = None
        pre_select_1 = self.windowCombo1.currentText()
        pre_select_2 = self.windowCombo2.currentText()
        pre_select_3 = self.windowCombo3.currentText()
        pre_select_4 = self.windowCombo4.currentText()
    
        # Get all subject's volumes
        volumes = []
        for moda in self.modalities:
            rdi2 = ReadDiskItem(moda, 'aims readable volume formats', requiredAttributes={'subject':patient, 'center':self.currentProtocol})
            volumes.extend(list(rdi2._findValues({}, None, False)))
        # Load all volumes
        dictionnaire_list_images = dict()
        objAtlas = []
        refT1pre = None
        self.vol_destrieux = []
        self.vol_dkt = []
        self.vol_hcp = []
        for t in volumes:
            if "skull_stripped" in t.fullName():
                continue
            # Progres bar
            if thread is not None:
                thread.progress_text.emit("Loading volume: " + t.attributes()['modality'] + "...")
            # Keep attributes of any image from this subject
            self.brainvisaPatientAttributes = t.attributes()
            na = None
            # Add elements in the display list
            if (t.attributes()['modality'] == 't1mri') and ('FreesurferAtlaspre' in t.attributes()['acquisition']):
                if ('MRI FreeSurfer' in dictionnaire_list_images.keys()):
                    errMsg += ["Multiple FreesurferAtlaspre images found"]
                dictionnaire_list_images.update({'MRI FreeSurfer':['FreesurferT1pre', 'electrodes']})
                na = 'FreesurferT1pre'
            elif (t.attributes()['modality'] == 't1mri') and ('T1pre' in t.attributes()['acquisition']) and (t.attributes()['normalized'] == 'no'):
                if ('MRI pre' in dictionnaire_list_images.keys()):
                    errMsg += ["Multiple T1pre images found"]
                dictionnaire_list_images.update({'MRI pre':['T1pre', 'electrodes']})
            elif (t.attributes()['modality'] == 't1mri') and ('T1pre' in t.attributes()['acquisition']) and (t.attributes()['normalized'] == 'yes') :
                if ('Normalized MRI pre' in dictionnaire_list_images.keys()):
                    errMsg += ["Multiple T1pre-Norm images found"]
                dictionnaire_list_images.update({'Normalized MRI pre':['preNorm', 'electrodes']})
                na = 'preNorm'         
            elif (t.attributes()['modality'] == 't1mri') and ('postOp' in t.attributes()['acquisition']) and (t.attributes()['normalized'] == 'no'):
                if ('MRI post-op' in dictionnaire_list_images.keys()):
                    errMsg += ["Multiple T1postOp images found"]
                dictionnaire_list_images.update({'MRI post-op':['T1postOp', 'electrodes']})
            elif (t.attributes()['modality'] == 't1mri') and ('postOp' in t.attributes()['acquisition']) and (t.attributes()['normalized'] == 'yes') :
                if ('Normalized MRI post-op' in dictionnaire_list_images.keys()):
                    errMsg += ["Multiple T1postOp-Norm images found"]
                dictionnaire_list_images.update({'Normalized MRI post-op':['postOpNorm', 'electrodes']})
                na = 'postOpNorm'   
            elif (t.attributes()['modality'] == 't1mri') and ('post' in t.attributes()['acquisition']) and not ('postOp' in t.attributes()['acquisition']):
                if ('MRI post' in dictionnaire_list_images.keys()):
                    errMsg += ["Multiple T1post images found"]
                dictionnaire_list_images.update({'MRI post':['T1post', 'electrodes']})
                if not pre_select_1:
                    pre_select_1 = 'MRI post'
                if not pre_select_2:
                    pre_select_2 = 'MRI post'
                if not pre_select_3:
                    pre_select_3 = 'MRI post'
            elif (t.attributes()['modality'] == 't2mri') and ('pre' in t.attributes()['acquisition']):
                if ('MRI pre T2' in dictionnaire_list_images.keys()):
                    errMsg += ["Multiple T2pre images found"]
                dictionnaire_list_images.update({'MRI pre T2':['T2pre', 'electrodes']})
            elif (t.attributes()['modality'] == 't2mri') and ('postOp' in t.attributes()['acquisition']):
                if ('MRI post-op T2' in dictionnaire_list_images.keys()):
                    errMsg += ["Multiple T2postOp images found"]
                dictionnaire_list_images.update({'MRI post-op T2':['T2postOp', 'electrodes']})
            elif (t.attributes()['modality'] == 't2mri') and ('post' in t.attributes()['acquisition']):
                if ('MRI post T2' in dictionnaire_list_images.keys()):
                    errMsg += ["Multiple T2post images found"]
                dictionnaire_list_images.update({'MRI post T2':['T2post', 'electrodes']})
            elif (t.attributes()['modality'] == 'ct') and ('post' in t.attributes()['acquisition']) and not ('postOp' in t.attributes()['acquisition']):
                if ('CT post' in dictionnaire_list_images.keys()):
                    errMsg += ["Multiple CTpost images found"]
                dictionnaire_list_images.update({'CT post':['CTpost', 'electrodes']})
                if not pre_select_1:
                    pre_select_1 = 'CT post'
                if not pre_select_2:
                    pre_select_2 = 'CT post'
                if not pre_select_3:
                    pre_select_3 = 'CT post'
            elif (t.attributes()['modality'] == 'ct') and ('pre' in t.attributes()['acquisition']):
                if ('CT pre' in dictionnaire_list_images.keys()):
                    errMsg += ["Multiple CTpre images found"]
                dictionnaire_list_images.update({'CT pre':['CTpre', 'electrodes']})
            elif (t.attributes()['modality'] == 'ct') and ('postOp' in t.attributes()['acquisition']):
                if ('CT post-op' in dictionnaire_list_images.keys()):
                    errMsg += ["Multiple CTpostOp images found"]
                dictionnaire_list_images.update({'CT post-op':['CTpostOp', 'electrodes']})
            elif (t.attributes()['modality'] == 'pet') and ('pre' in t.attributes()['acquisition']):
                if ('PET pre' in dictionnaire_list_images.keys()):
                    errMsg += ["Multiple PETpre images found"]
                dictionnaire_list_images.update({'PET pre':['PETpre', 'electrodes']})
            elif (t.attributes()['modality'] == 'flair') and ('pre' in t.attributes()['acquisition']):
                if ('FLAIR pre' in dictionnaire_list_images.keys()):
                    errMsg += ["Multiple FLAIRpre images found"]
                dictionnaire_list_images.update({'FLAIR pre':['FLAIRpre', 'electrodes']})
            elif (t.attributes()['modality'] == 'fgatir') and ('pre' in t.attributes()['acquisition']):
                if ('FGATIR pre' in dictionnaire_list_images.keys()):
                    errMsg += ["Multiple FGATIRpre images found"]
                dictionnaire_list_images.update({'FGATIR pre':['FGATIRpre', 'electrodes']})
            elif (t.attributes()['modality'] == 'fmri_epile') and ('pre' in t.attributes()['acquisition']):
                dictionnaire_list_images.update({'fMRI pre' + ' - ' + t.attributes()['subacquisition']:['fMRIpre', 'electrodes']})  # mettre le nom de la subacquisition
            elif t.attributes()['modality'] == 'statistic_data' and ('pre' in t.attributes()['acquisition']):
                dictionnaire_list_images.update({'Statistic Data' + ' - ' + t.attributes()['subacquisition']:['Statisticspre' + t.attributes()['subacquisition'], 'electrodes']})  # mettre le nom de la subacquisition
            elif t.attributes()['modality'] == 'statistic_data' and ('post' in t.attributes()['acquisition']) and not ('postOp' in t.attributes()['acquisition']):
                dictionnaire_list_images.update({'Statistic Data' + ' - ' + t.attributes()['subacquisition']:['Statisticspost' + t.attributes()['subacquisition'], 'electrodes']})  # mettre le nom de la subacquisition
            elif (t.attributes()['modality'] == 'resection'):
                dictionnaire_list_images.update({'Resection':['Resection', 'electrodes']})
            elif (t.attributes()['modality'] == 'freesurfer_atlas') and ('FreesurferAtlaspre' in t.attributes()['acquisition']):
                dictionnaire_list_images.update({'FreeSurfer Atlas (Destrieux)':['FreesurferAtlaspre', 'electrodes']})
                self.vol_destrieux = aims.read(str(t))
            elif (t.attributes()['modality'] == 'freesurfer_atlas') and ('DKT' in t.attributes()['acquisition']):
                dictionnaire_list_images.update({'DKT Atlas':['DKT', 'electrodes']})
                self.vol_dkt = aims.read(str(t))
                na = 'DKT'
            elif (t.attributes()['modality'] == 'freesurfer_atlas') and ('VEP' in t.attributes()['acquisition']):
                dictionnaire_list_images.update({'VEP Atlas':['VEP', 'electrodes']})
                self.vol_dkt = aims.read(str(t))
                na = 'VEP'
            elif (t.attributes()['modality'] == 'freesurfer_atlas') and ('HCP-MMP1' in t.attributes()['acquisition']):
                dictionnaire_list_images.update({'HCP-MMP1 Atlas':['HCP-MMP1', 'electrodes']})
                self.vol_hcp = aims.read(str(t))
                na = 'HCP-MMP1'
            elif (t.attributes()['modality'] == 'freesurfer_atlas') and ('Lausanne2008-33' in t.attributes()['acquisition']):
                dictionnaire_list_images.update({'Lausanne2008-33 Atlas':['Lausanne2008-33', 'electrodes']})
                na = 'Lausanne2008-33'
            elif (t.attributes()['modality'] == 'freesurfer_atlas') and ('Lausanne2008-60' in t.attributes()['acquisition']):
                dictionnaire_list_images.update({'Lausanne2008-60 Atlas':['Lausanne2008-60', 'electrodes']})
                na = 'Lausanne2008-60'
            elif (t.attributes()['modality'] == 'freesurfer_atlas') and ('Lausanne2008-125' in t.attributes()['acquisition']):
                dictionnaire_list_images.update({'Lausanne2008-125 Atlas':['Lausanne2008-125', 'electrodes']})
                na = 'Lausanne2008-125'
            elif (t.attributes()['modality'] == 'freesurfer_atlas') and ('Lausanne2008-250' in t.attributes()['acquisition']):
                dictionnaire_list_images.update({'Lausanne2008-250 Atlas':['Lausanne2008-250', 'electrodes']})
                na = 'Lausanne2008-250'
            elif (t.attributes()['modality'] == 'freesurfer_atlas') and ('Lausanne2008-500' in t.attributes()['acquisition']):
                dictionnaire_list_images.update({'Lausanne2008-500 Atlas':['Lausanne2008-500', 'electrodes']})
                na = 'Lausanne2008-500'

            # Simiplified acquisition name
            nameAcq = t.attributes()['acquisition']
            if not na:
                try:
                    # We try to get the acquisition name without the date (if there is one) : T1pre_2000-01-01 -> T1pre
                    if 'Statistics' in nameAcq:
                        na = nameAcq.split('_')[0] + t.attributes()['subacquisition']
                    else:
                        na = nameAcq.split('_')[0]
                except:
                    if moda == 'Electrode Implantation Coronal Image':
                        na = 'ImplantationCoro'
                    elif moda == 'Electrode Implantation Sagittal Image':
                        na = 'ImplantationSag'
                    else:
                        print "CANNOT find a nameAcq for ", repr(t)
                        na = 'unknown'
            
            # Try setting referential if not defined in the .minf
            if (na == 'T1pre') and (not 'referential' in t.attributes().keys()):
                print "*** DATABASE FIX: Adding the referential "
                # Look for referential file for this volume
                rdi = ReadDiskItem('Referential of Raw T1 MRI', 'Referential', exactType=True, requiredAttributes={'subject':t['subject'], 'center':t['center'], 'acquisition':t['acquisition']})
                volReferential = list(rdi.findValues({}, None, False))
                # If there is a referential available, add its UUID to the volume object attibutes
                if (len(volReferential) == 1):
                    t.setMinf('referential', volReferential[0].uuid())
                    neuroHierarchy.databases.insertDiskItem(t, update=True)
            if ((na == 'preNorm') or (na == 'postOpNorm')) and (not 'referential' in t.attributes().keys()):
                print "*** DATABASE FIX: Adding the referential "
                # Look for referential file for this volume
                rdi = ReadDiskItem('Referential of Raw T1 MRI', 'Referential', exactType=True, requiredAttributes={'subject':t['subject'], 'center':t['center'], 'acquisition':t['acquisition']})
                volReferential = list(rdi.findValues({}, None, False))
                # If there is a referential available, add its UUID to the volume object attibutes
                if (len(volReferential) == 1):
                    t.setMinf('referential', volReferential[0].uuid())
                    neuroHierarchy.databases.insertDiskItem(t, update=True)
            # Load volume in anatomist
            obj = self.loadAndDisplayObject(t, na)
            
            # === REFERENTIALS ===
            if (na == 'T1pre'):
                strVol = 'MRI pre'
                refT1pre = obj.getReferential()
                # Delete existing transformations, otherwise we can't match the T1pre and FreeSurferT1 scanner-based exactly
                self.deleteNormalizedTransf(obj)
                # Save volume center (in mm)
                if t.get('brainCenter'):
                    self.t1preCenter = t.get('brainCenter')
                elif t.get('volume_dimension') and t.get('voxel_size'):
                    volSize = t.get('volume_dimension')
                    voxSize = t.get('voxel_size')
                    self.t1preCenter = [volSize[0]*voxSize[0]/2, volSize[1]*voxSize[1]/2, volSize[2]*voxSize[2]/2]
                else:
                    self.t1preCenter = [128, 128, 128]
                try:
                    tr2sb = self.t1pre2ScannerBased()
                    if tr2sb is not None:
                        self.refConv.setAnatomistTransform("Scanner-based", tr2sb, toRef=True)
                        # Add the AC-centered Scanner-Based (for PTS importation using AC-centered Olivier David method
                        if self.refConv.isRefAvailable('AC-PC'):
                            acInScannerBased = self.refConv.anyRef2AnyRef([0.0, 0.0, 0.0], 'AC-PC', 'Scanner-based')
                            inf = tr2sb.getInfos()
                            rot = inf['rotation_matrix']
                            trans = [inf['translation'][0] - acInScannerBased[0], inf['translation'][1] - acInScannerBased[1], inf['translation'][2] - acInScannerBased[2]]
                            m = aims.Motion(rot[:3] + [trans[0]] + rot[3:6] + [trans[1]] + rot[6:] + [trans[2]] + [0, 0, 0, 1])
                            self.refConv.setTransformMatrix('AC-centered Scanner-Based', m.inverse(), m)
                except Exception, e:
                    print "Cannot load Scanner-based referential from T1 pre MRI : " + repr(e)

                # Load standard transformations (AC-PC, T1pre Scanner-based, BrainVisa Talairach)
                try:
                    self.refConv.loadACPC(t)
                except Exception, e:
                    print "Cannot load AC-PC referential from T1 pre MRI : " + repr(e)
                try:
                    self.refConv.loadTalairach(t)
                except Exception, e:
                    print "Cannot load Talairach referential from T1 pre MRI : " + repr(e)
            
            elif (na == 'FreesurferT1pre'):
                if (na == 'FreesurferT1pre'):
                    strVol = 'MRI FreeSurfer'
                # Delete existing transformations, otherwise we can't match the T1pre and FreeSurferT1 scanner-based exactly
                self.deleteNormalizedTransf(obj)
                # Load all related transformations
                self.loadVolTransformations(t)
                
            elif (na == 'FreesurferAtlaspre') or ('Lausanne2008' in na) or ('DKT' in na) or ('HCP-MMP1' in na) or ('VEP' in na):
                objAtlas.append(obj)
                # Create palette adapted to the volume
                if (na == 'FreesurferAtlaspre'):
                    colors = freesurfer_colormap
                elif (na == 'DKT'):
                    colors = freesurfer_colormap
                elif (na == 'VEP'):
                    colors = vep_colormap
                elif (na == 'HCP-MMP1'):
                    colors = hcp_colormap
                elif (na == 'Lausanne2008-33'):
                    colors = lausanne33_colormap
                elif (na == 'Lausanne2008-60'):
                    colors = lausanne60_colormap
                elif (na == 'Lausanne2008-125'):
                    colors = lausanne125_colormap
                elif (na == 'Lausanne2008-250'):
                    colors = lausanne250_colormap
                elif (na == 'Lausanne2008-500'):
                    colors = lausanne500_colormap
                # Create custom palette
                customPalette = self.a.createPalette(na)
                customPalette.setColors(colors=colors, color_mode='RGB')
                obj.setPalette(customPalette, minVal=0, maxVal=len(colors)/3-1, absoluteMode=True)

            # ===== SURFACES =====
            # For T1 MRI: Load surfaces and MarsAtlas
            if (na == 'T1pre') or (na == 'FreesurferT1pre'):
                # Progress bar
                if thread is not None:
                    thread.progress_text.emit("Loading cortex meshes...")
                    
                # === CORTEX MESHES ===
                # Get the hemisphere meshes for the acquisition : name = na + filename base : for example, if the acquisition is T1pre_2000-01-01 and the file head.gii, we want T1pre-head
                rdi = ReadDiskItem('Hemisphere Mesh', 'Anatomist mesh formats', requiredAttributes={'subject':patient, 'acquisition':nameAcq, 'center':self.currentProtocol})
                hemis = list(rdi._findValues({}, None, False))
                for hh in hemis:
                    self.loadAndDisplayObject(hh, na + '-' + hh.attributes()['side'] + 'Hemi', color=[0.8, 0.7, 0.4, 0.7])
                    dictionnaire_list_images.update({strVol + ' + ' + hh.attributes()['side'] + ' cortex':[na, na + '-' + hh.attributes()['side'] + 'Hemi', 'electrodes']})
                # Add display with both hemispheres
                if (len(hemis) >= 2):
                    dictionnaire_list_images.update({strVol + ' + cortex':[na, na + '-rightHemi', na + '-leftHemi', 'electrodes']})
                    dictionnaire_list_images.update({'Cortex':[na + '-rightHemi', na + '-leftHemi', 'electrodes']})
                    if not pre_select_4:
                        pre_select_4 = 'Cortex'
                    
                # === MARS ATLAS ===
                # Get MarsAtlas textures
                atlas_di = ReadDiskItem('hemisphere marsAtlas parcellation texture', 'aims Texture formats', requiredAttributes={ 'regularized': 'false', 'subject':patient, 'center':self.currentProtocol, 'acquisition':nameAcq })
                atlas_di_list = list(atlas_di._findValues({}, None, False))
                # If MarsAtlas found
                if len(atlas_di_list) > 0:
                    # Progress bar
                    if thread is not None:
                        thread.progress_text.emit("Loading MarsAtlas parcels...")
                    # Find white meshes
                    wm_di = ReadDiskItem('Hemisphere White Mesh', 'aims mesh formats', requiredAttributes={'subject':patient, 'center':self.currentProtocol })
                    # For each MarsAtlas texture
                    for atl in atlas_di_list:
                        # Look for corresponding white mesh
                        wm_side = wm_di.findValue(atl)
                        # Load object
                        self.loadAndDisplayObject(wm_side, na + '-' + atl.attributes()['side'] + 'MARSATLAS', texture_item=atl, palette='MarsAtlas', color=[0.8, 0.7, 0.4, 0.7])
                        dictionnaire_list_images.update({strVol + ' + ' + atl.attributes()['side'] + ' MarsAtlas':[na, na + '-' + atl.attributes()['side'] + 'MARSATLAS', 'electrodes']})

                # === HEAD MESH ===
                # Get head mesh for the acquisition
                rdi = ReadDiskItem('Head Mesh', 'Anatomist mesh formats', requiredAttributes={'subject':patient, 'acquisition':nameAcq, 'center':self.currentProtocol})
                head = list(rdi._findValues({}, None, False))
                # Only if the hemispheres and the head are available
                if (len(hemis) >= 2) and (len(head) > 0):
                    if thread is not None:
                        thread.progress_text.emit("Loading head mesh...")
                    self.loadAndDisplayObject(head[0], na + '-' + 'head', color=[0.0, 0.0, 0.8, 0.3])
                    dictionnaire_list_images.update({strVol + ' + cortex + head':[na, na + '-rightHemi', na + '-leftHemi', na + '-head', 'electrodes']})
                    dictionnaire_list_images.update({'Cortex + head':[na + '-rightHemi', na + '-leftHemi', na + '-head', 'electrodes']})
                    
                # === AMYGDALA + HIPPOCAMPUS ===
                # Get the amygdala+hippo meshes for the acquisition
                rdi = ReadDiskItem('Mesh', 'Anatomist mesh formats', requiredAttributes={'subject':patient, 'acquisition':nameAcq, 'center':self.currentProtocol})
                allMeshes = list(rdi._findValues({}, None, False))
                amygHippoMesh = []
                for hh in allMeshes:
                    hFilename = hh.fileNames()[0]
                    if ("Amygdala" in hFilename): 
                        self.loadAndDisplayObject(hh, hFilename, color=[1, 0, 0, 0.7])
                        amygHippoMesh.append(hFilename)
                    elif ("anteroHippocampus" in hFilename):
                        self.loadAndDisplayObject(hh, hFilename, color=[1, 1, 0, 0.7])
                        amygHippoMesh.append(hFilename)
                    elif ("posteroHippocampus" in hFilename):
                        self.loadAndDisplayObject(hh, hFilename, color=[1, 0.7, 0, 0.7])
                        amygHippoMesh.append(hFilename)
                # Add display with all the structures
                if amygHippoMesh:
                    dictionnaire_list_images.update({strVol + ' + hippocampus + amygdala':[na, 'electrodes'] + amygHippoMesh})

        # ===== FREESURFER TRANSFORMATION =====
        # FreeSurfer atlas was resliced using T1pre: force T1pre referential
        if objAtlas and refT1pre:
            for obj in objAtlas:
                obj.assignReferential(refT1pre)
        # Update interface
        if isGui:
            # Update list of available items in the combo boxes
            self.windowContent = dictionnaire_list_images
            self.updateComboboxes(pre_select_1, pre_select_2, pre_select_3, pre_select_4)
            # Display referential informations
            self.setWindowsReferential()
            if thread is not None:
                thread.progress_text.emit("Loading electrodes...")
        # Load electrodes
        self.loadElectrodes(patient, self.currentProtocol, isGui)
        # Update display
        if isGui:
            self.refreshReferentials()
            # Center view on AC
            self.centerCursor()
        return errMsg


    # Chargement d'un objet (MRI, mesh...) dans Anatomist et mise à jour de l'affichage
    def loadAndDisplayObject(self, diskitem, name = None, color=None, palette=None, texture_item = None):

        if name is None:
          return
        
        #Already exists ! Remove it.
        if name in self.dispObj:
            self.a.removeObjects([self.dispObj[name],], self.wins) # Remove from windows
            self.a.deleteObjects(self.dispObj[name]) # CURRENT
            del self.dispObj[name]
            del self.diskItems[name]

        obj = self.a.loadObject(diskitem)

        if 'ColorPalette' in diskitem.attributes():
            obj.setPalette(palette = diskitem.attributes()['ColorPalette'])
        elif palette is not None and texture_item is None:
            obj.setPalette(palette = palette)
        if texture_item is not None:
            texture = self.a.loadObject(texture_item)
            if palette is not None:
                texture.setPalette(palette = palette)
            textured_mesh = self.a.fusionObjects((obj,texture),method = 'FusionTexSurfMethod')
            #we need to keep the texture and the mesh as well as the fusion of both
            self.objtokeep[name + '_mesh'] = obj
            self.objtokeep[name + '_texture'] = texture
            obj = textured_mesh
        
        # Store the object
        self.dispObj[name] = obj
        self.diskItems[name] = diskitem
        # If this is a volume, smooth it :
        try:
            self.a.execute('TexturingParams', objects=[obj], filtering='linear')
        except:
            pass
        if color is not None:
            self.a.setMaterial(obj, diffuse=color)
        return obj


    def setWindowsReferential(self, ref=None):
        """ Get all available referentials from anatomist and tries to match identical referentials from SPM"""
        # If the T1pre image is already loaded
        if ref is None:
            if self.preReferential():
                # Assign the T1pre native referential to the windows
                self.currentWindowRef = self.preReferential()
                self.a.assignReferential(self.currentWindowRef, self.wins)
        else:
            self.currentWindowRef = ref
            self.a.assignReferential(ref, self.wins)


    # Get the click events
    def clickHandler(self, eventName, params):
        # If a real mouse click: save position
        if params:
            self.lastClickedPos = params['position']
            self.lastClickedRef = params['window'].getReferential()
        # No T1pre
        if not 'T1pre' in self.dispObj:
            return
        # Get coordinates
        pT1Pre = self.positionPreRef()
        if not pT1Pre:
            return
        # Labels
        if (self.coordsDisplayRef == "Destrieux") or (self.coordsDisplayRef == "DKT") or (self.coordsDisplayRef == "HCP-MMP1"):
            info_image = self.diskItems['T1pre'].attributes()
            pos_SB = [round(pT1Pre[i]/info_image['voxel_size'][i]) for i in range(3)]
            pos_fs = pos_SB
            if (self.coordsDisplayRef == "Destrieux"):
                fsIndex = self.vol_destrieux.value(pos_fs[0],pos_fs[1],pos_fs[2])
                self.positionLabel.setText(labels['Destrieux'][fsIndex])
            elif (self.coordsDisplayRef == "DKT"):
                fsIndex = self.vol_dkt.value(pos_fs[0],pos_fs[1],pos_fs[2])
                self.positionLabel.setText(labels['DKT'][fsIndex])
            elif (self.coordsDisplayRef == "HCP-MMP1"):
                fsIndex = self.vol_hcp.value(pos_fs[0],pos_fs[1],pos_fs[2])
                self.positionLabel.setText(labels['HCP-MMP1'][fsIndex])
        # Coordinates
        else:
            if self.coordsDisplayRef == 'Natif':
                coords = pT1Pre
            elif self.coordsDisplayRef == 'Scanner-based':
                infos = self.t1pre2ScannerBased().getInfos()
                rot = infos['rotation_matrix']
                trans = infos['translation']
                m = aims.Motion(rot[:3]+[trans[0]]+rot[3:6]+[trans[1]]+rot[6:]+[trans[2]]+[0,0,0,1])
                coords = m.transform(pT1Pre)
            else:
                try:
                    coords = self.refConv.real2AnyRef(pT1Pre, self.coordsDisplayRef)
                except:
                    coords = [0.0,0.0,0.0]
            if coords is None:
                coords = [0.0,0.0,0.0]
            self.positionLabel.setText("%.2f, %.2f, %.2f" % tuple(coords))


    def preReferential(self):
        if 'T1pre' in self.dispObj:
            return self.dispObj['T1pre'].getReferential()
        else:
            return None
  
  
    def positionPreRef(self):
        # MAY 2018: Not using linkCursorLastClickedPosition anymore because of a bug in BrainVISA 4.6
        # that causes linkCursorLastClickedPosition() to return coordinates in the referential currently
        # selected in the window, instead of the central referential of anatomist; and when called with
        # an argument, the function returns something completely random.

        #return list(self.a.linkCursorLastClickedPosition(self.preReferential()).items())

        # Get last click position and referentials
        pos = self.lastClickedPos
        refSrc = self.lastClickedRef
        refDest = self.preReferential()
        # Convert referential to T1pre_native
        if pos and refSrc and refDest and (refSrc != refDest):
            trm = self.a.getTransformation(refSrc, refDest)
            if trm:
                pos = trm.motion().transform([pos[0], pos[1], pos[2]])
        if pos:
            return [pos[0], pos[1], pos[2]]
        else:
            return None


    def t1pre2ScannerBased(self):
        """ Returns a Transformation object that transforms T1pre referential to T1pre Scanner-Based referential """
        if self.t1pre2ScannerBasedTransform is not None:
            if "dead" not in self.t1pre2ScannerBasedTransform.getInfos():
                return self.t1pre2ScannerBasedTransform
        rdi = ReadDiskItem('Transformation to Scanner Based Referential', 'Transformation matrix', exactType=True,\
                           requiredAttributes={'modality':'t1mri', 'subject':self.brainvisaPatientAttributes['subject'], 'center':self.brainvisaPatientAttributes['center']})
        allTransf = list (rdi._findValues( {}, None, False ) )
        for trsf in allTransf:
            if trsf.attributes()['acquisition'].startswith(u'T1pre'):
                srcrDiskItem = self.transfoManager.referential( trsf.attributes()['source_referential'] )
                srcr = self.a.createReferential(srcrDiskItem)
                dstrDiskItem = self.transfoManager.referential(trsf.attributes()['destination_referential'])
                self.t1pre2ScannerBasedId = trsf.attributes()['destination_referential']
                dstr = self.a.createReferential(dstrDiskItem)
                self.t1pre2ScannerBasedTransform = self.a.loadTransformation(trsf.fullPath(), srcr, dstr)
                return self.t1pre2ScannerBasedTransform
        return None
    

    def loadVolTransformations(self, t):
        allTransf = []
        # Find transformations
        rdi = ReadDiskItem('Transformation to Scanner Based Referential', 'Transformation matrix', exactType=True, requiredAttributes={'subject':t['subject'], 'center':t['center'], 'acquisition':t['acquisition']})
        allTransf += list(rdi.findValues({}, None, False))
        rdi = ReadDiskItem('Transform Raw T1 MRI to another image', 'Transformation matrix', exactType=True, requiredAttributes={'subject':t['subject'], 'center':t['center'], 'acquisition':t['acquisition']})
        allTransf += list(rdi.findValues({}, None, False))
        # Load all transformations
        loadedTrm = []
        for trm in allTransf:
            # Skip incomplete transformations
            if ('source_referential' not in trm.attributes().keys()) or ('destination_referential' not in trm.attributes().keys()):
                continue
            # Get source referential
            srcrDiskItem = self.transfoManager.referential(trm.attributes()['source_referential'])
            srcr = self.a.createReferential(srcrDiskItem)
            # Get destination referential
            dstrDiskItem = self.transfoManager.referential(trm.attributes()['destination_referential'])
            dstr = self.a.createReferential(dstrDiskItem)
            # Load transformation
            loadedTrm += [self.a.loadTransformation(trm.fullPath(), srcr, dstr)]
        return loadedTrm
    
    def deleteNormalizedTransf(self, obj):
        # Delete existing normalized transformations, otherwise we can't match the T1pre and FreeSurferT1 scanner-based exactly
        for trm in self.a.getTransformations():
            if trm.getInfos()\
            and ('source_referential' in trm.getInfos().keys()) \
            and (trm.getInfos()['source_referential'] == obj.getReferential().getInfo()['uuid']):
                # Remove only normalized transformations
                for ref in self.a.getReferentials():
                    if ('destination_referential' in trm.getInfos().keys()) \
                    and (ref.getInfo()['uuid'] == trm.getInfos()['destination_referential']) \
                    and ref.getInfo() and ('name' in ref.getInfo().keys()) \
                    and ((ref.getInfo()['name'] == 'Talairach-MNI template-SPM') or (ref.getInfo()['name'] == 'Talairach-AC/PC-Anatomist')):
                        self.a.deleteElements(trm)

    def mniReferentialId(self):
        return aims.StandardReferentials.mniTemplateReferentialID()

    def mniReferential(self):
        return self.a.mniTemplateRef

    def refreshReferentials(self):
        curr = str(self.referentialCombo.currentText())
        self.referentialCombo.clear()
        # Add items from the referential converter
        refs = self.refConv.availableReferentials().keys() + ['Natif',]
        self.referentialCombo.addItems(refs)
        # Add available atlases
        if ('FreesurferAtlaspre' in self.dispObj.keys()):
            self.referentialCombo.addItems(["Destrieux"])
            self.referentialCombo.setCurrentIndex( self.referentialCombo.count() - 1 )
        if ('DKT' in self.dispObj.keys()):
            self.referentialCombo.addItems(["DKT"])
        if ('HCP-MMP1' in self.dispObj.keys()):
            self.referentialCombo.addItems(["HCP-MMP1"])

    def updateCoordsDisplay(self, text):
        self.coordsDisplayRef = str(text)
        self.clickHandler(None, None)

    def updateDispMode(self, index):
        """ Update the display mode of all electrodes """
        mode = 'sphere'
        params = {}
        isbipole = False
        if index == 0:
            mode = 'real'
            self.colorConfigButton.setEnabled(False)
        elif index == 1:
            params = {'diameter':1.0}
            self.colorConfigButton.setEnabled(False)
        elif index == 2:
            params = {'diameter':2.0}
            self.colorConfigButton.setEnabled(False)
        elif index == 3:
            params = {'diameter':5.0}
            self.colorConfigButton.setEnabled(False)
        elif index == 4:
            mode = 'off'
            self.colorConfigButton.setEnabled(False)
        elif index == 5:
            mode = 'bipole'
            isbipole = True
            #is there a json file about the seeg results in the database.
            rdi_seeglabel = ReadDiskItem('Electrodes SEEG Labels','Electrode sEEG Label Format',requiredAttributes={'subject':self.brainvisaPatientAttributes['subject'], 'center':self.brainvisaPatientAttributes['center']})
            di_seeglabel = list(rdi_seeglabel.findValues({},None,False))
            #if not, ask for one
            if len(di_seeglabel) == 0:
                load_new_file = True               
            else:
                #ask if you want to replace the loaded data
                rep = QtGui.QMessageBox.warning(self, u'Use database or Import?', u"Use sEEG stim result from the database ? (if no, it will ask for a file)", QtGui.QMessageBox.Yes | QtGui.QMessageBox.No | QtGui.QMessageBox.Cancel, QtGui.QMessageBox.Cancel)
                if rep == QtGui.QMessageBox.Yes:
                    load_new_file = False 
                elif rep == QtGui.QMessageBox.No:    
                    load_new_file = True
                elif rep == QtGui.QMessageBox.Cancel:
                    return
              
            if load_new_file:
                wdi_seeg = WriteDiskItem('Electrodes SEEG Labels','Electrode sEEG Label Format')
                di_seeg = wdi_seeg.findValue({'subject':self.brainvisaPatientAttributes['subject'], 'center':self.brainvisaPatientAttributes['center']} )
                
                fichierseegLabel =  QtGui.QFileDialog.getOpenFileName(self, "Select a file containing seeg labels: ", "", "(*.xlsx *.csv *.json)")
                if os.path.basename(str(fichierseegLabel)).split('.')[-1] == 'json':
                    fin = open(str(fichierseegLabel),'rb')
                    new_label = json.loads(fin.read())
                    fin.close()
                    
                    try:
                        os.mkdir(os.path.dirname(str(di_seeg)))
                    except:
                        pass
                    fout = open(str(di_seeg),'w')
                    fout.write(json.dumps({'title':new_label['title'],'contacts':new_label['contacts']}))
                    fout.close()                
                    neuroHierarchy.databases.insertDiskItem(di_seeg, update=True )
                    
                elif os.path.basename(str(fichierseegLabel)).split('.')[-1] == 'xlsx':
                    contact_label_class = generate_contact_colors()
                    inter_label = contact_label_class.from_excel_files(str(fichierseegLabel))
                    #write the json and include it in the database
                    try:
                        os.mkdir(os.path.dirname(str(di_seeg)))
                    except:
                        pass
                    fout = open(str(di_seeg),'w')
                    fout.write(json.dumps({'title':inter_label[0],'contacts':inter_label[1]}))
                    fout.close() 
                    neuroHierarchy.databases.insertDiskItem(di_seeg, update=True )
                    new_label = {'title':inter_label[0],'contacts':inter_label[1]}
              
            else:
                fin = open(str(di_seeglabel[0]),'rb')
                new_label = json.loads(fin.read())
                fin.close()
    
            bipole_label_sorted = sorted(new_label['contacts'].keys(),key=natural_keys)
            plotsT1preRef = self.getPlotsT1preRef()
            info_plotsT1Ref= []
            for k,v in plotsT1preRef.iteritems():
                plot_name_split = k.split('-$&_&$-')
                info_plotsT1Ref.append((plot_name_split[0]+plot_name_split[1][4:].zfill(2),v.list()))
                
            plotsT1Ref_sorted = dict(sorted(info_plotsT1Ref, key=lambda plot_number: plot_number[0]))
              
            #remplir self.bipole
            for i_bipole in bipole_label_sorted:
                #get the name of both contact from the bipole name
                try:
                    pos_bipol = (numpy.array(plotsT1Ref_sorted[i_bipole.split()[0].title()]) + numpy.array(plotsT1Ref_sorted[i_bipole.split()[2].title()]))/2
                except:
                    raise Exception("problem plotsT1Ref")
                entry_bipole = numpy.array(plotsT1Ref_sorted[i_bipole.split()[0].title()])
                #il faut un orient vector
                self.addBipole(i_bipole,pos_bipol,self.preReferential().uuid(),entry_bipole) #on rajoute pour finir un vecteur

            self.colorConfigButton.setEnabled(True)

            contact_labels_seeg = {}
            params.update({'contact_colors':contact_labels_seeg})
            
            self.bipoleLabels = new_label
            self.configureColors()
    
        if mode == self.dispMode and params == self.dispParams:
            return
        self.updateElectrodeMeshes(clear=True)
        for elec in self.electrodes:
            elec['elecModel'].setDisplayMode(mode, params)
            elec['elecModel'].updateDisplay()
            elec['elecModel'].setDisplayReferential(elec['ref'])
        self.dispMode = mode
        self.dispParams = params
        self.updateElectrodeMeshes(bipole = isbipole)
        self.updateAllWindows()
        # Update the contact list meshes of the current electrode (so they can be selected)
        self.electrodeSelect(self.electrodeList.currentRow())
    

    def updateElectrodeMeshes(self, clear=False, bipole=False):
        if clear:
            if 'electrodes' not in self.dispObj.keys():
                return
            self.a.removeObjects(self.dispObj['electrodes'], self.wins)
            self.dispObj['electrodes'] = []
            return
        if not bipole:    
            self.dispObj['electrodes'] = [mesh for elec in self.electrodes for mesh in elec['elecModel'].getAnatomistObjects() if mesh is not None]
            self.setElectrodeMeshesNames()          
        elif bipole:
            self.dispObj['electrodes'] = [mesh for elec in self.bipoles for mesh in elec['elecModel'].getAnatomistObjects() if mesh is not None]
            self.setBipoleMeshesNames()


    # Add an electrode from a template
    def addElectrode(self, name=None, model=None, target=None, entry=None, refId=None, isUpdate=True, isGui=True):
        if name is None:
            name = str(self.nameEdit.text())
            self.isModified = True
        if model is None:
            model = str(self.typeComboBox.currentText())
        isTransparent = False
        if target is None:
            if self.t1preCenter:
                target = self.t1preCenter
            else:
                target = [0,0,0]
            isTransparent = True
        if entry is None:
            if self.t1preCenter:
                entry = [self.t1preCenter[0], self.t1preCenter[1], self.t1preCenter[2] - 0.1]
            else:
                entry = [0,0,-0.1]
            isTransparent = True
        if refId is not None:
            if self.preReferential().uuid() != refId:
                print "ERROR : the electrode is not defined in reference to the T1 pre image (%s) but in reference to %s !\nGoing on anyway..."%(self.preReferential().uuid(), refId)
        # Create electrode meshes
        for el in self.electrodes:
            if name == el['name']:
                name = self.findFreeElectrodeName()
        if str(model) not in self.elecModelListByName:
            print("ERROR: Could not find electrode model named " + str(model) + ". Ignoring electrode. Please check that the model is installed in IntrAnat database !")
            #QtGui.QMessageBox.warning(self, u'Error', u"Could not load electrode model %s. Please check it is installed in IntrAnat database !"%str(model))
            return
        (newRef, transf, elecModel) = createElectrode(target, entry, self.preReferential(), self.a,\
                                                      model = self.elecModelListByName[str(model)].fullPath(), dispMode = self.dispMode, dispParams = self.dispParams, isTransparent = isTransparent)
        # Reference electrodes in application
        self.electrodes.append({'name': name, 'ref':newRef, 'transf':transf, 'elecModel':elecModel,\
                                'target':target, 'entry':entry, 'model':model})
        # Add electrode name in GUI
        self.electrodeList.addItem(name)
        self.electrodeList.setCurrentRow(self.electrodeList.count() - 1)
        # Redraw electrodes
        if isUpdate and isGui:
            self.updateElectrodeMeshes()
            self.updateAllWindows(False)
        

    def addBipole(self, name=None, positionbip=[0,0,0], refId = None,entry_bipole = None):
        if name is None:
            print("error, bipole must have a name")
            return
        if refId is None:
            print("error, bipole has to be assigned to a referential")
            return
        if refId is not None:
            if self.preReferential().uuid() != refId:
                print "ERROR : the electrode is not defined in reference to the T1 pre image (%s) but in reference to %s !\nGoing on anyway..."%(self.preReferential().uuid(), refId)     
        #check is name is already taken
        if len(self.bipoles)>0:
            for bip in self.bipoles:
                if name == bip['name']:
                    print("error, bipole name already taken")
             
        #je met un objet elecmodel ici ou juste le mesh ? un objet elecmodel de 1 contact ?
        rdiEM = ReadDiskItem('Electrode Model', 'Electrode Model format', requiredAttributes={'center':self.currentProtocol})
        listEM = list(rdiEM.findValues({},None,False))
        matches = filter((lambda x: u"bipole" in str(x)), listEM)
            
        (newRef, transf, elecModel) = createBipole(positionbip.tolist(), entry_bipole.tolist(), self.preReferential(), self.a, model = matches[0], dispMode = 'bipole', dispParams = None)
        self.bipoles.append({'name': name, 'ref':newRef, 'transf':transf, 'target':positionbip.tolist(), 'entry': entry_bipole.tolist(), 'elecModel':elecModel}) #to change target and entry
        self.isModified = True
        

    # Setting names on meshes to get a nice tooltip for each mesh
    def setElectrodeMeshesNames(self, electrode = None):
        if electrode is None:
            electrodes = self.electrodes
        else:
            electrodes = [electrode,]
        for el in electrodes:
            for name, element in el['elecModel'].getDisplayed().iteritems():
                if element['mesh'] is not None:
                    if element['type'] == 'Plot':
                        element['mesh'].setName(name.replace('Plot', el['name']))
                    else:
                        element['mesh'].setName(el['name'])

          
    def setBipoleMeshesNames(self, bipole = None):
        if bipole is None:
            bipoles = self.bipoles
        else:
            bipoles = [bipole,]
        for bp in bipoles:
            for name, element in bp['elecModel'].getDisplayed().iteritems():
                element['mesh'].setName(bp['name'])  

    def findFreeElectrodeName(self):
        if len(self.electrodes) == 0:
            return 'A'
        n = [el['name'] for el in self.electrodes]
        newName = n[-1]
        firstletter = newName[::-1][-1]
        if firstletter in string.ascii_uppercase:
            newName = string.uppercase[(string.uppercase.find(firstletter)+1)%len(string.uppercase)] + newName[1:]
        if newName in n:
            while newName in n:
                newName = newName+'_'
        return newName

    def currentElectrode(self):
        idx = self.electrodeList.currentRow()
        if idx < 0:
          return None
        return self.electrodes[idx]

    def removeElectrode(self):
        """Remove an electrode (and all contacts)"""
        elec = self.currentElectrode()
        if elec is None:
            return
        idx = self.electrodes.index(elec)
        # Remove meshes from Anatomist
        self.updateElectrodeMeshes(clear=True)
        elec['elecModel'].clearDisplay()
        item = self.electrodeList.takeItem(idx)
        del item
        del self.electrodes[idx]
        self.updateElectrodeMeshes()
        self.updateAllWindows(False)
        self.isModified = True

    def updateElectrodeModel(self, model):
        # Get current electrode
        elec = self.currentElectrode()
        # Reload list of electrode models
        if str(model) == '[Reload electrode models]':
            # Wait cursor
            QtGui.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
            # Reload BrainVISA share database
            for database in neuroHierarchy.databases.iterDatabases():
                if '/share/brainvisa-share-4.' in database.name:
                    database.clear()
                    database.update()
            # Block callbacks
            self.typeComboBox.blockSignals(True)
            # Reload electrode models
            self.loadElectrodeModels()
            # Select again the previous model
            if elec:
                self.typeComboBox.setCurrentIndex(self.typeComboBox.findText(elec['model']))
            # Restore callbacks
            self.typeComboBox.blockSignals(False)
            # Restore cursor
            QtGui.QApplication.restoreOverrideCursor()
            return
        
        # No electrode selected: exit
        if elec is None:
            return
        if str(elec['model']) == str(model):
            return
        
        # Update electrode model for the selected electrode
        self.updateElectrodeMeshes(clear=True)
        elec['elecModel'].clearDisplay()
        del elec['elecModel']
        (newRef, transf, elecModel) = createElectrode(elec['target'], elec['entry'], self.preReferential(), self.a,\
                                                      model = self.elecModelListByName[str(model)].fullPath(), dispMode = self.dispMode, dispParams = self.dispParams)
        elec['elecModel'] = elecModel
        elec['model']=str(model)
        self.electrodeSelect(self.electrodeList.currentRow())
        self.updateElectrodeMeshes()
        self.updateAllWindows()
        self.isModified = True


    def updateEntry(self, e=None):
        """ Updates the current electrode entry point from the cursor position"""
        el = self.currentElectrode()
        pos = self.positionPreRef()
        # Just ignore if new entry is identical to target
        if pos == el['target']:
            return
        # Update electrode meshes
        meshes = el['elecModel'].getAnatomistObjects()
        (newRef, transf) = moveElectrode(el['target'], pos, self.preReferential(), el['ref'], self.a, meshes)
        # Save modifications
        el['entry'] = pos
        el['transf'] = transf
        self.isModified = True
        # If target was not set and now electrode is fully visible: make meshes visible
        if self.t1preCenter and (el['target'] != self.t1preCenter) and (el['entry'] != [self.t1preCenter[0], self.t1preCenter[1], self.t1preCenter[2] - 0.1]):
            for m in meshes:
                mat = m.getInternalRep().GetMaterial().genericDescription()['diffuse']
                m.setMaterial(diffuse=[mat[0], mat[1], mat[2], 1])
                
                
    def updateTarget(self, t=None):
        """ Updates the current electrode target point from the cursor position"""
        el = self.currentElectrode()
        if el is None:
            return
        pos = self.positionPreRef()
        # Just ignore if new target is identical to entry
        if pos == el['entry']:
            return
        # Update electrode meshes
        meshes = el['elecModel'].getAnatomistObjects()
        (newRef, transf) = moveElectrode(pos, el['entry'], self.preReferential(), el['ref'], self.a, meshes)
        # Save modifications
        el['target'] = pos
        el['transf'] = transf
        self.isModified = True
        # If target was not set and now electrode is fully visible: make meshes visible
        if self.t1preCenter and (el['target'] != self.t1preCenter) and (el['entry'] != [self.t1preCenter[0], self.t1preCenter[1], self.t1preCenter[2] - 0.1]):
            for m in meshes:
                mat = m.getInfos()["material"]["diffuse"]
                m.setMaterial(diffuse=[mat[0], mat[1], mat[2], 1])
                
            
    def editElectrodeName(self):
        """Update the electrode name of the selected contact"""
        name = str(self.nameEdit.text())
        idx = self.electrodeList.currentRow()
        # No electrode selected
        if idx == -1:
            return

        sameNameItems = self.electrodeList.findItems(name, QtCore.Qt.MatchFixedString)
        if len(sameNameItems) != 0:
            if sameNameItems[0] == self.electrodeList.item(idx): # Name was not changed
                return
            else:
                QtGui.QMessageBox.warning(self, u'Error', u"The name %s is already used by another electrode. Choose an other one !"%name)
                self.nameEdit.setText(self.electrodeList.item(idx).text())
        self.electrodes[idx]['name'] = name
        self.setElectrodeMeshesNames(self.electrodes[idx])
        self.electrodeList.currentItem().setText(name)
        self.electrodeSelect(idx) # updates the contacts' display
        self.isModified = True


    def electrodeSelect(self, idx):
        """Electrode/contact selection changed"""
        el = self.electrodes[idx]
        self.nameEdit.setText(el['name'])
        self.typeComboBox.setCurrentIndex(self.typeComboBox.findText(el['model']))
        # Select the contacts in anatomist
        g = self.a.getDefaultWindowsGroup()
        g.setSelection(el['elecModel'].plotMeshes())
        self.currentContacts = dict((name.replace('Plot', el['name']), element['mesh']) for name, element in el['elecModel'].getDisplayed().iteritems() if element['type'] == 'Plot' and element['mesh'] is not None)
        self.currentElectrodes = [el,]
        #Update contact list
        self.contactList.clear()
        self.contactList.addItems(sorted(self.currentContacts.keys(), key=natural_keys))
        # Update electrode referential (if selected)
        if (self.electrodeRefCheck.checkState() == QtCore.Qt.Checked):
            self.updateElectrodeView()

    def electrodeGo(self, idx = None, electrode = None):
        if idx is not None:
            if type(idx) == int:
                electrode = self.electrodes[idx]
            else:
                electrode = self.currentElectrodes[0]
        elif electrode is None:
            return
        self.contactGo(el = electrode)

    def contactSelect(self, item):
        g = self.a.getDefaultWindowsGroup()
        g.setSelection([self.currentContacts[str(item.text())],])

    def contactGo(self, item = None, el=None):
        # Put the cursor on the contact in anatomist
        try:
            if el is None:
                el = self.currentElectrodes[0]
            if item is not None:
                name = str(item.text()).replace(el['name'], 'Plot')
            else:
                name = 'Plot1'
            xyz = getPlotsCenters(el['elecModel'])[name]
            # setLinkedCursor uses window referential : must apply transform before setting the position
            # if self.currentWindowRef == self.preReferential():
            if self.wins[0].getReferential() == self.preReferential():
                xyz = el['transf'].transform(xyz)
                #self.a.execute('LinkedCursor', window=self.wins[0], position=xyz)
            elif  self.wins[0].getReferential() == el['ref']:
                # Nothing to change
                pass
            else:
                transf = self.a.getTransformation(el['ref'], self.wins[0].getReferential())
                xyz = transf.transform(xyz)
            self.wins[0].moveLinkedCursor(xyz)
        except:
            print "Error moving the cursor to the contact"

    # Center the cursor on all the anatomist views
    def centerCursor(self):
        if (self.t1preCenter) and (self.currentWindowRef == self.preReferential()):
            self.wins[0].moveLinkedCursor(self.t1preCenter)

    def getElectrodeTemplates(self):
        """Returns a list of Electrode objects, one for each available electrode template model available in the DB"""
        if self.electrodeTemplateStubs == []:
            self.electrodeTemplateStubs = dict([(n, ElectrodeModel(modelPath=model.fullPath(), dispMode='off')) for n, model in self.elecModelListByName.iteritems()])
        return self.electrodeTemplateStubs


    def fitElecModelToLength(self, target, entry):
        """ Tries to find a match between the length of the electrode and a model, but prefering uniform electrodes (no variable spacing between contacts"""
        length = linalg.norm(array(entry) - array(target))
    
        models = self.getElectrodeTemplates()
        uniformModelsLength = {}
        for n, model in models.iteritems():
            interPlotsM = []
            plotsM = sorted([[int(''.join(e for e in namePlot if e.isdigit())), ] + content['center'] for namePlot, content in model.getPlots().iteritems()], key=lambda p:p[0])
    
            for p in range(len(plotsM) - 1):
                interPlotsM.append(math.sqrt((plotsM[p][1] - plotsM[p + 1][1]) ** 2 + (plotsM[p][2] - plotsM[p + 1][2]) ** 2 + (plotsM[p][3] - plotsM[p + 1][3]) ** 2))
            lengthM = sum(interPlotsM)
            if len(set(interPlotsM)) == 1:  # All intervals are identical
                uniformModelsLength[n] = lengthM
        if length >= max(uniformModels.values()):
            return [m for m, l in uniformModels.iteritems() if l == max(uniformModels.values())][0]
        else:
            largerModels = dict([(m, l) for m, l in uniformModels.iteritems() if l >= length])
            return [m for m, l in largerModels.iteritems() if l == min(largerModels.values())][0]


    def fitElecModelToPlots(self, plots):
        """ Tries to find a match between a list of plots [[numPlot, x, y, z],[...], ...] and available electrode templates.
          Return None if it fail, [ModelName, [targetX, targetY, targetZ], [entryX, entryY, entryZ]] if it works
        """
        plots = sorted(plots, key=lambda p:int(p[0]))
        # Compute inter-contact distances
        interPlots=[]
        for p in range(len(plots)-1):
            interPlots.append(math.sqrt((plots[p][1]-plots[p+1][1])**2 + (plots[p][2]-plots[p+1][2])**2 + (plots[p][3]-plots[p+1][3])**2))
        # Find electrode models with a similar number of plots
        models = dict([(n,tpl) for n,tpl in self.getElectrodeTemplates().iteritems() if tpl.countPlots() == len(plots)]) # identical number
        if len(models) == 0:
            print "Cannot find a template with the right number of contacts (%i), trying longer ones"%len(plots)
            nbPlotsModels = dict([(n, tpl.countPlots()) for n,tpl in self.getElectrodeTemplates().iteritems() if tpl.countPlots() > len(plots)])
            if len(nbPlotsModels) == 0:
                print "Cannot find a template with enough contacts for this electrode ! (available : %s)\nWill match with the largest available electrodes !\n THIS WILL LOSE SOME PLOTS"%repr(nbPlotsModels.values())
                nbPlotsModels = dict([(n, tpl.countPlots()) for n,tpl in self.getElectrodeTemplates().iteritems()])
                models = dict([(n,tpl) for n,tpl in self.getElectrodeTemplates().iteritems() if tpl.countPlots() == max(nbPlotsModels.values())])
            else:
                models = dict([(n,tpl) for n,tpl in self.getElectrodeTemplates().iteritems() if tpl.countPlots() == min(nbPlotsModels.values())])

        # Now, check the interPlots distance in the model and compare with the template
        distanceFromModel = {}
        for n,model in models.iteritems():
            interPlotsM=[]
            plotsM = sorted([[int(''.join(e for e in namePlot if e.isdigit())), ] + content['center'] for namePlot, content in model.getPlots().iteritems()], key=lambda p:p[0])
            
            for p in range(len(plotsM)-1):
                interPlotsM.append(math.sqrt((plotsM[p][1]-plotsM[p+1][1])**2 + (plotsM[p][2]-plotsM[p+1][2])**2 + (plotsM[p][3]-plotsM[p+1][3])**2))
            distanceFromModel[n] = sum([(interPlotsM[i]-interPlots[i])**2 for i in range(min([len(interPlotsM), len(interPlots)]))])

        # Choose the model with the smallest distance, or the first one if a few are equally good, reject if sum is more than 1mm per contact
        minDist = min(distanceFromModel.values())
        if minDist > 1.0*len(plots):
            print "Cannot match a template : minimum distance found is %f mm, should be %f or less !"%(minDist, 1.0*len(plots))
            return None
        # We have a winner !
        goodModel = [m for m in distanceFromModel if distanceFromModel[m] == minDist][0]
        print "Found model %s, match distance = %f"%(goodModel, minDist)
        entry = plots[-1][1:]
        # Target is not the center of the first contact, but its end ! -> coords - lengthOfPlot/2*vector(plot0->plotLast)
        target = array(plots[0][1:]) - (array(plotsarr[-1][1:])-array(plots[0][1:]))/linalg.norm(array(plots[-1][1:])-array(plots[0][1:])) * 0.5 * self.getElectrodeTemplates()[goodModel].getCylinder('Plot1')['length']
        return [goodModel, target.tolist(), entry]


    def comboMessageBox(self, text, choices):
        """ Displays a message box with a choice (combobox), Ok and Cancel button
            Returns the selected value or None if it was cancelled"""
        if choices is None or choices == []:
            return None
        msgBox = QtGui.QMessageBox()
        msgBox.setText(text)
        combo = QtGui.QComboBox()
        combo.addItems(choices)
        msgBox.layout().addWidget(combo, 1, 0)
        msgBox.addButton(QtGui.QPushButton('Ok'), QtGui.QMessageBox.AcceptRole)
        msgBox.addButton(QtGui.QPushButton('Annuler'), QtGui.QMessageBox.RejectRole)
        ret = msgBox.exec_()
        if ret == QtGui.QMessageBox.Cancel:
            return None
        return str(combo.currentText())


    def loadPTS(self, path):
        """Load a PTS file (tries to find a suitable electrode model in the available templates)  """
        refId = self.preReferential().uuid()
        els = []
        elecs = {}
        f = open(path, 'r')
        lines = f.readlines()
        f.close()
    
        # The coordinates in the PTS are expressed in which referential ?
        refOfPts = self.comboMessageBox(u'Importation d\'un fichier PTS. Choisissez le référentiel utilisé (Scanner-based...)', sorted(self.refConv.availableReferentials().keys()))
        if refOfPts is None:
            print "User cancelled PTS importation, or no valid referential found"
            return
        print "PTS import referential : %s"%repr(refOfPts)
        lines.reverse()
        if not lines.pop().startswith('ptsfile'):
            print 'This is not a valid PTS file !'
            return (refId, [])
        lines.pop() # Useless line 1 1 1
        nb = int(lines.pop()) # Number of contacts
        for i in range(nb): # All lines (= all plots)
            l = lines.pop().rstrip().split() # name x y z 0 0 0 2 2.0 (last ones may be contact length/diameter ?)
            name = ''.join(e for e in l[0] if not e.isdigit())
            
            print l[0]
            plot = ''.join(e for e in l[0] if e.isdigit())
            if plot == '':
                continue
            plot = int(plot)
            
            if name not in elecs:
                elecs[name]=[]
            
            coords = list(self.refConv.anyRef2Real([float(l[1]), float(l[2]),float(l[3])], refOfPts))
            nameplot = l[0]
            elecs[name].append ( [plot,] + coords)
        # Iterate over all electrodes
        for k,l in elecs.iteritems():
            # Try to get a real model from l [[numPlot, x,y,z], [...], ...]
            res = self.fitElecModelToPlots(l)
            if res is not None:
                els.append( {'name':k.lower(), 'model':res[0], 'target': res[1], 'entry':res[2]} )
            else:
                print "Could not find model matching electrode %s, adding individual contacts."%k
                for pl in l:
                    coords = pl[1:]
                    els.append( {'name':k+str(pl[0]), 'model':'plot2mmD1mmCentered', 'target': coords, 'entry':coords[0:2]+[coords[2]+1.0,]} )
    
        return (refId, els)


    def loadElectrodeTXT(self, path):
        """Load an electrode.txt file with triplets of lines : nameOfElectrode\n X1 Y1 Z1\n X2 Y2 Z2 such as used by ImaGIN"""
        els=[]
        f = open(path, 'r')
        lines = f.readlines()
        f.close()
        lines.reverse()
    
        # The coordinates in the TXT are expressed in which referential ?
        refOfTxt = self.comboMessageBox(u'Importation d\'un fichier electrode TXT. Choisissez le référentiel utilisé (Scanner-based a priori)', sorted(self.refConv.availableReferentials().keys()))
        if refOfTxt is None:
            print "User cancelled electrode TXT importation, or no valid referential found"
            return
        print "TXT electrode import referential : %s"%repr(refOfTxt)
    
        while (len(lines) >= 3):
            try:
                name = lines.pop()
                targ = map (float,lines.pop().replace(',','.').split())
                entr = map (float,lines.pop().replace(',','.').split())
                if len(targ) == 3 and len(entr) == 3 and len(name)>0:
                    targ = list(self.refConv.anyRef2Real(targ, refOfTxt))
                    entr = list(self.refConv.anyRef2Real(entr, refOfTxt))
                    els.append( {'name':name.lower(), 'model':fitElecModelToLength(targ, entr), 'target': targ, 'entry':entr} )
                else:
                    print "Invalid lines in electrode txt file : %s, %s, %s"%(repr(name), repr(targ), repr(entr))
            except:
                pass
    
        refId = self.preReferential().uuid()
        return (refId, els)


    def loadElectrodes(self, patientName=None, patientCenter=None, isGui=True):
        """Load electrode implantation (if already entered) from BrainVisa or from a file"""
        path = None
        
        if patientCenter is None:
            path = str(QtGui.QFileDialog.getOpenFileName(self, "Open electrode implantation", "", "All implantations(*.elecimplant *.pts *.txt);;Electrode implantation (*.elecimplant);;PTS file (*.pts);;Electrode.txt (*.txt)"))
        else:
            rdi = ReadDiskItem( 'Electrode implantation', 'Electrode Implantation format', requiredAttributes={'subject':patientName, 'center':patientCenter} )
            elecs = list(rdi._findValues( {}, None, False ) )
            if len(elecs) == 1:
                path = elecs[0].fileName()
            elif len(elecs) > 1:
                print "CAREFUL : more than one electrode implantation are available, strange -> load the first found one" # TODO Dialogue de choix
                path = elecs[0].fileName()
            else: # No implantation available
                print "no electrode implantation found"
                return
    
        if not path:
            return
        # Check if we have a PTS/TXT/elecimplant file
        extension = os.path.splitext(path)[1].lower()
        els = []
        refId = None
        if extension == '.elecimplant':
            filein = open(path, 'rb')
            try:
                dic = json.loads(filein.read())
            except:
                filein.close()
                filein = open(path, 'rU')
                dic = pickle.load(filein)
            
            filein.close()
            els = dic['electrodes']
            refId = dic['ReferentialUuid']
    
        elif extension == '.txt':
            (refId, els) = self.loadElectrodeTXT(path) # Verifier que c'est un electrode.txt, si c'est un electrode_Pos ou electrode_Name on peut zapper ?
        elif extension == '.pts':
            (refId, els) = self.loadPTS(path)
          # Demander le référentiel, ou deviner si le nom contient _MNI ?
        else:
            print "file format unknown : %s !"%extension
            QtGui.QMessageBox.warning(self, u'Erreur', u"the file format has not been recognized : %s"%extension)
    
        if refId != self.preReferential().uuid():
            print "CAREFUL: electrodes load are defined in an other referential that the one of the T1 pre, problem possible !"
        for e in els:
            self.addElectrode(e['name'], e['model'], e['target'], e['entry'], refId, False, isGui)
        # Update display
        if isGui:
            self.updateElectrodeMeshes()
            # Update display if the function was called from a button click
            if not patientName:
                self.updateAllWindows(True)


    def saveElectrodes(self):
        """Save electrode implantation in BrainVisa Database"""
        # Saving : electrode model, name, target and entry point
        els = [dict((k,el[k]) for k in ['name', 'target', 'entry', 'model']) for el in self.electrodes]
        # Save Referential UID which is the base of electrode coordinates
        refId = self.preReferential().uuid()
        path = None
        di = None
        if self.brainvisaPatientAttributes is not None:
            wdi = WriteDiskItem( 'Electrode implantation', 'Electrode Implantation format' )
            di = wdi.findValue({'subject':self.brainvisaPatientAttributes['subject'], 'center':self.brainvisaPatientAttributes['center']} )
            if not di:
                print "Cannot find a valid path to store the electrode implantation for the current patient !"
                QtGui.QMessageBox.warning(self, u'Erreur', u"Cannot find a valid path to store the electrode implantation for the current patient !")
            else:
                path = di.fileName()
                createItemDirs(di)
        if not path:
            path = str(QtGui.QFileDialog.getSaveFileName(self, "Save electrode implantation", "", "Electrode implantation (*.elecimplant)"))
        if not path:
            return
        # If there is no extension, add the standard one !
        if os.path.splitext(path)[1] == '':
            path = path+'.elecimplant'
    
        plotsT1preRef = self.getPlotsT1preRef()
        info_plotsT1Ref= []
        for k,v in plotsT1preRef.iteritems():
            plot_name_split = k.split('-$&_&$-')
            info_plotsT1Ref.append((plot_name_split[0]+plot_name_split[1][4:].zfill(2),v.list()))
            #plots_label[k]=(label,label_name)
    
        plotsT1Ref_sorted = sorted(info_plotsT1Ref, key=lambda plot_number: plot_number[0])
        #previous_data.update({'plotsMNI':info_plotMNI})
        
        fileout = open(path, 'wb')
        fileout.write(json.dumps({'electrodes':els, 'ReferentialUuid':refId, '2mni':None, 'timestamp':time.time(),'plotsT1Nat':plotsT1Ref_sorted}))
        #pickle.dump({'electrodes':els, 'ReferentialUuid':refId, '2mni':None, 'timestamp':time.time()}, fileout)
        fileout.close()
        if di is not None:
            neuroHierarchy.databases.insertDiskItem( di, update=True )
        QtGui.QMessageBox.information(self, u'Implantation saved', u"Implantation has been saved in database.\nThe MNI coordinates of the contacts must be computed before the CSV generation.")
        self.isModified = False
        
    
    # EXPORT ALL INFO (INTERACTIVE)
    def exportAll(self):
        # Check if electrodes were modified
        if self.isModified:
            reply = QtGui.QMessageBox.question(self, 'Message', "Save modifications before exporting?",
                    QtGui.QMessageBox.Yes | QtGui.QMessageBox.No, QtGui.QMessageBox.No)
            if (reply == QtGui.QMessageBox.Yes):
                self.saveElectrodes()
        # If no subject is selected: stop
        selPatients = self.patientList.selectedItems()
        if not selPatients:
            return
        # Ask which options needed to be executed before the export
        dialog = DialogCheckbox([\
            "Compute MNI coordinates for all contacts",\
            "Compute parcels for all contacts",\
            "Compute MarsAtlas resection position",\
            "Compute parcel metrics",\
            "Save contact coordinates (.pts/.txt files)",\
            "Save contact info (CSV file)",\
            "Save contact info (BIDS .tsv)", \
            "Save screenshots",\
            "Save video (MP4)"],\
            "Export", "Select options to run:",\
            [True, True, False, False, True, True, self.bidspath != None, False, False])
        selOptions = dialog.exec_()
        # If user cancelled the selection
        if selOptions is None:
            return
        # If no patient is not loaded: load all the selected patients in a loop
        pSuccess = []
        pFailed = []
        isLoad = self.patientList.isEnabled()
        for p in selPatients:
            # Get patient name
            patient = p.text()
            # Load patient
            if isLoad:
                self.loadPatient(patient)
                
            # Run export with a progress bar
            res = ProgressDialog.call(lambda thr:self.exportAllWorker(selOptions, thr), True, self, "Processing...", "Export: " + patient)
            # res = self.exportAllWorker(selOptions)
            
            # Unload patient
            if isLoad:
                self.changePatient()

            # Display errors and new files (if processing only one subject)
            if (len(selPatients) == 1) and res:
                if res[1]:    # errMsg
                    QtGui.QMessageBox.critical(self, u'Export error', u"Errors occured during the export: \n\n" + u"\n".join(res[1]))
                if res[0]:    # newFiles
                    QtGui.QMessageBox.information(self, u'Export done', u"New files saved in the database: \n\n" + u"\n".join(res[0]))
            # If processing multiple subjects: log success/failure
            elif res:
                if res[1]:
                    pFailed.append(patient)
                else:
                    pSuccess.append(patient)
            else:
                pFailed.append(patient)
                
        # Display final report after a loop of export of multiple subjects
        if pFailed or pSuccess:
            strSummary = u"%d patients exported successfully\n%d patients failed\n\nCheck the console for more details."%(len(pSuccess), len(pFailed))
            strDetails = "Success: \n" + "\n".join(pSuccess) + "\n\nFailed:\n" + "\n".join(pFailed)
            # QtGui.QMessageBox.information(self, u'Export done', u"%d patients exported successfully\n%d patients failed\n\nCheck the console for more details."%(len(pSuccess), len(pFailed)))
            print("\n\n" + strDetails)
            msgBox = QtGui.QMessageBox()
            msgBox.setText(strSummary)
            msgBox.setWindowTitle(u'Export done')
            msgBox.setDetailedText(strDetails)
            msgBox.exec_()


    # EXPORT ALL INFO (WORKER)
    def exportAllWorker(self, selOptions, thread=None):
        """ Normalize and export the contact coordinates """
        
        # List of new files that are added to the database
        newFiles = []
        errMsg = []
        # Get coordinates of SEEG contact centers
        plots = self.getPlotsT1preScannerBasedRef()
        if not plots:
            errMsg += ["No electrodes available."]
            return [newFiles, errMsg]

        # Get options from selection
        isComputeMni = selOptions[0]
        isComputeParcels = selOptions[1]
        isMarsatlasResection = selOptions[2]
        isParcelMetrics = selOptions[3]
        isSavePts = selOptions[4]
        isSaveCsv = selOptions[5]
        isSaveTsv = selOptions[6]
        isScreenshot = selOptions[7]
        isVideo = selOptions[8]

        # ===== MNI COORDINATES =====
        plots_MNI = None
        if isComputeMni:
            if (self.getT1preMniTransform() is not None):
                if thread:
                    thread.progress.emit(5)
                    thread.progress_text.emit("Computing MNI coordinates...")
                [plots_MNI, mniFiles] = self.computeMniCoord()
                newFiles += mniFiles
            else:
                errMsg += ["Cannot compute MNI coordinates for the contacts."]
            
        # ===== Save TXT/PTS files ======
        if isSavePts:
            if thread:
                thread.progress.emit(20)
                thread.progress_text.emit("Exporting PTS files...")
            # ===== T1-pre scanner-based referential ======
            # Get current time stamp
            timestamp = time.time()
            # Get output files from the database
            wdipts = WriteDiskItem('Electrode Implantation PTS', 'PTS format', requiredAttributes={'no_ref_name':'True'}).findValue(self.diskItems['T1pre'])#,_debug= sys.stdout requiredAttributes={'ref_name':'default'}
            wditxtpos = WriteDiskItem('Electrode Implantation Position TXT', 'Text file', requiredAttributes={'no_ref_name':'True'}).findValue(self.diskItems['T1pre'])
            wditxtname = WriteDiskItem('Electrode Implantation Name TXT', 'Text file', requiredAttributes={'no_ref_name':'True'}).findValue(self.diskItems['T1pre'])
            # Save PTS and TXT files
            if (wdipts is not None) and (wditxtpos is not None) and (wditxtname is not None):
                # Save PTS
                self.savePTS(path = wdipts.fullPath(), contacts=plots)
                wdipts.setMinf('referential', self.t1pre2ScannerBasedId)
                wdipts.setMinf('timestamp', timestamp)
                neuroHierarchy.databases.insertDiskItem(wdipts, update=True )
                # Save TXT files
                self.saveTXT(pathPos=wditxtpos.fullPath(), pathName=wditxtname.fullPath(), contacts=plots)
                neuroHierarchy.databases.insertDiskItem(wditxtname, update=True )
                wditxtpos.setMinf('referential', self.t1pre2ScannerBasedId)
                wditxtpos.setMinf('timestamp', timestamp)
                neuroHierarchy.databases.insertDiskItem(wditxtpos, update=True )
                # List of new registered files
                newFiles += [wdipts.fullPath(), wditxtpos.fullPath(), wditxtname.fullPath()]
            else:
                errMsg += ["Cannot find a path to save T1pre PTS or TXT files in the database."]
                
            # ===== MNI referential =====
            # Get MNI coordinates if computation was not enforced previously
            if not isComputeMni:
                plots_MNI = self.getPlotsMNIRef(False)
                if plots_MNI is None:
                    errMsg += ["Cannot compute MNI coordinates for the contacts."]
            if plots_MNI is not None:
                # Get output files from the database
                wdiptsmni = WriteDiskItem('Electrode Implantation PTS', 'PTS format', requiredAttributes={'ref_name':'MNI'}).findValue(self.diskItems['T1pre'])
                wditxtmnipos = WriteDiskItem('Electrode Implantation Position TXT', 'Text file', requiredAttributes={'ref_name':'MNI'}).findValue(self.diskItems['T1pre'])
                wditxtmniname = WriteDiskItem('Electrode Implantation Name TXT', 'Text file', requiredAttributes={'ref_name':'MNI'}).findValue(self.diskItems['T1pre'])
                # Save PTS and TXT files
                if (wdiptsmni is not None) and (wditxtmnipos is not None) and (wditxtmniname is not None):
                    # Save PTS
                    self.savePTS(path = wdiptsmni.fullPath(), contacts = plots_MNI)
                    wdiptsmni.setMinf('referential', self.mniReferentialId())
                    wdiptsmni.setMinf('timestamp', timestamp)
                    neuroHierarchy.databases.insertDiskItem(wdiptsmni, update=True )
                    # Save TXT files
                    self.saveTXT(pathPos=wditxtmnipos.fullPath(), pathName=wditxtmniname.fullPath(), contacts = plots_MNI)
                    neuroHierarchy.databases.insertDiskItem(wditxtmniname, update=True )
                    wditxtmnipos.setMinf('referential', self.mniReferentialId())
                    wditxtmnipos.setMinf('timestamp', timestamp)
                    neuroHierarchy.databases.insertDiskItem(wditxtmnipos, update=True )
                    # List of new registered files
                    newFiles += [wdiptsmni.fullPath(), wditxtmnipos.fullPath(), wditxtmniname.fullPath()]
                else:
                    errMsg += ["Cannot find a path to save MNI PTS or TXT files in the database."]

            # ===== AC-PC referential =====
            if self.refConv.isRefAvailable('AC-PC'):
                # Get output files from the database
                wdiptsacpc = WriteDiskItem('Electrode Implantation PTS', 'PTS format', requiredAttributes={'ref_name':'ACPC'}).findValue(self.diskItems['T1pre'])
                # Save PTS
                if (wdiptsacpc is not None):
                    plotsACPC = self.getPlotsAnyReferential('AC-PC')
                    self.savePTS(path = wdiptsacpc.fullPath(), contacts=plotsACPC)
                    wdiptsacpc.setMinf('referential', 'AC-PC')
                    wdiptsacpc.setMinf('timestamp', timestamp)
                    neuroHierarchy.databases.insertDiskItem(wdiptsacpc, update=True )
                    # List of new registered files
                    newFiles += [wdiptsacpc.fullPath()]
                else:    
                    errMsg += ["Cannot find a path to save AC-PC PTS or TXT files in the database."]
                
                
        # ===== MARS ATLAS =====
        # Compute parcels 
        if isComputeParcels:
            if thread:
                thread.progress.emit(25)
                thread.progress_text.emit("Computing contact parcels...")
            res = self.computeParcels(self.diskItems['T1pre'])
            newFiles += res[0]
            errMsg += res[1]
        # Compute MarsAtlas resection
        if isMarsatlasResection:
            if thread:
                thread.progress.emit(40)
                thread.progress_text.emit("Computing resection parcels...")
            newFiles += self.computeMarsatlasResection()
        # Compute parcel metrics
        if isParcelMetrics:
            if thread:
                thread.progress.emit(50)
                thread.progress_text.emit("Computing parcel metrics...")
            newFiles += self.computeParcelMetrics()
            
        # ===== OTHER EXPORTS =====
        # Save CSV
        if isSaveCsv:
            if thread:
                thread.progress.emit(60)
                thread.progress_text.emit("Generating CSV files...")
            res = self.saveCSV(self.diskItems['T1pre'])
            newFiles += res[0]
            errMsg += res[1]
        # Save TSV
        if isSaveTsv:
            if thread:
                thread.progress.emit(70)
                thread.progress_text.emit("Generating BIDS TSV files...")
            res = self.saveBidsTsv(self.diskItems['T1pre'])
            newFiles += res[0]
            errMsg += res[1]
        # Save screenshots
        if isScreenshot:
            if thread:
                thread.progress.emit(80)
                thread.progress_text.emit("Generating screenshots...")
            newFiles += self.saveScreenshot()
        # Save video
        if isVideo:
            if thread:
                thread.progress.emit(90)
                thread.progress_text.emit("Generating videos...")
            newFiles += self.saveMP4()

        return [newFiles, errMsg]
    

    def saveScreenshot(self):
        #Check if all the volumes are available.
        #self.verifParcel allows us to know if the verification has already been made, in that case we won't do it again
        Mask_left = ReadDiskItem('Left Gyri Volume', 'Aims writable volume formats',requiredAttributes={'subject':self.brainvisaPatientAttributes['subject'], 'center':self.currentProtocol })
        diMaskleft = Mask_left.findValue(self.diskItems['T1pre'])
        
        Mask_right = ReadDiskItem('Right Gyri Volume', 'Aims writable volume formats',requiredAttributes={'subject':self.brainvisaPatientAttributes['subject'], 'center':self.currentProtocol })
        diMaskright = Mask_right.findValue(self.diskItems['T1pre'])
        
        MaskGW_right = ReadDiskItem('Right Grey White Mask','Aims writable volume formats',requiredAttributes={'subject':self.brainvisaPatientAttributes['subject'], 'center':self.currentProtocol })
        diMaskGW_right = MaskGW_right.findValue(self.diskItems['T1pre'])
        
        MaskGW_left = ReadDiskItem('Left Grey White Mask','Aims writable volume formats',requiredAttributes={'subject':self.brainvisaPatientAttributes['subject'], 'center':self.currentProtocol })
        diMaskGW_left = MaskGW_left.findValue(self.diskItems['T1pre'])
        
        if diMaskleft is None or diMaskright is None or diMaskGW_right is None or diMaskGW_left is None :
            #if a parcellation is missing we won't be able to make the screenshot, so we exit the function
            print "No parcellisation found."  
            return []
        else:
            #we take out the cursor
            self.a.config()[ 'linkedCursor' ] = 0
            #and replace the electrodes realistic model with 2mm spheres, which are more visible
            self.updateDispMode(2)
            
            #Right hemisphere
            #opens an Anatomist window and loads the parcellation and the electrodes
            w2=self.a.createWindow("Sagittal")
            contWin=self.windowContent['MRI pre + right MarsAtlas']
            content=[x for x in contWin if "rightMARSATLAS" in x]
            ele=[x for x in contWin if "electrode" in x]
            elec=self.dispObj[ele[0]] 
            for el in elec:
                w2.addObjects(el)
            image = self.dispObj[content[0]]
            w2.addObjects(image)
            time.sleep(1)
            w2.windowConfig(snapshot=os.path.join(getTmpDir(),self.brainvisaPatientAttributes['subject']+"RightHemiSagittal.png"))
            #rotates the view by 180°
            w2.camera(view_quaternion=[-0.5, 0.5, 0.5, -0.5])     
            #waits a second so that the window has the tame to be updated before the next screenshot, they are saved in the tmp repertory
            time.sleep(1)
            w2.windowConfig(snapshot=os.path.join(getTmpDir(),self.brainvisaPatientAttributes['subject']+"RightHemiSagittalRot.png"))
            
            w2.close()
            
            #Left hemisphere          
            w21=self.a.createWindow("Sagittal")
            contWin=self.windowContent['MRI pre + left MarsAtlas']
            content=[x for x in contWin if "leftMARSATLAS" in x]
            ele=[x for x in contWin if "electrode" in x]
            elec=self.dispObj[ele[0]] 
            for el in elec:
                w21.addObjects(el)
            image = self.dispObj[content[0]]
            w21.addObjects(image)
            time.sleep(1)
            w21.windowConfig(snapshot=os.path.join(getTmpDir(),self.brainvisaPatientAttributes['subject']+"LeftHemiSagittal.png"))
            w21.camera(view_quaternion=[-0.5, 0.5, 0.5, -0.5])          
            time.sleep(1)
            w21.windowConfig(snapshot=os.path.join(getTmpDir(),self.brainvisaPatientAttributes['subject']+"LeftHemiSagittalRot.png"))
            w21.close()
        
            #puts the cursor back
            self.a.config()[ 'linkedCursor' ] = 1
            #opens the pictures of the left and right brain which were taken previously
            images1 = map(Image.open, [os.path.join(getTmpDir(),self.brainvisaPatientAttributes['subject']+"RightHemiSagittal.png"),\
                                       os.path.join(getTmpDir(),self.brainvisaPatientAttributes['subject']+"RightHemiSagittalRot.png")])
            images2 = map(Image.open, [os.path.join(getTmpDir(),self.brainvisaPatientAttributes['subject']+"LeftHemiSagittal.png"),\
                                       os.path.join(getTmpDir(),self.brainvisaPatientAttributes['subject']+"LeftHemiSagittalRot.png")])
            
            #calculates the width and height of the new picture we are going to make
            widths1, heights1 = zip(*(i.size for i in images1))
            widths2, heights2 = zip(*(i.size for i in images2))
            
            total_width1 = sum(widths1)
            max_height1 = max(heights1)
            total_width2 = sum(widths2)
            max_height2 = max(heights2)
            
            #opens new empty images of the size calculated before
            new_im1 = Image.new('RGB', (total_width1, max_height1))
            new_im2 = Image.new('RGB', (total_width2, max_height2))
          
            #we are then going to paste the two images from right and left hemispheres together by the x axis
            x_offset = 0
            for im in images1:
                new_im1.paste(im, (x_offset,0))
                x_offset += im.size[0]
                
            new_im1.save(os.path.join(getTmpDir(), 'associated1.png'))
            
            x_offset = 0
            for im in images2:
                new_im2.paste(im, (x_offset,0))
                x_offset += im.size[0]   
            
            new_im2.save(os.path.join(getTmpDir(),'associated2.png'))
            
            #pastes the two images obtained previously by the y axis
            images = map(Image.open, [os.path.join(getTmpDir(),"associated1.png"), os.path.join(getTmpDir(),"associated2.png")])    
            new_im=Image.new('RGB', (total_width1, max_height1*2))
            
            y_offset = 0
            for im in images:
                new_im.paste(im, (0,y_offset))
                y_offset += im.size[1]  
            
            new_im.save(os.path.join(getTmpDir(), 'associated.png'))
            #verification that the path is creatable
            wdi = WriteDiskItem('Screenshot of Mars Atlas','PNG image')
            di=wdi.findValue(self.diskItems['T1pre'])
        
            if di is None:
                print "Can't generate file"
                return []
            else:
                try:
                    os.mkdir(os.path.dirname(di.fullPath()))                                                          
                    #line0 = runCmd(cmd0) 
                except:
                    line0=1
                cmd1 = ['mv', os.path.join(getTmpDir(),'associated.png'), str(di.fullPath())]
                line1 = runCmd(cmd1)
                #updates the database with the image of the 2 views of the 2 parcellation
                neuroHierarchy.databases.insertDiskItem(di, update=True )
        #puts back the realistic electrode model        
        self.updateDispMode(0)
        return [di.fullPath()]
          
          
    def saveMP4(self):
        from brainvisa import quaternion
        
        newFiles = []
        #takes the voxel size of the T1
        volume = aims.read(str(self.diskItems['T1pre']))
        sizeT1=volume.getVoxelSize()
        
        #verification that a CT or MRI is present, if both are, we take the CT
        T1post=None
        diCT = ReadDiskItem( 'CT', 'BrainVISA volume formats', requiredAttributes={'center':self.brainvisaPatientAttributes['center'], 'subject':self.brainvisaPatientAttributes['subject'] } )
        path = list(diCT.findValues({}, None, False ))
        #check if it is a CT post

        idCTpost = [x for x in range(len(path)) if 'CTpost' in str(path[x]) and not 'CTpostOp' in str(path[x])]

        try:
            path=path[idCTpost[0]].fullPath()
            volCT = aims.read(path)
            npCT = volCT.arraydata()
            PresCT = True
        except: 
            diMRI = ReadDiskItem( 'Raw T1 MRI', 'BrainVISA volume formats', requiredAttributes={'center':self.brainvisaPatientAttributes['center'], 'subject':self.brainvisaPatientAttributes['subject'] } )
            pathMRI = list(diMRI.findValues({}, None, False ))         
            id_post = [x for x in range(len(pathMRI)) if 'post' in str(pathMRI[x]) and not 'postOp' in str(pathMRI[x])]
            pathMRI=pathMRI[0].fullPath()
            volCT = aims.read(pathMRI)
            npCT = volCT.arraydata()
            PresCT=False
            if id_post!=[]:    
                T1post = pathMRI[id_post[0]]
        
        #if both CT and MRI are absent, we won't do the computation

        if T1post is None and idCTpost==[]:
            warning = QtGui.QMessageBox(self)
            warning.setText("No CT post or MRI post found")
            warning.setWindowTitle("Warning")
            return []
        
        #finding the brainMask if there is one
        try:
            brainMask = ReadDiskItem('Brain Mask', 'aims readable volume formats',requiredAttributes={'subject':self.brainvisaPatientAttributes['subject'], 'center':self.currentProtocol })
            diBrain = list(brainMask.findValues({}, None, False ))
            diBrain0=diBrain[0].fullPath()
            volBrainMask = aims.read(diBrain0)
            brainMaskArray = volBrainMask.arraydata() 
        except:
            brainMaskArray=None 
        
        #if the brainMask is present, the gif will be made taking only the brain tissues
        if brainMaskArray is not None:
            z=0
            while z<brainMaskArray.shape[1]:
                maxi=brainMaskArray[0,z,:,:].max()
                if maxi != 0:
                    i=round(z*sizeT1[2])
                    z=brainMaskArray.shape[1]
                z+=1
            z-=2
            while z>0:
                maxi=brainMaskArray[0,z,:,:].max()
                if maxi != 0:
                    limit=round(z*sizeT1[2])
                    z=0
                z-=1 
        else:
            i=1  
            limit=npCT.shape[1]
            
        #Erase all previous files
        os.system('find ' + os.path.join(getTmpDir(),'gif') + '*.png -type f -exec rm -f {} \;')
        
        #Take out the cursor
        self.a.config()[ 'linkedCursor' ] = 0
        #put the realistic electrode model back
        self.updateDispMode(2)
        
        #Electrodes GIF 
        #load the window
        w21=self.a.createWindow("Axial")
        w21.windowConfig(fullscreen=1)
        #if it is a CT
        if PresCT==True:
            contWin=self.windowContent['CT post']
            content=[x for x in contWin if "CTpost" in x]
            ele=[x for x in contWin if "electrode" in x]
        #if it is an MRI    
        else:
            contWin=self.windowContent['MRI post']
            content=[x for x in contWin if "T1post" in x]
            ele=[x for x in contWin if "electrode" in x]
        #add the image and electrodes    
        image = self.dispObj[content[0]]
        elec=self.dispObj[ele[0]] 
        w21.addObjects(image)
        for el in elec:
            w21.addObjects(el)
        #get all referentials    
        refs=self.a.getReferentials()
        #choose the T1 pre referential...
        for el in refs:
            info=str(el.getInfos())
            if 'T1pre' in info and 'native' in info:
                r1=el    
        #...and assign it        
        w21.assignReferential(r1)
        #activate the clipping
        self.a.execute('WindowConfig',windows = [w21],clipping=2,clip_distance=5.)
        if brainMaskArray is None:
            limit=npCT.shape[1]
        #take pictures of slices following the z axis    
        snapShots=[]
        while i<limit:
            self.a.execute( 'LinkedCursor', window=w21, position=( 0, 0, i ) )   
            w21.windowConfig(snapshot=os.path.join(getTmpDir(), "gifCT%03i.png"%(i)))
            snapShots.append(os.path.join(getTmpDir(), "gifCT%03i.png"%(i)))
            i+=1    
        w21.close()
        self.brainvisaContext.runProcess('mpegEncode_avconv', animation=os.path.join(getTmpDir(),"animationElec.mp4"),images = snapShots[:-1],encoding='h264',quality=75,passes=1)

        wdi = WriteDiskItem('MP4 of Electrodes','MP4 film')
        if PresCT==True:
            di=wdi.findValue(self.diskItems['CTpost'])
        else:
            di=wdi.findValue(self.diskItems['T1post'])
              
        if di is None:
            print "Can't generate file"
            return []
        else:
            try:
               os.mkdir(os.path.dirname(di.fullPath()))
            except:
               pass
               
            cmd1 = ['mv', os.path.join(getTmpDir(),"animationElec.mp4"), str(di.fullPath())]
            line1 = runCmd(cmd1)
            neuroHierarchy.databases.insertDiskItem(di, update=True )
            newFiles.append(di.fullPath())

        #Mars Atlas GIF
        #Check if all the volumes are available.
        Mask_left = ReadDiskItem('Left Gyri Volume', 'Aims writable volume formats',requiredAttributes={'subject':self.brainvisaPatientAttributes['subject'], 'center':self.currentProtocol })
        diMaskleft = Mask_left.findValue(self.diskItems['T1pre'])

        Mask_right = ReadDiskItem('Right Gyri Volume', 'Aims writable volume formats',requiredAttributes={'subject':self.brainvisaPatientAttributes['subject'], 'center':self.currentProtocol })
        diMaskright = Mask_right.findValue(self.diskItems['T1pre'])

        MaskGW_right = ReadDiskItem('Right Grey White Mask','Aims writable volume formats',requiredAttributes={'subject':self.brainvisaPatientAttributes['subject'], 'center':self.currentProtocol })
        diMaskGW_right = MaskGW_right.findValue(self.diskItems['T1pre'])

        MaskGW_left = ReadDiskItem('Left Grey White Mask','Aims writable volume formats',requiredAttributes={'subject':self.brainvisaPatientAttributes['subject'], 'center':self.currentProtocol })
        diMaskGW_left = MaskGW_left.findValue(self.diskItems['T1pre'])
            
        
        if diMaskleft is None or diMaskright is None or diMaskGW_right is None or diMaskGW_left is None :
            print "No parcellisation found."          
        else:
            #Right hemisphere
            w2=self.a.createWindow("Sagittal")
            w2.windowConfig(fullscreen=1)
            contWin=self.windowContent['MRI pre + right MarsAtlas']
            #contWin = self.windowContent['MRI pre + hemispheres']
            content=[x for x in contWin if "rightMARSATLAS" in x]
            #content = [x for x in contWin if "rightHemi" in x]
            ele=[x for x in contWin if "electrode" in x]
            elec=self.dispObj[ele[0]] 
            for el in elec:
                w2.addObjects(el)
            image = self.dispObj[content[0]]
            w2.addObjects(image)
            time.sleep(1)
            w2.windowConfig(snapshot=os.path.join(getTmpDir(),"gifR000.png"))
            i=1
            v0=[0.5, 0.5, 0.5, 0.5]
            v1=[-0.5,0.5,0.5,-0.5] 
            incr=float(1)/174
            while i<176:
                quat=quaternion.Quaternion((v0[0]*(1-incr*i)+v1[0]*incr*i,v0[1]*(1-incr*i)+v1[1]*incr*i,v0[2]*(1-incr*i)+v1[2]*incr*i,v0[3]*(1-incr*i)+v1[3]*incr*i)).normalized().vector()
                w2.camera(view_quaternion=quat)
                w2.windowConfig(snapshot=os.path.join(getTmpDir(),"gifR%03i.png"%(i)))
                i+=1  
            i=0    
            v1=[0.5, 0.5, 0.5, 0.5]
            v0=[-0.5, 0.5, 0.5, -0.5] 
            while i<175:
                quat=quaternion.Quaternion((v0[0]*(1-incr*i)-v1[0]*incr*i,v0[1]*(1-incr*i)-v1[1]*incr*i,v0[2]*(1-incr*i)-v1[2]*incr*i,v0[3]*(1-incr*i)-v1[3]*incr*i)).normalized().vector()
                w2.camera(view_quaternion=quat)                                                               
                w2.windowConfig(snapshot=os.path.join(getTmpDir(),"gifR%03i.png"%(i+175)))
                i+=1    
            w2.close()  
        
            #Left hemisphere
            w2=self.a.createWindow("Sagittal")
            w2.windowConfig(fullscreen=1)
            contWin=self.windowContent['MRI pre + left MarsAtlas']
            content=[x for x in contWin if "leftMARSATLAS" in x]
            ele=[x for x in contWin if "electrode" in x]
            elec=self.dispObj[ele[0]] 
            for el in elec:
                w2.addObjects(el)
            image = self.dispObj[content[0]]
            w2.addObjects(image)
            time.sleep(1)
            w2.windowConfig(snapshot=os.path.join(getTmpDir(),"gifL000.png"))
            i=1
            v0=[0.5, 0.5, 0.5, 0.5]
            v1=[-0.5,0.5,0.5,-0.5] 
            incr=float(1)/174
            while i<175:
                quat=quaternion.Quaternion((v0[0]*(1-incr*i)+v1[0]*incr*i,v0[1]*(1-incr*i)+v1[1]*incr*i,v0[2]*(1-incr*i)+v1[2]*incr*i,v0[3]*(1-incr*i)+v1[3]*incr*i)).normalized().vector()
                w2.camera(view_quaternion=quat)    
                w2.windowConfig(snapshot=os.path.join(getTmpDir(),"gifL%03i.png"%(i)))
                i+=1  
            i=0    
            v1=[0.5, 0.5, 0.5, 0.5]
            v0=[-0.5,0.5,0.5,-0.5] 
            while i<175:
                quat=quaternion.Quaternion((v0[0]*(1-incr*i)-v1[0]*incr*i,v0[1]*(1-incr*i)-v1[1]*incr*i,v0[2]*(1-incr*i)-v1[2]*incr*i,v0[3]*(1-incr*i)-v1[3]*incr*i)).normalized().vector()
                w2.camera(view_quaternion=quat)                                                               
                w2.windowConfig(snapshot=os.path.join(getTmpDir(),"gifL%03i.png"%(i+175)))
                i+=1    
            w2.close()  
        
            new_list_image = []
            j=0
            while j<350:
                images = map(Image.open,[os.path.join(getTmpDir(),"gifL%03i.png"%(j)), os.path.join(getTmpDir(),"gifR%03i.png"%(j))])
                widths, heights = zip(*(i.size for i in images))
                total_width = sum(widths)
                max_height = max(heights)
                new_im = Image.new('RGB', (total_width, max_height))
                x_offset = 0
                for im in images:
                    new_im.paste(im, (x_offset,0))
                    x_offset += im.size[0]
                new_im.save(os.path.join(getTmpDir(),'gifF%03i.jpg'%(j)), quality=70)
                new_list_image.append(os.path.join(getTmpDir(), 'gifF%03i.jpg'%(j)))
                j+=1    
            
            self.brainvisaContext.runProcess('mpegEncode_avconv', animation=os.path.join(getTmpDir(),'animationMA.mp4'),images = new_list_image,encoding='h264',quality=50,passes=1)
            wdi = WriteDiskItem('MP4 of Mars Atlas','MP4 film')
            di=wdi.findValue(self.diskItems['T1pre'])

            if di is None:
                print "Can't generate file"
                return []
            else:
                try:
                  os.mkdir(os.path.dirname(di.fullPath()))
                except:
                  pass

                cmd1 = ['mv', os.path.join(getTmpDir(),'animationMA.mp4'), str(di.fullPath())]
                line1 = runCmd(cmd1)
                neuroHierarchy.databases.insertDiskItem(di, update=True )
                newFiles.append(di.fullPath())
                
        self.updateDispMode(0)
        self.a.config()[ 'linkedCursor' ] = 1     
        self.a.execute('WindowConfig',windows = [self.wins[0]],clipping=0)
        
        print("MP4 done")
        return newFiles
    

    def fitvolumebyellipse(self,volumeposition):
        # from https://github.com/minillinim/ellipsoid/blob/master/ellipsoid.py
        tolerance = 0.025
        if len(volumeposition)>16000:
            volumeposition = [volumeposition[x] for x in range(len(volumeposition)) if x & 1 ==0]
        volumeposition=numpy.array(volumeposition)
        (N, dd) = numpy.shape(volumeposition)
        dd = float(dd)
        # Q will be our working array
        Q = numpy.vstack([numpy.copy(volumeposition.T), numpy.ones(N)])
        QT = Q.T
        # initializations
        err = 1.0 + tolerance
        u = (1.0 / N) * numpy.ones(N)
        # Khachiyan Algorithm
        while err > tolerance:
            V = numpy.dot(Q, numpy.dot(numpy.diag(u), QT))
            M = numpy.diag(numpy.dot(QT , numpy.dot(numpy.linalg.inv(V), Q)))    # M the diagonal vector of an NxN matrix
            j = numpy.argmax(M)
            maximum = M[j]
            step_size = (maximum - dd - 1.0) / ((dd + 1.0) * (maximum - 1.0))
            new_u = (1.0 - step_size) * u
            new_u[j] += step_size
            err = numpy.linalg.norm(new_u - u)
            u = new_u
        # center of the ellipse
        center = numpy.dot(volumeposition.T, u)
        # the A matrix for the ellipse
        AA = numpy.linalg.inv(numpy.dot(volumeposition.T, numpy.dot(numpy.diag(u), volumeposition)) - numpy.array([[a * b for b in center] for a in center])) / dd
        # Get the values we'd like to return
        UU, ss, rotation = numpy.linalg.svd(AA)
        radii = 1.0/numpy.sqrt(ss)
        return radii
   
          
    def computeParcelMetrics(self):
        import copy
        
        infos={}
        parcels={}
        tailleG=0
        
        #take the parcellations
        Mask_left = ReadDiskItem('Left Gyri Volume', 'Aims writable volume formats',requiredAttributes={'subject':self.brainvisaPatientAttributes['subject'], 'center':self.currentProtocol })
        diMaskleft = Mask_left.findValue(self.diskItems['T1pre'])
        
        Mask_right = ReadDiskItem('Right Gyri Volume', 'Aims writable volume formats',requiredAttributes={'subject':self.brainvisaPatientAttributes['subject'], 'center':self.currentProtocol })
        diMaskright = Mask_right.findValue(self.diskItems['T1pre'])
        
        #take the parcels names in one list       
        try:
            #transform the images to matrices
            Maskleft = aims.read(diMaskleft.fullPath())
            npleftMarsAtlas = Maskleft.arraydata()
        
            Maskright = aims.read(diMaskright.fullPath())
            nprightMarsAtlas = Maskright.arraydata()
        except:
            print "No parcellation found"
            return []
        #take voxel size of T1 pre
        volume = aims.read(str(self.diskItems['T1pre']))
        sizeT1=volume.getVoxelSize()
        
        #calculate the total volume of the parcels
        for key, value in labels['MarsAtlas'].items():
            if value[0].lower()=='l':
                parcel1 = numpy.where(npleftMarsAtlas == key)
            else:
                parcel1 = numpy.where(nprightMarsAtlas == key)
            coord=zip(parcel1[1],parcel1[2],parcel1[3])
            tailleG+=len(coord)
                
        print 'tailleG Done'     
        
        t=1
        #we do the computations for each parcel
        for key, value in labels['MarsAtlas'].items():
            print 'doing parcel ' ,value ,t
            t+=1
            par=0
            #we locate the parcel in the right or left hemispheres (if they begin with 'L':left, 'R':right)
            if value[0].lower()=='l':
                parcel1 = numpy.where(npleftMarsAtlas == key)
                try:
                    #we take only the parcel 
                    parcel=copy.deepcopy(npleftMarsAtlas[0,parcel1[1].min():parcel1[1].max()+1,parcel1[2].min():parcel1[2].max()+1,parcel1[3].min():parcel1[3].max()+1])
                except:
                    #if the parcel is absent we update it with the value 'not found'
                    print 'no parcel found'
                    parcels.update({value:"not found"})
                    par=None
            else:
                parcel1 = numpy.where(nprightMarsAtlas == key)
                try:
                    parcel=copy.deepcopy(nprightMarsAtlas[0,parcel1[1].min():parcel1[1].max()+1,parcel1[2].min():parcel1[2].max()+1,parcel1[3].min():parcel1[3].max()+1])
                except:
                    print 'no parcel found'
                    parcels.update({value:"not found"})
                    par=None
            #if the parcel has been computed        
            if par is not None:       
                #take the coordinates of the points
                coord=zip(parcel1[1],parcel1[2],parcel1[3])
                #volume of the parcel
                taille=len(coord)
                infos.update({"volume":int(taille)})
                parcelPositionN=[]
                parcelPositionNT=[]
                parcelPosition = [[round(parcel1[1][i]*sizeT1[2]),round(parcel1[2][i]*sizeT1[1]),round(parcel1[3][i]*sizeT1[0])] for i in range(len(parcel1[1]))]
                i=0
                #we take one over 3 values (downsampling by 3)
                while i<len(parcelPosition):
                    parcelPositionN.append(parcelPosition[i])
                    i+=3
                d=0
                if len(parcelPositionN)>10000:
                    downSamp=int(math.ceil(float(len(parcelPositionN))/10000))
                    while d<len(parcelPositionN):
                        parcelPositionNT.append(parcelPositionN[d])
                        d+=downSamp
                    print "Downsampling more"  
                    print 'downSamp',downSamp
                    parcelPositionN=parcelPositionNT
                #approximation by an ellips which will give the length width and height 
                fit=self.fitvolumebyellipse(parcelPositionN)
                print fit
                width=fit[0]
                height=fit[1]
                length=fit[2]
                infos.update({"width":width})
                infos.update({"length":length})
                infos.update({"height":height})
                parcel[parcel!=key]=0
                #center of mass of the parcel
                cM=ndimage.measurements.center_of_mass(parcel) 
                centreMasse=(cM[0]+parcel1[1].min()*sizeT1[2],cM[1]+parcel1[2].min()*sizeT1[1],cM[2]+parcel1[3].min()*sizeT1[0])
                infos.update({"center of mass":centreMasse})
                rapport=float(taille)/tailleG
                infos.update({"ratio":rapport})
                parcels.update({value:infos})
                infos={} 
        
        #transformation to the MNI referential
        points=[]
        for element in parcels:
            points.append(parcels[element]['center of mass'])
        
        field = self.getT1preMniTransform()
        if field is None:
            print "MNI deformation field not found"
            return []
        tmpOutput = os.path.join(getTmpDir(),'test.csv') #pour tester
        arr = numpy.asarray(points) #tous tes centres de masses pour toutes les parcels tel quel ([ [1,2,3], [4,5,6], [7,8,9] ])
        numpy.savetxt(tmpOutput, arr, delimiter=",")
        errMsg = matlabRun(spm12_normalizePoints % ("'"+self.spmpath+"'",field, tmpOutput, tmpOutput))
        if errMsg:
            print errMsg
            return []
        out = numpy.loadtxt(tmpOutput, delimiter=",")
        os.remove(tmpOutput)
        for element in out:
            for parcel in parcels:
                if all([points[[list(x) for x in tuple(out)].index(list(element))][i] == parcels[parcel]['center of mass'][i] for i in range(len(parcels[parcel]['center of mass']))]):
                #if parcels[parcel]['center of mass']==points[[list(x) for x in tuple(out)].index(list(element))]:
                    parcels[parcel]['center of mass']=tuple(element)
        valPatient={"Mars Atlas":parcels} 
        wdi = WriteDiskItem('Parcels Infos','Atlas metrics')
        di=wdi.findValue(self.diskItems['T1pre'])     
        
        #writes the file
        fout = open(di.fullPath(),'w')
        fout.write(json.dumps(valPatient))
        fout.close()
        
        #updates the database
        neuroHierarchy.databases.insertDiskItem(di, update=True )
        
        print "parcel metrics done"
        return [di.fullPath()]
    

    # ===== EXPORT CSV FILES =====
    def saveCSV(self, diT1pre, plots_MNI=None):
        # === GET DATABASE FILES === 
        # eleclabel
        rdi_eleclabel = ReadDiskItem('Electrodes Labels','Electrode Label Format',requiredAttributes={'subject':diT1pre['subject'], 'center':diT1pre['center']})
        di_eleclabel = rdi_eleclabel.findValue(diT1pre)
        if not di_eleclabel:
            return [[], ["File eleclabel not found."]]
        fileEleclabel = di_eleclabel.fullPath()
        # reseclabel
        rdi_reseclabel = ReadDiskItem('Resection Description','Resection json',requiredAttributes={'subject':diT1pre['subject'], 'center':diT1pre['center']})
        di_reseclabel = rdi_reseclabel.findValue(diT1pre)
        if di_reseclabel:
            fileResecLabel = di_reseclabel.fullPath()
        else:
            fileResecLabel = None
        # Exported CSV
        wdi = WriteDiskItem('Final Export Dictionaries','CSV file')
        di = wdi.findValue(diT1pre)
        if not di:
            return [[], ["Could not find where to save .csv file."]]
        fileCsv = di.fullPath()
        
        # ===== GET COORDINATES =====
        # Get scanner-based coordinates
        if not plots_MNI:
            plots_SB = self.getPlotsT1preScannerBasedRef()
            plots_SB = self.sortContacts(plots_SB)
            # Bipolar montage: Scanner-based coordinates
            bip_SB = []
            for pindex in range(1,len(plots_SB)):
                prev_elec = "".join([i for i in plots_SB[pindex-1][0] if not i.isdigit()])
                cur_elec = "".join([i for i in plots_SB[pindex][0] if not i.isdigit()])
                prev_ind = int("".join([i for i in plots_SB[pindex-1][0] if i.isdigit()]))
                cur_ind = int("".join([i for i in plots_SB[pindex][0] if i.isdigit()]))
                if (prev_elec == cur_elec) and (prev_ind == cur_ind - 1):
                    bip_SB.append((plots_SB[pindex-1][0] + '-' + plots_SB[pindex][0],(numpy.array(plots_SB[pindex][1])+numpy.array(plots_SB[pindex-1][1]))/2 ))
            plots_SB = dict(plots_SB)
            bip_SB = dict(bip_SB)
        else:
            plots_SB = {}
            bip_SB = {}
            
        # MNI coordinates
        if not plots_MNI:
            plots_MNI = self.getPlotsMNIRef(False)
            if plots_MNI is None:
                return [[], ["MNI coordinates are not available."]]
            # Sort contacts by name and index
            plots_MNI = self.sortContacts(plots_MNI)
        else:
            # Convert from dict to list
            plots_MNI = [(pname, plots_MNI[pname]) for pname in sorted(plots_MNI.keys())]
        # Bipolar montage: MNI coordinates
        bip_MNI = []
        for pindex in range(1,len(plots_MNI)):
            prev_elec = "".join([i for i in plots_MNI[pindex-1][0] if not i.isdigit()])
            cur_elec = "".join([i for i in plots_MNI[pindex][0] if not i.isdigit()])
            prev_ind = int("".join([i for i in plots_MNI[pindex-1][0] if i.isdigit()]))
            cur_ind = int("".join([i for i in plots_MNI[pindex][0] if i.isdigit()]))
            if (prev_elec == cur_elec) and (prev_ind == cur_ind - 1):
                bip_MNI.append((plots_MNI[pindex-1][0] + '-' + plots_MNI[pindex][0],(numpy.array(plots_MNI[pindex][1])+numpy.array(plots_MNI[pindex-1][1]))/2 ))
        plots_MNI = dict(plots_MNI)
        bip_MNI = dict(bip_MNI)
        
        
        # === GENERATE CSV ===
        # Read eleclabel file            
        fin = open(fileEleclabel, 'r')
        eleclabel = json.loads(fin.read())
        fin.close()
        # Write file
        with open(fileCsv, 'w') as fout:
            # Get new field InitialSegmentation
            if 'InitialSegmentation' in eleclabel['Template'].keys():
                InitialSegmentation = eleclabel['Template']['InitialSegmentation']
            else:
                InitialSegmentation = 'missing_info'
            # Write file header
            writer = csv.writer(fout, delimiter='\t')
            writer.writerow([u'Contacts Positions'])
            writer.writerow([u'Use of MNI Template', \
                'MarsAtlas', eleclabel['Template']['MarsAtlas'], \
                'Freesurfer', eleclabel['Template']['Freesurfer'], \
                'InitialSegmentation', InitialSegmentation])
            
            # List of parcel names to print
            parcelNames = [\
                'MarsAtlas', \
                'Destrieux', \
                'DKT', \
                'VEP', \
                'HCP-MMP1', \
                'Lausanne2008-33', \
                'Lausanne2008-60', \
                'Lausanne2008-125', \
                'Lausanne2008-250', \
                'Lausanne2008-500', \
                'GreyWhite', \
                'IntrAnat-MarsAtlas', \
                'MNI-MarsAtlas', \
                'MNI-Destrieux', \
                'MNI-DKT', \
                'MNI-Lausanne2008-33', \
                'MNI-Lausanne2008-60', \
                'MNI-Lausanne2008-125', \
                'MNI-Lausanne2008-250', \
                'MNI-Lausanne2008-500', \
                'MNI-AAL', \
                'MNI-AALDilate', \
                'MNI-AAL1_2018', \
                'MNI-AAL2', \
                'MNI-AAL3', \
                'MNI-Brodmann', \
                'MNI-BrodmannDilate', \
                'MNI-Hammers', \
                'MNI-HCP-MMP1', \
                'MNI-AICHA', \
                'MNI-JulichBrain', \
                'Resection rate']
            # Add list of column names
            colNames = [u'contact', 'MNI', 'T1pre Scanner Based', 'MarsAtlasFull'] + parcelNames
            # Print column names
            writer.writerow(colNames)
            
            # Print all contacts
            dict_sorted_tmp = OrderedDict(sorted(eleclabel['plots_label'].items()))
            for kk,vv in dict_sorted_tmp.iteritems():
                # Upper case for all the letters except from "p" that stand for ' (prime)
                contact_label = list(kk)
                for i in range(len(contact_label)):
                    # if not kk[i].isdigit() and ((i == 0) or (kk[i] != 'p')):
                    # We cannot consider that "Tp" should be kept unchanged, otherwise, "t'" is also converted to "Tp"
                    # Convention: In IntrAnat, all chars are upper case and electrode names can include "'", converted to 
                    # "p" in the .csv files
                    if contact_label[i].isalpha():
                        contact_label[i] = contact_label[i].upper()
                    elif (contact_label[i] == "'"):
                        contact_label[i] = 'p'
                contact_label = "".join(contact_label)
                listwrite = [contact_label]
                listwrite.append([float(format(plots_MNI[kk][i],'.3f')) for i in range(3)])
                if plots_SB:
                    listwrite.append([float(format(plots_SB[kk][i],'.3f')) for i in range(3)])
                else:
                    listwrite.append("N/A")
                if 'MarsAtlasFull' in vv.keys():
                    listwrite.append(vv['MarsAtlasFull'])
                else:
                    listwrite.append("N/A")
                for p in parcelNames:
                    if p in vv.keys():
                        listwrite.append(vv[p][1])
                    else:
                        listwrite.append("N/A")
                writer.writerow(listwrite)
            writer.writerow([])
            writer.writerow([])
            
            # Print all dipoles contacts
            dict_sorted_tmp = OrderedDict(sorted(eleclabel['plots_label_bipolar'].items()))
            for kk,vv in dict_sorted_tmp.iteritems():
                # Upper case for all the letters except from "p" that stand for ' (prime)
                contact_label = list(kk)
                for i in range(len(contact_label)):
                    # if not kk[i].isdigit() and ((i == 0) or (kk[i] != 'p')):
                    # We cannot consider that "Tp" should be kept unchanged, otherwise, "t'" is also converted to "Tp"
                    # Convetion: In IntrAnat, all chars are upper case and electrode names can include "'", converted to 
                    # "p" in the .csv files
                    if contact_label[i].isalpha():
                        contact_label[i] = contact_label[i].upper()
                    elif (contact_label[i] == "'"):
                        contact_label[i] = 'p'
                contact_label = "".join(contact_label)
                listwrite = [contact_label]
                listwrite.append([float(format(bip_MNI[kk][i],'.3f')) for i in range(3)])
                if bip_SB:
                    listwrite.append([float(format(bip_SB[kk][i],'.3f')) for i in range(3)])
                else:
                    listwrite.append("N/A")
                if 'MarsAtlasFull' in vv.keys():
                    listwrite.append(vv['MarsAtlasFull'])
                else:
                    listwrite.append("N/A")
                for p in parcelNames:
                    if p in vv.keys():
                        listwrite.append(vv[p][1])
                    else:
                        listwrite.append("N/A")
                writer.writerow(listwrite)
            writer.writerow([])
            writer.writerow([])

        # resecLabel
        if fileResecLabel:
            fin = open(fileResecLabel, 'r')
            info_label_resec = json.loads(fin.read())
            fin.close()
            with open(fileCsv, 'a') as fout:
                writer = csv.writer(fout, delimiter='\t')
                writer.writerow([u'Resection Information'])
                for kk,vv in info_label_resec.iteritems():
                    #writer.writerow([kk])
                    if type(vv) == type(float()):
                        listwrite = [kk,format(vv,'.1f')]
                        writer.writerow(listwrite)
                    else:
                        writer.writerow([kk])
                        for ll,bb in vv.iteritems():
                            listwrite = [ll, format(float(bb[1]),'.1f')]
                            writer.writerow(listwrite)

        # Reference csv file in the database
        neuroHierarchy.databases.insertDiskItem(di, update=True)
        print "Export csv done."
        return [[fileCsv], []]
    
    
    # ===== EXPORT CSV FILES =====
    def saveBidsTsv(self, diT1pre):
        # === GET BIDS SUBJECT ===
        # If BIDS database not set
        if not self.bidspath:
            return [[], ["BIDS database not set: Open ImageImport preferences and select BIDS folder."]]
        # Parse subject name: SIT_SESS_SUBj
        splitName = diT1pre['subject'].split('_')
        if (len(splitName) != 3):
            return [[], ["Subject '" + diT1pre['subject'] + "' not formatted as SIT_SESS_SUBj."]]
        bidsSes = 'ses-' + splitName[1].lower()
        bidsSub = 'sub-' + splitName[2].lower()
        # Check if folder exists
        sesPath = os.path.join(self.bidspath, bidsSub, bidsSes)
        if not os.path.exists(sesPath):
            return [[], ["Session does not exist in BIDS database: " + sesPath]]
        # Create ieeg folder if missing
        ieegPath = os.path.join(sesPath, 'ieeg')
        if not os.path.exists(ieegPath):
            try:
                os.mkdir(ieegPath)
            except:
                return [[], ["Could not create folder: " + ieegPath]]
        # Get T1pre file
        bidsT1w = 'N/A'
        anatPath = os.path.join(sesPath, 'anat')
        if os.path.exists(anatPath):
            anatDir = os.listdir(anatPath)
            anatNii = [x for x in anatDir if ('_T1w.nii' in x) and ('_acq-pre_' in x)]
            if len(anatNii) == 1:
                bidsT1w = os.path.join(bidsSub, bidsSes, 'anat', anatNii[0])
        # Target tsv spaces
        T1 = 'Other'
        MNI = 'IXI549Space'
          
        # === GET ELEC FILE === 
        # Find eleclabel
        rdi_eleclabel = ReadDiskItem('Electrodes Labels','Electrode Label Format',requiredAttributes={'subject':diT1pre['subject'], 'center':diT1pre['center']})
        di_eleclabel = rdi_eleclabel.findValue(diT1pre)
        if not di_eleclabel:
            return [[], ["File eleclabel not found."]]
        # Read eleclabel file            
        fin = open(di_eleclabel.fullPath(), 'r')
        eleclabel = json.loads(fin.read())
        fin.close()
        
        # === GET COORDINATES ===
        # Scanner-based coordinates
        plots = dict()
        plots[T1] = self.getPlotsT1preScannerBasedRef()
        plots[T1] = dict(self.sortContacts(plots[T1]))
        # MNI coordinates
        plots[MNI] = self.getPlotsMNIRef(False)
        if plots[MNI] is None:
            return [[], ["MNI coordinates are not available."]]
        plots[MNI] = dict(self.sortContacts(plots[MNI]))
        
        # === GENERATE ALL TSV ===
        outfiles = []
        for space in plots.keys():
            # Output TSV filename
            fileTsv = os.path.join(ieegPath, bidsSub + '_' + bidsSes + '_acq-intranat_space-' + space + '_electrodes.tsv')
            fileJson = os.path.join(ieegPath, bidsSub + '_' + bidsSes + '_acq-intranat_space-' + space + '_coordsystem.json')
            # List of parcel names to print
            if space == T1:
                parcelNames = [\
                    'MarsAtlas', \
                    'MarsAtlasFull', \
                    'Destrieux', \
                    'DKT', \
                    'HCP-MMP1', \
                    'VEP', \
                    'Lausanne2008-33', \
                    'Lausanne2008-60', \
                    'Lausanne2008-125', \
                    'Lausanne2008-250', \
                    'Lausanne2008-500', \
                    'GreyWhite']
                jsonCoord = {\
                        "iEEGCoordinateSystem": space, \
                        "iEEGCoordinateSystemDescription": "(x,y,z) in the scanner-based coordinate system of the pre-implantation T1w (.nii qform)", \
                        "iEEGCoordinateUnits": "mm", \
                        "iEEGCoordinateProcessingDescription": "SEEG contacts were positioned manually on the post-implantation T1w/CT with the IntrAnat software. Pre- and post-implantation images were coregistered using SPM12.", \
                        "IntendedFor": bidsT1w}
            elif space == MNI:
                parcelNames = [\
                    'IntrAnat-MarsAtlas', \
                    'MNI-MarsAtlas', \
                    'MNI-Destrieux', \
                    'MNI-DKT', \
                    'MNI-Lausanne2008-33', \
                    'MNI-Lausanne2008-60', \
                    'MNI-Lausanne2008-125', \
                    'MNI-Lausanne2008-250', \
                    'MNI-Lausanne2008-500', \
                    'MNI-AAL', \
                    'MNI-AALDilate', \
                    'MNI-AAL1_2018', \
                    'MNI-AAL2', \
                    'MNI-AAL3', \
                    'MNI-Brodmann', \
                    'MNI-BrodmannDilate', \
                    'MNI-Hammers', \
                    'MNI-HCP-MMP1', \
                    'MNI-AICHA', \
                    'MNI-JulichBrain']
                jsonCoord = {\
                        "iEEGCoordinateSystem": space, \
                        "iEEGCoordinateUnits": "mm", \
                        "iEEGCoordinateProcessingDescription": "SEEG contacts were positioned manually on the post-implantation T1w/CT with the IntrAnat software. Pre- and post-implantation images were coregistered using SPM12. Pre-implantation T1w image was normalized to IXI549Space using SPM12.", \
                        "IntendedFor": bidsT1w}
            # Write electrodes.tsv file
            with open(fileTsv, 'w') as fout:
                # Write file header
                writer = csv.writer(fout, delimiter='\t')
                # Add list of column names
                colNames = ['name', 'x', 'y', 'z', 'size', 'type'] + [p.replace('MNI-','') for p in parcelNames]
                # Print column names
                writer.writerow(colNames)
                # Print all contacts
                dict_sorted_tmp = OrderedDict(sorted(eleclabel['plots_label'].items()))
                for kk,vv in dict_sorted_tmp.iteritems():
                    # name
                    listwrite = [self.fixContactLabel(kk)]
                    # x y z
                    for i in range(3):
                        listwrite.append(format(plots[space][kk][i],'.3f'))
                    # size
                    listwrite.append('N/A')
                    # type
                    listwrite.append('N/A')
                    # List parcel names
                    for p in parcelNames:
                        if not vv[p]:
                            listwrite.append('N/A')
                        elif p == 'MarsAtlasFull':
                            listwrite.append(str(vv[p]))
                        elif isinstance(vv[p], list) and (len(vv[p]) == 2):
                            listwrite.append(vv[p][1])
                        else:
                            listwrite.append(str(vv[p]))
                    # Write line in TSV file
                    writer.writerow(listwrite)
                # Write coordsystem.json
                with open(fileJson, 'w') as fout:
                    json.dump(jsonCoord, fout, indent=2)
                # List of output files
                outfiles += [fileTsv]
                outfiles += [fileJson]
        # Return success
        print "Export BIDS TSV done."
        return [outfiles, []]
    
    
    # Convert between naming conventions (' or p for left/right electrodes)
    def fixContactLabel(self, contact_label):
        # Upper case for all the letters except from "p" that stand for ' (prime)
        contact_label = list(contact_label)
        for i in range(len(contact_label)):
            # if not kk[i].isdigit() and ((i == 0) or (kk[i] != 'p')):
            # We cannot consider that "Tp" should be kept unchanged, otherwise, "t'" is also converted to "Tp"
            # Convention: In IntrAnat, all chars are upper case and electrode names can include "'", converted to 
            # "p" in the .csv files
            if contact_label[i].isalpha():
                contact_label[i] = contact_label[i].upper()
            elif (contact_label[i] == "'"):
                contact_label[i] = 'p'
        return "".join(contact_label)
    

    def sortContacts(self, dict_contacts):
        contacts = []
        for plotName, plotCoords in dict_contacts.iteritems():
            plotName_split = plotName.split('-$&_&$-')
            ct = dict()
            ct['elecc'] = plotName_split[0]
            ct['ind'] = int(float(plotName_split[1][4:]))
            ct['name'] = plotName_split[0] + plotName_split[1][4:].zfill(2)
            ct['coord'] = plotCoords
            contacts.append(ct)
        contacts = sorted(contacts, key=lambda ct:(ct['elecc'], ct['ind']))
        contacts_sorted = [(ct['name'], ct['coord']) for ct in contacts] 
        return contacts_sorted
    

    def computeMarsatlasResection(self):
        wdi_resec = ReadDiskItem('Resection', 'NIFTI-1 image', requiredAttributes={'subject':self.brainvisaPatientAttributes['subject'], 'center':self.currentProtocol})
        di_resec = list(wdi_resec.findValues({}, None, False ))
        if len(di_resec)==0:
            print('Warning: No resection image found')
            return []
        print "TODO: computeMarsatlasResection: Use the FreeSurfer segmentation when available"
        # MarsAtlas left
        Mask_left = ReadDiskItem('Left Gyri Volume', 'Aims writable volume formats',requiredAttributes={'subject':self.brainvisaPatientAttributes['subject'], 'center':self.currentProtocol })
        diMaskleft = Mask_left.findValue(self.diskItems['T1pre'])
        if diMaskleft is None:
            print('left gyri volume failed, perform export mars atlas export contact first')
        else:
            vol_left = aims.read(diMaskleft.fileName())
        # MarsAltas right
        Mask_right = ReadDiskItem('Right Gyri Volume', 'Aims writable volume formats',requiredAttributes={'subject':self.brainvisaPatientAttributes['subject'], 'center':self.currentProtocol })
        diMaskright = Mask_right.findValue(self.diskItems['T1pre'])
        if diMaskright is None:
            print('right gyri volume failed, perform export mars atlas export contact first')
        else:
            vol_right = aims.read(diMaskright.fileName())
        # Get FreeSurfer atlas:  saved in subject/FreesurferAtlas/FreesurferAtlaspre<acq_name>
        diFreesurfer = ReadDiskItem('FreesurferAtlas', 'BrainVISA volume formats',requiredAttributes={'subject':self.brainvisaPatientAttributes['subject'], 'center':self.currentProtocol, 'modality': 'freesurfer_atlas'})
        rdiFreesurfer = list(diFreesurfer.findValues({}, None, False ))
        iFS = [i for i in range(len(rdiFreesurfer)) if ('FreesurferAtlaspre' in rdiFreesurfer[i].attributes()["acquisition"])]
        if iFS:
            rdiFreesurfer = [rdiFreesurfer[iFS[0]]]
        else:
            rdiFreesurfer = None
        # Load resection volume
        Vol_resec = aims.read(di_resec[0].fileName())
        Vox_size_resection = Vol_resec.getVoxelSize().arraydata()
        Vol_resection_mm = Vol_resec.arraydata().sum()*Vox_size_resection.prod() #choppe les trois premiers
        
        #intersection avec mars atlas label
        if diMaskleft is not None and diMaskright is not None:
            Vol_mask_tot = vol_left.arraydata()+vol_right.arraydata()
            Vol_resec_rsz = numpy.resize(Vol_resec.arraydata(), (len(Vol_mask_tot),len(Vol_mask_tot[0]),len(Vol_mask_tot[0][0]),len(Vol_mask_tot[0][0][0]))) #to make vol_resec of the same size as Vol_mask_tot
            inter_resec_mars_atlas = numpy.multiply(Vol_resec_rsz, Vol_mask_tot)
            label_resec_mars_atlas = numpy.histogram(inter_resec_mars_atlas,bins = 255, range = (0,255))
            total_label_mars_atlas = numpy.histogram(Vol_mask_tot,bins = 255)
            percent_resec_mars_atlas = numpy.divide(label_resec_mars_atlas[0],total_label_mars_atlas[0],dtype=float)*100
            non_zero_inter = numpy.add(numpy.nonzero(label_resec_mars_atlas[0][1:-1]),1)
        if rdiFreesurfer:
            vol_FS = aims.read(rdiFreesurfer[0].fileName())
            inter_resec_FS = numpy.multiply(Vol_resec.arraydata(),vol_FS.arraydata())
            label_resec_FS = numpy.histogram(inter_resec_FS,bins = 20000, range = (0,20000))
            total_label_FS = numpy.histogram(vol_FS.arraydata(),bins = 20000, range = (0,20000))
            percent_resec_FS = numpy.divide(label_resec_FS[0],total_label_FS[0],dtype=float)*100
            #we keep only value between 1 and 100 to not display the thousands freesurfer parcels...
            interesting_percent_resec_FS = numpy.where((percent_resec_FS >1) & (percent_resec_FS < 100))
        # Write into database
        wdi = WriteDiskItem('Resection Description','Resection json')
        di = wdi.findValue(self.diskItems['Resection'])
        if di is None:
            print("Can't generate files")
            return []
        fout = open(di.fullPath(),'w')
        fout.write(json.dumps({'Volume resection (mm3): ': Vol_resection_mm}))
        fout.close() 
        neuroHierarchy.databases.insertDiskItem(di, update=True )
        print "export resection info done"
        return [di.fullPath()]


    # ===== GET VOXELS WITHIN SPHERE =====
    def getSphVoxels(self, vol, pos, Nsph, sphere_size):
        if (not pos) or (len(pos) < 3) or \
           (pos[0]-Nsph[2] < 0) or (pos[0]+Nsph[2] >= vol.getSizeX()) or \
           (pos[1]-Nsph[2] < 0) or (pos[1]+Nsph[2] >= vol.getSizeY()) or \
           (pos[2]-Nsph[2] < 0) or (pos[2]+Nsph[2] >= vol.getSizeZ()):
            return []
        else:
            vox = [vol.value(pos[0]+vox_i, pos[1]+vox_j, pos[2]+vox_k) \
                for vox_k in range(-Nsph[2], Nsph[2]+1) \
                for vox_j in range(-Nsph[1], Nsph[1]+1) \
                for vox_i in range(-Nsph[0], Nsph[0]+1) \
                if math.sqrt(vox_i**2 + vox_j**2 + vox_k**2) < sphere_size]
            return [int(round(x)) for x in vox if not math.isnan(x)]


    # ===== COMPUTE PARCELS =====
    def computeParcels(self, diT1pre, plots_MNI_dict=None):
        errMsg = []
        # Get subject and protocol
        subject = diT1pre['subject']
        center = diT1pre['center']
        
        # ===== GET: NATIVE COORDINATES =====
        if plots_MNI_dict is None:
            # Get contact coordinates in T1pre coordinates
            plots = self.getPlotsT1preRef()
            if len(plots)==0:
                errMsg += ['No contact found']
                return [[], errMsg]
            # Sort by contact names
            info_plot = []
            for k,v in plots.iteritems():
                plot_name_split = k.split('-$&_&$-')
                info_plot.append((plot_name_split[0]+plot_name_split[1][4:].zfill(2),v))
            plots = sorted(info_plot, key=lambda plot_number: plot_number[0])
        else:
            plots = dict()
        # List of loaded volumes
        vol = {}
        # Get image voxel size
        info_image = diT1pre.attributes()
        # Default freesurfer coordinates
        plots_fs = copy.deepcopy(plots)
        info_fs = info_image
        
        # ===== GET: MNI COORDINATES =====
        if plots_MNI_dict is None:
            plots_MNI_dict = self.getPlotsMNIRef()
            if plots_MNI_dict is None:
                errMsg += ["No MNI coordinates available, please go back to ImageImport and compute the SPM normalization first."]
                return [[], errMsg]
            missingMni = list(set([i[0] for i in plots]) - set(plots_MNI_dict.keys()))
            if missingMni:
                errMsg += ["Recomputing MNI coordinates is required. Missing contacts: " + " ".join(missingMni)]
                return [[], errMsg]
        # Convert MNI coordinates to voxels in MNI atlas files 
        matrix_MNI_Nativ = numpy.matrix([[  -1.,    0.,    0.,   90.],[0.,   -1.,    0.,   91.],[0.,    0.,   -1.,  109.],[0.,    0.,    0.,    1.]])
        plots_MNI = {}
        for vv,kk in plots_MNI_dict.iteritems():
            inter_pos = [kk[0], kk[1], kk[2], 1]
            inter_pos = numpy.matrix(inter_pos).reshape([4,1])
            result_pos = numpy.dot(matrix_MNI_Nativ,inter_pos)
            plots_MNI.update({vv:[result_pos.tolist()[0][0],result_pos.tolist()[1][0],result_pos.tolist()[2][0]]})
        
        # ===== READ: MARSATLAS =====
        # Get surface MarsAtlas: Only to (re-)compute the volume version
        LeftGyri = ReadDiskItem('hemisphere parcellation texture','Aims texture formats',requiredAttributes={'side':'left' ,'subject':subject, 'center':center, 'parcellation_type':'marsAtlas' })
        diLeftGyri = list(LeftGyri.findValues({}, None, False))
        if diLeftGyri:
            # If there are multiple files, use the most recent one
            if (len(diLeftGyri) > 1):
                errMsg += ["Export CSV: Multiple MarsAtlas segmentations found: using the first one."]
            diLeftGyri = diLeftGyri[0]
            # Generate 3D version of the MarsAtlas surface-based atlas
            try:
                self.brainvisaContext.runProcess('2D Parcellation to 3D parcellation', Side = "Both", left_gyri = diLeftGyri)
            except Exception, e:
                errMsg += ["Could not interpolate MarsAtlas in 3D: " + repr(e)]
            # Get volume MarsAtlas
            Mask_left = ReadDiskItem('Left Gyri Volume', 'Aims writable volume formats', requiredAttributes={'subject':subject, 'center':center, 'acquisition':diLeftGyri['acquisition']})
            diMaskleft = Mask_left.findValue(diLeftGyri)
            Mask_right = ReadDiskItem('Right Gyri Volume', 'Aims writable volume formats',requiredAttributes={'subject':subject, 'center':center, 'acquisition':diLeftGyri['acquisition']})
            diMaskright = Mask_right.findValue(diLeftGyri)
            # Read volumes
            if diMaskleft and diMaskright:
                vol['MarsAtlas'] = aims.read(diMaskleft.fileName())\
                                 + aims.read(diMaskright.fileName())
            else:
                errMsg += ["MarsAtlas not found (volume)"]
        else:
            errMsg += ["MarsAtlas not found"]
        
        # ===== READ: GREY/WHITE =====
        # If MarsAtlas available: Find MarsAtlas grey-white segmentation
        if diLeftGyri:
            MaskGW_left = ReadDiskItem('Left Grey White Mask','Aims writable volume formats',requiredAttributes={'subject':subject, 'center':center, 'acquisition':diLeftGyri['acquisition']})
            MaskGW_right = ReadDiskItem('Right Grey White Mask','Aims writable volume formats',requiredAttributes={'subject':subject, 'center':center, 'acquisition':diLeftGyri['acquisition']})
            diMaskGW_left = MaskGW_left.findValue(diLeftGyri)
            diMaskGW_right = MaskGW_right.findValue(diLeftGyri)
        # Otherwise: Look for grey-white anywhere in the subject
        else:
            MaskGW_left = ReadDiskItem('Left Grey White Mask','Aims writable volume formats',requiredAttributes={'subject':subject, 'center':center})
            MaskGW_right = ReadDiskItem('Right Grey White Mask','Aims writable volume formats',requiredAttributes={'subject':subject, 'center':center})
            diMaskGW_left = list(MaskGW_left.findValues({}, None, False))
            diMaskGW_right = list(MaskGW_right.findValues({}, None, False))
            # If segmentations found: pick the first ones
            if diMaskGW_left and diMaskGW_right:
                diMaskGW_left = diMaskGW_left[0]
                diMaskGW_right = diMaskGW_right[0]
        # Read gray-white volume segmentation
        if diMaskGW_left and diMaskGW_right:
            vol['GW'] = aims.read(diMaskGW_left.fileName()) \
                      + aims.read(diMaskGW_right.fileName())
        else:
            errMsg += ['Grey-white mask not found']
            
        # ===== INITIAL SEGMENTATION =====        
        # Initial segmentation not identified
        if not diMaskGW_left:
            initSegmentation = 'N/A'
        # Initial segmentation: FreeSurfer =>  Compute FreeSurfer coordinates
        elif ('FreesurferAtlaspre' in diMaskGW_left['acquisition']):
            initSegmentation = "FreeSurfer"
            # Get FreeSurfer => Scanner-based transformation
            rdi = ReadDiskItem('Transformation to Scanner Based Referential', 'Transformation matrix', exactType=True, requiredAttributes={'modality':diMaskGW_left['modality'], 'subject':subject, 'center':center, 'acquisition':diMaskGW_left['acquisition']})
            rdiTransFS = list(rdi.findValues({}, None, False ))
            TransFS = aims.read(rdiTransFS[0].fullName())
            # Get T1pre => Scanner-based transformation
            rdi = ReadDiskItem('Transformation to Scanner Based Referential', 'Transformation matrix', exactType=True, requiredAttributes={'modality':diT1pre['modality'], 'subject':subject, 'center':center, 'acquisition':diT1pre['acquisition']})
            rdiTransT1 = list(rdi.findValues({}, None, False ))
            TransT1 = aims.read(rdiTransT1[0].fullName())
            # Compute transformation: T1pre=>FreeSurfer
            Transf = TransFS.inverse() * TransT1
            # Apply to contact coordinates
            for i in range(len(plots_fs)):
                plots_fs[i] = (plots_fs[i][0], Transf.transform(plots_fs[i][1]))
            # Get FS image information
            info_fs = diMaskGW_left.attributes()
        else:
            initSegmentation = "BrainVISA"

        # ===== READ: HIPPOCAMPUS =====
        # HIPPOCAMPUS: LEFT
        diHipLeft = ReadDiskItem('leftHippocampusNII', 'BrainVISA volume formats', requiredAttributes={'subject':subject, 'center':center})
        rdiHipLeft = list(diHipLeft.findValues({}, None, False ))
        if not rdiHipLeft:
            # If file is not properly registered in DB, try to search volume by filename directly
            lhippoFile = os.path.dirname(str(diT1pre)) + os.path.sep + 'default_analysis' + os.path.sep + 'segmentation' + os.path.sep + 'leftHippocampus' + os.path.basename(str(diT1pre))
            if os.path.isfile(lhippoFile):
                rdiHipLeft = [lhippoFile]
        # HIPPOCAMPUS: RIGHT
        diHipRight = ReadDiskItem('rightHippocampusNII', 'BrainVISA volume formats', requiredAttributes={'subject':subject, 'center':center})
        rdiHipRight = list(diHipRight.findValues({}, None, False ))
        if not rdiHipRight:
            # If file is not properly registered in DB, try to search volume by filename directly
            rhippoFile = os.path.dirname(str(diT1pre)) + os.path.sep + 'default_analysis' + os.path.sep + 'segmentation' + os.path.sep + 'rightHippocampus' + os.path.basename(str(diT1pre))
            if os.path.isfile(rhippoFile):
                rdiHipRight = [rhippoFile]
        # HIPPOCAMPUS: Merge left+right
        if rdiHipLeft and rdiHipRight:
            vol['hip'] = aims.read(str(rdiHipLeft[0])) \
                       + aims.read(str(rdiHipRight[0]))
            
        # ===== READ: LAUSANNE2008 =====
        # Get all the FreeSurfer atlases:
        diFs = ReadDiskItem('FreesurferAtlas', 'BrainVISA volume formats',requiredAttributes={'subject':subject, 'center':center})
        rdiFs = list(diFs.findValues({}, None, False))
        for name in ['FreesurferAtlaspre', 'DKT', 'VEP', 'HCP-MMP1', 'Lausanne2008-33', 'Lausanne2008-60', 'Lausanne2008-125', 'Lausanne2008-250', 'Lausanne2008-500']:
            iFile = [i for i in range(len(rdiFs)) if name in rdiFs[i].attributes()["acquisition"]]
            if iFile:
                if name == 'FreesurferAtlaspre':
                    vol['Destrieux'] = aims.read(rdiFs[iFile[0]].fileName())
                else:
                    vol[name] = aims.read(rdiFs[iFile[0]].fileName())
        if not 'Destrieux' in vol.keys():
            errMsg += ["Freesurfer atlas not found"]
            
        # ===== READ: RESECTION =====
        diResec = ReadDiskItem('Resection', 'BrainVISA volume formats', requiredAttributes={'subject':subject, 'center':center})
        rdiResec = list(diResec.findValues({}, None, False ))
        if rdiResec:
            vol['resec'] = aims.read(rdiResec[0].fileName())

        # ===== READ: MNI ATLASES =====
        # Define all atlases to generate
        files_MNI = {}
        files_MNI['IntrAnat-MarsAtlas']   = {'vol':None,                                                 'labels':labels['MarsAtlas']}
        files_MNI['MNI-MarsAtlas']        = {'vol':'MNI_Atlases/colin27_MNI_MarsAtlas_reslice.nii.gz',   'labels':labels['MarsAtlas']}
        files_MNI['MNI-AAL']              = {'vol':'MNI_Atlases/rAALSEEG12.nii.gz',                      'labels':'MNI_Atlases/rAALSEEG12_labels.txt'}
        files_MNI['MNI-AALDilate']        = {'vol':'MNI_Atlases/rAALSEEG12Dilate.nii.gz',                'labels':'MNI_Atlases/rAALSEEG12Dilate_labels.txt'}
        files_MNI['MNI-AAL1_2018']        = {'vol':'MNI_Atlases/AAL1_2018_reslice.nii.gz',               'labels':'MNI_Atlases/AAL1_2018_labels.txt'}
        files_MNI['MNI-AAL2']             = {'vol':'MNI_Atlases/AAL2_reslice.nii.gz',                    'labels':'MNI_Atlases/AAL2_labels.txt'}
        files_MNI['MNI-AAL3']             = {'vol':'MNI_Atlases/AAL3_reslice.nii.gz',                    'labels':'MNI_Atlases/AAL3_labels.txt'}
        files_MNI['MNI-BrodmannDilate']   = {'vol':'MNI_Atlases/rBrodmannSEEG3spm12.nii.gz',             'labels':None}
        files_MNI['MNI-Brodmann']         = {'vol':'MNI_Atlases/rbrodmann.nii.gz',                       'labels':None}
        files_MNI['MNI-Hammers']          = {'vol':'MNI_Atlases/rHammersSEEG12.nii.gz',                  'labels':'MNI_Atlases/rHammersSEEG12_labels.txt'}
        files_MNI['MNI-HCP-MMP1']         = {'vol':'MNI_Atlases/HCP-MMP1_on_MNI305_reslice.nii.gz',      'labels':'MNI_Atlases/HCP-MMP1_on_MNI305_labels.txt'}
        files_MNI['MNI-AICHA']            = {'vol':'MNI_Atlases/AICHA_reslice.nii.gz',                   'labels':'MNI_Atlases/AICHA_labels.txt'}
        files_MNI['MNI-JulichBrain']      = {'vol':'MNI_Atlases/JuBrain_Map_icbm_v25_lr_reslice.nii.gz', 'labels':'MNI_Atlases/JuBrain_Map_icbm_v25_lr_labels.txt'}
        files_MNI['MNI-Destrieux']        = {'vol':'MNI_Atlases/freesurfer_parcelisation_mni2.nii.gz',   'labels':labels['Destrieux']}
        files_MNI['MNI-DKT']              = {'vol':'MNI_Atlases/DKT40_reslice.nii.gz',                   'labels':labels['DKT']}
        files_MNI['MNI-Lausanne2008-33']  = {'vol':'MNI_Atlases/Lausanne2008-33.nii.gz',                 'labels':labels['Lausanne2008-33']}
        files_MNI['MNI-Lausanne2008-60']  = {'vol':'MNI_Atlases/Lausanne2008-60.nii.gz',                 'labels':labels['Lausanne2008-60']}
        files_MNI['MNI-Lausanne2008-125'] = {'vol':'MNI_Atlases/Lausanne2008-125.nii.gz',                'labels':labels['Lausanne2008-125']}
        files_MNI['MNI-Lausanne2008-250'] = {'vol':'MNI_Atlases/Lausanne2008-250.nii.gz',                'labels':labels['Lausanne2008-250']}
        files_MNI['MNI-Lausanne2008-500'] = {'vol':'MNI_Atlases/Lausanne2008-500.nii.gz',                'labels':labels['Lausanne2008-500']}
        # Load: all MNI volumes and atlases
        for atlas in files_MNI:
            if files_MNI[atlas]['vol']:
                vol[atlas] = aims.read(files_MNI[atlas]['vol'])
            if files_MNI[atlas]['labels']:
                if isinstance(files_MNI[atlas]['labels'], dict):
                    labels[atlas] = files_MNI[atlas]['labels']
                else:
                    labels[atlas] = readLabels(files_MNI[atlas]['labels'])
        # Load: IntrAnat-MarsAtlas
        vol['IntrAnat-MarsAtlas'] = aims.read('MNI_Brainvisa/Gre_2016_MNI1_L_gyriVolume.nii.gz') \
                                  + aims.read('MNI_Brainvisa/Gre_2016_MNI1_R_gyriVolume.nii.gz')
        # Load: MNI-Hippocampus
        vol['IntrAnat-hip'] = aims.read('MNI_Brainvisa/rightHippocampusGre_2016_MNI1.nii.gz') \
                            + aims.read('MNI_Brainvisa/leftHippocampusGre_2016_MNI1.nii.gz')

        # ===== PROCESS ALL CONTACTS =====
        # Define sphere size
        sphere_size = 3
        # All the computation happends in this subfunction
        plots_label, plots_name = self.computeParcelsSub(plots, plots_fs, plots_MNI, sphere_size, info_image, info_fs, vol, labels, files_MNI, initSegmentation)

        # ===== BIPOLAR MONTAGE =====
        # Scanner-based coordinates
        bip_SB = []
        for pindex in range(1,len(plots)):
            prev_elec = "".join([i for i in plots[pindex-1][0] if not i.isdigit()])
            cur_elec = "".join([i for i in plots[pindex][0] if not i.isdigit()])
            prev_ind = int("".join([i for i in plots[pindex-1][0] if i.isdigit()]))
            cur_ind = int("".join([i for i in plots[pindex][0] if  i.isdigit()]))
            if (prev_elec == cur_elec) and (prev_ind == cur_ind - 1):
                 bip_SB.append((plots[pindex-1][0] + '-' + plots[pindex][0],(plots[pindex][1]+plots[pindex-1][1])/2 ))
        # FreeSurfer coordinates
        bip_fs = []
        for pindex in range(1,len(plots_fs)):
            prev_elec = "".join([i for i in plots_fs[pindex-1][0] if not i.isdigit()])
            cur_elec = "".join([i for i in plots_fs[pindex][0] if not i.isdigit()])
            prev_ind = int("".join([i for i in plots_fs[pindex-1][0] if i.isdigit()]))
            cur_ind = int("".join([i for i in plots_fs[pindex][0] if  i.isdigit()]))
            if (prev_elec == cur_elec) and (prev_ind == cur_ind - 1):
                 bip_fs.append((plots_fs[pindex-1][0] + '-' + plots_fs[pindex][0],(plots_fs[pindex][1]+plots_fs[pindex-1][1])/2 ))
        # MNI coordinates
        bip_MNI = {}
        if not plots:  # Reconstruct missing list of contact names
            plots = []
            for k in sorted(plots_MNI.keys()):
                plots.append((k, plots_MNI[k]))
        for pindex in range(1,len(plots)):
            prev_elec = "".join([i for i in plots[pindex-1][0] if not i.isdigit()])
            cur_elec = "".join([i for i in plots[pindex][0] if not i.isdigit()])
            prev_ind = int("".join([i for i in plots[pindex-1][0] if i.isdigit()]))
            cur_ind = int("".join([i for i in plots[pindex][0] if  i.isdigit()]))
            if (prev_elec == cur_elec) and (prev_ind == cur_ind - 1):
                bip_MNI.update({plots[pindex-1][0] + '-' + plots[pindex][0]:(numpy.array(plots_MNI[plots[pindex][0]])+numpy.array(plots_MNI[plots[pindex-1][0]]))/2})
        # Search sphere size
        sphere_size = 5
        # Compute all parcels
        bip_label, bip_name = self.computeParcelsSub(bip_SB, bip_fs, bip_MNI, sphere_size, info_image, info_fs, vol, labels, files_MNI, initSegmentation)

        # ===== SAVE ELEC FILE =====
        # Get electrodes label file
        wdi = WriteDiskItem('Electrodes Labels','Electrode Label Format')
        di = wdi.findValue(diT1pre)
        if di is None:
            errMsg += ["Could not find where to save file .eleclabel"]
            return [[], errMsg]
        fileEleclabel = di.fullPath()
        # Create folder if it doesn't exist
        createItemDirs(di)
        # Write file
        fout = open(fileEleclabel,'w')
        fout.write(json.dumps({ \
            'Template' : {\
                'MarsAtlas' : str(not 'MarsAtlas' in vol.keys()), \
                'Freesurfer' : str(not 'Freesurfer' in vol.keys()), \
                'InitialSegmentation' : initSegmentation}, \
            'plots_label' : plots_label, \
            'plots_label_bipolar' : bip_label, \
            }))
        fout.close()
        # Reference file in database
        neuroHierarchy.databases.insertDiskItem(di, update=True )
        return [[fileEleclabel], errMsg]
        
    
    # ===== COMPUTE PARCES: SUBFUNCTION =====
    # Warning:  plots is a list
    #           plots_MNI is a dict
    def computeParcelsSub(self, plots, plots_fs, plots_MNI, sphere_size, info_image, info_fs, vol, labels, files_MNI, initSegmentation):
        # Size of search sphere in voxels
        if plots:
            Nsph = [int(round(sphere_size/info_image['voxel_size'][i])) for i in range(3)]
            Nsph_fs = [int(round(sphere_size/info_fs['voxel_size'][i])) for i in range(3)]
        if plots_MNI:
            Nsph_MNI = [sphere_size, sphere_size, sphere_size] #because mni has a 1 mm isotropic resolution
            
        # Process all the contacts (patient space and/or MNI space)
        plots_label = {}
        plots_name = []
        for pindex in range(max(len(plots), len(plots_MNI))):
            
            # === PATIENT SPACE ===
            if plots:
                # Contact name
                pname = plots[pindex][0]
                plots_label[pname] = {}
                # Contact positions in voxels in the different types of volumes (ScannerBased, FreeSurfer, MNI)
                pos_SB = [round(plots[pindex][1][i] / info_image['voxel_size'][i]) for i in range(3)]
                pos_fs = [round(plots_fs[pindex][1][i] / info_fs['voxel_size'][i]) for i in range(3)]
                if math.isnan(pos_SB[0]) or math.isnan(pos_fs[0]):
                    print("ERROR: Invalid coordinates for contact: " + plots[pindex][0])
                    continue
                # Positions depending on the segmentation: MarsAtlas, Grey-White
                if initSegmentation == "FreeSurfer":
                    pos_seg = pos_fs
                    Nsph_seg = Nsph_fs
                else:
                    pos_seg = pos_SB
                    Nsph_seg = Nsph
                    
                # === PROCESS: MARS ATLAS ===
                value = 0
                label = 'N/A'
                full_MA = 'N/A'
                if 'MarsAtlas' in vol.keys():
                    vox = self.getSphVoxels(vol['MarsAtlas'], pos_seg, Nsph_seg, sphere_size)
                    if vox:
                        vox = [x for x in vox if x != 0 and x !=255 and x != 100]
                        if vox:
                            value,N = Counter(vox).most_common(1)[0]
                            label = labels['MarsAtlas'][value]
                            full_count = Counter(vox).most_common()
                            full_MA = [(labels['MarsAtlas'][iLabel[0]],float(iLabel[1])/len(vox)*100) for iLabel in full_count]
                    else:
                        print pname + '-MarsAtlas: Invalid coordinates (%d,%d,%d)'%(pos_seg[0],pos_seg[1],pos_seg[2])
                plots_label[pname]['MarsAtlas'] = (value, label)
                plots_label[pname]['MarsAtlasFull'] = full_MA          
                
                # === PROCESS GREY/WHITE ===
                value = 255
                if ('GW' in vol.keys()):
                    vox = self.getSphVoxels(vol['GW'], pos_seg, Nsph_seg, sphere_size)
                    if vox:
                        vox = [x for x in vox if x !=255 and x !=0]
                        if vox:
                            most_common2,num_most_common2 = Counter(vox).most_common(1)[0]
                            value = most_common2
                        else:
                            value = vol['GW'].value(pos_seg[0],pos_seg[1],pos_seg[2])
                    else:
                        print pname + '-GW: Invalid coordinates (%d,%d,%d)'%(pos_seg[0],pos_seg[1],pos_seg[2])                
                GW_label_name = {0:'N/A',100:'GreyMatter',200:'WhiteMatter',255:'N/A'}[value]
                plots_label[pname]['GreyWhite'] = (value, GW_label_name)
                
                # === PROCESS: DESTRIEUX ===
                value = 0
                label = 'N/A'
                if ('Destrieux' in vol.keys()):
                    # Destrieux Atlas is reinterpolated in the t1pre space 
                    vox = self.getSphVoxels(vol['Destrieux'], pos_SB, Nsph, sphere_size)
                    if vox:
                        vox = [x for x in vox if x != 0 and x != 2 and x != 41] #et 2 et 41 ? left and right white cerebral matter
                        if vox:
                            value,N = Counter(vox).most_common(1)[0]
                            # Check if it's in the hippocampus
                            if (value == 53 or value == 17) and ('hip' in vol.keys()):
                                vox = self.getSphVoxels(vol['hip'], pos_SB, Nsph, sphere_size)
                                vox = [x for x in vox if x != 0 and x != 2 and x != 41]
                                try:
                                    value,N = Counter(vox).most_common(1)[0]
                                    if value == 3403:
                                        raise Exception("?????")
                                    label = labels['Destrieux'][value]
                                except:
                                    print("ERROR: Hippocampus/Amygdala surfaces not aligned with MRI")
                                    label = 'ERROR hipp/amyg surfaces not aligned with MRI'
                                    value = vol['Destrieux'].value(pos_SB[0],pos_SB[1],pos_SB[2])
                            else:
                                label = labels['Destrieux'][value]
                    else:
                        print pname + '-Destrieux: Invalid coordinates (%d,%d,%d)'%(pos_SB[0],pos_SB[1],pos_SB[2])
                plots_label[pname]['Destrieux'] = (value, label)
    
                # === PROCESS: OTHER FREESURFER ATLASES ===
                for name in ['DKT', 'VEP', 'HCP-MMP1', 'Lausanne2008-33', 'Lausanne2008-60', 'Lausanne2008-125', 'Lausanne2008-250', 'Lausanne2008-500']:
                    value = 0
                    label = 'N/A'
                    if (name in vol.keys()):
                        # These volumes seem to be saved in the T1pre/orig.mgz/Destrieux space (maybe a 1mm shift???)
                        vox = self.getSphVoxels(vol[name], pos_SB, Nsph, sphere_size)
                        if vox:
                            vox = [x for x in vox if x != 0]
                            if vox:
                                value,N = Counter(vox).most_common(1)[0]
                                try:
                                    label = labels[name][value]
                                except:
                                    label = 'N/A'
                                    print('labels ' + str(name) + ' -> ' + str(value) + ' NOT FOUND')
                                    #import pdb; pdb.set_trace()
                            else:
                                value = vol[name].value(pos_SB[0],pos_SB[1],pos_SB[2])
                        else:
                            print pname + '-' + name + ': Invalid coordinates (%d,%d,%d)'%(pos_SB[0],pos_SB[1],pos_SB[2])
                    plots_label[pname][name] = (value, label)
    
                # === PROCESS: RESECTION ===
                if ('resec' in vol.keys()):
                    vox = self.getSphVoxels(vol['resec'], pos_SB, Nsph, sphere_size)
                    valResec,N = Counter(vox).most_common(1)[0]
                    value = max(vox) 
                    per_mc = (float(N)/float(size(vox)))*100 #percentage of most_common value inside the sphere created for that contact
                    if value == 0:
                        per_mc = 100 - per_mc  #because the per_mc previously calculated was the percentage of voxels with value 0
                else:
                    value = 255
                    per_mc = 0
                label = {0:str(round(per_mc,2)), 1:str(round(per_mc,2)), 2:str(round(per_mc,2)), 255:'N/A'}[value]
                plots_label[pname]['Resection rate'] = (value, label)
            else:
                pname = None


            # === PROCESS: MNI ATLASES ===
            if plots_MNI:
                # Contact name
                if not pname:
                    pname = plots_MNI.keys()[pindex]
                    plots_label[pname] = {}
                # Get MNI coordinates
                pos_MNI = [round(plots_MNI[pname][i]) for i in range(3)]
                if math.isnan(pos_MNI[0]):
                    print("ERROR: Invalid MNI coordinates for contact: " + pname)
                    continue
                # Process all MNI volumes
                for atlas in files_MNI:
                    value = 0
                    label = 'N/A'
                    vox = self.getSphVoxels(vol[atlas], pos_MNI, Nsph_MNI, sphere_size)
                    if vox:
                        vox = [x for x in vox if x != 0]
                        if vox:
                            most_common,N = Counter(vox).most_common(1)[0]
                            value = int(round(most_common))
                            if atlas == 'MNI-Brodmann':
                                if pos_MNI[0]>90:
                                    label = "%d" % (value)
                                else:
                                    label = "%d" % (value+100)
                            elif atlas == 'MNI-BrodmannDilate':
                                if value>48:
                                    label = "%d" % (value-48+100)
                                else:
                                    label = "%d" % (value)
                            else:
                                if value in labels[atlas].keys():
                                    label = labels[atlas][value]
                                else:
                                    print atlas + '-' + pname + ": Missing label: " + str(value)
                                    label = 'N/A'
                    else:
                        print pname + '-' + atlas + ': Invalid coordinates (%d,%d,%d)'%(pos_MNI[0],pos_MNI[1],pos_MNI[2])
                    plots_label[pname][atlas] = (value, label)

            # Add to list of exported plot names
            plots_name += [pname]
        return plots_label, plots_name
        

    def getPlotsT1preRef(self):
        """Return a dictionary {'ElectrodeName-$&_&$-PlotName':[x,y,z], ...} where x,y,z is in the T1pre native referential"""
        return dict((el['name']+'-$&_&$-'+plotName, el['transf'].transform(plotCoords)) for el in self.electrodes for plotName, plotCoords in getPlotsCenters(el['elecModel']).iteritems())

    def getPlotsT1preScannerBasedRef(self):
        """Return a dictionary {'ElectrodeName-$&_&$-PlotName':[x,y,z], ...} where x,y,z is in the T1pre scanner-based referential"""
        transfo = self.t1pre2ScannerBased()
        if not transfo:
            return None
        return dict((key, transfo.transform(coords)) for key, coords in self.getPlotsT1preRef().iteritems())

    def getPlotsAnyReferential(self, referential):
        """Return a dictionary {'ElectrodeName-$&_&$-PlotName':[x,y,z], ...} where x,y,z is in the provided referential (must be available in referential converter self.refConv)"""
        if not self.refConv.isRefAvailable(referential):
            print "Trying to convert to unknown referential %s"%repr(referential)
            return None
        return dict((key, self.refConv.real2AnyRef(coords, referential)) for key, coords in self.getPlotsT1preRef().iteritems())

    def getPlotsMNIRef(self, isShortName=True):
        """Return a dictionary {'ElectrodeName-$&_&$-PlotName':[x,y,z], ...} where x,y,z is in the MNI referential"""

        # Get elecimplant file
        rdi = ReadDiskItem( 'Electrode implantation', 'Electrode Implantation format',requiredAttributes={'subject':self.brainvisaPatientAttributes['subject'], 'center':self.currentProtocol})
        impl = list(rdi.findValues({},None,False))
        if len(impl) == 0:
            #QtGui.QMessageBox.critical(self, "Error", "Can't find implantation, you have to add electrodes and save your modifications first.")
            print "ERROR: Can't find implantation, you have to add electrodes and save your modifications first."
            return None
        elif not os.path.exists(str(impl[0])):
            #QtGui.QMessageBox.critical(self, "Error", "Can't find implantation file. Update your BrainVISA database.")
            print "ERROR: Can't find implantation file. Update your BrainVISA database."
            return None
        # Read implantation file
        filein = open(str(impl[0]), 'rb')
        try:
            dic_impl = json.loads(filein.read())
        except:
            filein.close()
            filein = open(str(impl[0]), 'rU')
            dic_impl = pickle.load(filein)
            filein.close()
        # Return existing coordinates
        if 'plotsMNI' in dic_impl.keys():
            dict_fileMNI = dict(dic_impl['plotsMNI'])
            # Convert keys to long plot names: eg. "A01" => "A-$&_&$-Plot1"
            if not isShortName:
                plots_MNI = dict()
                for el in self.electrodes:
                    for plotName, plotCoords in getPlotsCenters(el['elecModel']).iteritems():
                        contactName = el['name'] + plotName[4:].zfill(2)
                        if contactName in dict_fileMNI:
                            plots_MNI[el['name']+'-$&_&$-'+plotName] = dict_fileMNI[contactName]
                return plots_MNI
            else:
                return dict_fileMNI
        # If not available: Compute MNI coordinates
        else:
            print "WARNING: MNI position not available. Compute SPM normalization first."
            return None
        
        
    def computeMniCoord(self):
        """Compute MNI coordinates for the all the SEEG contacts, save them in database"""

        # Get contacts coordinates
        plots = self.getPlotsT1preScannerBasedRef()
        coords = [plots[k] for k in sorted(plots.keys())]
        # Get MNI transformation field
        field = self.getT1preMniTransform()
        if field is None:
            print "ERROR: MNI deformation field not found."
            return [None,[]]
        
        # Save coordinates in a .csv file, to be passed to MATLAB/SPM
        tmpOutput = getTmpFilePath('csv')
        arr = numpy.asarray(coords)
        numpy.savetxt(tmpOutput, arr, delimiter=",")
        # Run SPM normalization
        errMsg = matlabRun(spm12_normalizePoints % ("'"+self.spmpath+"'", field, tmpOutput, tmpOutput))
        if errMsg:
            print errMsg
            return [None,[]]
        # Read MNI coordinates from output .csv file
        out = numpy.loadtxt(tmpOutput, delimiter=",")
        os.remove(tmpOutput)
        if numpy.array_equal(out, arr):
            print "ERROR: Points to MNI: Result is identical to input"
            return [None,[]]
        if (out.shape != arr.shape):
            print "ERROR: Points to MNI: Result (%s) has not the same number of elements as input (%s)"%(repr(out),repr(arr))
            return [None,[]]
        mniCoords = out.tolist()
        if mniCoords is None:
            return [None,[]]
        # Convert back to an electrode dictionnary
        plotsMNI = dict(zip(sorted(plots.keys()), mniCoords))

        # Get elecimplant file
        rdi = ReadDiskItem( 'Electrode implantation', 'Electrode Implantation format', requiredAttributes={'center':self.brainvisaPatientAttributes['center'],'subject':self.brainvisaPatientAttributes['subject']})
        di = rdi.findValues({},None,False)
        ldi = list(rdi.findValues({},None,False))
        # If there is already a elecimplant file: load the existing file
        previous_data = dict()
        if (len(ldi) > 0):
            if (os.path.exists(str(ldi[0]))):
                filein = open(str(ldi[0]), 'rb')
                try:
                    previous_data = json.loads(filein.read())
                except:
                    filein.close()
                    filein = open(str(ldi[0]), 'rU')
                    previous_data = pickle.load(filein)
                filein.close()
    
        # Add new coordinates to the existing json
        info_plotMNI= []
        for k,v in plotsMNI.iteritems():
            plot_name_split = k.split('-$&_&$-')
            info_plotMNI.append((plot_name_split[0]+plot_name_split[1][4:].zfill(2),v))
        plotMNI_sorted = sorted(info_plotMNI, key=lambda plot_number: plot_number[0])
        previous_data.update({'plotsMNI':info_plotMNI})
    
        # Resave as json file
        fout = open(str(ldi[0]),'w')
        fout.write(json.dumps(previous_data))
        fout.close()
        neuroHierarchy.databases.insertDiskItem([x for x in di][0], update=True )
        print ".elecimplant saved with MNI"

        return [plotsMNI, [ldi[0].fullPath() + u"\xa0(MNI)"]]


    def saveTXT(self, contacts=None, path=None, pathPos=None, pathName=None):
        """ Saves two txt files electrode_Name.txt and electrode_Pos.txt. Path should be supplied as /myPath/electrode.txt
           or as an alternative, both file path should be supplied as pathPos and pathName
           contacts should be a dictionary {'ElectrodeName-$&_&$-PlotName:[x,y,z],...} in the appropriate referential
         """
        # Get a valid path
        if not path and not (pathPos is not None and pathName is not None):
            path = str(QtGui.QFileDialog.getSaveFileName(self, "Save TXT files", "", "Electrode implantation TXT files (*.txt)"))
            if not path:
                return None
        if not contacts:
            return None
        if path is not None:
            path = os.path.splitext(path)[0]
            pathName = path+'_Name.txt'
            pathPos = path+"_Pos.txt"
        fileName = open(pathName, 'wb')
        filePos = open(pathPos, 'wb')
        # Iteration over electrodes
        for contactName in sorted(contacts.keys(), key=natural_keys):
            (elName, plotName) = contactName.split('-$&_&$-')
            # Name of each contact is name of electrode (p for prime) + number of the plot (electrode A' contact 5 is "Ap5")
            fileName.write(elName.replace('\'', 'p') + plotName.replace('Plot','') + "\n")
            filePos.write("%f %f %f\n"%tuple(contacts[contactName]))
        fileName.close()
        filePos.close()
        return pathPos

    def savePTS(self, contacts=None, path=None):
        """ Save a PTS file with all contacts coordinates.
           contacts parameter should be a dictionary {'ElectrodeName-$&_&$-PlotName:[x,y,z],...} in the appropriate referential
        """
        # Get a valid path
        if not path:
            path = str(QtGui.QFileDialog.getSaveFileName(self, "Save PTS file (%s referential)"%referential, "", "PTS Electrode implantation (*.pts)"))
            if not path:
                return None
        if contacts is None:
            return None
    
        plots=[]
        for contactName in sorted(contacts.keys(), key=natural_keys):
            coords = contacts[contactName]
            (elName, plotName) = contactName.split('-$&_&$-')
            # Name of each contact is name of electrode (prime ' replaced by the letter p) + number of the plot (electrode A' contact 5 is "Ap5")
            plots.append("%s\t%.2f\t%.2f\t%.2f\t0\t0\t0\t2\t2.0\n"%(elName.replace('\'', 'p') + plotName.replace('Plot',''), coords[0], coords[1], coords[2]))
    
        fileout = open(path, 'wb')
        fileout.write("ptsfile\n1\t1\t1\n%s\n"%str(len(plots)))
        for p in plots:
            fileout.write(p)
    
        fileout.close()
        return path

    # Open configuration dialog
    def configureColors(self):
        self.bipoleSEEGColors=bipoleSEEGColors(self)
        self.bipoleSEEGColors.show()
    
    
    # Set number of windows
    def setWindowNumber(self, nWin):
        if nWin == 1:
            self.windowCombo2.setVisible(False)
            self.windowContainer2.setVisible(False)
            self.windowCombo3.setVisible(False)
            self.windowContainer3.setVisible(False)
            self.windowCombo4.setVisible(False)
            self.windowContainer4.setVisible(False)
            self.verticalLayout_wright.parent().setStretchFactor(self.verticalLayout_wright, 0)
        elif nWin == 2:
            self.windowCombo2.setVisible(True)
            self.windowContainer2.setVisible(True)
            self.windowCombo3.setVisible(False)
            self.windowContainer3.setVisible(False)
            self.windowCombo4.setVisible(False)
            self.windowContainer4.setVisible(False)
            self.verticalLayout_wright.parent().setStretchFactor(self.verticalLayout_wright, 1)
        elif nWin == 4:
            self.windowCombo2.setVisible(True)
            self.windowContainer2.setVisible(True)
            self.windowCombo3.setVisible(True)
            self.windowContainer3.setVisible(True)
            self.windowCombo4.setVisible(True)
            self.windowContainer4.setVisible(True)
            self.verticalLayout_wright.parent().setStretchFactor(self.verticalLayout_wright, 1)
        # Update window orientation (to force redraw)
        for i in range(nWin):
            self.updateWindowOrient(self.wins[i])
            
            
    # Update all anatomist windows
    def updateAllWindows(self, isUpdateOrient=False):
        self.updateWindow(0, self.windowCombo1.currentText(), isUpdateOrient)
        self.updateWindow(1, self.windowCombo2.currentText(), isUpdateOrient)
        self.updateWindow(2, self.windowCombo3.currentText(), isUpdateOrient)
        self.updateWindow(3, self.windowCombo4.currentText(), isUpdateOrient)

    # Update a window content
    def updateWindow(self, winId, key, isUpdateOrient=False):
        key = str(key) # Can be a QString !
        if key not in self.windowContent:
            return #when QT interface but, there was a variable generated in "frame" and we were not able to delete it, then we were not able to load another patient, there was leftover from previous patient
        w = self.wins[winId]
        # Loop on ll the options that should be displayed
        for obj in self.dispObj:
            # Is object already displayed
            isInWindow = (isinstance(self.dispObj[obj], list) and self.dispObj[obj] and (w not in self.dispObj[obj][0].getWindows())) or \
                         (not isinstance(self.dispObj[obj], list) and (w not in self.dispObj[obj].getWindows()))
            
            if (obj in self.windowContent[key]):
                w.addObjects(self.dispObj[obj])
            else:
                w.removeObjects(self.dispObj[obj])#CURRENT
        # Redraw figure
        if isUpdateOrient:
            self.updateWindowOrient(w)
   
    # Update window orient
    def updateWindowOrient(self, w):
            viewType = w.getInternalRep().viewType()
            if (viewType == 0):
                w.getInternalRep().muteAxial()
                w.getInternalRep().muteOblique()
            elif (viewType == 1):
                w.getInternalRep().muteCoronal()
                w.getInternalRep().muteAxial()
            elif (viewType == 2):
                w.getInternalRep().muteCoronal()
                w.getInternalRep().muteSagittal()
            elif (viewType == 3):
                w.getInternalRep().muteAxial()
                w.getInternalRep().muteCoronal()
            elif (viewType == 4):
                w.getInternalRep().muteAxial()
                w.getInternalRep().mute3D()
            
             
    # Update combo boxes
    def updateComboboxes(self, default1=None, default2=None, default3=None, default4=None):
        # Disable combobox callbacks
        isAlreadyBlocked = self.windowCombo1.signalsBlocked()
        if not isAlreadyBlocked:
            self.windowCombo1.blockSignals(True)
            self.windowCombo2.blockSignals(True)
            self.windowCombo3.blockSignals(True)
            self.windowCombo4.blockSignals(True)
        # Empty lists
        self.windowCombo1.clear()
        self.windowCombo2.clear()
        self.windowCombo3.clear()
        self.windowCombo4.clear()
        # Get the list of items
        items = sorted(self.windowContent.keys())
        if not items:
            return
        # Move "Cortex" and "Cortex + ..." at the end
        iCortex = [i for i in range(len(items)) if items[i] == 'Cortex']
        if iCortex:
            items.append(items[iCortex[0]])
            del items[iCortex[0]]
        iCortex = [i for i in range(len(items)) if items[i] == 'Cortex + head']
        if iCortex:
            items.append(items[iCortex[0]])
            del items[iCortex[0]]
        # Move "FreeSurfer Atlas" at the end
        iFreeSurfer = [i for i in range(len(items)) if (('FreeSurfer Atlas' in items[i]) or ('Lausanne' in items[i]) or ('DKT' in items[i]) or ('VEP' in items[i]) or ('HCP-MMP1' in items[i]))]
        for i in iFreeSurfer:
            items.append(items[i])
        for i in reversed(iFreeSurfer):
            del items[i]
        # Set new list of items
        self.windowCombo1.addItems(items)
        self.windowCombo2.addItems(items)
        self.windowCombo3.addItems(items)
        self.windowCombo4.addItems(items)
        # Set defaults
        if default1:
            self.windowCombo1.setCurrentIndex(max(self.windowCombo1.findText(default1),0))
        if default2:
            self.windowCombo2.setCurrentIndex(max(self.windowCombo2.findText(default2),0))
        if default3:
            self.windowCombo3.setCurrentIndex(max(self.windowCombo3.findText(default3),0))
        if default4:
            self.windowCombo4.setCurrentIndex(max(self.windowCombo4.findText(default4),0))   
        # Enable combobox callbacks
        if not isAlreadyBlocked:
            self.windowCombo1.blockSignals(False)
            self.windowCombo2.blockSignals(False)
            self.windowCombo3.blockSignals(False)
            self.windowCombo4.blockSignals(False)


    def updateElectrodeView(self, checkStatus=None):
        """Sets the view to electrode referential or back to native T1 referential
        checkStatus is Qt.CheckState"""
        if checkStatus is None:
            checkStatus = self.electrodeRefCheck.checkState()
            isRestoreNative = False
        else:
            isRestoreNative = True
        if checkStatus == QtCore.Qt.Checked:
            if self.electrodeList.count() > 0 and self.electrodeList.currentRow() <= len(self.electrodes):
                #Change to electrode referential
                el = self.currentElectrode()
                if el is None:
                    return
                if 'ref' in el:
                    self.setWindowsReferential(el['ref'])
                    self.electrodeGo(electrode = el)
        else:
            # Change back to T1pre native ref
            if isRestoreNative:
                self.setWindowsReferential()
            # Jump to position of the electrode
            if self.electrodeList.count() > 0 and self.electrodeList.currentRow() <= len(self.electrodes):
                el = self.currentElectrode()
                self.electrodeGo(electrode = el)
        # Update window orientation (to force redraw)
        for i in range(4):
            self.updateWindowOrient(self.wins[i])
            

    def clippingUpdate(self):
        for i in range(3):
            if self.Clipping_checkbox.isChecked():
                self.a.execute('WindowConfig',windows = [self.wins[i]],clipping=2,clip_distance=5.)
            else:
                self.a.execute('WindowConfig',windows = [self.wins[i]],clipping=0)

  
    def makeFusion(self):
        #get objects
        Text_win1 = self.windowCombo1.currentText()
        Text_win2 = self.windowCombo2.currentText()
        
        for obj in self.dispObj:
            if obj in self.windowContent[str(Text_win1)][0]:
                obj1 = obj
        for obj in self.dispObj:
            if obj in self.windowContent[str(Text_win2)][0]:
                obj2 = obj
        
        if 'obj1' in locals() and 'obj2' in locals():
            fusion_obj = self.a.fusionObjects((self.dispObj[obj1], self.dispObj[obj2]), method='Fusion2DMethod')
            self.a.addObjects(fusion_obj, self.wins[1])
            # Add the fusion in the disObj and the windowCombo
            self.dispObj[obj1+'+'+obj2] = fusion_obj
            # Update list of available items in the combo boxes
            self.windowContent.update({obj1+'+'+obj2:[obj1+'+'+obj2,'electrodes']})
            self.updateComboboxes(Text_win1, obj1+'+'+obj2)
            # Update anatomist windows
            self.updateAllWindows()
        else:
            print "one of the image is not recognized"
            return


    def deleteMarsAtlasFiles(self):
        # No subject loaded
        if not self.brainvisaPatientAttributes or not self.brainvisaPatientAttributes['subject']:
            return
        # Ask for confirmation
        rep = QtGui.QMessageBox.warning(self, u'Confirmation', u"<font color='red'><b>ATTENTION</b><br/>You are gonna delete MarsAtlas files, are you sure?<br/><b>DELETE MARSATLAS ?</b></font>", QtGui.QMessageBox.Yes | QtGui.QMessageBox.No, QtGui.QMessageBox.No)
        if not (rep == QtGui.QMessageBox.Yes):
            return
        # Get files to delete
        atlas_di = ReadDiskItem('hemisphere marsAtlas parcellation texture', 'aims Texture formats', requiredAttributes={ 'regularized': 'false','subject':self.brainvisaPatientAttributes['subject'], 'center':self.currentProtocol })
        atlas_di_list = list(atlas_di._findValues({}, None, False ))
        Mask_left = ReadDiskItem('Left Gyri Volume', 'Aims writable volume formats',requiredAttributes={'subject':self.brainvisaPatientAttributes['subject'], 'center':self.currentProtocol })
        diMaskleft = list(Mask_left.findValues({}, None, False ))
        Mask_right = ReadDiskItem('Right Gyri Volume', 'Aims writable volume formats',requiredAttributes={'subject':self.brainvisaPatientAttributes['subject'], 'center':self.currentProtocol })
        diMaskright = list(Mask_right.findValues({}, None, False ))
        # Delete them
        if len(atlas_di_list)>0:
            for i,infoi in enumerate(atlas_di_list):
                removeFromDB(infoi.fullPath(), neuroHierarchy.databases.database(infoi.get("_database")))
        if len(diMaskleft)>0:
            removeFromDB(diMaskleft[0].fullPath(), neuroHierarchy.databases.database(diMaskleft[0].get("_database")))
        if len(diMaskright)>0:
            removeFromDB(diMaskright[0].fullPath(), neuroHierarchy.databases.database(diMaskright[0].get("_database")))
        print("MarsAtlas files deleted.")


#     def generateResection(self):
# 
#         #look for T1 pre and T1 postop
#         T1 = ReadDiskItem('Raw T1 MRI', 'aims readable volume formats',requiredAttributes={'subject':self.brainvisaPatientAttributes['subject'], 'center':self.currentProtocol, 'normalized':'no' })
#         diT1 = list(T1.findValues({}, None, False ))
#         for t in diT1:
#             if 'T1pre' in t.attributes()['acquisition']:
#                 T1pre = t.fullPath()
#                 try:
#                     SB_transf = t.attributes()['SB_Transform']
#                 except:
#                     SB_transf = t.attributes()['transformations'][0]
#         
#             elif 'postOp' in t.attributes()['acquisition']:
#                 T1postop = t.fullPath()
#         
#         if not T1pre:
#             print('don t find the t1pre')
#             return
#         
#         if 'T1postop' not in locals():
#             print('don t find the t1postop')
#             print('look for a CTpostop')
#             CTs = ReadDiskItem('CT', 'aims readable volume formats',requiredAttributes={'subject':self.brainvisaPatientAttributes['subject'], 'center':self.currentProtocol })
#             diCTs = list(CTs.findValues({},None,False))
#             for t in diCTs:
#                 if 'postOp' in t.attributes()['acquisition']:
#                     CTpostop = t.fullPath()
#                     method = 'CT'
#             if 'CTpostop' not in locals():
#                 print("can't find a CTpostop either")
#                 return
#         else:
#             method = 'T1'
#         
#         T1pre = None
#         T1postop = None
#         
#         id_pre = [x for x in range(len(diT1)) if 'T1pre' in str(diT1[x])]
#         T1pre = diT1[id_pre[0]]
#         
#         if method == 'T1':          
#             id_postop = [x for x in range(len(diT1)) if 'postOp' in str(diT1[x])]
#             T1postop = diT1[id_postop[0]]
#             
#         elif method == 'CT': 
#             id_ctpostop = [x for x in range(len(diCTs)) if 'postOp' in str(diCTs[x])]
#             CTpostop = diCTs[id_ctpostop[0]]
#         
#         if method == 'T1':
#             wdiTransform = ReadDiskItem('Transform Raw T1 MRI to another image', 'Transformation matrix', exactType=True, requiredAttributes = {'subject':self.brainvisaPatientAttributes['subject'], 'center':self.currentProtocol }) #pourquoi la suite marche pas ?, requiredAttributes = {'modalityTarget':T1pre.attributes()['modality'], 'acquisitionTarget':T1pre.attributes()['acquisition']}
#             diTransform = list(wdiTransform.findValues({}, None, False ))
#             
#             for t in diTransform:
#                 if t.attributes()['modality'] == 't1mri' and 'postOp' in t.attributes()['acquisition']:
#                     trmpostop_to_pre = t
#             
#             wdiTransform2 = ReadDiskItem('Transformation to Scanner Based Referential', 'Transformation matrix', exactType=True, requiredAttributes = {'subject':self.brainvisaPatientAttributes['subject'], 'center':self.currentProtocol })
#             diTransform2 = wdiTransform2.findValues({}, None, False )
#             
#             for t in diTransform2:
#                 if t.attributes()['modality'] == 't1mri' and 'postOp' in t.attributes()['acquisition']:
#                     trmpostop_to_SB = t
#                 if t.attributes()['modality'] == 't1mri' and 'T1pre' in t.attributes()['acquisition']:
#                     trmpre_to_SB = t
# 
#         if method == 'CT':
#             wdiTransform = ReadDiskItem('Transform CT to another image', 'Transformation matrix', exactType=True, requiredAttributes = {'subject':self.brainvisaPatientAttributes['subject'], 'center':self.currentProtocol }) #pourquoi la suite marche pas ?, requiredAttributes = {'modalityTarget':T1pre.attributes()['modality'], 'acquisitionTarget':T1pre.attributes()['acquisition']}
#             diTransform = list(wdiTransform.findValues({}, None, False ))
#             
#             for t in diTransform:
#                 if t.attributes()['modality'] == 'ct' and 'postOp' in t.attributes()['acquisition']:
#                     trmpostop_to_pre = t
#                  
#             wdiTransform2 = ReadDiskItem('Transformation to Scanner Based Referential', 'Transformation matrix', exactType=True, requiredAttributes = {'subject':self.brainvisaPatientAttributes['subject'], 'center':self.currentProtocol })
#             diTransform2 = wdiTransform2.findValues({}, None, False )
#             
#             for t in diTransform2:
#                 if t.attributes()['modality'] == 'ct' and 'postOp' in t.attributes()['acquisition']:
#                     trmpostop_to_SB = t
#                 if t.attributes()['modality'] == 't1mri' and 'T1pre' in t.attributes()['acquisition']:
#                     trmpre_to_SB = t
# 
#         import copy
#         trmpre_to_SBinvpath = trmpre_to_SB.fullPath().split('/')
#         trmpre_to_SBinvpath[-1] = 'inv'+trmpre_to_SBinvpath[-1]
#         trmpostop_to_pre_path = copy.deepcopy(trmpre_to_SBinvpath)
#         trmpostop_to_pre_path[-1] = 'postop_to_pre.trm'
#         trmpre_to_postop_path = copy.deepcopy(trmpre_to_SBinvpath)
#         trmpre_to_postop_path[-1] = 'pre_to_postop.trm'
#         trmpre_to_SBinvpath = '/'.join(trmpre_to_SBinvpath)
#         trmpostop_to_pre_path = '/'.join(trmpostop_to_pre_path)
#         trmpre_to_postop_path = '/'.join(trmpre_to_postop_path)
# 
#         ret = subprocess.call(['AimsInvertTransformation','-i',trmpre_to_SB.fullPath(),'-o', trmpre_to_SBinvpath])
#         ret = subprocess.call(['AimsComposeTransformation', '-o',trmpostop_to_pre_path, trmpre_to_SBinvpath, trmpostop_to_pre.fullPath(), trmpostop_to_SB.fullPath()])
#         ret = subprocess.call(['AimsInvertTransformation','-i',trmpostop_to_pre_path,'-o',trmpre_to_postop_path])
#         
#         print 'select the center of the resection'
#         center_seg = QtGui.QMessageBox(self)
#         center_seg.setText("Enter the center (approximately) of the resection")
#         center_seg.setWindowTitle("Resection center")
#         center_seg.setWindowModality(QtCore.Qt.NonModal)
#         center_seg.show()
#         center_seg.exec_()
#         
#         #convert brainvisa voxel position into spm_voxel position
#         ResecCenterCoord = self.positionPreRef()
#         print ResecCenterCoord
#         
#         if method == 'T1':
#             transfo_pre_to_postop = aims.read(trmpre_to_postop_path).toMatrix()
#             
#             brainMask = ReadDiskItem('Brain Mask', 'aims readable volume formats',requiredAttributes={'subject':self.brainvisaPatientAttributes['subject'], 'center':self.currentProtocol })
#             diBrain = list(brainMask.findValues({}, None, False ))
#             
#             id_pre = [x for x in range(len(diBrain)) if 'T1pre' in str(diBrain[x])]
#             id_postop = [x for x in range(len(diBrain)) if 'postOp' in str(diBrain[x])]
#             
#             #NOT NECESSARY FOR 'ResectionStart' AND TAKES TIME
#             morphologist = getProcessInstance('morphologist')
#             #morphologist.executionNode().PrepareSubject.setSelected(True)
#             morphologist.executionNode().BiasCorrection.setSelected(True)
#             morphologist.executionNode().PrepareSubject.StandardACPC.setSelected(False)
#             morphologist.executionNode().PrepareSubject.Normalization.setSelected(True)
#             morphologist.executionNode().HistoAnalysis.setSelected(True)
#             morphologist.executionNode().BrainSegmentation.setSelected(True)
#             morphologist.executionNode().Renorm.setSelected(False)
#             morphologist.executionNode().SplitBrain.setSelected(False)
#             morphologist.executionNode().TalairachTransformation.setSelected(False)
#             morphologist.executionNode().HeadMesh.setSelected(False)
#             morphologist.executionNode().HemispheresProcessing.setSelected(False)
#             morphologist.executionNode().SulcalMorphometry.setSelected(False)
#             
#             
#             if not id_pre:             
#                 self.brainvisaContext.runProcess(morphologist, t1mri = T1pre, perform_normalization = True, anterior_commissure = '',\
#                         posterior_commissure = '', interhemispheric_point = '', left_hemisphere_point = '', perform_sulci_recognition = False)
#                 
#             if not id_postop:
#                 self.brainvisaContext.runProcess(morphologist, t1mri = T1postop, perform_normalization = True, anterior_commissure = '',\
#                     posterior_commissure = '', interhemispheric_point = '', left_hemisphere_point = '', perform_sulci_recognition = False)
#                 
#             self.resectionStart(trmpostop_to_pre_path, ResecCenterCoord, method = 'T1')
#             
#             #self.brainvisaContext.runInteractiveProcess(lambda x='',trm=trmpostop_to_pre_path,resec_coord=ResecCenterCoord, methodo=method:self.resectionStart(trm,resec_coord,methodo), 
#             #                                            morphologist, t1mri = T1postop, perform_normalization = False)
#                                 #anterior_commissure = Ac_vector_postop[0:3].T.tolist()[0],\
#                                 #posterior_commissure = Pc_vector_postop[0:3].T.tolist()[0], interhemispheric_point = Ih_vector_postop[0:3].T.tolist()[0], left_hemisphere_point = Lh_postop.tolist()[0], perform_sulci_recognition = False)
#             #self.resectionStart(trmpostop_to_pre_path, ResecCenterCoord, method = 'T1')
#         if method == 'CT':
#             self.resectionStart(trmpostop_to_pre_path, ResecCenterCoord,method = 'CT')

  
#     def resectionStart(self,trm_postop_to_pre, resec_coord,method = 'T1'):
#         wdi_resec = WriteDiskItem('Resection', 'NIFTI-1 image')
#         di_resec = wdi_resec.findValue({'subject':self.brainvisaPatientAttributes['subject'], 'center':self.currentProtocol, 'acquisition':'Resection'})
#         
#         #if the path doesn't exist, create it
#         if not os.path.exists(os.path.dirname(di_resec.fullPath())):
#             os.makedirs( os.path.dirname(di_resec.fullPath()))
#         
#         
#         if method == 'T1':
#             brainMask = ReadDiskItem('Brain Mask', 'aims readable volume formats',requiredAttributes={'subject':self.brainvisaPatientAttributes['subject'], 'center':self.currentProtocol })
#             diBrain = list(brainMask.findValues({}, None, False ))
#             
#             id_pre = [x for x in range(len(diBrain)) if 'T1pre' in str(diBrain[x])]
#             id_postop = [x for x in range(len(diBrain)) if 'postOp' in str(diBrain[x])]
#             
#             fullname_postop = diBrain[id_postop[0]].fullPath()
#             fullpost_split = fullname_postop.split('/')
#             fullpost_split[-1] = 'brainpostop_on_pre.nii'
#             fullpost = '/'.join(fullpost_split)
#             
#             ret = subprocess.call(['AimsResample', '-i', diBrain[id_postop[0]].fullPath(), '-m', trm_postop_to_pre, '-o', fullpost, '-t', 'n', '-r', diBrain[id_pre[0]].fullPath()])
#             ret = subprocess.call(['AimsMorphoMath', '-m', 'dil', '-i', fullpost, '-o', fullpost, '-r', '1.2'])
#             ret = subprocess.call(['AimsMorphoMath', '-m', 'dil', '-i', diBrain[id_pre[0]].fullPath(), '-o', diBrain[id_pre[0]].fullPath(), '-r', '0.75'])
#             ret = subprocess.call(['AimsMask', '-i',  fullpost, '-m', diBrain[id_pre[0]].fullPath(), '-o', fullpost])
#             ret = subprocess.call(['AimsLinearComb', '-i', diBrain[id_pre[0]].fullPath(), '-j', fullpost, '-c', '-1', '-o', di_resec.fullPath()])
#             ret = subprocess.call(['AimsThreshold', '-i', di_resec.fullPath(), '-m', 'ge', '-t', '250' , '-o', di_resec.fullPath()])
#             ret = subprocess.call(['AimsMorphoMath', '-m', 'ero', '-i', di_resec.fullPath(), '-o', di_resec.fullPath(), '-r', '1']) #to remove small structures out of interest
#             ret = subprocess.call(['AimsConnectComp', '-i', di_resec.fullPath(), '-o', di_resec.fullPath(), '-c', '6'])
#                 
#         if method == 'CT':
#             brainMask = ReadDiskItem('Brain Mask', 'aims readable volume formats',requiredAttributes={'subject':self.brainvisaPatientAttributes['subject'], 'center':self.currentProtocol })
#             diBrain = list(brainMask.findValues({}, None, False ))
#             
#             id_pre = [x for x in range(len(diBrain)) if 'T1pre' in str(diBrain[x])]
#             
#             CTs = ReadDiskItem('CT', 'aims readable volume formats',requiredAttributes={'subject':self.brainvisaPatientAttributes['subject'], 'center':self.currentProtocol })
#             diCTs = list(CTs.findValues({},None,False))
#             for t in diCTs:
#                 if 'postOp' in t.attributes()['acquisition']:
#                     CTpostop = t.fullPath()
# 
#             ret = subprocess.call(['AimsResample', '-i', str(CTpostop), '-m', trm_postop_to_pre, '-o',  di_resec.fullPath(), '-t', 'c', '-r',diBrain[id_pre[0]].fullPath()])
#             
#             ret = subprocess.call(['AimsThreshold', '-i',di_resec.fullPath(),'-m', 'be','-t','-20','-u','25','-b','-o', di_resec.fullPath()])
#             ret = subprocess.call(['cartoLinearComb.py','-i',str(di_resec.fullPath()),'-i',str(diBrain[id_pre[0]]),'-o',str(di_resec.fullPath()),'-f','(I1/32767*I2/255)'])
#             
#             #apply brainmask
#             ret = subprocess.call(['AimsMorphoMath', '-m','ope', '-i', di_resec.fullPath(), '-o', di_resec.fullPath(), '-r', '2'])
# 
#             ret = subprocess.call(['AimsConnectComp', '-i',di_resec.fullPath(),'-o',di_resec.fullPath(),'-c','6',])
#         
#         vol_connectcomp = aims.read(di_resec.fullPath())
#         #a small sphere here as for plots:
#         sphere_size = 4
#         Nsph = [int(round(sphere_size/vol_connectcomp.getVoxelSize()[i])) for i in range(0,3)]
#         voxel_within_sphere = [vol_connectcomp.value(resec_coord[0]/vol_connectcomp.getVoxelSize()[0]+vox_i,resec_coord[1]/vol_connectcomp.getVoxelSize()[1]+vox_j,resec_coord[2]/vol_connectcomp.getVoxelSize()[2]+vox_k) for vox_k in range(-Nsph[2],Nsph[2]+1) for vox_j in range(-Nsph[1],Nsph[1]+1) for vox_i in range(-Nsph[0],Nsph[0]+1) if math.sqrt(vox_i**2+vox_j**2+vox_k**2) < sphere_size]
#         
#         voxel_to_keep = [x for x in voxel_within_sphere if x != 0]
#         from collections import Counter
#         most_common,num_most_common = Counter(voxel_to_keep).most_common(1)[0]
#         
#         ret = subprocess.call(['AimsThreshold', '-i',di_resec.fullPath(),'-m', 'eq','-t',str(most_common) , '-o', di_resec.fullPath()])
#         ret = subprocess.call(['AimsMorphoMath', '-m', 'dil', '-i', di_resec.fullPath(), '-o', di_resec.fullPath(), '-r', '1.75'])
#         if ret < 0:
#             print "Importation error"
#             QtGui.QMessageBox.warning(self, "Error", "Brainvisa Importation error/ AimsFileConvert")
#             return
#         
#         neuroHierarchy.databases.insertDiskItem( di_resec, update=True )
#         self.transfoManager.setReferentialTo(di_resec, self.diskItems['T1pre'].attributes()['referential'] )
#         
#         wdi_resec_roi = WriteDiskItem( 'ROI IntrAnat', 'Graph' )
#         di_resec_roi = wdi_resec_roi.findValue({'subject':self.brainvisaPatientAttributes['subject'], 'center':self.currentProtocol, 'acquisition':'Resection'})
# 
#         self.brainvisaContext.runInteractiveProcess(lambda x='',di_roi=di_resec_roi,di_res=di_resec:self.roiconversionDone(di_roi,di_resec),'Volume To ROI Graph Converter', read = di_resec, write = di_resec_roi)  #, sulcus_identification ='label')


#     def roiconversionDone(self,di_resec_roi,di_resec):
#         neuroHierarchy.databases.insertDiskItem( di_resec_roi, update=True )
#         self.transfoManager.setReferentialTo(di_resec_roi, self.diskItems['T1pre'].attributes()['referential'] )
#         obj_roi = self.a.loadObject(di_resec_roi)
#         
#         Text_win1 = self.windowCombo1.currentText()
#         obj = self.a.loadObject(di_resec)
#         self.diskItems['Resection'] = di_resec
#         self.dispObj['Resection']=obj
#         # Update list of available items in the combo boxes
#         self.windowContent.update({'Resection':['Resection','electrodes']})
#         self.updateComboboxes(Text_win1, 'Resection')
#         # Update anatomist windows
#         self.updateAllWindows()
        
#     def ROIResectiontoNiftiResection(self):
#         wdi_resec = WriteDiskItem('Resection', 'NIFTI-1 image')
#         di_resec = wdi_resec.findValue({'subject':self.brainvisaPatientAttributes['subject'], 'center':self.currentProtocol, 'acquisition':'Resection'})
#         
#         wdi_resec_roi = ReadDiskItem( 'ROI Intranat', 'Graph', requiredAttributes={'subject':self.brainvisaPatientAttributes['subject'], 'center':self.currentProtocol} )
#         di_resec_roi = list(wdi_resec_roi.findValues({}, None, False))
# 
#         self.brainvisaContext.runInteractiveProcess(lambda x='',di_roi=di_resec_roi[0],di_res=di_resec:self.roiconversionDone(di_roi,di_resec),'Graph To Volume Converter', read = di_resec_roi[0], write = di_resec) #removeSource, False, extract_contours, 'No'


# =============================================================================
# MAIN: Main function that starts the interface
# =============================================================================
def main(noapp=0):
    # Create application
    app = None
    if noapp == 0:
        if hasattr(QtCore.Qt, 'AA_X11InitThreads'):
            QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_X11InitThreads)
        app = QtWidgets.QApplication(sys.argv)
        axon.initializeProcesses()
    # Show main window
    w = LocateElectrodes(app = app)
    #w.setWindowFlags(QtCore.Qt.Window)
    w.show()
    # Run the application
    if noapp == 0 and app is not None:
        app.exec_()

if __name__ == "__main__":
    main()


