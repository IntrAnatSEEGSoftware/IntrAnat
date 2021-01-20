# -*- coding: utf-8 -*-
#
# Importation of data in the database + image registration and normalization
#(c) Inserm U836 2012-2021 - Manik Bhattacharjee, Pierre Deman, Francois Tadel 
#
# License GNU GPL v3

import os, subprocess, pickle, shutil, tempfile, json, numpy
from scipy import ndimage
from shutil import copyfile

from soma import aims
from brainvisa import axon, anatomist
from brainvisa.configuration import neuroConfig
neuroConfig.gui = True
from brainvisa.configuration.qt4gui import neuroConfigGUI
from brainvisa.data import neuroHierarchy
import brainvisa.registration as registration
from brainvisa.processes import *
from soma.qt_gui.qt_backend import uic, QtGui, QtCore
from brainvisa.data.readdiskitem import ReadDiskItem
from brainvisa.data.writediskitem import WriteDiskItem
import brainvisa.data.neuroDiskItems
from soma.wip.application.api import Application
from freesurfer.brainvisaFreesurfer import *

from externalprocesses import *
import patientinfo, pathologypatientinfo
from TimerMessageBox import *
from progressbar import ProgressDialog
import DialogCheckbox


#  Matlab code : coregister file1 to file2
def matlab_cellstr(listOfStrings):
    """ Converts a list of filenames as a matlab-suitable cell table of strings and adds ',1' at the end of each string (for spm file lists) """
    return "{'"+",1' '".join([str(stri) for stri in listOfStrings])+",1'}"

# Coregister two files and saves the matrix in an ASCII file (space-separated)
spm_coregister = """try
    addpath(genpath(%s));
    VF=spm_vol(%s);
    VG=spm_vol(%s);
    centermatF = eye(4);
    center_F = %s;
    F_orient = %s;
    F_orient = reshape(F_orient,4,4)';
    centermatF(:,4) = F_orient*[center_F 1]';
    centermatG = eye(4);
    center_G = %s;
    G_orient = %s;
    G_orient = reshape(G_orient,4,4)';
    centermatG(:,4) = G_orient*[center_G 1]';
    centeredmatF = inv(centermatF)*VF.private.mat;
    centeredmatG = inv(centermatG)*VG.private.mat;
    VF.mat = centeredmatF;
    VF.private.mat = centeredmatF;
    VF.private.mat0 = centeredmatF;
    VG.mat = centeredmatG;
    VG.private.mat = centeredmatG;
    VG.private.mat0 = centeredmatG;
    x = spm_coreg(VF,VG);
    matCoreg = spm_matrix(x(:)');
    trm = centermatG*matCoreg*inv(centermatF);
    trm = [trm(1:3,4)';trm(1:3,1:3)];
    dlmwrite(%s,trm, 'delimiter',' ','precision',16); \n
catch
    disp 'AN ERROR OCCURED'; 
end
quit;"""


# SPM coregistration onto file1 ({'/home/manik/data/epilepsie/IRM-testOD/Pre/3DT1.img,1'}) of file 2 ({'/home/manik/data/epilepsie/IRM-testOD/Post/SagT2.img,1'}) with reslicing
spm_coregisterReslice = """try
    addpath(genpath(%s));
    spm('defaults', 'FMRI');spm_jobman('initcfg');
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref == %s;
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = %s;
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 1;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
    spm_jobman('run',matlabbatch);
catch
    disp 'AN ERROR OCCURED';
end
quit;"""

spm8_normalise = """try
    addpath(genpath(%s)); 
    spm('defaults', 'FMRI');
    spm_jobman('initcfg');
    matlabbatch{1}.spm.spatial.normalise.estwrite.subj.source = %s;
    matlabbatch{1}.spm.spatial.normalise.estwrite.subj.wtsrc = '';
    matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = %s;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.template = %s;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.weight = '';
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.smosrc = 8;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.smoref = 0;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.regtype = 'mni';
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.cutoff = 25;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.nits = 16;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = 1;
    matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.preserve = 0;
    matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.bb = [-78 -112 -50
                                                                 78 76 85];
    matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.vox = [1 1 1];
    matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.interp = 1;
    matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.prefix = 'w';
    spm_jobman('run',matlabbatch);
catch
    disp 'AN ERROR OCCURED';
end
quit;"""

spm12_normalise = """try
    addpath(genpath(%s)); 
    spm('defaults', 'FMRI');
    spm_jobman('initcfg');
    matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = %s;
    matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = %s;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm = {%s};
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
                                                                 78 76 85];
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = [1 1 1];
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.interp = 4;
    spm_jobman('run',matlabbatch);
catch
    disp 'AN ERROR OCCURED';
end
quit;"""

matlab_removeGado = """try
    addpath(genpath(%s));
    spm('defaults', 'FMRI');
    spm_jobman('initcfg');
    matlabbatch{1}.spm.spatial.preproc.channel.vols = {%s};
    matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
    matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {%s};
    matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {%s};
    matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {%s};
    matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {%s};
    matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {%s};
    matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {%s};
    matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{1}.spm.spatial.preproc.warp.samp = 2;
    matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];
    spm_jobman('run',matlabbatch);
    c1 = spm_vol(%s);
    c2 = spm_vol(%s);
    c3 = spm_vol(%s);
    c4 = spm_vol(%s);
    Yc1 = spm_read_vols(c1);
    Yc2 = spm_read_vols(c2);
    Yc3 = spm_read_vols(c3);
    Yc4 = spm_read_vols(c4);
    keepC1 = find(Yc1 ~= 0);
    keepC2 = find(Yc2 ~= 0);
    keepC3 = find(Yc3 ~= 0);
    keepC4 = find(Yc4 ~= 0);
    all_keep = unique([keepC1; keepC2]);
    fullImage = spm_vol(%s);
    YfullImage = spm_read_vols(fullImage);
    all_remove = [1:1:size(YfullImage(:))];
    all_remove(all_keep)=[];
    YfullImage(all_remove)=0;
    fullImage.fname=%s;
    spm_write_vol(fullImage,YfullImage);
catch
    disp 'AN ERROR OCCURED';
end
quit;"""

spm_inverse_y_field12 = """try
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

spm_MNItoScannerBased = """try
    addpath(genpath(%s));
    spm('defaults', 'FMRI');
    spm_jobman('initcfg');
    clear matlabbatch;
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {%s};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {%s};
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-150 -150 -150
                                                              150 150 150];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'tmpMNItoScannerBased';
    matlabbatch{2}.spm.spatial.realign.write.data = {
                                                     %s
                                                     %s
                                                     };
    matlabbatch{2}.spm.spatial.realign.write.roptions.which = [2 1];
    matlabbatch{2}.spm.spatial.realign.write.roptions.interp = 4;
    matlabbatch{2}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
    matlabbatch{2}.spm.spatial.realign.write.roptions.mask = 1;
    matlabbatch{2}.spm.spatial.realign.write.roptions.prefix = '';
    spm_jobman('run',matlabbatch);
catch
    disp 'AN ERROR OCCURED';
end
quit;"""


class ImageImport (QtGui.QDialog):
    """ImageImport is the main dialog class of the Image Importer software"""

    def __init__ (self, app=None):
        QtGui.QWidget.__init__(self)
        self.ui = uic.loadUi("ImageImport.ui", self)
        self.setWindowTitle('Image Import - NOT FOR MEDICAL USE')
        self.app = app
    
        self.seriesUIDbyName = {}
        self.studiesUIDbyName = {}
        self.currentProtocol = None
        self.currentSubject = None
        self.AcPc = {}
        self.brainCenter = None
        self.defaultAcqDate = QtCore.QDate(1900,1,1)
        self.ui.acqDate.setDate(self.defaultAcqDate)
        self.modas = {'t1mri':'Raw T1 MRI', 't2mri':'T2 MRI', 'ct': 'CT', 'pet': 'PET','fmri_epile':'fMRI-epile', 'statistic_data':'Statistic-Data', 'flair': 'FLAIR', 'freesurfer_atlas':'FreesurferAtlas', 'fgatir':'FGATIR'}

        # Allow calling brainvisa
        self.brainvisaContext = defaultContext()
        # Get Transformation Manager
        self.transfoManager = registration.getTransformationManager()
        # Load anatomist
        self.a = anatomist.Anatomist('-b')
        # Force all the scanner-based referentials to be linked 
        self.a.config()['setAutomaticReferential'] = 1
        self.a.config()['commonScannerBasedReferential'] = 1

        layout = QtGui.QVBoxLayout( self.ui.windowContainer1 )
        layout2 = QtGui.QHBoxLayout( )
        layout3 = QtGui.QHBoxLayout( )
        self.axWindow = self.a.createWindow( 'Axial' )#, no_decoration=True )
        self.axWindow2 = self.a.createWindow( 'Sagittal' )#, no_decoration=True )
        self.axWindow3 = self.a.createWindow( 'Axial' )#, no_decoration=True )
        self.axWindow4 = self.a.createWindow( 'Sagittal' )#, no_decoration=True )
        self.wins = [self.axWindow, self.axWindow2, self.axWindow3, self.axWindow4]
        self.axWindow.setParent(self.ui.windowContainer1)
        self.axWindow2.setParent(self.ui.windowContainer1)
        self.axWindow3.setParent(self.ui.windowContainer1)
        self.axWindow4.setParent(self.ui.windowContainer1)
        layout2.addWidget( self.axWindow.getInternalRep() )
        layout2.addWidget( self.axWindow2.getInternalRep() )
        layout3.addWidget( self.axWindow3.getInternalRep() )
        layout3.addWidget( self.axWindow4.getInternalRep() )
        layout.addLayout(layout2)
        layout.addLayout(layout3)
        
        # Hide control window
        self.a.getControlWindow().setVisible(False)
        self.dispObj = []
        self.initialize()

        self.threads = []
        self.waitingForThreads = {}
    
        # Initialize patient informations tab
        self.patientInfo = patientinfo.buildUI(self.ui.patGroupBox)
        cdate = self.patientInfo['currentDate'].date()
        bdate = self.patientInfo['patientBirthday'].date()
        page = bdate.daysTo(cdate)/365
        self.patientInfo['patientAge'].setText(str(page))
        self.patientInfo['comoraucune'].stateChanged.connect(self.checkbox_comor)
        self.patientInfo['patientBirthday'].dateChanged.connect(self.update_patientAge)
    
        # TABS: Reload at each tab change
        self.connect(self.ui.tabWidget,  QtCore.SIGNAL('currentChanged(int)'), self.analyseBrainvisaDB)
        self.connect(self.ui.tabWidget,  QtCore.SIGNAL('currentChanged(int)'), self.setTabGuide)
        # TAB1: Database
        self.connect(self.ui.bvProtocolCombo, QtCore.SIGNAL('currentIndexChanged(QString)'), lambda s:self.selectProtocol(s,self.ui.bvSubjectCombo))
        self.connect(self.ui.bvSubjectCombo, QtCore.SIGNAL('currentIndexChanged(QString)'), self.selectBvSubject)
        self.connect(self.ui.bvSubjectCombo, QtCore.SIGNAL('activated(QString)'), self.setCurrentSubject)
        self.connect(self.ui.bvImageList, QtCore.SIGNAL('itemDoubleClicked(QListWidgetItem*)'), self.selectBvImage)
        self.connect(self.ui.bvDeleteImageButton, QtCore.SIGNAL('clicked()'), self.deleteBvImage)
        self.connect(self.ui.bvDeleteSubjectButton, QtCore.SIGNAL('clicked()'), self.deleteBvSubject)
        self.connect(self.ui.bvEditPref, QtCore.SIGNAL('clicked()'), self.editBvPref)
        self.connect(self.ui.bvUpdateDb, QtCore.SIGNAL('clicked()'), self.updateBvDb)
        self.connect(self.ui.bvImportBids, QtCore.SIGNAL('clicked()'), self.importFromBids)
        # TAB2: Add subject
        self.connect(self.ui.subjectSiteCombo, QtCore.SIGNAL('activated(QString)'), self.updatePatientCode)
        self.connect(self.ui.subjectSiteCombo, QtCore.SIGNAL('editTextChanged(QString)'), self.updatePatientCode)
        self.connect(self.ui.subjectYearSpinbox, QtCore.SIGNAL('valueChanged(int)'), self.updatePatientCode)
        self.connect(self.ui.subjectPatientName, QtCore.SIGNAL('textChanged(QString)'), self.updatePatientCode)
        self.connect(self.ui.subjectPatientFirstName, QtCore.SIGNAL('textChanged(QString)'), self.updatePatientCode)
        self.connect(self.ui.subjectAddSubjectButton, QtCore.SIGNAL('clicked()'), self.storePatientInDB)
        # TAB3: Patient info
        self.connect(self.ui.patProtocolCombo, QtCore.SIGNAL('currentIndexChanged(QString)'), lambda s:self.selectProtocol(s,self.ui.patSubjectCombo))
        self.connect(self.ui.patProtocolCombo, QtCore.SIGNAL('activated(QString)'),lambda :pathologypatientinfo.buildUI(self.ui.patpatGroupBox,self.currentProtocol)) #lambda :pathologypatientinfo.buildUI(self.ui.patpatGroupBox,self.currentProtocol)) self._pathologyinfo(self.ui.patpatGroupBox,self.currentProtocol)
        self.connect(self.ui.patSubjectCombo, QtCore.SIGNAL('currentIndexChanged(QString)'), self.selectPatSubject)
        self.connect(self.ui.patSubjectCombo, QtCore.SIGNAL('activated(QString)'), self.setCurrentSubject)
        self.connect(self.ui.patInfoValidate,QtCore.SIGNAL('clicked()'),self.ValidatePatInfo)
        # TAB4: Import to BrainVisa
        self.connect(self.ui.niftiImportButton, QtCore.SIGNAL('clicked()'), self.importNifti)
        self.connect(self.ui.niftiSetBraincenterButton, QtCore.SIGNAL('clicked()'), self.setBrainCenter)
        self.connect(self.ui.ImportFSoutputspushButton, QtCore.SIGNAL('clicked()'),self.importFSoutput)
        self.connect(self.ui.buttonImportLausanne, QtCore.SIGNAL('clicked()'),self.importLausanne2008)
        self.connect(self.ui.niftiProtocolCombo, QtCore.SIGNAL('currentIndexChanged(QString)'), lambda s:self.selectProtocol(s,self.ui.niftiSubjectCombo))
        self.connect(self.ui.niftiSubjectCombo, QtCore.SIGNAL('activated(QString)'), self.setCurrentSubject)
        self.connect(self.ui.niftiSeqType, QtCore.SIGNAL('currentIndexChanged(QString)'),self.enable_disable_gadooption)
        # TAB5: Coregistration
        self.connect(self.ui.regProtocolCombo, QtCore.SIGNAL('currentIndexChanged(QString)'), lambda s:self.selectProtocol(s,self.ui.regSubjectCombo))
        self.connect(self.ui.regSubjectCombo, QtCore.SIGNAL('currentIndexChanged(QString)'), self.selectRegSubject)
        self.connect(self.ui.regSubjectCombo, QtCore.SIGNAL('activated(QString)'), self.setCurrentSubject)
        self.connect(self.ui.regImageList, QtCore.SIGNAL('itemDoubleClicked(QListWidgetItem*)'), self.selectRegImage)
        self.connect(self.ui.regImageList2, QtCore.SIGNAL('itemDoubleClicked(QListWidgetItem*)'), self.selectRegImage2)
        self.connect(self.ui.registerNormalizeSubjectButton, QtCore.SIGNAL('clicked()'), self.registerNormalizeSubject)
        self.connect(self.ui.segmentationHIPHOPbutton,QtCore.SIGNAL('clicked()'),self.runPipelineBV)
        self.connect(self.ui.runMarsAtlasFreesurferButton,QtCore.SIGNAL('clicked()'),self.runPipelineFS)
        self.connect(self.ui.runHiphopOnly, QtCore.SIGNAL('clicked()'), self.runProcessHiphop)
        # TAB6: Preferences
        self.connect(self.ui.prefSpmTemplateButton, QtCore.SIGNAL('clicked()'), self.setSpmTemplatePath)
        self.connect(self.ui.prefANTsButton, QtCore.SIGNAL('clicked()'), self.setANTsPath)
        self.connect(self.ui.prefFreesurferButton, QtCore.SIGNAL('clicked()'), self.setFreesurferPath)
        self.connect(self.ui.prefBidsButton, QtCore.SIGNAL('clicked()'), self.setBidsPath)
        self.connect(self.ui.prefANTScheckbox, QtCore.SIGNAL('clicked()'), lambda: self.setPrefCoregister('ANTS'))
        self.connect(self.ui.prefSPMcheckbox, QtCore.SIGNAL('clicked()'), lambda: self.setPrefCoregister('SPM'))
        self.ui.radioProjFtract.toggled.connect(self.updatePatientCode)
        self.ui.radioProjFtract.toggled.connect(self.switchProjectButton)
        self.ui.radioOther.toggled.connect(self.switchProjectButton)
        self.ui.radioOther.toggled.connect(self.updatePatientCode)
        self.ui.radioClassic.toggled.connect(self.switchProjectButton)
        self.ui.radioClassic.toggled.connect(self.updatePatientCode)
        self.ui.radioProjNeuro.toggled.connect(self.switchProjectButton)
        self.ui.radioProjNeuro.toggled.connect(self.updatePatientCode)
        self.ui.SavePreferencespushButton.clicked.connect(self.savePreferences)
    
        self.warningMEDIC()

        # Rest of the interface
        def setDictValue(d, k, v,button):
            d[k]=v
            button.setStyleSheet("background-color: rgb(90, 255, 95);")
        self.connect(self.ui.regAcButton, QtCore.SIGNAL('clicked()'), lambda:setDictValue(self.AcPc,'AC',list(self.a.linkCursorLastClickedPosition()),self.ui.regAcButton))
        self.connect(self.ui.regPcButton, QtCore.SIGNAL('clicked()'), lambda:setDictValue(self.AcPc,'PC',list(self.a.linkCursorLastClickedPosition()),self.ui.regPcButton))
        self.connect(self.ui.regIhButton, QtCore.SIGNAL('clicked()'), lambda:setDictValue(self.AcPc,'IH',list(self.a.linkCursorLastClickedPosition()),self.ui.regIhButton))
        self.connect(self.ui.regLhButton, QtCore.SIGNAL('clicked()'), lambda:setDictValue(self.AcPc,'LH',list(self.a.linkCursorLastClickedPosition()),self.ui.regLhButton))
        self.connect(self.ui.regAcPcValidateButton, QtCore.SIGNAL('clicked()'), self.validateAcPc)
        
        # Finds the available protocols in brainvisa database and fills the comboboxes
        self.analyseBrainvisaDB()
        self.setGuide("Choose on of the tabs to import data in the BrainVisa database")
        self.enableACPCButtons(False)
        # Show anatomist
        self.showAnatomist.setIcon(QtGui.QIcon('logoAnatomist.png'))
        self.showAnatomist.setIconSize(QtGui.QSize(24, 24))
        self.connect(self.showAnatomist, QtCore.SIGNAL('clicked()'), self.toggleAnatomistWindow)


    def __del__ (self):
        self.ui = None

    def closeEvent(self, event):
        QtGui.QFrame.closeEvent(self, event)
        self.savePreferences()
        # Kill all hat may still be running (stuck) or better, wait for them to complete
        maxTimeToCompletion = 1000*60*20 # 20min in ms -> maximum duration for any thread (it will be killed after that)
        qthr = [t for t in self.threads if isinstance(t, QtCore.QThread)]
        thr = [t for t in self.threads if not isinstance(t, QtCore.QThread)]
        print "Closing all qthreads"
        unfinished = [t for t in qthr if not t.wait(maxTimeToCompletion)]
        print "Killing %i QThreads still running"%len(unfinished)
        [t.terminate() for t in unfinished]
        def checkKill(p, timeout):
            """ This function checks if the thread is still running.
            If so, it waits until timeout (in ms), then kills it if it is still running"""
            try:
                print "Trying to kill "+repr(p)
                t = 0
                while(p.poll() is None and t < timeout):
                    time.sleep(5)
                    t = t + 5000
                if p.poll() is None:
                    p.kill()
            except:
                pass

        print "Closing all threads"
        finished = [checkKill(t, maxTimeToCompletion) for t in thr]
        print "Quitting BrainVisa"
        axon.cleanup()
        print "Quit !"
        self.app.quit()


    def setStatus(self, text):
        """ Sets the status text displayed at the bottom of the dialog"""
        self.ui.statusCombo.insertItem(0,QtCore.QTime.currentTime().toString("H:m:s") + ': ' + text)
        self.ui.statusCombo.setCurrentIndex(0)

    def setGuide(self, text):
        """ Set the text of the User Guide on the top of the window"""
        self.ui.guideLabel.setText(text)

    def setTabGuide(self, tabId):
        helpTab = [u'Look for images in BrainVisa database. Double click to load an images in anatomist windows on the side',\
                   u'Add a new patient in the Database (perform an anonymization as well)',\
                   u'Precise patient information (general and pathology specific)',\
                   u'Import images (format: Nifti/GIS/Analyse/(mgz in some case). Select the file and fill information needed',\
                   u'Coregister all images to the 3DT1 before electrode implantation. Normalize this 3DTI to the MNI, perform grey/white matter segmentation, mesh of the cortex, and MarsAtlas parcellation',\
                   u'Set software preferences',\
                   u'Import SEEG data' ,\
                   u'Importez des images depuis un dossier DICOM : Sélectionnez une image, importez et convertissez en Nifti, puis importez-là dans la base',\
                   u'Importez des images depuis un PACS. Sélectionnez une image, importez et convertissez en Nifti, puis importez-là dans la base',\
                   u'Importez les schémas d\'implantation scannés',\
                   ]
        self.setGuide(helpTab[tabId])
        

    def initialize(self):
        """Reads the preference file, checks the available patients and sequences"""
        #Localisation du fichier de preference
        prefpath = os.path.join(os.path.expanduser('~'), '.imageimport')

        # Chargement du fichier et stockage dans self.prefs
        self.prefs = {'sites':['Gre',],
                      'coregisterMethod':'spm'}
        try:
            if (os.path.exists(prefpath)):
                filein = open(prefpath, 'rU')
                try:
                    self.prefs = json.loads(filein.read())
                except:
                    filein.close()
                    filein = open(prefpath, 'rU')
                    self.prefs = pickle.load(filein)
    
                filein.close()
        except:
            print("There is a problem with %s opening apparently"%prefpath)

        if 'currentProtocol' in self.prefs:
            self.currentProtocol = self.prefs['currentProtocol']
        if 'currentSubject' in self.prefs:
            self.currentSubject = self.prefs['currentSubject']
        # Grenoble/Lyon/Other sites...
        if 'sites' in self.prefs:
            self.ui.subjectSiteCombo.clear()
            self.ui.subjectSiteCombo.addItems(self.prefs['sites'])

        #check what spm path has been precised in BrainVisa
        configuration = Application().configuration
        brainvisa_spm12_path = configuration.SPM.spm12_path
        brainvisa_freesurfer_home_path = configuration.freesurfer._freesurfer_home_path

        if 'spm' in self.prefs:
            if brainvisa_spm12_path != self.prefs['spm']:
                QtGui.QMessageBox.warning(self, u"SPM", u"SPM path different between IntrAnat and BrainVisa, strange, you should check that, by default keep the one precised in IntrAnat")
            if os.path.isfile(self.prefs['spm']+os.sep+'Contents.m'):
                op_cont = open(self.prefs['spm']+os.sep+'Contents.m')
                line = op_cont.readline()
                line = op_cont.readline()
                spm_version = line.split(' ')[3]
                if spm_version == '(SPM12)':
                    self.setSpmTemplatePath(self.prefs['spm']) #+os.sep+'toolbox/OldNorm/T1.nii')
                elif spm_version == '(SPM8)':
                    QtGui.QMessageBox.warning(self, u"SPM", u"SPM version not supported anymore")
                    #self.setSpmTemplatePath(self.prefs['spm']) #+os.sep+'templates/T1.nii')
                else:
                    QtGui.QMessageBox.warning(self, u"SPM", u"SPM version unknown")
            else:
                QtGui.QMessageBox.warning(self, u"SPM", u"SPM path is not defined in tab 'Préférences'\n Normalization and SPM coregistration won't work !")
        else:
            if len(configuration.SPM.spm12_path)>0:
                print("used brainvisa spm path")
                self.setSpmTemplatePath(brainvisa_spm12_path)
                self.prefs['spm']=brainvisa_spm12_path
            else:
                QtGui.QMessageBox.warning(self, u"SPM", u"SPM path is not defined (or wrong) in tab 'Préférences'\n Normalization and SPM coregistration won't work !")

        if 'freesurfer' in self.prefs:
            if brainvisa_freesurfer_home_path != self.prefs['freesurfer']:
                QtGui.QMessageBox.warning(self, u"Freesurfer", u"Freesurfer path different between IntrAnat and BrainVisa, strange, you should check that, by default keep the one precised in IntrAnat")
            self.setFreesurferPath(self.prefs['freesurfer'])
        else:
            self.setFreesurferPath(brainvisa_freesurfer_home_path)

        if 'ants' in self.prefs:
            self.setANTsPath(self.prefs['ants'])
            
        if 'bids' in self.prefs:
            self.setBidsPath(self.prefs['bids'])
            
        if 'coregisterMethod' in self.prefs:
            self.coregMethod = self.prefs['coregisterMethod']
            if self.coregMethod == 'ANTs':
                self.ui.prefSPMcheckbox.setCheckState(False)
                self.ui.prefANTScheckbox.setCheckState(True)
            if self.coregMethod == 'spm':
                self.ui.prefSPMcheckbox.setCheckState(True)
                self.ui.prefANTScheckbox.setCheckState(False)

        if 'projectSelected' in self.prefs:
            if self.prefs['projectSelected'] != []:
                projects= [self.ui.radioProjFtract, self.ui.radioProjNeuro, self.ui.radioClassic, self.ui.radioOther]
                projects[self.prefs['projectSelected'][0]].setChecked(True)
        else:
            self.prefs.update({'projectSelected':[2]})
            projects= [self.ui.radioProjFtract, self.ui.radioProjNeuro, self.ui.radioClassic, self.ui.radioOther]
            projects[self.prefs['projectSelected'][0]].setChecked(True)

    def warningMEDIC(self):
        shortwarning = TimerMessageBox(5,self)
        shortwarning.exec_()

    def toggleAnatomistWindow(self):
        self.a.getControlWindow().setVisible(self.showAnatomist.isChecked())
        
        
    def savePreferences(self):
        prefpath = os.path.join(os.path.expanduser('~'), '.imageimport')

        self.prefs['currentProtocol'] = self.currentProtocol
        self.prefs['currentSubject'] = self.currentSubject
        self.prefs['sites'] = [str(self.ui.subjectSiteCombo.itemText(item)) for item in range(self.ui.subjectSiteCombo.count())]
        # SPM path
        if len(str(self.ui.prefSpmTemplateEdit.text()))>0:
            self.prefs['spm'] = str(self.ui.prefSpmTemplateEdit.text())
        else:
            self.prefs.pop('spm',None)
        # ANTs path
        if len(str(self.ui.prefANTsEdit.text()))>0:
            self.prefs['ants']=str(self.ui.prefANTsEdit.text())
        else:
            self.prefs.pop('ants',None)
        # FreeSurfer path
        if len(str(self.ui.prefFreesurferEdit.text()))>0:
            self.prefs['freesurfer'] = str(self.ui.prefFreesurferEdit.text())
        else:
            self.prefs.pop('freesurfer',None)
        # BIDS database
        if len(str(self.ui.prefBidsEdit.text()))>0:
            self.prefs['bids'] = str(self.ui.prefBidsEdit.text())
        else:
            self.prefs.pop('bids',None)
            
        if self.ui.prefANTScheckbox.isChecked():
            self.prefs['coregisterMethod'] = 'ANTs'
        elif self.ui.prefSPMcheckbox.isChecked():
            self.prefs['coregisterMethod'] = 'spm'
        else:
            self.prefs['coregisterMethod'] = 'spm'

        self.coregMethod = self.prefs['coregisterMethod']
        projects= [self.ui.radioProjFtract.isChecked(),self.ui.radioProjNeuro.isChecked(), self.ui.radioClassic.isChecked(),self.ui.radioOther.isChecked()]
        project_selected = [x for x in range(len(projects)) if projects[x]==True]

        self.prefs['projectSelected'] = project_selected
        # Save file
        fileout = open(prefpath, 'wb')
        pickle.dump(self.prefs, fileout)
        fileout.close()

    def createItemDirs(self, item):
        """ Create the directories containing the provided WriteDiskItem and insert them in the BrainVisa database """
        # Copied from brainvisa.in ExecutionContext._processExecution()
        dirname = os.path.dirname( item.fullPath() )
        dir=dirname
        dirs = []
        while not os.path.exists( dir ):
            dirs.append(dir)
            dir=os.path.dirname(dir)
        if dirs:
            try:
                os.makedirs( dirname )
            except OSError, e:
                if not e.errno == errno.EEXIST:
                    # filter out 'File exists' exception, if the same dir has
                    # been created concurrently by another instance of BrainVisa
                    # or another thread
                    raise
            for d in dirs:
                dirItem=neuroHierarchy.databases.createDiskItemFromFileName(d, None)


    def storeImageReferentialsAndTransforms(self, image):
        """Stores referential and transformation files for an image from its headers -> native, scanner-based, coregistered to"""

        # Notes : on prend le header de l'image, on crée un objet referentiel pour le référentiel natif,
        #  puis un pour chaque transformation définie dans le header
        # On stocke les transfos entre le natif et chacun de ces référentiels à partir des infos du header.
        # Comment alors faire pour lier le 'coregistered to blabla' au 'scanner based' d'une autre image ? Stocker une transfo identité ou utiliser copyReferential ?
        # Et est-ce que copyReferential met à jour toutes les transfos stockées qui utilisaient le référentiel écrasé ?
        # D'autre part, comment repérer le référentiel Scanner-based d'une image en particulier ? Pour créer la transfo t1post vers t1pre à partir du recalage SPM...

        att = image.attributes()
        refName = att['modality']+'_'+att['acquisition']
        print "Storing refs/transfos for "+str(refName)
    
        # Create a Referential object in the database to represent the native referential of the image
        moda = self.modas[att['modality']]
        nativeRef = self.transfoManager.referential(image)
        if nativeRef is None:
            print '...Creating native ref'
            nativeRef = self.transfoManager.createNewReferentialFor(image, name = refName+'native', description = 'Native Anatomist Referential for ' + refName + ' image', referentialType='Referential of '+moda)
            try:
                neuroHierarchy.databases.createDiskItemFromFileName(os.path.dirname( nativeRef.fullPath() ))
            except:
                pass

        # Create the scanner-based referential
        refScan = self.transfoManager.findOrCreateReferential('Scanner Based Referential', image, name = refName+'_Scanner-Based', description = 'Scanner-Based referential for ' + refName + ' image')
        try:
            neuroHierarchy.databases.createDiskItemFromFileName(os.path.dirname( refScan.fullPath() ))
        except:
            pass

        if refScan is None:
            print '...Creating scanner-based ref'
            refScan = self.transfoManager.createNewReferentialFor(wdiScan, name = refName+'_Scanner-Based', description = 'Scanner-Based referential for ' + refName + ' image', referentialType='Scanner Based Referential')
            try:
                neuroHierarchy.databases.createDiskItemFromFileName(os.path.dirname( refScan.fullPath() ))
            except:
                pass

        # Find where to store the scanner-base transformation
        wdiTransformScan = WriteDiskItem( 'Transformation to Scanner Based Referential', 'Transformation matrix', exactType=True )
        transformScan = wdiTransformScan.findValue(image)

        # Read the header of the image and read all transformations available
        f = aims.Finder()
        f.check(image.fileName())
        trs = list(f.header()['transformations'])
        refs = list(f.header()['referentials'])
        trm_to_scannerBased = None

        if ("Scanner-based anatomical coordinates" in refs ):

            trm_to_scannerBased = trs[refs.index("Scanner-based anatomical coordinates")]

        if trm_to_scannerBased is None:
            print "No Scanner-based transfo in image header, using first available transfo instead."
            try:
                trm_to_scannerBased = trs[0]
            except:
                pass

        if transformScan is None:
            print "Cannot find path for Scanner-based transform in database -> transformation not stored !"
        elif trm_to_scannerBased is None:
            print "Cannot find Scanner-based transformation in header !"
        else:
            print '...storing scanner-based transfo'
            #Store information into the trm file
            mot = aims.Motion( trm_to_scannerBased )
            image.setMinf('SB_Transform',str(trm_to_scannerBased))
            aims.write( mot, transformScan.fullPath() )
            try:
                neuroHierarchy.databases.createDiskItemFromFileName(os.path.dirname( transformScan.fullPath() ))
            except:
                pass
            #set and update database
            print ".... setting transfo info : from %s to %s for %s"%(repr(nativeRef), repr(refScan), repr(transformScan))
            self.transfoManager.setNewTransformationInfo( transformScan, source_referential=nativeRef, destination_referential= refScan)
        

    def analyseBrainvisaDB(self, tab=None):
        """ Analyze the BrainVisa Database to fill the data for the brainvisa and registration tab (provide the tab number to update only the chosen one, 0 or 4) """
        rdi = ReadDiskItem( 'Center', 'Directory' )#, requiredAttributes={'center':'Epilepsy'} )
        protocols = list( rdi._findValues( {}, None, False ) )
        protocols = sorted([str(p.attributes()['center']) for p in protocols])
        if len(protocols) == 0:
            rep = QtGui.QMessageBox.warning(self, "Center (former protocol)", u"No center/protocole defined in BrainVisa Database !")
            text, ok = QtGui.QInputDialog.getText(self, 'Input Dialog', 'Enter the name of the center/protocole:')
    
            if (ok & bool(text)):
                wdi = WriteDiskItem('Center', 'Directory' )    #'gz compressed NIFTI-1 image' )
                di = wdi.findValue({'center' : str(text)})
                os.mkdir(di.fullPath())
                neuroHierarchy.databases.insertDiskItem( di, update=True )
    
        if self.currentProtocol in protocols:
            currentIndex = protocols.index(self.currentProtocol)
        else:
            currentIndex = 0

        def selectSubj(subjCombo):
            " Subfunction to select the current subject "
            if self.currentSubject is None:
                return
            idx = subjCombo.findText(self.currentSubject)
            # print "search current subj"
            if idx != -1:
                # print "found at "+repr(idx)
                subjCombo.setCurrentIndex(idx)

        if tab == 0 or tab is None: # Tab 0 is BrainVisa database tab
            self.ui.bvProtocolCombo.clear()
            self.ui.bvProtocolCombo.addItems(protocols)
            self.ui.bvProtocolCombo.setCurrentIndex(currentIndex)
            selectSubj(self.ui.bvSubjectCombo)
        if tab == 1 or tab is None: # Tab 1 is Create Subject tab
            self.ui.subjectProtocolCombo.clear()
            self.ui.subjectProtocolCombo.addItems(protocols)
            self.ui.subjectProtocolCombo.setCurrentIndex(currentIndex)
        if tab == 2 or tab is None: # Tab 2 is Patient Information tab
            self.ui.patProtocolCombo.clear()
            self.ui.patProtocolCombo.addItems(protocols)
            self.ui.patProtocolCombo.setCurrentIndex(currentIndex)
            self.pathologypatientInfo = pathologypatientinfo.buildUI(self.ui.patpatGroupBox,self.currentProtocol)
            selectSubj(self.ui.patSubjectCombo)
        if tab == 3 or tab is None: # Tab 3 is Import Nifti tab
            self.ui.niftiProtocolCombo.clear()
            self.ui.niftiProtocolCombo.addItems(protocols)
            self.ui.niftiProtocolCombo.setCurrentIndex(currentIndex)
            selectSubj(self.ui.niftiSubjectCombo)
        if tab == 4 or tab is None: # Tab 4 is Registration tab
            self.ui.regProtocolCombo.clear()
            self.ui.regProtocolCombo.addItems(protocols)
            self.ui.regProtocolCombo.setCurrentIndex(currentIndex)
            selectSubj(self.ui.regSubjectCombo)

    def selectProtocol(self, proto, comboToFillWithSubjects):
        """ A BrainVisa protocol was selected : query the database to get the list of subjects and put them in the provided QComboBox"""
        self.currentProtocol = str(proto)
        rdi = ReadDiskItem( 'Subject', 'Directory', requiredAttributes={'center':str(proto)} )
        subjects = list( rdi._findValues( {}, None, False ) )
        comboToFillWithSubjects.clear()#self.ui.bvSubjectCombo.clear()
        comboToFillWithSubjects.addItems(sorted([s.attributes()['subject'] for s in subjects]))

    def setCurrentSubject(self, subj):
        if not subj: # ignore empty values
            return
        print "current subj : "+str(subj)
        self.currentSubject = str(subj)
        self.clearAnatomist()
        

    def findAllImagesForSubject(self, protocol, subj):
        """ Returns all images for a subject from the brainvisa database.
            @return an array of ReadDiskItem
        """
        rdi = ReadDiskItem( 'Raw T1 MRI', 'BrainVISA volume formats', requiredAttributes={'center':str(protocol), 'subject':str(subj), 'normalized':'no' } )
        images = list( rdi._findValues( {}, None, False ) )
        images = [x for x in images if "skull_stripped" not in x.fullName()]
        rdi = ReadDiskItem( 'T2 MRI', 'BrainVISA volume formats', requiredAttributes={'center':str(protocol), 'subject':str(subj) } )
        images += list( rdi._findValues( {}, None, False ) )
        rdi = ReadDiskItem( 'CT', 'BrainVISA volume formats', requiredAttributes={'center':str(protocol), 'subject':str(subj) } )
        images += list( rdi._findValues( {}, None, False ) )
        rdi = ReadDiskItem( 'PET', 'BrainVISA volume formats', requiredAttributes={'center':str(protocol), 'subject':str(subj) } )
        images += list( rdi._findValues( {}, None, False ) )
        rdi = ReadDiskItem( 'fMRI-epile', 'BrainVISA volume formats', requiredAttributes={'center':str(protocol), 'subject':str(subj) } )
        images += list( rdi._findValues( {}, None, False ) )
        rdi = ReadDiskItem( 'Statistic-Data', 'BrainVISA volume formats', requiredAttributes={'center':str(protocol), 'subject':str(subj) } )
        images += list( rdi._findValues( {}, None, False ) )
        rdi = ReadDiskItem( 'FLAIR', 'BrainVISA volume formats', requiredAttributes={'center':str(protocol), 'subject':str(subj) } )
        images += list( rdi._findValues( {}, None, False ) )
        rdi = ReadDiskItem( 'FreesurferAtlas', 'BrainVISA volume formats', requiredAttributes={'center':str(protocol), 'subject':str(subj) } )
        images += list( rdi._findValues( {}, None, False ) )
        rdi = ReadDiskItem( 'FGATIR', 'BrainVISA volume formats', requiredAttributes={'center':str(protocol), 'subject':str(subj) } )
        images += list( rdi._findValues( {}, None, False ) )
        return images


    def selectBvSubject(self, subj):
        """ A BrainVisa subject was selected : query the database to get the list of images"""
        # Change subject
        self.clearAnatomist()
        # Display "Date : XXXX-XX-XX - Seq: T1 - Acq : T1Pre
        self.ui.bvImageList.clear()
        images = self.findAllImagesForSubject(self.ui.bvProtocolCombo.currentText(), subj)
        self.ui.bvImageList.addItems(sorted([i.attributes()['modality'] + ' - '+ i.attributes()['acquisition'] + ' - ' + i.attributes()['subacquisition'] if 'subacquisition' in i.attributes().keys()\
                else i.attributes()['modality'] + ' - '+ i.attributes()['acquisition'] for i in images ]))

        self.bvImages = {}
        for i in images:
            if 'subacquisition' in i.attributes().keys():
                self.bvImages.update({i.attributes()['modality'] + ' - ' + i.attributes()['acquisition'] + ' - ' + i.attributes()['subacquisition']:i})
            else:
                self.bvImages.update({i.attributes()['modality'] + ' - ' + i.attributes()['acquisition']:i})
               

    def selectBvImage(self, item):
        """ A BrainVisa image was double-clicked : display it !"""
        # Get image to display
        image = self.bvImages[str(item.text())]
        # Display images
        self.displayImage(image, self.wins)
            

    def removeFromDB(self, file, db=None):
        """
        If the file is a directory, recursive call to remove all its content before removing the directory.
        Corresponding diskitem is removed from the database if it exists.
        Taken from brainvisa-4.3.0/python/brainvisa/data/qt4gui/hierarchyBrowser.py
        """
        if db is None:
            try:
                db=neuroHierarchy.databases.database(neuroHierarchy.databases.getDiskItemFromFileName(file).get("_database"))
            except:
                pass
    
        if os.path.isdir(file):
            for f in os.listdir(file):
                self.removeFromDB(os.path.join(file, f), db)
            os.rmdir(file)
        else:
            os.remove(file)
        if db:
            diskItem=db.getDiskItemFromFileName(file, None)
            if diskItem:
                db.removeDiskItem(diskItem)

    def deleteBvSubject(self):
        protocol = str(self.ui.bvProtocolCombo.currentText())
        subj = str(self.ui.bvSubjectCombo.currentText())
        rep = QtGui.QMessageBox.warning(self, u'Confirmation', u"<font color='red'><b>WARNING</b><br/>You are about to erase ALL the patient data.<br/><br/><b>DELETE PATIENT \"%s\" ?</b></font>"%subj, QtGui.QMessageBox.Yes | QtGui.QMessageBox.No, QtGui.QMessageBox.No)
        if rep == QtGui.QMessageBox.Yes:
            print "Deleting subject %s"%subj
            rdi = ReadDiskItem( 'Subject', 'Directory', requiredAttributes={'center':str(protocol), 'subject':str(subj) } )
            di = list( rdi._findValues( {}, None, False ) )
            if len(di) != 1:
                print "%s subject(s) found ! Cannot delete if there is more than 1 !"%repr(len(di))
                return
            self.removeFromDB(di[0].fullPath(), neuroHierarchy.databases.database(di[0].get("_database")))
            self.setStatus(u"Subject %s deleted from the database"%subj)
            # Reset the list of subjects :
            self.selectProtocol(str(self.ui.bvProtocolCombo.currentText()),self.ui.bvSubjectCombo)


    def deleteBvImage(self):
        #Delete the current acquisition (including transforms and segmentations : the full acquisition directory
        imageName = None
        try:
            imageName = str(self.ui.bvImageList.currentItem().text())
        except: # Probably no image available or no image selected
            return
        
        rep = QtGui.QMessageBox.warning(self, u'Confirmation', u"<font color='red'><b>WARNING</b><br/>You are about to delete the selected image and all linked data.<br/><br/><b>DELETE IMAGE \"%s\" ?</b></font>"%imageName, QtGui.QMessageBox.Yes | QtGui.QMessageBox.No, QtGui.QMessageBox.No)
        if rep == QtGui.QMessageBox.Yes:
            image = self.bvImages[imageName] # Gives the path of the image
            di = neuroHierarchy.databases.getDiskItemFromFileName(image.fileName())
            acqPath = os.path.dirname(image.fileName())
            if str(os.path.basename(acqPath)) != str(di.attributes()['acquisition']):
                print "CANNOT REMOVE IMAGE : acquisition path does not match acquisition name !"
                return
        
            if "subacquisition" in di.attributes().keys():
                self.removeDiskItems(di,eraseFiles=True)
            else:
                self.removeFromDB(acqPath)
            self.setStatus(u"Image %s deleted from the database"%imageName)
            self.selectBvSubject(str(self.ui.bvSubjectCombo.currentText()))


    # ===== EDIT BRAINVISA PREFERENCES =====
    def editBvPref(self):
        """ Open a dialog window to edit the BrainVISA preferences """
        # Edit preferences
        neuroConfigGUI.editConfiguration()
        # Reset the list of subjects
        self.analyseBrainvisaDB()

    # ===== UPDATE BRAINVISA DATABASES =====
    def updateBvDb(self):
        """ Open a dialog window to update BrainVISA databases """
        from brainvisa.processing.qt4gui import neuroProcessesGUI
        # Create dialog to update the databases
        d = QtGui.QDialog()
        mainLayout = QtGui.QVBoxLayout()
        # Display update database process
        updateDialog = neuroProcessesGUI.ProcessView( brainvisa.processes.getProcessInstance( 'updateDatabases' ), d )
        mainLayout.addWidget(updateDialog)
        d.setLayout(mainLayout)
        d.move(100, 100)
        d.resize(600, 801)
        d.exec_()
        # Reset the list of subjects
        self.analyseBrainvisaDB()
        
    # ===== IMPORT FROM BIDS =====
    def importFromBids(self):
        """ Open a dialog window to import subjects from a BIDS database """
        if self.prefs['bids'] is None:
            QtGui.QMessageBox.warning(self, u"BIDS", u"Set the path to the BIDS database in the Prefs tab.")
            return
        # Parse subfolders to list subjects
        subBids = [f[4:] for f in os.listdir(self.prefs['bids']) if ((len(f) > 4) and (f[:4] == 'sub-') and os.path.isdir(os.path.join(self.prefs['bids'], f)) )]
        if not subBids:
            return
        # Get list of subjects in database
        rdi = ReadDiskItem( 'Subject', 'Directory', requiredAttributes={'center':str(self.ui.bvProtocolCombo.currentText())} )
        subBvFiles = list(rdi._findValues( {}, None, False))
        if not subBvFiles:
            subBv = [os.path.split(s.attributes()['subject'])[1] for s in subBvFiles]
        # Difference between lists
        subMissing = list(set(subBids) - set(subBv))
        # No subjects to import
        if not subMissing:
            QtGui.QMessageBox.warning(self, "Error", u"All BIDS subjects are already imported in the IntrAnat database.")
        # Ask user to select which ones to import
        dialog = DialogCheckbox.DialogCheckbox(subMissing, "Import BIDS", "Select the patients to import:")
        selSub = dialog.exec_()
        
        # Loop on patients to import
        for i in range(length(subMissing)):
            # Skip ignore subjects
            if not selSub[i]:
                continue
            # Create patient in database
            print("Import subject: " + subMissing(i))
            
    
    # ===== ANATOMIST =====    
    def clearAnatomist(self, windows=None):
        """ If "windows" is provided, just empties the provided windows.
            If not, all loaded objects are closed and all windows are emptied """
        # Delete all graphical objects
        if windows is not None:
            self.a.removeObjects(self.dispObj, windows)
            return
        self.a.removeObjects(self.dispObj, self.wins)
        self.a.deleteObjects(self.dispObj)
        self.dispObj = []
        # Remove ununsed referentials
        referentials = self.a.getReferentials()
        for element in referentials:
            if element.getInfos().get('name') not in ('Talairach-MNI template-SPM', 'Talairach-AC/PC-Anatomist'):
                self.a.deleteElements(element)


    def displayAnatomist(self, win1Path=None, win2Path=None):
        self.clearAnatomist()
        if win1Path:
            im1 = self.a.loadObject(win1Path)
            self.a.addObjects(im1, self.wins[0])
            self.dispObj.append(im1)
        if win2Path:
            im2 = self.a.loadObject(win2Path)
            self.a.addObjects(im2, self.wins[1])
            self.dispObj.append(im2)

    # ===== TAB4: IMPORT =====
    def setBrainCenter(self):
        brainCenterCoord = list(self.a.linkCursorLastClickedPosition())
        self.ui.brainCenterLabel.setText("%.1f, %.1f, %.1f"%tuple([float(x) for x in brainCenterCoord][:3]))
        self.brainCenter = brainCenterCoord

    def importNifti(self, path = None, patient = None, proto=None, modality = None, acq = None):
        """ Imports a Nifti file in the BrainVisa database.
        path is the input image, patient is e.g. Dupont_Jean, acq is the acquisition name e.g. 'T1pre_2000-01-01'
        If the parameters are empty, the data is taken from the Nifti Tab in the UI"""
    
        if self.brainCenter == None:
            print "brain center not set"
            QtGui.QMessageBox.warning(self, "Error",u"You haven't selected the BrainCenter !")
            return
        if self.ui.acqDate.date() == self.defaultAcqDate:
            QtGui.QMessageBox.warning(self, "Error",u"Acquisition date is not valid!")
            return
        if path is None:
            path = str(self.ui.niftiFileLabel.text())
        if not path:
            QtGui.QMessageBox.warning(self, "Error", u"Choose a file to import !")
            return
        filename = os.path.split(path)[1]
        if patient is None:
            patient = str(self.ui.niftiSubjectCombo.currentText())
        if proto is None:
            proto = str(self.ui.niftiProtocolCombo.currentText())
        if acq is None:
            d = self.ui.acqDate.date()
            acq = str(self.ui.niftiSeqType.currentText() + self.ui.niftiAcqType.currentText() + '_'+str(d.year()) + '-' + str(d.month()) + '-' + str(d.day()))
        if modality is None:
            # Determines the BrainVisa file type (depends on the imaging modality)
            modality = str(self.ui.niftiSeqType.currentText()) # T1 / T2 / CT / PET /fMRI
    
        filetype = {'T1':'Raw T1 MRI', 'T2':'T2 MRI', 'CT':'CT', 'TEP':'PET', 'PET':'PET', 'fMRI':'fMRI-epile', 'Statistics':'Statistic-Data', 'FLAIR':'FLAIR', 'FreesurferAtlas':'FreesurferAtlas', 'FGATIR':'FGATIR'}[modality]
    
        write_filters = { 'center': proto, 'acquisition': str(acq), 'subject' : str(patient) }
    
        if filetype == 'fMRI-epile' or filetype == 'Statistic-Data':
            text, ok = QtGui.QInputDialog.getText(self, 'Input Dialog', 'Enter the name of the test:')
            if (ok & bool(text)):
                write_filters.update({'subacquisition':str(text)})
    
    
        wdi = WriteDiskItem(filetype, 'NIFTI-1 image' )#'gz compressed NIFTI-1 image' )
        print "Finding BV path for acquisition %s and patient %s"%(acq, patient)
        di = wdi.findValue(write_filters)
    
        if di is None:
            QtGui.QMessageBox.warning(self, "Error", u"Impossible to find a valid path to import the Image in Brainvisa database (%s, %s)"%(patient, acq))
            return
        # Copy the file : lancer l'importation standard de T1 MRI pour la conversion de format et
        if os.path.exists(di.fileName()):
            reply = QtGui.QMessageBox.question(self, 'Overwrite', u"This image already exists.\nOverwrite the file ?", QtGui.QMessageBox.Yes | QtGui.QMessageBox.No, QtGui.QMessageBox.No)
            if reply == QtGui.QMessageBox.No:
                return
        print "Importing file as "+di.fileName()
    
        #current images loaded for this patient:
        self.selectBvSubject(str(patient))
    
        if filetype != 'fMRI-epile' and filetype != 'Statistic-Data' and filetype != 'FreesurferAtlas':
            if len(self.bvImages) > 0:
                ImAlreadyHere = [i for i in range(len(self.bvImages)) if str(self.ui.niftiSeqType.currentText() + self.ui.niftiAcqType.currentText()+'_') in self.bvImages.keys()[i]]
                if len(ImAlreadyHere):
                    QtGui.QMessageBox.warning(self, 'WARNING', u"There is already a %s image, delete it before importing a new one."%str(self.ui.niftiSeqType.currentText() + self.ui.niftiAcqType.currentText()))
                    self.setStatus(u"Sequence %s importation not performed (already one equivalent)"%acq)
                    return
    
        # Call import in a separate function with a progress bar
        res = ProgressDialog.call(lambda thr:self.importNiftiWorker(di, path, filetype, thr), True, self, "Processing...", "Import image", False)
        # self.importNiftiWorker(di, path, filetype, proto, patient)
        
        di.setMinf('brainCenter', self.brainCenter)
        if filetype == 'Statistic-Data':
            textMNI, okMNI = QtGui.QInputDialog.getItem(self,'MNI or not','Is the statistic image generated in the MNI',['True','False'])
            if (okMNI & bool(textMNI)):
                di.setMinf('MNI',str(textMNI))
            di.setMinf('ColorPalette','Yellow-Red-White-Blue-Green')
        neuroHierarchy.databases.insertDiskItem( di, update=True )

        if filetype == 'Statistic-Data':
            try:
                self.StatisticDataMNItoScannerBased(proto, patient, acq)
            except:
                pass

        self.setStatus(u"Sequence %s importation done"%acq)
        
        
        
    def importNiftiWorker(self, di, path, filetype, thr=None):
        # Create directories that do not exist yet
        self.createItemDirs(di)
    
        # We check that there is a Scanner based transformation, if not and if both transformations are equal, we rename one as scanner based transformation
        temp_nii = aims.read(path)
        refs_nii = temp_nii.header()['referentials']
        if len(refs_nii)>1:
            isScannerBased = 0
            sum_transfo = numpy.zeros(16,dtype='float32')
            for ii in range(len(refs_nii)):
                if refs_nii[ii] == 'Scanner-based anatomical coordinates':
                    isScannerBased = 1
                sum_transfo = sum_transfo + temp_nii.header()['transformations'][ii]
            mean_transfo = sum_transfo/len(refs_nii)
            if isScannerBased == 0:
                diff = mean_transfo - temp_nii.header()['transformations'][0]
                if numpy.absolute(diff).round(2).sum() <= 0.03:
                    temp_nii.header()['referentials'][0] = 'Scanner-based anatomical coordinates'
        else:
            if refs_nii[0] == 'Coordinates aligned to another file or to anatomical truth':
                temp_nii.header()['referentials'][0] = 'Scanner-based anatomical coordinates'
    
        destination = di.fileName()
        aims.write(temp_nii, destination)
    
        if (filetype != 'fMRI-epile') & (filetype != 'PET') & (filetype != 'Statistic-Data'):  # or filetype != 'MTT'
            ret = subprocess.call(['AimsFileConvert', '-i', str(destination), '-o', str(destination), '-t', 'S16'])
        elif filetype == 'PET':
            print "conversion PET float to S16"
            ret = subprocess.call(['AimsFileConvert', '-i', str(destination), '-o', str(destination), '-t', 'S16', '-r', 'True'])
            di.setMinf('ColorPalette','Blue-Red-fusion')
        else:
            print "no conversion to grayscale"
            ret = subprocess.call(['AimsFileConvert', '-i', str(destination), '-o', str(destination)])
        if ret < 0:
            print "Importation error: BrainVisa / AimsFileConvert"
            return

        if str(self.ui.niftiSeqType.currentText()) == 'T1':
            if self.ui.radioButtonGado.isChecked():
                di.setMinf('Gado',True)
            if self.ui.radioButtonNoGado.isChecked():
                di.setMinf('Gado',False)
            neuroHierarchy.databases.insertDiskItem( di, update=True )

        
    def importFSoutput(self,subject=None):
        # Get current subject
        if not subject:
            subject = self.currentSubject
        # Find T1pre of the subject
        rT1BV = ReadDiskItem('Raw T1 MRI', 'BrainVISA volume formats',requiredAttributes={'subject':subject})
        allT1 = list(rT1BV.findValues({},None,False))
        if not allT1:
            QtGui.QMessageBox.warning(self, "Error", u"No T1 MRI found for patient: " + subject)
            return
        idxT1pre = [i for i in range(len(allT1)) if 'T1pre' in str(allT1[i])]
        if not idxT1pre:
            QtGui.QMessageBox.warning(self, "Error", u"No T1pre found for this patient: " + subject)
            return
        diT1pre = allT1[idxT1pre[0]]

        # Find FreeSurfer database
        FsSubjDir = None
        for db in neuroHierarchy.databases.iterDatabases():
            if db.directory.lower().find("freesurfer") != -1:
                FsSubjDir = db.directory
                print("FreeSurfer database folder: " + FsSubjDir)
                break
        if not FsSubjDir:
            QtGui.QMessageBox.warning(self, "Error", u"No local FreeSurfer database found.")
            return
        # Ask file name
        importDir = str(QtGui.QFileDialog.getExistingDirectory(self, "Select FreeSurfer folder", "", QtGui.QFileDialog.ShowDirsOnly))
        if not importDir:
            return
        
        # Get all the filenames expected in the FreeSurfer output folder
        allFiles = dict()
        allFiles['anat'] =           {'side':None,    'type':'RawFreesurferAnat',          'format':'FreesurferMGZ',              'file':importDir + '/mri/orig/001.mgz'}
        allFiles['T1_orig'] =        {'side':None,    'type':'T1 FreesurferAnat',          'format':'FreesurferMGZ',              'file':importDir + '/mri/orig.mgz'}
        allFiles['nu'] =             {'side':None,    'type':'Nu FreesurferAnat',          'format':'FreesurferMGZ',              'file':importDir + '/mri/nu.mgz'}
        allFiles['ribbon'] =         {'side':None,    'type':'Ribbon Freesurfer',          'format':'FreesurferMGZ',              'file':importDir + '/mri/ribbon.mgz'}
        allFiles['destrieux'] =      {'side':None,    'type':'Freesurfer Cortical Parcellation using Destrieux Atlas',  'format':'FreesurferMGZ',  'file':importDir + '/mri/aparc.a2009s+aseg.mgz'}
        allFiles['aseg'] =           {'side':None,    'type':'Freesurfer aseg',            'format':'FreesurferMGZ',              'file':importDir + '/mri/aseg.mgz'}
        allFiles['lausanne33'] =     {'side':None,    'type':None,                          'format':None,                        'file':importDir + '/parcellation_Lausanne2008/ROIv_HR_th_scale33.nii.gz'}
        allFiles['lausanne60'] =     {'side':None,    'type':None,                          'format':None,                        'file':importDir + '/parcellation_Lausanne2008/ROIv_HR_th_scale60.nii.gz'}
        allFiles['lausanne125'] =    {'side':None,    'type':None,                          'format':None,                        'file':importDir + '/parcellation_Lausanne2008/ROIv_HR_th_scale125.nii.gz'}
        allFiles['lausanne250'] =    {'side':None,    'type':None,                          'format':None,                        'file':importDir + '/parcellation_Lausanne2008/ROIv_HR_th_scale250.nii.gz'}
        allFiles['lausanne500'] =    {'side':None,    'type':None,                          'format':None,                        'file':importDir + '/parcellation_Lausanne2008/ROIv_HR_th_scale500.nii.gz'}
        allFiles['xfm'] =            {'side':None,    'type':'Talairach Auto Freesurfer',  'format':'MINC transformation matrix', 'file':importDir + '/mri/transforms/talairach.auto.xfm'}
        allFiles['leftPial'] =       {'side':'left',  'type':'BaseFreesurferType',         'format':'FreesurferPial',             'file':importDir + '/surf/lh.pial'}
        allFiles['leftWhite'] =      {'side':'left',  'type':'BaseFreesurferType',         'format':'FreesurferWhite',            'file':importDir + '/surf/lh.white'}
        allFiles['leftSphereReg'] =  {'side':'left',  'type':'BaseFreesurferType',         'format':'FreesurferSphereReg',        'file':importDir + '/surf/lh.sphere.reg'}
        allFiles['leftThickness'] =  {'side':'left',  'type':'BaseFreesurferType',         'format':'FreesurferThickness',        'file':importDir + '/surf/lh.thickness'}
        allFiles['leftCurv'] =       {'side':'left',  'type':'BaseFreesurferType',         'format':'FreesurferCurv',             'file':importDir + '/surf/lh.curv'}
        allFiles['leftAvgCurv'] =    {'side':'left',  'type':'BaseFreesurferType',         'format':'FreesurferAvgCurv',          'file':importDir + '/surf/lh.avg_curv'}
        allFiles['leftCurvPial'] =   {'side':'left',  'type':'BaseFreesurferType',         'format':'FreesurferCurvPial',         'file':importDir + '/surf/lh.curv.pial'}
        allFiles['leftGyri'] =       {'side':'left',  'type':'FreesurferGyriTexture',      'format':'FreesurferParcellation',     'file':importDir + '/label/lh.aparc.annot'}
        allFiles['leftSulciGyri'] =  {'side':'left',  'type':'FreesurferSulciGyriTexture', 'format':'FreesurferParcellation',     'file':importDir + '/label/lh.aparc.a2009s.annot'}
        allFiles['rightPial'] =      {'side':'right', 'type':'BaseFreesurferType',         'format':'FreesurferPial',             'file':importDir + '/surf/rh.pial'}
        allFiles['rightWhite'] =     {'side':'right', 'type':'BaseFreesurferType',         'format':'FreesurferWhite',            'file':importDir + '/surf/rh.white'}
        allFiles['rightSphereReg'] = {'side':'right', 'type':'BaseFreesurferType',         'format':'FreesurferSphereReg',        'file':importDir + '/surf/rh.sphere.reg'}
        allFiles['rightThickness'] = {'side':'right', 'type':'BaseFreesurferType',         'format':'FreesurferThickness',        'file':importDir + '/surf/rh.thickness'}
        allFiles['rightCurv'] =      {'side':'right', 'type':'BaseFreesurferType',         'format':'FreesurferCurv',             'file':importDir + '/surf/rh.curv'}
        allFiles['rightAvgCurv'] =   {'side':'right', 'type':'BaseFreesurferType',         'format':'FreesurferAvgCurv',          'file':importDir + '/surf/rh.avg_curv'}
        allFiles['rightCurvPial'] =  {'side':'right', 'type':'BaseFreesurferType',         'format':'FreesurferCurvPial',         'file':importDir + '/surf/rh.curv.pial'}
        allFiles['rightGyri'] =      {'side':'right', 'type':'FreesurferGyriTexture',      'format':'FreesurferParcellation',     'file':importDir + '/label/rh.aparc.annot'}
        allFiles['rightSulciGyri'] = {'side':'right', 'type':'FreesurferSulciGyriTexture', 'format':'FreesurferParcellation',     'file':importDir + '/label/rh.aparc.a2009s.annot'}
        # Check that all the files exist (skip the lausanne files)
        for key in allFiles:
            if (not 'lausanne' in key) and (not os.path.isfile(allFiles[key]['file'])):
                QtGui.QMessageBox.warning(self, "Error", u"FreeSurfer file not found:\n" + allFiles[key]['file'])
                return
        
        # Run copy and conversion in a separate thread
        res = ProgressDialog.call(lambda thr:self.importFSoutputWorker(FsSubjDir, subject, allFiles, diT1pre, thr), True, self, "Processing...", "Import FreeSurfer output")
        #res = self.importFSoutputWorker(FsSubjDir, subject, allFiles, diT1pre)
        

    def importFSoutputWorker(self, FsSubjDir, subject, allFiles, diT1pre, thread=None):
        # Where to copy the new files
        acq =  str(diT1pre.attributes()['acquisition']).replace('T1','FreesurferAtlas')
        write_filters = { 'center': str(self.ui.niftiProtocolCombo.currentText()), 'acquisition': str(acq), 'subject' : subject }
        # Progress bar
        if thread:
            thread.emit(QtCore.SIGNAL("PROGRESS_TEXT"), "Copying files to freesurfer database...")
        # Add all the files to the local freesurfer database
        for key in allFiles:
            # Skip files that do not have a format or placeholder in the freesurfer ontology
            if not allFiles[key]['type']:
                continue
            # Get target into the local FreeSurfer database
            if allFiles[key]['side']:
                wdi = WriteDiskItem(allFiles[key]['type'], allFiles[key]['format'], requiredAttributes={'subject':subject,'_ontology':'freesurfer', 'side':allFiles[key]['side']})
            else:
                wdi = WriteDiskItem(allFiles[key]['type'], allFiles[key]['format'], requiredAttributes={'subject':subject,'_ontology':'freesurfer'})
            # Type 'T1 FreesurferAnat' returns both nu.mgz and orig.mgz: need to filter the list
            if (key == 'T1_orig'):
                di = wdi.findValues(write_filters, None, True)
                allT1 = list(di)
                idxT1orig = [i for i in range(len(allT1)) if 'orig.mgz' in str(allT1[i])]
                di = allT1[idxT1orig[0]]
            else:
                di = wdi.findValue(write_filters)
            # If there is an error while finding where to save the file
            if not di:
                print("Error: Cannot import file: " + allFiles[key]['file'])
                return
            # Create target folder
            self.createItemDirs(di)
            # Copy file into the local FS datbase
            print("Copy: " + allFiles[key]['file'] + " => " + di.fullPath())
            copyfile(allFiles[key]['file'], di.fullPath())
            # Add reference in the database (creates .minf)
            neuroHierarchy.databases.insertDiskItem(di, update=True)
       
        # Progress bar
        if thread:
            thread.emit(QtCore.SIGNAL("PROGRESS_TEXT"), "Converting Destrieux atlas...")
        # Importing Destrieux atlas to BrainVISA database
        wdi = WriteDiskItem('FreesurferAtlas', 'NIFTI-1 image')
        diFSDestrieux = wdi.findValue(write_filters)
        # Create the folder if doesn't exist
        self.createItemDirs(diFSDestrieux)
        # Reslice volume to match T1pre
        launchFreesurferCommand(self.brainvisaContext, None, 'mri_convert', '-i',str(allFiles['destrieux']['file']),'-o',str(diFSDestrieux.fullPath()),'-rl',str(diT1pre.fullPath()),'-rt','nearest','-nc')
        # Convert to AIMS
        ret = subprocess.call(['AimsFileConvert', '-i', str(diFSDestrieux.fullPath()), '-o', str(diFSDestrieux.fullPath()), '-t', 'S16'])
        # Add reference in the database (creates .minf)
        neuroHierarchy.databases.insertDiskItem(diFSDestrieux, update=True)

        # Importing all the Lausanne2008 parcellations to BrainVISA database
        for n in [33, 60, 125, 250, 500]:
            # File reference
            key = 'lausanne{}'.format(n)
            # Skip if file doesn't exist
            if (not os.path.isfile(allFiles[key]['file'])):
                continue
            # Progress bar
            if thread:
                thread.emit(QtCore.SIGNAL("PROGRESS_TEXT"), "Importing Lausanne2008-{} parcellation...".format(n))
            # Where to copy the new files
            acqLaus = str(diT1pre.attributes()['acquisition']).replace('T1','Lausanne2008-{}-'.format(n))
            # Importing Destrieux atlas to BrainVISA database
            wdi = WriteDiskItem('FreesurferAtlas', 'NIFTI-1 image')
            diLaus = wdi.findValue({'center' : write_filters['center'], 'acquisition' : acqLaus, 'subject' : subject})
            # Create the folder if doesn't exist
            self.createItemDirs(diLaus)

            # NEW VERSION: 25-Oct-2018
            # Lausanne2008 atlases are not resliced on the original T1 anymore, and stay in FreeSurfer format (256x256x256)
            # They are resliced on the T1pre at the moment of their importation in the IntrAnat/BrainVISA database,
            # just like the FreeSurfer/Destrieux atlas (previous section of this function)
            
            # Reslice volume to match T1pre
            launchFreesurferCommand(self.brainvisaContext, None, 'mri_convert', '-i',str(allFiles[key]['file']),'-o',str(diLaus.fullPath()),'-rl',str(diT1pre.fullPath()),'-rt','nearest','-nc')
            # Convert to AIMS
            ret = subprocess.call(['AimsFileConvert', '-i', str(diLaus.fullPath()), '-o', str(diLaus.fullPath()), '-t', 'S16'])
            # Add reference in the database (creates .minf)
            neuroHierarchy.databases.insertDiskItem(diLaus, update=True)
        
        # Generate amygdala and hippocampus meshes
        if thread:
            thread.emit(QtCore.SIGNAL("PROGRESS_TEXT"), "Creating hippcampus and amygdala meshes...")
        self.generateAmygdalaHippoMesh(str(self.ui.niftiProtocolCombo.currentText()), subject, acq, diFSDestrieux)


    # Sub-selection of the operations done when importing the full FreeSurfer segmentation
    # Added only for O.David to add the Lausanne 2008 parcellation to patients for which FreeSurfer was already imported
    def importLausanne2008(self):
        # Get current subject
        subject = self.currentSubject
        # Find T1pre of the subject
        rT1BV = ReadDiskItem('Raw T1 MRI', 'BrainVISA volume formats',requiredAttributes={'subject':subject})
        allT1 = list(rT1BV.findValues({},None,False))
        idxT1pre = [i for i in range(len(allT1)) if 'T1pre' in str(allT1[i])]
        if not idxT1pre:
            QtGui.QMessageBox.warning(self, "Error", u"No T1pre found for this patient: " + subject)
            return
        diT1pre = allT1[idxT1pre[0]]

        # Ask input folder
        importDir = str(QtGui.QFileDialog.getExistingDirectory(self, "Select Lausanne2008 folder", "", QtGui.QFileDialog.ShowDirsOnly))
        if not importDir:
            return
        # Get files of interest in the input folder
        allFiles = {'lausanne33'  : importDir + '/ROIv_HR_th_scale33.nii.gz',\
                    'lausanne60'  : importDir + '/ROIv_HR_th_scale60.nii.gz',\
                    'lausanne125' : importDir + '/ROIv_HR_th_scale125.nii.gz',\
                    'lausanne250' : importDir + '/ROIv_HR_th_scale250.nii.gz',\
                    'lausanne500' : importDir + '/ROIv_HR_th_scale500.nii.gz'}
        # Check that all the files exist (skip the lausanne files)
        for key in allFiles:
            if not os.path.isfile(allFiles[key]):
                QtGui.QMessageBox.warning(self, "Error", u"Lausanne2008 file not found:\n" + allFiles[key])
                return
        # Run copy and conversion in a separate thread
        res = ProgressDialog.call(lambda thr:self.importLausanne2008Worker(subject, allFiles, diT1pre, thr), True, self, "Processing...", "Import FreeSurfer output")

        
    def importLausanne2008Worker(self, subject, allFiles, diT1pre, thread=None):
        # Importing all the Lausanne2008 parcellations to BrainVISA database
        for n in [33, 60, 125, 250, 500]:
            # File reference
            key = 'lausanne{}'.format(n)
            # Skip if file doesn't exist
            if (not os.path.isfile(allFiles[key])):
                continue
            # Progress bar
            if thread:
                thread.emit(QtCore.SIGNAL("PROGRESS_TEXT"), "Importing Lausanne2008-{} parcellation...".format(n))
            # Where to copy the new files
            acqLaus = str(diT1pre.attributes()['acquisition']).replace('T1','Lausanne2008-{}-'.format(n))
            # Importing Destrieux atlas to BrainVISA database
            wdi = WriteDiskItem('FreesurferAtlas', 'NIFTI-1 image')
            diLaus = wdi.findValue({'center' : str(self.ui.niftiProtocolCombo.currentText()), 'acquisition' : acqLaus, 'subject' : subject})
            # Create the folder if doesn't exist
            self.createItemDirs(diLaus)
            # Reslice volume to match T1pre
            launchFreesurferCommand(self.brainvisaContext, None, 'mri_convert', '-i',str(allFiles[key]),'-o',str(diLaus.fullPath()),'-rl',str(diT1pre.fullPath()),'-rt','nearest','-nc')
            # Convert to AIMS
            ret = subprocess.call(['AimsFileConvert', '-i', str(diLaus.fullPath()), '-o', str(diLaus.fullPath()), '-t', 'S16'])
            # Add reference in the database (creates .minf)
            neuroHierarchy.databases.insertDiskItem(diLaus, update=True)


    #******************************** Add New Subject
    # ANONYMIZATION TAKES PLACE HERE
    def updatePatientCode(self):
        # <site>_<year>_<three letters from name><one letter from first name><number if duplicate entries, but should be entered manually anyway>
        if self.prefs['projectSelected'][0] == 1 or self.prefs['projectSelected'][0] == 2:
            print "neuropsynov or classique"
            self.ui.subjectPatientCode.setText(self.ui.subjectSiteCombo.currentText() + '_' + str(self.ui.subjectYearSpinbox.value()) + '_'  + str(self.ui.subjectPatientName.text())[:3].upper() + str(self.ui.subjectPatientFirstName.text())[:1].lower())
        if self.prefs['projectSelected'][0] == 0:
            print "f-Tract"
            self.ui.subjectPatientCode.setText(str(self.ui.subjectPatientName.text()).upper() + str(self.ui.subjectSiteCombo.currentText()).upper())
        if self.prefs['projectSelected'][0] == 3:
            print "other"
            self.ui.subjectPatientCode.setText(self.ui.subjectSiteCombo.currentText() + '_' + str(self.ui.subjectYearSpinbox.value()) + '_'  + str(self.ui.subjectPatientName.text())[:3].upper() + str(self.ui.subjectPatientFirstName.text())[:1].lower())

    def storePatientInDB(self):
        anon = str(self.ui.subjectPatientCode.text())
        if self.prefs['projectSelected'][0] == 1 or self.prefs['projectSelected'][0] == 2:
            if len(anon) < 10:
                QtGui.QMessageBox.warning(self, u"Subject adding", u"Impossible to add the subject, subject identifier is too short and doesn't match the subject identifier model chosen in preferences")
                return
        if self.prefs['projectSelected'][0] == 0:
            if len(anon) < 16:
                QtGui.QMessageBox.warning(self, u"Subject adding", u"Impossible to add the subject, subject identifier is too short and doesn't match the subject identifier model chosen in preferences")
                return
            elif len(anon) > 16:
                QtGui.QMessageBox.warning(self, u"Subject adding", u"Impossible to add the subject, subject identifier is too long and doesn't match the subject identifier model chosen in preferences")
                return
            else:
                if not (anon[0:4].isdigit() and not anon[4:8].isdigit() and anon[8:].isdigit() and anon[7] == '_'):
                    QtGui.QMessageBox.warning(self, u"Subject adding", u"Impossible to add the subject, subject identifier doesn't match the subject identifier model chosen in preferences: ex: 0001GRE_yyyymmdd")
                    return
                else:
                    if not int(anon[8:12]) > 1950 and not int(anon[8:12]) < 2100 and not int(anon[12:14]) >0 and not int(anon[12:14]) < 13 and not int(anon[14:16]) > 0 and not int(anon[14:16]) < 32:
                        QtGui.QMessageBox.warning(self, u"Subject adding", u"Impossible to add the subject, subject identifier doesn't match the subject identifier model chosen in preferences: ex: 0001GRE_yyyymmdd")
                        return
    
    
        wdi = WriteDiskItem( 'Subject', 'Directory' )
        di = wdi.findValue( { 'center': str(self.ui.subjectProtocolCombo.currentText()), 'subject' : anon } )
        if di is None:
          QtGui.QMessageBox.warning(self, u"Error", u"Impossible to find a valid path for the patient")
          self.setStatus(u"Subject add failed : the path is not valid")
          return
        # Create directories that do not exist yet
        self.createItemDirs(di)
        try:
          os.mkdir(di.fileName())
          neuroHierarchy.databases.insertDiskItem( di, update=True )
          self.setStatus(u"Subject %s added to the database"%anon)
          self.currentSubject = anon
        except:
          self.setStatus(u"Subject adding failed : impossible to create subject folder")
    
        self.ui.tabWidget.setCurrentIndex(2)

    # ******************************* Registration

    def selectRegSubject(self, subj):
        """ A BrainVisa subject was selected : query the database to get the list of images"""
        # Display "Date : XXXX-XX-XX - Seq: T1 - Acq : T1Pre
        self.ui.regImageList.clear()
        self.ui.regImageList2.clear()
        images = self.findAllImagesForSubject(self.ui.regProtocolCombo.currentText(), subj)
        #name = modality + acquisition + subacquisition if exist subacquisition key else modality + acquisition
    
        dict_temporaire1 = {}
        for i in images:
            if 'subacquisition' in i.attributes().keys():
                dict_temporaire1.update({i.attributes()['modality'] + ' - ' + i.attributes()['acquisition'] + ' - ' + i.attributes()['subacquisition']:i.fileName()})
            else:
                dict_temporaire1.update({i.attributes()['modality'] + ' - ' + i.attributes()['acquisition']:i.fileName()})
    
        dict_temporaire2 = {}
        for i in images:
            if 'subacquisition' in i.attributes().keys():
                dict_temporaire2.update({i.attributes()['modality'] + ' - ' + i.attributes()['acquisition'] + ' - ' + i.attributes()['subacquisition']:i})
            else:
                dict_temporaire2.update({i.attributes()['modality'] + ' - ' + i.attributes()['acquisition']:i})
    
        self.ui.regImageList.addItems(sorted([i.attributes()['modality'] + ' - '+ i.attributes()['acquisition'] + ' - ' + i.attributes()['subacquisition'] if 'subacquisition' in i.attributes().keys()\
                else i.attributes()['modality'] + ' - '+ i.attributes()['acquisition'] for i in images ]))
        self.ui.regImageList2.addItems(sorted([i.attributes()['modality'] + ' - '+ i.attributes()['acquisition'] + ' - ' + i.attributes()['subacquisition'] if 'subacquisition' in i.attributes().keys()\
                else i.attributes()['modality'] + ' - '+ i.attributes()['acquisition'] for i in images ]))
    
    
        self.regImagePaths = dict_temporaire1
        self.regImages = dict_temporaire2


    def selectRegImage(self, item):
        """ A BrainVisa image was double-clicked : display it !"""
        image = self.regImages[str(item.text())]
        self.displayImage(image, self.wins[:2])

    def selectRegImage2(self, item):
        """ A BrainVisa image was double-clicked in the second list: display it !"""
        image = self.regImages[str(item.text())]
        self.displayImage(image, self.wins[2:])

    
    def displayImage(self, image, wins):
        # Clear existing windows
        self.clearAnatomist(wins)
        # Load image
        mri = self.a.loadObject(image)
        # Force reading Scanner-based referential from .nii file
        print "WARNING: Using referential from .nii header, ignoring transformations from BrainVISA database."
        self.a.execute('LoadReferentialFromHeader', objects=[mri])
        # Add to anatomist windows
        self.a.assignReferential(mri.getReferential(), wins)
        self.a.addObjects(mri, wins)
        self.dispObj.append(mri)
        
        # Guess center of image
        attr = mri.getInternalRep().attributed()
        if (attr['volume_dimension'] and attr['voxel_size']):
            volSize = attr['volume_dimension']
            voxSize = attr['voxel_size']
            center = [volSize[0]*voxSize[0]/2, volSize[1]*voxSize[1]/2, volSize[2]*voxSize[2]/2];
        else:
            center = [128, 128, 128]
        # Center view on center of image
        wins[0].moveLinkedCursor(center)


    def registerNormalizeSubject(self, progressThread=None):
        """ Registers all images of the subject with SPM, then store the transforms in BrainVisa database, and launch T1pre morphologist analysis (brain segmentation) """

        proto = str(self.ui.regProtocolCombo.currentText())
        subj = str(self.ui.regSubjectCombo.currentText())
        images = self.findAllImagesForSubject(proto, subj)
        
        # Call process in a different thread
        errMsg = ProgressDialog.call(lambda thr:self.registerNormalizeSubjectWorker(proto, subj, images, thr), True, self, "Processing...", "Coregister and normalize")
        #errMsg = self.registerNormalizeSubjectWorker(proto, subj, images)
        
        # Display error messages
        if errMsg:
            QtGui.QMessageBox.critical(self, u'Regitration error', u"Errors occured during the normalization or registration: \n\n" + u"\n".join(errMsg))

    
    def registerNormalizeSubjectWorker(self, proto, subj, images, progressThread=None):
        # NORMALIZE: Find all images, normalize the T1s, find the T1pre if there is one, read image headers, store their referentials and transformations in the DB
        t1preImage = None
        for image in images:
            patient = image.attributes()['subject']
            acq = image.attributes()['acquisition']
            # Store Scanner-based referential and referential files in the DB
            if not image.attributes()['modality'] == 'statistic_data' and not image.attributes()['modality'] == 'freesurfer_atlas':   # and not image.attributes()['modality'] == 'hippofreesurfer_atlas':
                print "do the coregistration"
                self.storeImageReferentialsAndTransforms(image)
            # If this is a T1, normalize it to MNI (pre or post). If the pre is badly normalized, "using the post" should be stored in the DB
            # if acq.startswith('T1'):
            # FT 9-Nov-2018: Logic change: normalize only the T1pre volumes, the other ones are useless
            if acq.startswith('T1pre'):
                self.setStatus(u"SPM normalization %s..."%acq)
                if progressThread:
                    progressThread.emit(QtCore.SIGNAL("PROGRESS_TEXT"), "SPM normalization: " + subj + "/" + image.attributes()['modality'] + "...")
                errMsg = self.spmNormalize(image.fileName(), proto, patient, acq)
                if errMsg:
                    return errMsg
                self.taskfinished(u"SPM normalization done")
                # If there is a T1pre, remember the image
                #if acq.find('T1pre') == 0:
                t1preImage = image
        # No T1pre : nothing else to do
        if t1preImage is None:
            return
    
        # COREGISTER
        self.setStatus(u"Coregistration of all images to T1pre...")
        for image in images:
            acq = image.attributes()['acquisition']
            patient = image.attributes()['subject']
            # A T1pre is there, coregister all images to it
            if image == t1preImage:
                continue
            # FreeSurfer MRI: Skip this step, transformation was saved at the import time
            elif ('FreesurferAtlas' in image.attributes()['acquisition']) or ('FreesurferAtlas' in acq):
                print("Coregistration: Skipping " + image.attributes()['acquisition'] + "...")
                continue
            # Statistics: Same referential as t1pre
            elif image.attributes()['modality'] == 'statistic_data' or image.attributes()['modality'] == 'freesurfer_atlas':  # or image.attributes()['modality'] == 'hippofreesurfer_atlas':
                print "Coregistration: Attribute T1pre referential to this modality {}".format(image.attributes()['modality'])
                self.transfoManager.setReferentialTo(image, t1preImage.attributes()['referential'] )
                continue
    
            # ===== ANTS =====
            if self.coregMethod == 'ANTs':
                print("Coregistration method: ANTs")
                if progressThread:
                    progressThread.emit(QtCore.SIGNAL("PROGRESS_TEXT"), "ANTS coregistration: " + subj + "/" + image.attributes()['modality'] + "...")
                temp_folder_ants = tempfile.mkdtemp('ANTs_IntrAnat') +'/'
                ants_call = 'antsRegistrationSyN.sh -d 3 -f {} -m {} -o {} -t r'.format(str(t1preImage.fullPath()),str(image.fullPath()),temp_folder_ants)
                # print("ANTs call: " + ants_call)
                # Run ANTs system command
                runCmd(ants_call.split())
                # Register transformation in the database
                self.setANTstrm_database(image, temp_folder_ants)
                
            # ===== SPM =====
            elif self.coregMethod == 'spm':
                print("Coregistration method: SPM")
                # SPM coregister
                if self.ui.regResampleCheck.isChecked():
                    if progressThread:
                        progressThread.emit(QtCore.SIGNAL("PROGRESS_TEXT"), "SPM coregistration+resample: " + subj + "/" + image.attributes()['modality'] + "...")
                    # Call SPM coregister+resample
                    call = spm_coregisterReslice%("'"+str(self.prefs['spm'])+"'","{'"+str(t1preImage)+",1'}", "{'"+str(image.fileName())+",1'}")
                    spl = os.path.split(image.fileName())
                    registeredPath = os.path.join(spl[0], 'r'+spl[1])
                    errMsg = matlabRun(call)
                    if errMsg:
                        return errMsg
                    # Update resampled volume
                    self.setResampledToT1pre(image, registeredPath)
                    
                # SPM coregister
                else:
                    if progressThread:
                        progressThread.emit(QtCore.SIGNAL("PROGRESS_TEXT"), "SPM coregistration: " + subj + "/" + image.attributes()['modality'] + "...")
                    self.setStatus('Coregistration of the image: '+ acq)

                    # Temporary txt file to store the trm transformation
                    tmpOutput = getTmpFilePath('txt')
                    imageFileName = image.fileName()
                    if ('data_type' in image.attributes().keys()) and (image.attributes()['data_type'] == 'RGB'):
                        print "it is RGB"
                        imageFileName = getTmpFilePath('nii')
                        ret = subprocess.call(['AimsFileConvert', '-i', str(image.fileName()), '-o', str(imageFileName), '-t', 'S16'])
                        if ret < 0:
                            print "Conversion to S16 error: "+repr(registeredPath) #terminal
                            return
                    # Run registration
                    if 'brainCenter' not in image.attributes() or 'brainCenter' not in t1preImage.attributes():
                        call = spm_coregister%("'"+str(self.prefs['spm'])+"'","'"+str(imageFileName)+",1'", "'"+str(t1preImage)+",1'", str([0, 0 ,0]), str([1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1]), str([0, 0, 0]), str([1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1]),"'"+tmpOutput+"'")
                    else:
                        call = spm_coregister%("'"+str(self.prefs['spm'])+"'","'"+str(imageFileName)+",1'", "'"+str(t1preImage)+",1'", str(image.attributes()['brainCenter']), str(image.attributes()['SB_Transform']), str(t1preImage.attributes()['brainCenter']), str(t1preImage.attributes()['SB_Transform']), "'"+tmpOutput+"'")
                    errMsg = matlabRun(call)
                    if errMsg:
                        return errMsg
                    # Register transformation in the database
                    self.insertTransformationToT1pre(tmpOutput, image)
                    if ('data_type' in image.attributes().keys()) and (image.attributes()['data_type'] == 'RGB'):
                        os.remove(imageFileName)
                        os.remove(imageFileName+'.minf')

        self.taskfinished(u"Coregistration done")
        # Clear all the views
        self.clearAnatomist()


    def deleteExistingMeshes(self, rdi, acq):
        delFiles = []
        subject = rdi.attributes()['subject']
        protocol = rdi.attributes()['center']
        # Get reference to current database
        db = neuroHierarchy.databases.database(rdi.get("_database"))
        # Delete folder default_analysis/segmentation/mesh/surface_analysis                
        delPath = os.path.dirname(os.path.abspath(rdi.fullPath())) + '/surface_analysis'
        if os.path.isdir(delPath):
            removeFromDB(delPath, db)
        # Get pial meshes
        rdi = ReadDiskItem('Hemisphere Mesh', 'Anatomist mesh formats', requiredAttributes={'subject':subject, 'center':protocol})
        delFiles += list(rdi._findValues({}, None, False))
        # Get white meshes
        rdi = ReadDiskItem('Hemisphere White Mesh', 'Anatomist mesh formats', requiredAttributes={'subject':subject, 'center':protocol})
        delFiles += list(rdi._findValues({}, None, False))
        # Get head meshe
        rdi = ReadDiskItem('Head Mesh', 'Anatomist mesh formats', requiredAttributes={'subject':subject, 'center':protocol})
        delFiles += list(rdi._findValues({}, None, False))
        # Delete all files
        for h in delFiles:
            # Delete only the ones in T1pri
            if acq in h.attributes()["acquisition"]:
                removeFromDB(h.fullPath(), db)
                # Delete associated .minf
                if os.path.isfile(h.fullPath() + '.minf'):
                    os.remove(h.fullPath() + '.minf')


    def deleteExistingNormalization(self, rdi, acq):
        delFiles = []
        subject = rdi.attributes()['subject']
        protocol = rdi.attributes()['center']
        # Get reference to current database
        db = neuroHierarchy.databases.database(rdi.get("_database"))
        # Get normalized volume
        rdi = ReadDiskItem('T1 SPM resampled in MNI', 'NIFTI-1 image', requiredAttributes={'subject':subject, 'center':protocol, 'acquisition':acq})
        delFiles += list(rdi._findValues({}, None, False))
        # Get white meshes
        rdi = ReadDiskItem('SPM normalization deformation field', 'NIFTI-1 image', requiredAttributes={'subject':subject, 'center':protocol, 'acquisition':acq})
        delFiles += list(rdi._findValues({}, None, False))
        # Get head meshe
        rdi = ReadDiskItem('SPM normalization inverse deformation field', 'NIFTI-1 image', requiredAttributes={'subject':subject, 'center':protocol, 'acquisition':acq})
        delFiles += list(rdi._findValues({}, None, False))
        # Delete all files
        for h in delFiles:
            removeFromDB(h.fullPath(), db)
            # Delete associated .minf
            if os.path.isfile(h.fullPath() + '.minf'):
                os.remove(h.fullPath() + '.minf')


    def deleteLockFiles(self, subjPath):
        # Search .lock files in the output folder
        lockFiles = []
        for root, dirs, files in os.walk(subjPath):
            for file in files:
                if file.endswith(".lock"):
                    lockFiles += [os.path.join(root, file)]
        # Offer to delete lock files
        if lockFiles:
            rep = QtGui.QMessageBox.warning(self, u'Confirmation', \
                    u"<b>WARNING:</b><br/>This subject can\'t be processed because it contains .lock files.<br/><br/>" + \
                    u"Option #1: Someone else is processing the same subject right now. In this case you have to wait until the other computation is done.<br/><br/>" + \
                    u"Option #2: The previous execution crashed or was closed before the end and BrainVISA could not delete these .lock files properly, they must be deleted manually.<br/><br/>" + \
                    u"Delete all the lock files now?", QtGui.QMessageBox.Yes | QtGui.QMessageBox.No, QtGui.QMessageBox.No)
            if (rep == QtGui.QMessageBox.Yes):
                for file in lockFiles:
                    os.remove(file)


    def runPipelineBV(self):

        proto = str(self.ui.regProtocolCombo.currentText())
        subj = str(self.ui.regSubjectCombo.currentText())
        images = self.findAllImagesForSubject(proto, subj)
        # Find all images, normalize the T1s, find the T1pre if there is one, read image headers, store their referentials and transformations in the DB
        t1preImage = None
        for image in images:
            patient = image.attributes()['subject']
            acq = image.attributes()['acquisition']
            # Store Scanner-based referential and referential files in the DB
            if acq.startswith('T1'):
                # If there is a T1pre, remember the image
                if acq.find('T1pre') == 0:
                    t1preImage = image
        # No T1pre : nothing else to do
        if t1preImage is None:
            return
        # Check that there is no FreeSurfer segmentation available: otherwise delete it
        rdi = ReadDiskItem('Hemisphere Mesh', 'Anatomist mesh formats', requiredAttributes={'subject':t1preImage['subject'], 'center':t1preImage['center'], "modality":"t1mri"})
        hemis = list(rdi._findValues({}, None, False))
        iFS = [i for i in range(len(hemis)) if 'FreesurferAtlaspre' in hemis[i].attributes()["acquisition"]]
        # Existing segmentation: Ask for confirmation before deleting
        if iFS:
            rep = QtGui.QMessageBox.warning(self, u'Confirmation', u"<font color='red'><b>WARNING:</b> There is an existing FreeSurfer+Morphologist segmentation for this subject. Running this pipeline will delete the existing meshes and MarsAtlas parcels.<br/><br/>Delete existing meshes and MarsAtlas parcels?</font>", QtGui.QMessageBox.Yes | QtGui.QMessageBox.No, QtGui.QMessageBox.No)
            if (rep == QtGui.QMessageBox.Yes):
                self.deleteExistingMeshes(hemis[iFS[0]], 'FreesurferAtlaspre')
            else:
                return

        # Delete .lock files
        subjPath = os.path.normpath(os.path.join(t1preImage.fullName(), os.pardir, os.pardir, os.pardir))
        self.deleteLockFiles(subjPath)

        # Run Morphologist + HipHop
        self.mriAcPc = t1preImage
        self.runMorphologistBV(t1preImage)


    def runPipelineFS(self):

        # Get protocol and subject
        protocol = str(self.ui.regProtocolCombo.currentText())
        subject = str(self.ui.regSubjectCombo.currentText())
        # Find FreeSurfer orig.mgz (FreeSurfer DB)
        rT1 = ReadDiskItem('T1 FreesurferAnat', 'FreesurferMGZ', requiredAttributes={'subject':subject, '_ontology':'freesurfer'})
        allT1 = list(rT1.findValues({},None,False))
        iOrig = [i for i in range(len(allT1)) if 'orig.mgz' in str(allT1[i])]
        if not iOrig:
            QtGui.QMessageBox.warning(self, 'Error', "You must import the output of the FreeSurfer segmentation first.\n(Tab \"Import\")")
            return
        diOrig = allT1[iOrig[0]]

        # Find T1pre of the subject
        rT1BV = ReadDiskItem('Raw T1 MRI', 'BrainVISA volume formats',requiredAttributes={'subject':subject})
        allT1 = list(rT1BV.findValues({},None,False))
        idxT1pre = [i for i in range(len(allT1)) if 'T1pre' in str(allT1[i])]
        diT1pre = allT1[idxT1pre[0]]
        
        # Get output acquisition name
        acq = str(diT1pre.attributes()['acquisition']).replace('T1','FreesurferAtlas')
        # Get output T1 (BrainVISA DB)
        wdi = WriteDiskItem("Raw T1 MRI", "NIFTI-1 image", requiredAttributes={'center':protocol, 'subject':subject})
        diOut = wdi.findValue({'center':protocol, 'subject':subject, 'acquisition':acq, 'modality':'freesurfer_atlas'})
        T1_output = str(diOut.fullPath())
        
        # Check that there is no BrainVISA segmentation available: otherwise delete it
        rdi = ReadDiskItem('Hemisphere Mesh', 'Anatomist mesh formats', requiredAttributes={'subject':subject, 'center':protocol})
        hemis = list(rdi._findValues({}, None, False))
        idxT1pre = [i for i in range(len(hemis)) if 'T1pre' in hemis[i].attributes()["acquisition"]]
        # Existing segmentation: Ask for confirmation before deleting
        if idxT1pre:
            rep = QtGui.QMessageBox.warning(self, u'Confirmation', u"<font color='red'><b>WARNING:</b> There is an existing BrainVISA/Morphologist segmentation for this subject. Running this pipeline will delete the existing meshes and MarsAtlas parcels.<br/><br/>Delete existing meshes and MarsAtlas parcels?</font>", QtGui.QMessageBox.Yes | QtGui.QMessageBox.No, QtGui.QMessageBox.No)
            if (rep == QtGui.QMessageBox.Yes):
                self.deleteExistingMeshes(hemis[idxT1pre[0]], 'T1pre')
            else:
                return

        # Delete .lock files
        subjPath = os.path.normpath(os.path.join(diT1pre.fullName(), os.pardir, os.pardir, os.pardir))
        self.deleteLockFiles(subjPath)

        # Run Morphologist+Hiphop on FreeSurfer mesh (in a different thread)
        ProgressDialog.call(lambda thr:self.runPipelineFSWorker(diOrig, diOut, thr), True, self, "Running Morphologist on FreeSurfer output...", "FreeSufer")
        #self.runPipelineFSWorker(diOrig, diOut)
        # Update list of images
        self.selectRegSubject(subject)
        
        
    def runPipelineFSWorker(self, diOrig, diOut, thread=None):

        # ==== STEP 1: IMPORT =====
        # Run import process
        if thread:
            thread.emit(QtCore.SIGNAL("PROGRESS"), 10)
        proc = getProcessInstance('Import_FROM_FreeSurfer_TO_Morpho')
        proc.T1_orig = str(diOrig.fullPath())
        proc.T1_output = str(diOut.fullPath())
        self.brainvisaContext.runProcess(proc)
        
        # ===== STEP 2: REGISTER =====
        # Add identity transform between "FreeSurfer MRI scanner-based" and "T1pre MRI scanner-based"
        self.insertFreesurferTransformation(diOut)
        
        # ===== STEP 3: MARS ATLAS =====
        # Start hip-hop
        if thread:
            thread.emit(QtCore.SIGNAL("PROGRESS_TEXT"), "Running Hip-Hop...")
            thread.emit(QtCore.SIGNAL("PROGRESS"), 50)
        self.hiphopStart(diOut.attributes()['center'], diOut.attributes()['subject'], diOut.attributes()['acquisition'])
        
        
    def runProcessHiphop(self):
        # Get current protocol/subject
        center = str(self.ui.regProtocolCombo.currentText())
        subject = str(self.ui.regSubjectCombo.currentText())
        # Get possible inputs
        Lrdi = ReadDiskItem('Labelled Cortical folds graph', 'Graph and data', requiredAttributes={'side': 'left', 'subject':subject, 'center':center})
        Lrdi = list( Lrdi._findValues( {}, None, False ) )
        # Display list of available files
        if not Lrdi:
            QtGui.QMessageBox.warning(self, "Labelled Cortical folds graph", "No segmentation available.", QtGui.QMessageBox.Ok, QtGui.QMessageBox.Ok)
            return
        elif (len(Lrdi) > 1):
            iSel = [i for i in range(len(Lrdi)) if 'FreesurferAtlaspre' in str(Lrdi[i])]
            if not iSel:
                iSel = 0
            else:
                iSel = iSel[0]
            strFiles = ""
            for r in Lrdi:
                strFiles += r.fullPath() + "\n"
            QtGui.QMessageBox.information(self, "Labelled Cortical folds graph", "Available files:\n" + strFiles + " \nSelecting:\n" + Lrdi[iSel].fullPath(), QtGui.QMessageBox.Ok, QtGui.QMessageBox.Ok)
        else:
            iSel = 0
        # Running hiphop
        self.hiphopStart(center, subject, Lrdi[iSel].attributes()['acquisition'])


    def setANTstrm_database(self,im, tmp_folder):
        file_to_open = tmp_folder + '0GenericAffine.mat'
        tmp_trm_path = getTmpFilePath('txt')
        # ants_call ='{}/ConvertTransformFile 3 {} {} --RAS --homogeneousMatrix'.format(self.prefs['ants'],file_to_open,tmp_trm_path)
        ants_call ='ConvertTransformFile 3 {} {} --RAS --homogeneousMatrix'.format(file_to_open,tmp_trm_path)
        runCmd(ants_call.split())
        info_mat = numpy.loadtxt(tmp_trm_path)
        info_mat = numpy.linalg.inv(info_mat)
        rot_mat = info_mat[:3,:3]
        tr_mat = info_mat[:3,3]
        result_mat = numpy.vstack((tr_mat,rot_mat))
    
        numpy.savetxt(tmp_trm_path,result_mat,delimiter =' ',fmt='%8.8f')
        self.insertTransformationToT1pre(tmp_trm_path,im)
        if 'ANTs_IntrAnat' in tmp_folder:
            shutil.rmtree(tmp_folder)


    def enableACPCButtons(self, trueorfalse):
        """ Enables or disables all buttons for entering AC/PC reference points """
        self.ui.regAcButton.setEnabled(trueorfalse)
        self.ui.regPcButton.setEnabled(trueorfalse)
        self.ui.regIhButton.setEnabled(trueorfalse)
        self.ui.regLhButton.setEnabled(trueorfalse)
        self.ui.regACPCLabel.setEnabled(trueorfalse)
        self.ui.regAcPcValidateButton.setEnabled(trueorfalse)

    def runMorphologistBV(self, image):
        """ Runs the morphologist analysis to get segmented brain and hemispheres.
        Just activates the buttons to enter AC/PC. The validate button will run the
        analysis itself"""

        # Display images
        self.displayImage(image.fileName(), self.wins)
        # Enable buttons
        self.enableACPCButtons(True)
        # Ask the user to provide AC/PC information to launch MRI segmentation
        QtGui.QMessageBox.information(self, "AC/PC", "Enter AC, PC, IH, LH then validate to segment the brain",   QtGui.QMessageBox.Ok, QtGui.QMessageBox.Ok)
    

    def validateAcPc(self, ):
        if not all(k in self.AcPc for k in ('AC','PC','IH','LH')):
            # Some values are missing
            QtGui.QMessageBox.warning(self, 'Missing points', u"You have entered the following points %s, You have to enter AC, PC, IH and LH"%repr(self.AcPc.keys()))
            return

        # Close the display, disable the buttons and forget AcPc coordinates (they could be mistakenly used for another subject !)
        self.clearAnatomist()
        self.enableACPCButtons(False)
        # Reset buttons color
        self.ui.regAcButton.setStyleSheet("")
        self.ui.regPcButton.setStyleSheet("")
        self.ui.regIhButton.setStyleSheet("")
        self.ui.regLhButton.setStyleSheet("")
        self.setStatus(u"Starting segmentation BrainVisa/Morphologist2015 of the image T1 pre-implantation")
        # Start computation in a separate thread
        errMsg = ProgressDialog.call(lambda thr:self.validateAcPcWorker(thr), True, self, "Processing...", "BrainVISA segmentation")
        # Display error messages
        if errMsg:
            QtGui.QMessageBox.critical(self, u'Segmentation error', u"Errors occured during the segmentation: \n\n" + u"\n".join(errMsg))
        
        
    def validateAcPcWorker(self, thread=None):
        # No gado
        if (not 'Gado' in self.mriAcPc.attributes().keys()) or (not self.mriAcPc.attributes()['Gado']):
            if thread:
                thread.emit(QtCore.SIGNAL("PROGRESS_TEXT"), "Running Morphologist 2015...")
                thread.emit(QtCore.SIGNAL("PROGRESS"), 10)
            self.brainvisaContext.runProcess('Morphologist 2015', t1mri = self.mriAcPc, perform_normalization = False, anterior_commissure = self.AcPc['AC'],\
                    posterior_commissure = self.AcPc['PC'], interhemispheric_point = self.AcPc['IH'], left_hemisphere_point = self.AcPc['LH'], perform_sulci_recognition = True)
        # Gado
        else:
            # Morphologist step #1: Prepare subject
            if thread:
                thread.emit(QtCore.SIGNAL("PROGRESS_TEXT"), "Preparing subject...")
                thread.emit(QtCore.SIGNAL("PROGRESS"), 5)
            processPrepare = getProcessInstance('preparesubject')
            processPrepare.T1mri = self.mriAcPc
            processPrepare.Normalised = "No"
            processPrepare.Anterior_Commissure = self.AcPc['AC']
            processPrepare.Posterior_Commissure = self.AcPc['PC']
            processPrepare.Interhemispheric_Point = self.AcPc['IH']
            processPrepare.Left_Hemisphere_Point = self.AcPc['LH']
            processPrepare.allow_flip_initial_MRI = True
            self.brainvisaContext.runProcess(processPrepare)
            
            # Morphologist step #2: Bias correction
            if thread:
                thread.emit(QtCore.SIGNAL("PROGRESS_TEXT"), "T1 Bias correction...")
                thread.emit(QtCore.SIGNAL("PROGRESS"), 10)
            processBias = getProcessInstance('t1biascorrection')
            processBias.t1mri = self.mriAcPc
            self.brainvisaContext.runProcess(processBias)
            
            # Remove GADO with SPM
            if thread:
                thread.emit(QtCore.SIGNAL("PROGRESS_TEXT"), "SPM image segmentation...")
                thread.emit(QtCore.SIGNAL("PROGRESS"), 15)
            print self.currentSubject + ": Running segmentation to remove Gado on T1 no bias..."
            nobiasRDI = ReadDiskItem("T1 MRI Bias Corrected", 'BrainVISA volume formats',requiredAttributes={"center":self.currentProtocol,"subject":self.currentSubject})
            nobiasimages = list( nobiasRDI._findValues( {}, None, False ) )
            id_pre = [x for x in range(len(nobiasimages)) if 'pre' in str(nobiasimages[x])]
            nobiasPre = str(nobiasimages[id_pre[0]])
            # Run SPM segmentation
            pathTPMseg = os.path.join(str(self.prefs['spm']),'tpm','TPM.nii')
            import copy
            splittedName = nobiasPre.split('/')
            c1Name = copy.deepcopy(splittedName)
            c2Name = copy.deepcopy(splittedName)
            c3Name = copy.deepcopy(splittedName)
            c4Name = copy.deepcopy(splittedName)
            c1Name[-1] = str("c1")+c1Name[-1]
            c1Name = '/'.join(c1Name)
            c2Name[-1] = str("c2")+c2Name[-1]
            c2Name = '/'.join(c2Name)
            c3Name[-1] = str("c3")+c3Name[-1]
            c3Name = '/'.join(c3Name)
            c4Name[-1] = str("c4")+c4Name[-1]
            c4Name = '/'.join(c4Name)
            splittedName[-1]=str("WithoutGado.nii")
            nogadoPre = str('/'.join(splittedName))
            # Run Matlab
            call = matlab_removeGado%("'"+str(self.prefs['spm'])+"'","'"+nobiasPre+",1'","'"+str(pathTPMseg)+",1'","'"+str(pathTPMseg)+",2'","'"+str(pathTPMseg)+",3'","'"+str(pathTPMseg)+",4'","'"+str(pathTPMseg)+",5'","'"+str(pathTPMseg)+",6'",\
                   "'"+c1Name+"'","'"+c2Name+"'","'"+c3Name+"'","'"+c4Name+"'","'"+nobiasPre+"'","'"+nogadoPre+"'")
            errMsg = matlabRun(call)
            if errMsg:
                return errMsg
            print self.currentSubject + ": Segmentation gado done."
            
            # Replace segmented nobias image with segmented image
            if thread:
                thread.emit(QtCore.SIGNAL("PROGRESS_TEXT"), "Saving new segmented image...")
                thread.emit(QtCore.SIGNAL("PROGRESS"), 20)
            print self.currentSubject + ": Replacing nobias.nii with segmented image..."
            nobiasBak = os.path.join(getTmpDir(),self.currentSubject + 'backup.nii')
            cmd1 = ['mv', nobiasPre, nobiasBak]
            cmd2 = ['cp', nogadoPre, nobiasPre]
            line1 = runCmd(cmd1)
            line2 = runCmd(cmd2)
            # Force SPM volume to be saved in S16 (otherwise brainvisa4.6 crashes)
            ret = subprocess.call(['AimsFileConvert', '-i', nobiasPre, '-o', nobiasPre, '-t', 'S16'])
            if ret < 0:
                #QtGui.QMessageBox.warning(self, "Error", u"Error in the conversion of the segmented image to int16.")
                print "Error: Error in the conversion of the segmented image to int16." 
                return
            
            # Execute the rest of the Morphologist pipeline
            if thread:
                thread.emit(QtCore.SIGNAL("PROGRESS_TEXT"), "Calling Morpologist 2015... (GADO)")
                thread.emit(QtCore.SIGNAL("PROGRESS"), 30)
            morphologist = getProcessInstance('morphologist')
            morphologist.executionNode().PrepareSubject.setSelected(False)
            morphologist.executionNode().BiasCorrection.setSelected(False)
            morphologist.executionNode().HistoAnalysis.setSelected(True)
            morphologist.executionNode().BrainSegmentation.setSelected(True)
            morphologist.executionNode().Renorm.setSelected(True)
            morphologist.executionNode().SplitBrain.setSelected(True)
            morphologist.executionNode().TalairachTransformation.setSelected(False)
            morphologist.executionNode().HeadMesh.setSelected(True)
            morphologist.executionNode().HemispheresProcessing.setSelected(True)
            morphologist.executionNode().SulcalMorphometry.setSelected(True)
            # Synchronous computation
            self.brainvisaContext.runProcess(morphologist, t1mri = self.mriAcPc, perform_normalization = False, anterior_commissure = self.AcPc['AC'],\
                    posterior_commissure = self.AcPc['PC'], interhemispheric_point = self.AcPc['IH'], left_hemisphere_point = self.AcPc['LH'], perform_sulci_recognition = True)

            # Task finishesd
            self.taskfinished(self.currentSubject + u': BrainVISA segmentation and meshes generation')
            # Restore initial nobias image
            print self.currentSubject + ": Restoring original nobias.nii..."
            cmd = ['mv', nobiasBak, nobiasPre]
            line1 = runCmd(cmd)
            
        # Start hip-hop
        if thread:
            thread.emit(QtCore.SIGNAL("PROGRESS_TEXT"), "Running Hip-Hop...")
            thread.emit(QtCore.SIGNAL("PROGRESS"), 50)
        self.hiphopStart(self.mriAcPc.attributes()['center'], self.mriAcPc.attributes()['subject'], self.mriAcPc.attributes()['acquisition'])
        

    def hiphopStart(self, center, subject, acq):
        self.taskfinished(self.currentSubject + u': Morphologist segmentation and meshes generation')
        self.setStatus(self.currentSubject + u": Starting Hip-Hop")

        Lrdi = ReadDiskItem('Labelled Cortical folds graph', 'Graph and data', requiredAttributes={ 'side': 'left', 'subject':subject, 'center':center, 'acquisition':acq})
        Rrdi = ReadDiskItem('Labelled Cortical folds graph', 'Graph and data', requiredAttributes={ 'side': 'right', 'subject':subject, 'center':center, 'acquisition':acq})
        Lrdi = list( Lrdi._findValues( {}, None, False ) )
        if len(Lrdi) == 0:
            print('no left sulci label found, CANNOT RUN HIP HOP')
            return
    
        Rrdi = list( Rrdi._findValues( {}, None, False ) )
        if len(Rrdi) == 0:
            print('no right sulci label found, CANNOT RUN HIP HOP')
            return
    
        self.brainvisaContext.runProcess('Hip-Hop Cortical Parameterization', Lgraph = Lrdi[0], Rgraph = Rrdi[0], sulcus_identification ='label')
        self.taskfinished(u'Hip-Hop done')

    def spmNormalize(self, image, protocol, patient, acq):
        """ Normalize one image (filepath of nifti file) to the T1 template MNI referential""" 
        # Delete existing normalization
        rdi_y = ReadDiskItem('SPM normalization deformation field','NIFTI-1 image',requiredAttributes={'center':str(protocol), 'subject':str(patient),'acquisition':str(acq) })
        rdi_y = list(rdi_y.findValues({}, None, False))
        if rdi_y:
            print("Deleting existing SPM normalization...")
            self.deleteExistingNormalization(rdi_y[0], acq)
        # SPM call
        call = spm12_normalise%("'"+str(self.prefs['spm'])+"'","{'"+str(image)+",1'}", "{'"+str(image)+",1'}", "{'"+str(self.ui.prefSpmTemplateEdit.text())+os.sep+'tpm/TPM.nii'+",1'}")
        # Call SPM normalization
        errMsg = matlabRun(call)
        if errMsg:
            return errMsg
        # Register new files
        self.insertSPMdeformationFile(protocol, patient, acq)
        self.StatisticDataMNItoScannerBased(protocol, patient, acq)


    def getScannerBasedRef(self, imageDiskItem):
        rdi = ReadDiskItem('Scanner Based Referential', 'Referential')
        return rdi.findValue(imageDiskItem)

    def getT1preScannerBasedRef(self, protocol, patient):
        rdi = ReadDiskItem('Scanner Based Referential', 'Referential', requiredAttributes={'modality':'t1mri', 'subject':patient, 'center':protocol})
        allT1s = list (rdi._findValues( {}, None, False ) )
        for t1 in allT1s:
            if t1.attributes()['acquisition'].startswith(u'T1pre'):
                return t1
        return None

    def getT1preNativeRef(self, protocol, patient):
        rdi = ReadDiskItem('Referential of Raw T1 MRI', 'Referential', requiredAttributes={'subject':patient, 'center':protocol})
        allT1s = list (rdi._findValues( {}, None, False ) )
        for t1 in allT1s:
            if t1.attributes()['acquisition'].startswith(u'T1pre'):
                return t1
        return None

    def setResampledToT1pre(self, image, registeredPath):
        ret = subprocess.call(['AimsFileConvert', '-i', str(registeredPath), '-o', str(image.fileName()), '-t', 'S16'])
        if ret < 0:
            print "ERROR: Cannot import the image resampled by SPM : "+repr(registeredPath)
            #QtGui.QMessageBox.warning(self, "Error", u"the image has not been resampled by SPM !")
            return
        # TODO # Should destroy existing referentials for this image
        # The referential of this image is the same as the T1 pre
        self.transfoManager.setReferentialTo(image, self.getT1preNativeRef(image.attributes()['center'], image.attributes()['subject'])) #replace self.getT1preNativeRef(image.attributes()['center'], image.attributes()['subject']) by talairachMNIReferentialId
        os.remove(registeredPath)
        

    def insertFreesurferTransformation(self, diOut):
        # Create identity transform
        trmpath = getTmpFilePath('trm')
        id_trm = numpy.matrix('0 0 0; 1 0 0; 0 1 0; 0 0 1')
        numpy.savetxt(trmpath, id_trm, delimiter =' ', fmt='%d')
        # Re-generate scanner-based transformations
        self.storeImageReferentialsAndTransforms(diOut)
        # Add identity transform from "Freesurfer T1 scanner-based" to "BrainVISA T1 scanner-based"
        transformT1 = self.insertTransformationToT1pre(trmpath, diOut)
        if not transformT1:
            return
        

    def insertTransformationToT1pre(self, trmpath, image):
        """Inserts a TRM file in the DB as the transformation between image and the T1pre scanner-based ref from the same subject"""
        # Check that transformation file exists
        if not os.path.exists(trmpath):
            print "Error : file %s does not exist : cannot insert it as a TRM for image %s"%(trmpath, image.fileName())
        # Get the source (image) and destination (T1pre) referentials
        t1preRef = self.getT1preScannerBasedRef(image.attributes()['center'], image.attributes()['subject'])
        imageRef = self.getScannerBasedRef(image)
        if t1preRef is None or imageRef is None:
            print "Error: Cannot find referentials to declare coregister transform to T1pre."
            return
        # Get new transformation file
        wdiTransformT1 = WriteDiskItem( 'Transform '+ self.modas[image.attributes()['modality']] +' to another image', 'Transformation matrix', exactType=True, requiredAttributes = {'modalityTarget':t1preRef.attributes()['modality'], 'acquisitionTarget':t1preRef.attributes()['acquisition']} )
        transformT1 = wdiTransformT1.findValue(image)
        if transformT1 is None:
            print "Error: Cannot find path for T1pre transform in database for %s-> transformation NOT stored !"%(image.fileName())
            return
        # Remove the existing .trm file first
        if os.path.exists(transformT1.fullPath()):
            os.remove(transformT1.fullPath())
        # Move the trm file to the database (not with shutil because it crashes when using mounted CIFS shares)
        runCmd(['mv', trmpath, transformT1.fullPath()])

        # Create new folder
        try:
            neuroHierarchy.databases.createDiskItemFromFileName(os.path.dirname( transformT1.fullPath() ))
        except:
            pass
        # Set and update database
        self.transfoManager.setNewTransformationInfo( transformT1, source_referential=imageRef, destination_referential=t1preRef)
        print "Added transformation: from %s to %s for %s"%(repr(imageRef), repr(t1preRef), repr(transformT1))
        return transformT1

    def insertSPMdeformationFile(self, protocol, patient, acq):
        """ Should be called when a _sn file (SPM normalization transformation file) was created to update it in the database """
        wdi = WriteDiskItem('SPM normalization deformation field', 'NIFTI-1 image' )
        di = wdi.findValue( { 'center': protocol, 'subject' : patient, 'acquisition':acq } )
        if di is None:
            #QtGui.QMessageBox.warning(self, "Error", "Impossible to find a valid path to import SPM normalization into BrainVisa")
            print("Error: Impossible to find a valid path to import SPM normalization into BrainVisa")
            return
    
        wdi_write = WriteDiskItem('T1 SPM resampled in MNI','NIFTI-1 image')
        di_write = wdi_write.findValue({'center': protocol, 'subject' : patient, 'acquisition':acq})
        if not os.path.isfile(di_write.fileName()):
            #QtGui.QMessageBox.warning(self, "Error", "Impossible dto find a valid path for the T1 coregistered into the MNI")
            print("Error: Impossible to find a valid path for the T1 coregistered into the MNI")
        else:
            print "Declaring T1 registered MNI in BrainVisa DB : " + di_write.fileName()
            neuroHierarchy.databases.insertDiskItem(di_write, update = True)
            self.setNormalizedToT1pre(di_write,di_write.fileName())
    
        # The file should already be there : if it is not, abort with an error, otherwise declare it in the DB
        if os.path.isfile(di.fileName()):
            print "Declaring SPM normalization in BrainVisa DB : " + di.fileName()
            neuroHierarchy.databases.insertDiskItem( di, update=True )
            # Compute deformation fields
            #wdi = WriteDiskItem('SPM normalization inverse deformation field','aims readable volume formats')
            #di2 = wdi.findValue(di)
            #neuroHierarchy.databases.insertDiskItem( di2, update=True )
    
    
        else:
            print "No SPM normalization file found ! The normalization probably failed..."
            #QtGui.QMessageBox.warning(self, "Error", u"SPM normalization file %s unfound ! \n Normalization probably failed !"%di.fileName())
            return

    def StatisticDataMNItoScannerBased(self, protocol, patient, acq):

        #is there any statistic data to convert to the Scanner-Based space
        rdi = ReadDiskItem('Statistic-Data', 'BrainVISA volume formats', requiredAttributes={'center':str(protocol), 'subject':str(patient) })
        di=list(rdi.findValues({}, None, False))
        StatToConvert = [i for i in range(len(di)) if eval(di[i].attributes()['MNI'])]

        if len(StatToConvert)>0:

            diT1 = ReadDiskItem( 'Raw T1 MRI', 'BrainVISA volume formats', requiredAttributes={'center':protocol, 'subject':patient, 'normalized':'no' } )
            allT1 = list(diT1.findValues({},None,False))
            idxT1pre = [i for i in range(len(allT1)) if 'T1pre' in str(allT1[i])]
            self.mriAcPc = allT1[idxT1pre[0]]

            for i_stat2conv in range(len(StatToConvert)):

                self.setStatus(u"Start MNI to ScannerBased conversion of %s "%str(di[StatToConvert[i_stat2conv]].fileName()))

                #check if inverse deformation field exist.
                rdi_defField_read = ReadDiskItem('SPM normalization inverse deformation field','NIFTI-1 image',requiredAttributes={'center':str(protocol), 'subject':str(patient),'acquisition':self.mriAcPc.attributes()['acquisition']})
                di_defField_read=list(rdi_defField_read.findValues({}, None, False))


                tempNameMNI2SB = str(di[StatToConvert[i_stat2conv]].fileName()).split('/')
                tempNameMNI2SB[-1]= "tmpMNItoScannerBased"+tempNameMNI2SB[-1]
                tempNameMNI2SB = "/".join(tempNameMNI2SB)

                if len(di_defField_read) > 0:
                    print "inverse deformation field found and used"
                    matlabRun(spm_MNItoScannerBased%("'"+str(self.prefs['spm'])+"'","'"+str(di_defField_read[0].fileName())+"'","'"+str(di[StatToConvert[i_stat2conv]].fileName())+",1'","'"+str(self.mriAcPc)+",1'","'"+tempNameMNI2SB+",1'"))
    
                    cmd1 = ['mv', str(di[StatToConvert[i_stat2conv]].fileName()), di[StatToConvert[i_stat2conv]].fileName()[:-4]+"_MNI.nii.backup"]
                    line1 = runCmd(cmd1)
    
                    rname = tempNameMNI2SB.split('/')
                    rname[-1]="r"+rname[-1]
                    rname = "/".join(rname)
                    cmd2 = ['cp',rname,  str(di[StatToConvert[i_stat2conv]].fileName())]
                    line2 = runCmd(cmd2)
    
                    os.remove(tempNameMNI2SB)
                    os.remove(rname)
                    os.remove(str(di[StatToConvert[i_stat2conv]].fileName())+".minf")
    
                    write_filters = { 'center': protocol, 'acquisition': str(di[StatToConvert[i_stat2conv]].attributes()['acquisition']), 'subject' : str(patient) }
                    wdi_new = WriteDiskItem('Statistic-Data', 'NIFTI-1 image' )#'gz compressed NIFTI-1 image' )
                    write_filters.update({'subacquisition':di[StatToConvert[0]].attributes()['subacquisition']})
    
                    di_new = wdi_new.findValue(write_filters)
    
                    ret = subprocess.call(['AimsFileConvert', '-i', str(di[StatToConvert[i_stat2conv]].fileName()), '-o', str(di[StatToConvert[i_stat2conv]].fileName())])
                    di[StatToConvert[i_stat2conv]].setMinf('MNI','False')
                    di[StatToConvert[i_stat2conv]].setMinf('ColorPalette','Yellow-Red-White-Blue-Green')
                    neuroHierarchy.databases.insertDiskItem( di_new, update=True )
                    self.transfoManager.setReferentialTo(di_new, self.mriAcPc)

                else:
                    #look for a y_file second
                    rdi_y = ReadDiskItem('SPM normalization deformation field','NIFTI-1 image',requiredAttributes={'center':str(protocol), 'subject':str(patient),'acquisition':self.mriAcPc.attributes()['acquisition'] })
                    di_y = list(rdi_y.findValues({}, None, False))
    
                    if len(di_y) == 0:
                        print "No deformation field found in database"
                        self.setStatus(u"MNI to ScannerBased conversion of %s not performed, no deformation field found"%str(di[StatToConvert[i_stat2conv]].fileName()))
                        return
                    else:
                        print "deformation field found and used"
                        wdi_inverse = WriteDiskItem('SPM normalization inverse deformation field','NIFTI-1 image')
                        dir_yinv_split = str(di_y[0].fileName()).split('/')
                        name_yinverse = dir_yinv_split.pop()[2:]
                        #name_yinverse.replace('.nii','_inverse.nii')
                        dir_yinverse = "/".join(dir_yinv_split)
                        di_inverse = wdi_inverse.findValue(di_y[0])
                        #on fait l'inversion de la deformation
                        #pour le moment ce bout de code ne marche qu'avec spm12
                        matlabRun(spm_inverse_y_field12%("'"+str(self.prefs['spm'])+"'","'"+str(di_y[0].fileName())+"'","'"+str(self.mriAcPc)+"'","'"+name_yinverse.replace('.nii','_inverse.nii')+"'","'"+dir_yinverse+"'"))
                        neuroHierarchy.databases.insertDiskItem( di_inverse, update=True )
                        matlabRun(spm_MNItoScannerBased%("'"+str(self.prefs['spm'])+"'","'"+str(di_inverse)+"'","'"+str(di[StatToConvert[i_stat2conv]].fileName())+",1'","'"+str(self.mriAcPc)+",1'","'"+tempNameMNI2SB+",1'"))
        
                        cmd1 = ['mv', str(di[StatToConvert[i_stat2conv]].fileName()), di[StatToConvert[i_stat2conv]].fileName()[:-4]+"_MNI.nii"]
                        line1 = runCmd(cmd1)
        
                        rname = tempNameMNI2SB.split('/')
                        rname[-1]="r"+rname[-1]
                        rname = "/".join(rname)
                        cmd2 = ['cp',rname,  str(di[StatToConvert[i_stat2conv]].fileName())]
                        line2 = runCmd(cmd2)
        
                        os.remove(tempNameMNI2SB)
                        os.remove(rname)
                        os.remove(str(di[StatToConvert[i_stat2conv]].fileName())+".minf")
        
                        write_filters = { 'center': protocol, 'acquisition': str(di[StatToConvert[i_stat2conv]].attributes()['acquisition']), 'subject' : str(patient) }
                        wdi_new = WriteDiskItem('Statistic-Data', 'NIFTI-1 image' )#'gz compressed NIFTI-1 image' )
                        write_filters.update({'subacquisition':di[StatToConvert[i_stat2conv]].attributes()['subacquisition']})
        
                        di_new = wdi_new.findValue(write_filters)
        
                        ret = subprocess.call(['AimsFileConvert', '-i', str(di[StatToConvert[i_stat2conv]].fileName()), '-o', str(di[StatToConvert[i_stat2conv]].fileName())])
                        di[StatToConvert[i_stat2conv]].setMinf('MNI','False')
                        di[StatToConvert[i_stat2conv]].setMinf('ColorPalette','Yellow-Red-White-Blue-Green')
                        neuroHierarchy.databases.insertDiskItem( di[StatToConvert[i_stat2conv]], update=True )
                        self.transfoManager.setReferentialTo(di[StatToConvert[i_stat2conv]], self.mriAcPc)
    
                self.setStatus(u"MNI to ScannerBased conversion of %s done"%str(di[StatToConvert[i_stat2conv]].fileName()))

        else:
            print "no statistic data to convert from MNI to Scanner-Based"



    def generateAmygdalaHippoMesh(self, protocol, patient, acq, diFS):

        print("start generation of amygdala and hippocampus meshes")
        #attention ça merde quand l'irm est trop haute résolution

        volDestrieux = aims.read(diFS.fullPath())
        npDestrieux = volDestrieux.arraydata()

        lefthippopx = numpy.where(npDestrieux==17)
        notlefthippopx = numpy.where(npDestrieux!=17)
        lefthippoposition = [[lefthippopx[1][i]*volDestrieux.getVoxelSize()[2],lefthippopx[2][i]*volDestrieux.getVoxelSize()[1],lefthippopx[3][i]*volDestrieux.getVoxelSize()[0]] for i in range(len(lefthippopx[1]))]
        righthippopx = numpy.where(npDestrieux==53)
        notrighthippopx = numpy.where(npDestrieux!=53)
        righthippoposition = [[righthippopx[1][i]*volDestrieux.getVoxelSize()[2],righthippopx[2][i]*volDestrieux.getVoxelSize()[1],righthippopx[3][i]*volDestrieux.getVoxelSize()[0]] for i in range(len(righthippopx[1]))]

        leftamygdalapx = numpy.where(npDestrieux == 18)
        notleftamygdalapx = numpy.where(npDestrieux != 18)
        rightamygdalapx = numpy.where(npDestrieux == 54)
        notrightamygdalapx = numpy.where(npDestrieux != 54)

        #add here a verification of number of pixel per parcels. if not enough, no point to generate a mesh
        AmygdalaLeft = True
        AmygdalaRight = True
        HippoLeft = True
        HippoRight = True

        if len(leftamygdalapx[0]) < 50:
            AmygdalaLeft = False
        if len(rightamygdalapx[0]) < 50:
            AmygdalaRight = False

        if len(lefthippopx[0]) < 50:
            HippoLeft = False
        if len(righthippopx[0]) < 50:
            HippoRight = False

        di = ReadDiskItem( 'Raw T1 MRI', 'BrainVISA volume formats', requiredAttributes={'center':protocol, 'subject':patient, 'normalized':'no' } )
        allT1 = list(di.findValues({},None,False))
        idxT1pre = [i for i in range(len(allT1)) if 'T1pre' in str(allT1[i])]
        T1pre = allT1[idxT1pre[0]]
        self.storeImageReferentialsAndTransforms(T1pre)
        constraints =  { 'center':protocol, 'subject':patient, 'acquisition':T1pre.attributes()['acquisition'] }

        #il faudrait mettre le referentiel à freesurferatlas la aussi ça serait fait comme ça.
        self.transfoManager.setReferentialTo(diFS, T1pre.attributes()['referential'] )

        if AmygdalaRight:
            # Save amigdala mask in a temporary .nii file (set all the non-amygdala values to zero)
            xyzAmyg = ([0] * notrightamygdalapx[3].size, notrightamygdalapx[1].tolist(), notrightamygdalapx[2].tolist(), notrightamygdalapx[3].tolist())
            iAmyg = numpy.ravel_multi_index(xyzAmyg, volDestrieux.arraydata().shape)
            volDestrieux.arraydata().put(iAmyg, 0)
            aims.write(volDestrieux, os.path.join(getTmpDir(),'rightamygdala.nii'))
            # Read FreeSurfer again atlas because it was modified
            volDestrieux = aims.read(diFS.fullPath())
            
            # Create right amygdala mesh
            wdirightamygdala = WriteDiskItem('rightAmygdala', 'GIFTI file' )
            dirightamygdala = wdirightamygdala.findValue(constraints)
            if dirightamygdala is None:
                print "mesh export : could not find valid path"
            else:
                if not os.path.isdir(os.path.dirname(dirightamygdala.fullPath())):
                    os.makedirs(os.path.dirname(dirightamygdala.fullPath()))
                ret = subprocess.call(['AimsMeshBrain', '-i', os.path.join(getTmpDir(),'rightamygdala.nii'), '-o', dirightamygdala.fullPath()])
                neuroHierarchy.databases.insertDiskItem(dirightamygdala, update = True)
                if 'referential' in T1pre.attributes().keys():
                    self.transfoManager.setReferentialTo(dirightamygdala, T1pre.attributes()['referential'] )

        if AmygdalaLeft:
            # Save amigdala mask in a temporary .nii file (set all the non-amygdala values to zero)
            xyzAmyg = ([0] * notleftamygdalapx[3].size, notleftamygdalapx[1].tolist(), notleftamygdalapx[2].tolist(), notleftamygdalapx[3].tolist())
            iAmyg = numpy.ravel_multi_index(xyzAmyg, volDestrieux.arraydata().shape)
            volDestrieux.arraydata().put(iAmyg, 0)
            aims.write(volDestrieux, os.path.join(getTmpDir(),'leftamygdala.nii'))
            # Read FreeSurfer again atlas because it was modified
            volDestrieux = aims.read(diFS.fullPath())

            # Create left amygdala mesh
            wdileftamygdala = WriteDiskItem('leftAmygdala', 'GIFTI file' )
            dileftamygdala = wdileftamygdala.findValue(constraints)
            if dileftamygdala is None:
                print "mesh export : could not find valid path"
            else:
                if not os.path.isdir(os.path.dirname(dileftamygdala.fullPath())):
                    os.makedirs(os.path.dirname(dileftamygdala.fullPath()))
                ret = subprocess.call(['AimsMeshBrain', '-i', os.path.join(getTmpDir(),'leftamygdala.nii'), '-o', dileftamygdala.fullPath()])
                neuroHierarchy.databases.insertDiskItem(dileftamygdala, update = True)
                if 'referential' in T1pre.attributes().keys():
                    self.transfoManager.setReferentialTo(dileftamygdala, T1pre.attributes()['referential'] )

        #pour pouvoir decouper les hippocampes en deux
        if isinstance(T1pre.attributes()['SB_Transform'],basestring):
            import ast
            m = aims.Motion(ast.literal_eval(T1pre.attributes()['SB_Transform']))
        else:
            m = aims.Motion(T1pre.attributes()['SB_Transform'])

        if HippoRight:
            # Save hippocampush mask in a temporary .nii file (set all the non-hippocampus values to zero)
            xyzHippo = ([0] * notrighthippopx[3].size, notrighthippopx[1].tolist(), notrighthippopx[2].tolist(), notrighthippopx[3].tolist())
            iHippo = numpy.ravel_multi_index(xyzHippo, volDestrieux.arraydata().shape)
            volDestrieux.arraydata().put(iHippo, 0)
            aims.write(volDestrieux, os.path.join(getTmpDir(),'righthippo.nii'))
            # Read FreeSurfer again atlas because it was modified
            volDestrieux = aims.read(diFS.fullPath())
                
            # Create right hippocampus mesh
            wdirighthippo =  WriteDiskItem('rightHippo', 'GIFTI file' )
            dirightHippo = wdirighthippo.findValue(constraints)
            if dirightHippo is None:
                print "mesh export : could not find valid path"
            else:
                if not os.path.isdir(os.path.dirname(dirightHippo.fullPath())):
                    os.makedirs(os.path.dirname(dirightHippo.fullPath()))
            ret = subprocess.call(['AimsMeshBrain', '-i', os.path.join(getTmpDir(),'righthippo.nii'), '-o', dirightHippo.fullPath()])
            neuroHierarchy.databases.insertDiskItem(dirightHippo, update = True)
            if 'referential' in T1pre.attributes().keys():
                self.transfoManager.setReferentialTo(dirightHippo, T1pre.attributes()['referential'] )

            wdirighthippoantero =  WriteDiskItem('rightanteroHippocampus', 'GIFTI file' )
            wdirighthippopostero =  WriteDiskItem('rightposteroHippocampus', 'GIFTI file' )
            dirighthippoantero = wdirighthippoantero.findValue(constraints)
            dirighthippopostero = wdirighthippopostero.findValue(constraints)

            UU, ss, rotation, center = self.fitvolumebyellipse(righthippoposition)

            #comment deviner automatiquement quel bout est antero et quel bout est postero...

            mesh = aims.read(str(dirightHippo.fullPath()))
            testhippocut = aims.AimsSurfaceTriangle()
            testhippocut2 = aims.AimsSurfaceTriangle()
            border = aims.AimsTimeSurface(2)
            border2 = aims.AimsTimeSurface(2)
            aims.SurfaceManip.cutMesh(mesh, [rotation[2,2], rotation[2,1], rotation[2,0], -numpy.inner(rotation[2,:],center)], testhippocut, border)
            aims.SurfaceManip.cutMesh(mesh, [-rotation[2,2], -rotation[2,1], -rotation[2,0], numpy.inner(rotation[2,:],center)], testhippocut2, border2)
            planmesh = aims.SurfaceManip.meshPlanarPolygon([rotation[2,2], rotation[2,1], rotation[2,0], -numpy.inner(rotation[2,:],center)], border)
            planmesh2 = aims.SurfaceManip.meshPlanarPolygon([-rotation[2,2], -rotation[2,1], -rotation[2,0], numpy.inner(rotation[2,:],center)], border2)

            aims.SurfaceManip.meshMerge(testhippocut,planmesh)
            aims.SurfaceManip.meshMerge(testhippocut2,planmesh2)
            aims.write(testhippocut,os.path.join(getTmpDir(),'testhippocut.gii'))
            aims.write(testhippocut2,os.path.join(getTmpDir(),'testhippocut2.gii'))
            
            #equivalent à : test = aims.read(-numpy.inner(rotation[2,:],center))  truc = aims.AimsTimeSurface_2()  bidule=aims.AimsTimeSurface()  aims.SurfaceManip.cutMesh(test1,[rotation[2,2],rotation[2,1],rotation[2,0],-numpy.inner(rotation[2,:],center)],bidule,truc) aims.write(bidule,os.path.join(getTmpDir(),'bidule.gii'))
            hippogii = aims.read(os.path.join(getTmpDir(),'testhippocut.gii'))
            hippovertex = numpy.array(hippogii.vertex().list())
            centerhippo = numpy.average(hippovertex,axis=0)

            hippogii2 = aims.read(os.path.join(getTmpDir(),'testhippocut2.gii'))
            hippovertex2 = numpy.array(hippogii2.vertex().list())
            centerhippo2 = numpy.average(hippovertex2,axis=0)

            coords1 = m.transform(centerhippo)
            coords2 = m.transform(centerhippo2)

            if coords2[1] > coords1[1]:
                #2 est l'antero; 1 est posterieur
                if dirighthippoantero is None or dirighthippopostero is None:
                    print "mesh export : could not find valid path"
                else:
                    cmd1 = ['mv', os.path.join(getTmpDir(),'testhippocut2.gii'), str(dirighthippoantero.fullPath())]
                    cmd2 = ['mv', os.path.join(getTmpDir(),'testhippocut.gii'), str(dirighthippopostero.fullPath())]
                    line1 = runCmd(cmd1)
                    line2 = runCmd(cmd2)
                    neuroHierarchy.databases.insertDiskItem(dirighthippoantero, update = True)
                    neuroHierarchy.databases.insertDiskItem(dirighthippopostero, update = True)
                    if 'referential' in T1pre.attributes().keys():
                        self.transfoManager.setReferentialTo(dirighthippoantero, T1pre.attributes()['referential'] )
                        self.transfoManager.setReferentialTo(dirighthippopostero, T1pre.attributes()['referential'] )
            elif coords2[1] < coords1[1]:
                #1 est l'antero; 2 est le posterieur
                if dirighthippoantero is None or dirighthippopostero is None:
                    print "mesh export : could not find valid path"
                else:
                    cmd1 = ['mv', os.path.join(getTmpDir(),'testhippocut.gii'), str(dirighthippoantero.fullPath())]
                    cmd2 = ['mv', os.path.join(getTmpDir(),'testhippocut2.gii'), str(dirighthippopostero.fullPath())]
                    line1 = runCmd(cmd1)
                    line2 = runCmd(cmd2)
                    neuroHierarchy.databases.insertDiskItem(dirighthippoantero, update = True)
                    neuroHierarchy.databases.insertDiskItem(dirighthippopostero, update = True)
                    if 'referential' in T1pre.attributes().keys():
                        self.transfoManager.setReferentialTo(dirighthippoantero, T1pre.attributes()['referential'] )
                        self.transfoManager.setReferentialTo(dirighthippopostero, T1pre.attributes()['referential'] )

            #  aims.TimeTexture() puis AimsMeshParcellation2VolumeParcellation
            if dirighthippoantero is not None and dirighthippopostero is not None:
                #we read it, estimate the number of vertex, and attribute the good value to the good number of vertex
                #nb vertex of right antero
                voxel_size_T1 = [T1pre.attributes()['voxel_size'][0], T1pre.attributes()['voxel_size'][1], T1pre.attributes()['voxel_size'][2], 1.0]

                hippo_vol_antero=aims.Volume(*T1pre.attributes()['volume_dimension'][:3],dtype='S16')
                hippo_vol_postero=aims.Volume(*T1pre.attributes()['volume_dimension'][:3],dtype='S16')
                hippo_vol_full=aims.Volume(*T1pre.attributes()['volume_dimension'][:3],dtype='S16')
                hippo_vol_antero.header()['voxel_size']= voxel_size_T1
                hippo_vol_postero.header()['voxel_size']= voxel_size_T1
                hippo_vol_full.header()['voxel_size']= voxel_size_T1

                wdirighthippoNII =  WriteDiskItem('rightHippocampusNII', 'NIFTI-1 image' )
                dirightHippoNII = wdirighthippoNII.findValue(constraints)

                meshantero = aims.read(dirighthippoantero.fullPath())
                meshpostero = aims.read(dirighthippopostero.fullPath())
                aims.SurfaceManip.rasterizeMesh(meshantero,hippo_vol_antero,1)
                aims.SurfaceManip.rasterizeMesh(meshpostero,hippo_vol_postero,1)
                #fill the insides voxel
                for zslices in range(hippo_vol_antero.arraydata().shape[1]):
                    hippo_vol_antero.arraydata()[0,zslices,:,:] = ndimage.morphology.binary_fill_holes(hippo_vol_antero.arraydata()[0,zslices,:,:]).astype(int)

                for zslices in range(hippo_vol_postero.arraydata().shape[1]):
                    hippo_vol_postero.arraydata()[0,zslices,:,:] = ndimage.morphology.binary_fill_holes(hippo_vol_postero.arraydata()[0,zslices,:,:]).astype(int)

                hippo_vol_antero *= 5301
                hippo_vol_postero *= 5302
                hippo_vol_full = hippo_vol_antero + hippo_vol_postero
                for z in xrange(hippo_vol_full.getSizeZ()):
                    for y in xrange(hippo_vol_full.getSizeY()):
                        for x in xrange(hippo_vol_full.getSizeX()):
                            if hippo_vol_full.value(x, y, z) >5302:
                                hippo_vol_full.setValue(5302, x, y, z)

                aims.write(hippo_vol_full,str(dirightHippoNII.fullPath()))
                neuroHierarchy.databases.insertDiskItem(dirightHippoNII, update = True)
                if 'referential' in T1pre.attributes().keys():
                    self.transfoManager.setReferentialTo(dirightHippoNII, T1pre.attributes()['referential'] )


        if HippoLeft:
            # Save hippocampush mask in a temporary .nii file (set all the non-hippocampus values to zero)
            xyzHippo = ([0] * notlefthippopx[3].size, notlefthippopx[1].tolist(), notlefthippopx[2].tolist(), notlefthippopx[3].tolist())
            iHippo = numpy.ravel_multi_index(xyzHippo, volDestrieux.arraydata().shape)
            volDestrieux.arraydata().put(iHippo, 0)
            aims.write(volDestrieux, os.path.join(getTmpDir(),'lefthippo.nii'))
            # No need to read again the modified FreeSurfer atlas: it's not used anymore
            
            # Create right hippocampus mesh
            wdilefthippo =  WriteDiskItem('leftHippo', 'GIFTI file' )
            dileftHippo = wdilefthippo.findValue(constraints)
            if dileftHippo is None:
                print "mesh export : could not find valid path"
            else:
                if not os.path.isdir(os.path.dirname(dileftHippo.fullPath())):
                    os.makedirs(os.path.dirname(dileftHippo.fullPath()))
                ret = subprocess.call(['AimsMeshBrain', '-i', os.path.join(getTmpDir(),'lefthippo.nii'), '-o', dileftHippo.fullPath()])
                neuroHierarchy.databases.insertDiskItem(dileftHippo, update = True)
                if 'referential' in T1pre.attributes().keys():
                    self.transfoManager.setReferentialTo(dileftHippo, T1pre.attributes()['referential'] )

                wdilefthippoantero =  WriteDiskItem('leftanteroHippocampus', 'GIFTI file' )
                wdilefthippopostero =  WriteDiskItem('leftposteroHippocampus', 'GIFTI file' )
                dilefthippoantero = wdilefthippoantero.findValue(constraints)
                dilefthippopostero = wdilefthippopostero.findValue(constraints)

                UU, ss, rotation, center = self.fitvolumebyellipse(lefthippoposition)

                #comment deviner automatiquement quel bout est antero et quel bout est postero...
                mesh = aims.read(str(dileftHippo.fullPath()))
                testhippocut = aims.AimsSurfaceTriangle()
                testhippocut2 = aims.AimsSurfaceTriangle()
                border = aims.AimsTimeSurface(2)
                border2 = aims.AimsTimeSurface(2)
                aims.SurfaceManip.cutMesh(mesh, [rotation[2,2], rotation[2,1], rotation[2,0], -numpy.inner(rotation[2,:],center)], testhippocut, border)
                aims.SurfaceManip.cutMesh(mesh, [-rotation[2,2], -rotation[2,1], -rotation[2,0], numpy.inner(rotation[2,:],center)], testhippocut2, border2)
                planmesh = aims.SurfaceManip.meshPlanarPolygon([rotation[2,2], rotation[2,1], rotation[2,0], -numpy.inner(rotation[2,:],center)], border)
                planmesh2 = aims.SurfaceManip.meshPlanarPolygon([-rotation[2,2], -rotation[2,1], -rotation[2,0], numpy.inner(rotation[2,:],center)], border2)

                aims.SurfaceManip.meshMerge(testhippocut,planmesh)
                aims.SurfaceManip.meshMerge(testhippocut2,planmesh2)
                aims.write(testhippocut, os.path.join(getTmpDir(),'testhippocut.gii'))
                aims.write(testhippocut2, os.path.join(getTmpDir(),'testhippocut2.gii'))

                hippogii = aims.read(os.path.join(getTmpDir(),'testhippocut.gii'))
                hippovertex = numpy.array(hippogii.vertex().list())
                centerhippo = numpy.average(hippovertex,axis=0)

                hippogii2 = aims.read(os.path.join(getTmpDir(),'testhippocut2.gii'))
                hippovertex2 = numpy.array(hippogii2.vertex().list())
                centerhippo2 = numpy.average(hippovertex2,axis=0)

                coords1 = m.transform(centerhippo)
                coords2 = m.transform(centerhippo2)

                if coords2[1] > coords1[1]:
                    #2 est l'antero; 1 est posterieur
                    if dilefthippoantero is None or dilefthippopostero is None:
                        print "mesh export : could not find valid path"
                    else:
                        cmd1 = ['mv', os.path.join(getTmpDir(),'testhippocut2.gii'), str(dilefthippoantero.fullPath())]
                        cmd2 = ['mv', os.path.join(getTmpDir(),'testhippocut.gii'), str(dilefthippopostero.fullPath())]
                        line1 = runCmd(cmd1)
                        line2 = runCmd(cmd2)
                        neuroHierarchy.databases.insertDiskItem(dilefthippoantero, update = True)
                        neuroHierarchy.databases.insertDiskItem(dilefthippopostero, update = True)
                        if 'referential' in T1pre.attributes().keys():
                            self.transfoManager.setReferentialTo(dilefthippoantero, T1pre.attributes()['referential'] )
                            self.transfoManager.setReferentialTo(dilefthippopostero, T1pre.attributes()['referential'] )
                elif coords2[1] < coords1[1]:
                    #1 est l'antero; 2 est le posterieur
                    if dilefthippoantero is None or dilefthippopostero is None:
                        print "mesh export : could not find valid path"
                    else:
                        cmd1 = ['mv', os.path.join(getTmpDir(),'testhippocut.gii'), str(dilefthippoantero.fullPath())]
                        cmd2 = ['mv', os.path.join(getTmpDir(),'testhippocut2.gii'), str(dilefthippopostero.fullPath())]
                        line1 = runCmd(cmd1)
                        line2 = runCmd(cmd2)
                        neuroHierarchy.databases.insertDiskItem(dilefthippoantero, update = True)
                        neuroHierarchy.databases.insertDiskItem(dilefthippopostero, update = True)
                        if 'referential' in T1pre.attributes().keys():
                            self.transfoManager.setReferentialTo(dilefthippoantero, T1pre.attributes()['referential'] )
                            self.transfoManager.setReferentialTo(dilefthippopostero, T1pre.attributes()['referential'] )

                if dilefthippoantero is not None and dilefthippopostero is not None:
                    #we read it, estimate the number of vertex, and attribute the good value to the good number of vertex
                    #nb vertex of right antero
                    voxel_size_T1 = [T1pre.attributes()['voxel_size'][0], T1pre.attributes()['voxel_size'][1], T1pre.attributes()['voxel_size'][2], 1.0]

                    hippo_vol_antero=aims.Volume(*T1pre.attributes()['volume_dimension'][:3],dtype='S16')
                    hippo_vol_postero=aims.Volume(*T1pre.attributes()['volume_dimension'][:3],dtype='S16')
                    hippo_vol_full=aims.Volume(*T1pre.attributes()['volume_dimension'][:3],dtype='S16')
                    hippo_vol_antero.header()['voxel_size']= voxel_size_T1
                    hippo_vol_postero.header()['voxel_size']= voxel_size_T1
                    hippo_vol_full.header()['voxel_size']= voxel_size_T1

                    wdilefthippoNII =  WriteDiskItem('leftHippocampusNII', 'NIFTI-1 image' )
                    dileftHippoNII = wdilefthippoNII.findValue(constraints)

                    meshantero = aims.read(dilefthippoantero.fullPath())
                    meshpostero = aims.read(dilefthippopostero.fullPath())
                    aims.SurfaceManip.rasterizeMesh(meshantero,hippo_vol_antero,1)
                    aims.SurfaceManip.rasterizeMesh(meshpostero,hippo_vol_postero,1)
                    #fill the insides voxel
                    for zslices in range(hippo_vol_antero.arraydata().shape[1]):
                        hippo_vol_antero.arraydata()[0,zslices,:,:] = ndimage.morphology.binary_fill_holes(hippo_vol_antero.arraydata()[0,zslices,:,:]).astype(int)

                    for zslices in range(hippo_vol_postero.arraydata().shape[1]):
                        hippo_vol_postero.arraydata()[0,zslices,:,:] = ndimage.morphology.binary_fill_holes(hippo_vol_postero.arraydata()[0,zslices,:,:]).astype(int)

                    hippo_vol_antero *= 1701
                    hippo_vol_postero *= 1702
                    hippo_vol_full = hippo_vol_antero + hippo_vol_postero
                    for z in xrange(hippo_vol_full.getSizeZ()):
                        for y in xrange(hippo_vol_full.getSizeY()):
                            for x in xrange(hippo_vol_full.getSizeX()):
                                if hippo_vol_full.value(x, y, z) >1702:
                                    hippo_vol_full.setValue(1702, x, y, z)

                    aims.write(hippo_vol_full,str(dileftHippoNII.fullPath()))
                    neuroHierarchy.databases.insertDiskItem(dileftHippoNII, update = True)
                    if 'referential' in T1pre.attributes().keys():
                        self.transfoManager.setReferentialTo(dileftHippoNII, T1pre.attributes()['referential'] )


        print "generation of amygdala and hippocamp meshes done"
        self.setStatus(u"generation of amygdala and hippocamp meshes done")


    def fitvolumebyellipse(self,volumeposition):

        #cut the hippocampus in two:
        # from https://github.com/minillinim/ellipsoid/blob/master/ellipsoid.py
        tolerance = 0.015
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

        return UU, ss, rotation, center

    def setNormalizedToT1pre(self, image, registeredPath):

        ret = subprocess.call(['AimsFileConvert', '-i', str(registeredPath), '-o', str(image.fileName()), '-t', 'S16'])
        if ret < 0:
            print "Error: Cannot import the SPM resampled image : "+repr(registeredPath)
            #QtGui.QMessageBox.warning(self, "Error", u"The image has not been resampled by SPM !")
            return
        # TODO # Should destroy existing referentials for this image
        # The referential of this image is the same as the T1 pre
        self.transfoManager.setReferentialTo(image, registration.talairachMNIReferentialId) #replace self.getT1preNativeRef(image.attributes()['center'], image.attributes()['subject']) by talairachMNIReferentialId
        #os.remove(registeredPath)
    

    def taskfinished(self, message, threadObj=None):
        """ When a task (thread) is finished, display the message and launch the waiting callbacks"""
        self.setStatus(u"Task done : "+message) #+' : '+repr(threadObj)
        if threadObj is None:
            return
        # Looking for callback functions waiting
        #print "****** WAITING FOR THREADS *******\n"+repr(self.waitingForThreads) # DEBUG
        toDelete = [] # Store the dictionary keys to delete (cannot be done while iterating on it)
        for callback in self.waitingForThreads:
            for thr in self.waitingForThreads[callback]: # All threads this callback is waiting for
                if thr not in self.threads: # This one does not exist, it was probably deleted before the callback was set, so just remove it from the list
                    self.waitingForThreads[callback].remove(thr)
            if threadObj in self.waitingForThreads[callback]: # This thread was in the list
                self.waitingForThreads[callback].remove(threadObj) # So remove it !
            if self.waitingForThreads[callback] == []: # If there is no thread left
                callback() # Call the function !
                toDelete.append(callback)
        for tD in toDelete:
            del self.waitingForThreads[tD]
    
        if threadObj in self.threads:
            idx = self.threads.remove(threadObj) # Now remove the thread

    # Preference tab : set the SPM T1 template path
    def setSpmTemplatePath(self, path=None):
        if path is None:
            path = QtGui.QFileDialog.getExistingDirectory(self, u"Select SPM path")
        if path is not None:
            self.ui.prefSpmTemplateEdit.setText(path)

    def setANTsPath(self, path=None):
        if path is None:
            path = QtGui.QFileDialog.getExistingDirectory(self, u"Select ANTs path")
        if path is not None:
            self.ui.prefANTsEdit.setText(path)

    def setFreesurferPath(self, path=None):
        if path is None:
            path = QtGui.QFileDialog.getExistingDirectory(self, u"Select FREESURFER path")
        if path is not None:
            self.ui.prefFreesurferEdit.setText(path)

    def setBidsPath(self, path=None):
        if path is None:
            path = QtGui.QFileDialog.getExistingDirectory(self, u"Select BIDS database")
        if path is not None:
            self.ui.prefBidsEdit.setText(path)
            
    def setPrefCoregister(self,key):
        if key == 'ANTS':
            if self.ui.prefANTScheckbox.isChecked():
                self.ui.prefSPMcheckbox.setCheckState(False)
        elif key == 'SPM':
            if self.ui.prefSPMcheckbox.isChecked():
                self.ui.prefANTScheckbox.setCheckState(False)
    
    def switchProjectButton(self):
        projects= [self.ui.radioProjFtract.isChecked(),self.ui.radioProjNeuro.isChecked(), self.ui.radioClassic.isChecked(),self.ui.radioOther.isChecked()]
        project_selected = [x for x in range(len(projects)) if projects[x]==True]
        self.prefs['projectSelected'] = project_selected


    ################################## Patient Information Tab #####################################
    def selectPatSubject(self, subj):
        pass

    def checkbox_comor(self):
        if self.patientInfo['comoraucune'].isChecked():
            self.patientInfo['comorpsy'].setEnabled(False)
            self.patientInfo['comorautre'].setEnabled(False)
            self.patientInfo['comorneuro'].setEnabled(False)
        else:
            self.patientInfo['comorpsy'].setEnabled(True)
            self.patientInfo['comorautre'].setEnabled(True)
            self.patientInfo['comorneuro'].setEnabled(True)

    def update_patientAge(self):
        patAge = self.patientInfo['patientBirthday'].date().daysTo( self.patientInfo['currentDate'].date())/365
        self.patientInfo['patientAge'].setText(str(patAge))

    def ValidatePatInfo(self):
        patListToCheck = ['previousHistory','previousFamilyHistory','causaleDisease','mriLesion','comorneuro','comorpsy','comorautre']
        patpatListToCheck = {'Epilepsy':['Aura','TraitementTried','TraitementNow']}

        #read choice2 list from intersubject database
        write_filters = { 'center': self.currentProtocol, 'subject' : self.currentSubject }
        wdi_global = ReadDiskItem('PatientInfoTemplate','Patient Template format')
        di_global = list(wdi_global.findValues({}, None, False))

        if len(di_global):
            if os.path.isfile(str(di_global[0])):
                print "read previous patienttemplate json"
                from codecs import open as opencodecs
                fin = opencodecs(str(di_global[0]),'r','latin1')
                info_dicti = json.loads(fin.read().decode('latin1'))
                fin.close()
    
                previous_lists_not_path = info_dicti['notPathoSpecific']
                previous_lists_path_full = info_dicti['PathoSpecific']
                previous_lists_path_protocol = info_dicti['PathoSpecific'][self.currentProtocol]
    
        else:
            previous_lists_not_path = {}
            previous_lists_path_full = {}
            previous_lists_path_protocol = {}
    
        #register non pathology specific patient Info
        PatInfo = {}
        for kk, vv in self.patientInfo.iteritems():
            if isinstance(vv,QtGui.QDateEdit):
                PatInfo.update({kk:unicode(vv.date().toPyDate())})
                if vv.date() == self.defaultAcqDate:
                    QtGui.QMessageBox.warning(self, "Error",u"The acquisition date is not valid !  " + kk)
                    return
            elif isinstance(vv,QtGui.QComboBox):
                #faire deux conditions lorsque c'est une liste checkable et lorsque non
                if vv.model().item(1).isCheckable():
                    item_checked = [unicode(vv.itemText(x)) for x in range(1,vv.count()) if vv.model().item(x).checkState()]
                    if len(item_checked) > 1:
                        if 'Inconnu' in item_checked or 'Aucun' in item_checked or 'Aucune' in item_checked:
                            QtGui.QMessageBox.warning(self, "Error",u"You can't have 'Inconnu' or 'Aucune' selected whith another choice, it doesn't make sens" + kk)
                            return
                    #c est une combo box editable du coup il faut verifier si il faut updater la liste des choix possibles
    
                else:
                    item_checked = unicode(vv.currentText())
                PatInfo.update({kk:item_checked})
            elif isinstance(vv,QtGui.QDialogButtonBox):
                item_selected = [unicode(vv.children()[i].text()) for i in range(1,len(vv.children())) if vv.children()[i].isChecked()]
                PatInfo.update({kk:item_selected})
            elif isinstance(vv,QtGui.QLineEdit):
                PatInfo.update({kk:unicode(vv.text())})
            elif isinstance(vv,QtGui.QCheckBox):
                PatInfo.update({kk:vv.isChecked()})
            else:
                print "qtgui format not recognized"
    
        patPatInfo = {}
        for kk, vv in self.pathologypatientInfo.iteritems():
            if isinstance(vv,QtGui.QDateEdit):
                patPatInfo.update({kk:unicode(vv.date().toPyDate())})
                if vv.date() == self.defaultAcqDate:
                    QtGui.QMessageBox.warning(self, "Error",u"The acquisition date is not valid !  " + kk)
                    return
            elif isinstance(vv,QtGui.QComboBox):
                #faire deux conditions lorsque c'est une liste checkable et lorsque non
                if vv.model().item(1).isCheckable():
                    item_checked = [unicode(vv.itemText(x)) for x in range(1,vv.count()) if vv.model().item(x).checkState()]
                    if len(item_checked) > 1:
                        if 'Inconnu' in item_checked or 'Aucun' in item_checked or 'Aucune' in item_checked:
                            QtGui.QMessageBox.warning(self, "Error",u"You can't have 'Inconnu' or 'Aucune' selected whith another choice, it doesn't make sens"+ kk)
                            return
                else:
                    item_checked = unicode(vv.currentText())
                patPatInfo.update({kk:item_checked})
            elif isinstance(vv,QtGui.QDialogButtonBox):
                item_selected = [unicode(vv.children()[i].text()) for i in range(1,len(vv.children())) if vv.children()[i].isChecked()]
                patPatInfo.update({kk:item_selected})
            elif isinstance(vv,QtGui.QLineEdit):
                patPatInfo.update({kk:unicode(vv.text())})
            elif isinstance(vv,QtGui.QCheckBox):
                patPatInfo.update({kk:vv.isChecked()})
            else:
                print "qtgui format not recognized"
    
        wdi = WriteDiskItem('SubjectInfo','Subject Information format')
        di = wdi.findValue(write_filters)
        if di is None:
            print('Can t generate files')
            return
    
        full_dictio = {'notPathoSpecific':PatInfo,'PathoSpecific':patPatInfo}
        fout = open(di.fullPath(),'w')
        fout.write(json.dumps(full_dictio, ensure_ascii=False))
        #fout.write(json.dumps({'PathoSpecific':patPatInfo}))
        fout.close()
    
        neuroHierarchy.databases.insertDiskItem(di, update=True )
    
        #update patientinfo list si une nouvelle entree a ete ajoutée dans les choices list
        wdi_global = WriteDiskItem('PatientInfoTemplate','Patient Template format')
        di_global = wdi_global.findValue(write_filters)
    
        for ii in range(len(patListToCheck)):
            if len(PatInfo[patListToCheck[ii]])>0:
                previous_lists_not_path.update({patListToCheck[ii]:PatInfo[patListToCheck[ii]]})
    
        if self.currentProtocol in patpatListToCheck:
            for jj in range(len(patpatListToCheck[self.currentProtocol])):
                if len(patPatInfo[patpatListToCheck[self.currentProtocol][jj]])>0:
                    previous_lists_path_protocol.update({patpatListToCheck[self.currentProtocol][jj]:patPatInfo[patpatListToCheck[self.currentProtocol][jj]]})
    
        previous_lists_path_full.update({self.currentProtocol:previous_lists_path_protocol})
        full_dictio_inter = {'notPathoSpecific':previous_lists_not_path,'PathoSpecific':previous_lists_path_full}
    
        fout = open(di_global.fullPath(),'w')
        fout.write(json.dumps(full_dictio_inter,ensure_ascii=False))
        fout.close()
        neuroHierarchy.databases.insertDiskItem(di_global, update=True )


    def enable_disable_gadooption(self):

        if str(self.ui.niftiSeqType.currentText()) == 'T1':
            self.ui.radioButtonGado.setEnabled(True)
            self.ui.radioButtonNoGado.setEnabled(True)
        else:
            self.ui.radioButtonGado.setEnabled(False)
            self.ui.radioButtonNoGado.setEnabled(False)


# =============================================================================
# MAIN: Main function that starts the interface
# =============================================================================
def main():
    
    # Create application
    if hasattr(QtCore.Qt, 'AA_X11InitThreads'):
        QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_X11InitThreads)
    app = QtGui.QApplication(sys.argv)
    axon.initializeProcesses()
    
    # Show main window
    w = ImageImport(app = app)
    w.setWindowFlags(QtCore.Qt.Window)
    w.show()
    
    # Kill application when windows are closed
    # QtGui.QObject.connect(app, QtCore.SIGNAL('lastWindowClosed()'), app, QtCore.SLOT('quit()'))
    # Debug -> evite un pb entre ipython, pdb et qt
    # pyqtRemoveInputHook()
    
    # Run the application
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()


