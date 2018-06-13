# -*- coding: utf-8 -*-
#
# Importation of data in the database + image registration and normalization
#(c) Inserm U836 2012-2014 - Manik Bhattacharjee
#
# License GNU GPL v3
#
#
#TODO
#
# Tester et finaliser l'import depuis un PACS
#
# Verifier qu'on ne convertit pas une image déjà récupérée (surtout DICOM-> Nifti !)
#
# Améliorer le choix de l'image pre/post dans BrainVisa : afficher des aperçus avec Anatomist après l'import DICOM ?
# Ou coder en C++ un créateur d'aperçu DICOM avec la lib dcmtk, ou utiliser http://sourceforge.net/projects/qdcws/
#
#Optimisation : on peut supposer que le nombre de slices
# (SliceNumberMR) quand on ouvre un fichier permet de
# sauter les SliceNumberMR-1 suivants car ce sont les
# autres slices de la même image
#
# Onglet registration : mettre en vert les images qui ont déjà un recalage vers la T1pre, en gras l'image normalisée vers le template SPM
# FIXME  La version suivante de BrainVisa (post 4.3) devrait corriger un bug dans l'insertion des items dans la base (en insérant les répertoires qui les contiennent)
# Il faudra alors virer la fonction createItemDirs que j'ai créé pour résoudre ce pb. D'autre part pour résoudre le pb pour les fichiers registration
# j'ai ajouté l'insertion des répertoires dans data/neuroDiskItems.py (fonction createParentDirectory) mais ça ne marche pas (import circulaires !)
# Donc à chaque ajout de référentiel ou de transfo, j'ai ajouté une ligne neuroHierarchy.databases.createDiskItemFromFileName A SUPPRIMER quand BrainVisa sera corrigé (après BV 4.3)
        #dirItem=neuroHierarchy.databases.createDiskItemFromFileName(d, None)
# ipython -q4thread  locateElectrodes.py # Pour avoir ipython et Qt actifs en même temps au cours du debug


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

from externalprocesses import *
from dicomutilities import *
import seegprocessing, patientinfo, pathologypatientinfo
from checkSpmVersion import *
from freesurfer.brainvisaFreesurfer import *
from TimerMessageBox import *
from progressbar import ProgressDialog


#  Matlab code : coregister file1 to file2
def matlab_cellstr(listOfStrings):
    """ Converts a list of filenames as a matlab-suitable cell table of strings and adds ',1' at the end of each string (for spm file lists) """
    return "{'"+",1' '".join([str(stri) for stri in listOfStrings])+",1'}"

# SPM coregistration onto file1 ({'/home/manik/data/epilepsie/IRM-testOD/Pre/3DT1.img,1'}) of file 2 ({'/home/manik/data/epilepsie/IRM-testOD/Post/SagT2.img,1'})
#spm_coregister = "try, spm('defaults', 'FMRI');spm_jobman('initcfg');\
#matlabbatch{1}.spm.spatial.coreg.estimate.ref = %s;\
#matlabbatch{1}.spm.spatial.coreg.estimate.source = %s;\
#matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};\
#matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';\
#matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];\
#matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];\
#matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];\
#spm_jobman('run',matlabbatch); catch, disp 'AN ERROR OCCURED'; end;quit;"

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


# SPM normalisation of file 1 ({'/home/manik/data/epilepsie/IRM-testOD/Pre/3DT1.img,1'}) onto template
#spm_template_t1 = "{'/home/manik/prog/implantation/matlab/spm8/templates/T1.nii,1'}"
#spm_template_t1 = "{'/matlab_prog/spm8_5236/templates/T1.nii,1'}"
#spm8_normalise = """try, addpath(genpath(%s)); spm('defaults', 'FMRI');spm_jobman('initcfg');
#matlabbatch{1}.spm.spatial.normalise.est.subj.source = %s;
#matlabbatch{1}.spm.spatial.normalise.est.subj.wtsrc = '';
#matlabbatch{1}.spm.spatial.normalise.est.eoptions.template = %s;
#matlabbatch{1}.spm.spatial.normalise.est.eoptions.weight = '';
#matlabbatch{1}.spm.spatial.normalise.est.eoptions.smosrc = 8;
#matlabbatch{1}.spm.spatial.normalise.est.eoptions.smoref = 0;
#matlabbatch{1}.spm.spatial.normalise.est.eoptions.regtype = 'mni';
#matlabbatch{1}.spm.spatial.normalise.est.eoptions.cutoff = 25;
#matlabbatch{1}.spm.spatial.normalise.est.eoptions.nits = 16;
#matlabbatch{1}.spm.spatial.normalise.est.eoptions.reg = 1;
#spm_jobman('run',matlabbatch); catch, disp 'AN ERROR OCCURED'; end;quit;"""

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
        #self.modas = {'t1mri':'Raw T1 MRI', 't2mri':'T2 MRI', 'ct': 'CT', 'pet': 'PET','fmri_epile':'fMRI-epile', 'statistic_data':'Statistic-Data', 'flair': 'FLAIR', 'freesurfer_atlas':'FreesurferAtlas', 'fgatir':'FGATIR', 'hippofreesurfer_atlas':'HippoFreesurferAtlas'}
        self.modas = {'t1mri':'Raw T1 MRI', 't2mri':'T2 MRI', 'ct': 'CT', 'pet': 'PET','fmri_epile':'fMRI-epile', 'statistic_data':'Statistic-Data', 'flair': 'FLAIR', 'freesurfer_atlas':'FreesurferAtlas', 'fgatir':'FGATIR'}

        # Allow calling brainvisa
        self.brainvisaContext = defaultContext()
        # Get Transformation Manager
        self.transfoManager = registration.getTransformationManager()

        # Load anatomist
        self.a = anatomist.Anatomist('-b')
        # Force all the scanner-based referentials to be linked 
        print("\nAnatomist is configured to ignore differences in Scanner-based referentials:\n   - setAutomaticReferential = 1\n   - commonScannerBasedReferential = 1\n")
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

        self.dispObj = []

        self.initialize()

        self.ui.dicomOutputPathEdit.setText(self.dicomOutputPath)
        self.ui.niftiOutputPathEdit.setText(self.niftiOutputPath)

        self.threads = []
        self.waitingForThreads={}
        # Init the patient informations tab
        self.patientInfo = patientinfo.buildUI(self.ui.patGroupBox)
        cdate = self.patientInfo['currentDate'].date()
        bdate = self.patientInfo['patientBirthday'].date()
        page = bdate.daysTo(cdate)/365
        self.patientInfo['patientAge'].setText(str(page))
    
        self.patientInfo['comoraucune'].stateChanged.connect(self.checkbox_comor)
    
        self.patientInfo['patientBirthday'].dateChanged.connect(self.update_patientAge)
        #self.pathologypatientInfo = pathologypatientinfo.buildUI(self.ui.patpatGroupBox,self.currentProtocol)
    
        # Recharge la BDD BrainVisa à chaque changement d'onglet
        self.connect(self.ui.tabWidget,  QtCore.SIGNAL('currentChanged(int)'), self.analyseBrainvisaDB)
        self.connect(self.ui.tabWidget,  QtCore.SIGNAL('currentChanged(int)'), self.setTabGuide)
        # Brainvisa database tab
        self.connect(self.ui.bvProtocolCombo, QtCore.SIGNAL('currentIndexChanged(QString)'), lambda s:self.selectProtocol(s,self.ui.bvSubjectCombo))
        self.connect(self.ui.bvSubjectCombo, QtCore.SIGNAL('currentIndexChanged(QString)'), self.selectBvSubject)
        self.connect(self.ui.bvSubjectCombo, QtCore.SIGNAL('activated(QString)'), self.setCurrentSubject)
        self.connect(self.ui.bvImageList, QtCore.SIGNAL('itemDoubleClicked(QListWidgetItem*)'), self.selectBvImage)
        self.connect(self.ui.bvDeleteImageButton, QtCore.SIGNAL('clicked()'), self.deleteBvImage)
        self.connect(self.ui.bvDeleteSubjectButton, QtCore.SIGNAL('clicked()'), self.deleteBvSubject)
        self.connect(self.ui.bvEditPref, QtCore.SIGNAL('clicked()'), self.editBvPref)
        self.connect(self.ui.bvUpdateDb, QtCore.SIGNAL('clicked()'), self.updateBvDb)

        # Add Subject tab
        self.connect(self.ui.subjectSiteCombo, QtCore.SIGNAL('activated(QString)'), self.updatePatientCode)
        self.connect(self.ui.subjectSiteCombo, QtCore.SIGNAL('editTextChanged(QString)'), self.updatePatientCode)
        self.connect(self.ui.subjectYearSpinbox, QtCore.SIGNAL('valueChanged(int)'), self.updatePatientCode)
        self.connect(self.ui.subjectPatientName, QtCore.SIGNAL('textChanged(QString)'), self.updatePatientCode)
        self.connect(self.ui.subjectPatientFirstName, QtCore.SIGNAL('textChanged(QString)'), self.updatePatientCode)
        self.connect(self.ui.subjectAddSubjectButton, QtCore.SIGNAL('clicked()'), self.storePatientInDB)
        # DICOM tab
        self.connect(self.ui.chooseDirButton, QtCore.SIGNAL('clicked()'), self.openDicomDir)
        self.connect(self.ui.dirImportButton, QtCore.SIGNAL('clicked()'), self.importDataDir)
        self.connect(self.ui.dicomOutputButton, QtCore.SIGNAL('clicked()'), lambda :self.selectOutput('Dicom'))
        self.connect(self.ui.niftiOutputButton, QtCore.SIGNAL('clicked()'), lambda :self.selectOutput('Nifti'))
        # PACS tab
        self.connect(self.ui.pacsSearchButton, QtCore.SIGNAL('clicked()'), self.searchPacs)
        self.connect(self.ui.dirPatientList, QtCore.SIGNAL('itemSelectionChanged()'), lambda :self.patientSelectionChanged('dir'))
        self.connect(self.ui.pacsPatientList, QtCore.SIGNAL('itemSelectionChanged()'), lambda :self.patientSelectionChanged('pacs'))
        self.connect(self.ui.dirStudiesList, QtCore.SIGNAL('itemSelectionChanged()'), lambda :self.studiesSelectionChanged('dir'))
        self.connect(self.ui.pacsStudiesList, QtCore.SIGNAL('itemSelectionChanged()'), lambda :self.studiesSelectionChanged('pacs'))
        self.connect(self.ui.chooseNiftiButton, QtCore.SIGNAL('clicked()'), self.chooseNifti)
        self.connect(self.ui.pacsImportButton, QtCore.SIGNAL('clicked()'), self.downloadPacs)
        # Import to BrainVisa for DICOM/NIFTI/PACS
        #NIFTI
        self.connect(self.ui.niftiImportButton, QtCore.SIGNAL('clicked()'), self.importNifti)
        self.connect(self.ui.niftiSetBraincenterButton, QtCore.SIGNAL('clicked()'), self.setBrainCenter)
        self.connect(self.ui.ImportFSoutputspushButton, QtCore.SIGNAL('clicked()'),self.importFSoutput)
        #PACS
        self.connect(self.ui.pacsImportBVButton, QtCore.SIGNAL('clicked()'), self.importPacs)
        #DICOM
        self.connect(self.ui.dirImportBVButton, QtCore.SIGNAL('clicked()'), self.importDir)
        # NIFTI
        self.connect(self.ui.niftiProtocolCombo, QtCore.SIGNAL('currentIndexChanged(QString)'), lambda s:self.selectProtocol(s,self.ui.niftiSubjectCombo))
        self.connect(self.ui.niftiSubjectCombo, QtCore.SIGNAL('activated(QString)'), self.setCurrentSubject)
        self.connect(self.ui.niftiSeqType, QtCore.SIGNAL('currentIndexChanged(QString)'),self.enable_disable_gadooption)
        # Registration tab
        self.connect(self.ui.regProtocolCombo, QtCore.SIGNAL('currentIndexChanged(QString)'), lambda s:self.selectProtocol(s,self.ui.regSubjectCombo))
        self.connect(self.ui.regSubjectCombo, QtCore.SIGNAL('currentIndexChanged(QString)'), self.selectRegSubject)
        self.connect(self.ui.regSubjectCombo, QtCore.SIGNAL('activated(QString)'), self.setCurrentSubject)
        self.connect(self.ui.regImageList, QtCore.SIGNAL('itemDoubleClicked(QListWidgetItem*)'), self.selectRegImage)
        self.connect(self.ui.regImageList2, QtCore.SIGNAL('itemDoubleClicked(QListWidgetItem*)'), self.selectRegImage2)
        #self.connect(self.ui.registerNormalizeSubjectButton, QtCore.SIGNAL('clicked()'), self.registerNormalizeSubject)
        self.connect(self.ui.registerNormalizeSubjectButton, QtCore.SIGNAL('clicked()'), lambda :ProgressDialog.call(self.registerNormalizeSubject, False, self, "Processing...", "Coregister and normalize"))
        self.connect(self.ui.segmentationHIPHOPbutton,QtCore.SIGNAL('clicked()'),self.runPipelineBV)
        self.connect(self.ui.FreeSurferReconAllpushButton,QtCore.SIGNAL('clicked()'),self.runFreesurferReconAll)
        self.ui.FreeSurferReconAllpushButton.setEnabled(False)
        self.connect(self.ui.runMarsAtlasFreesurferButton,QtCore.SIGNAL('clicked()'),self.runPipelineFS)
    
        self.warningMEDIC()

        def setDictValue(d, k, v,button):
            d[k]=v
            button.setStyleSheet("background-color: rgb(90, 255, 95);")
        self.connect(self.ui.regAcButton, QtCore.SIGNAL('clicked()'), lambda:setDictValue(self.AcPc,'AC',list(self.a.linkCursorLastClickedPosition()),self.ui.regAcButton))
        self.connect(self.ui.regPcButton, QtCore.SIGNAL('clicked()'), lambda:setDictValue(self.AcPc,'PC',list(self.a.linkCursorLastClickedPosition()),self.ui.regPcButton))
        self.connect(self.ui.regIhButton, QtCore.SIGNAL('clicked()'), lambda:setDictValue(self.AcPc,'IH',list(self.a.linkCursorLastClickedPosition()),self.ui.regIhButton))
        self.connect(self.ui.regLhButton, QtCore.SIGNAL('clicked()'), lambda:setDictValue(self.AcPc,'LH',list(self.a.linkCursorLastClickedPosition()),self.ui.regLhButton))
        self.connect(self.ui.regAcPcValidateButton, QtCore.SIGNAL('clicked()'), self.validateAcPc)
        # Implantation tab (surgery drawings)
        self.connect(self.ui.implProtocolCombo, QtCore.SIGNAL('currentIndexChanged(QString)'), lambda s:self.selectProtocol(s,self.ui.implSubjectCombo))
        self.connect(self.ui.implSubjectCombo, QtCore.SIGNAL('currentIndexChanged(QString)'), self.selectImplSubject)
        self.connect(self.ui.implSubjectCombo, QtCore.SIGNAL('activated(QString)'), self.setCurrentSubject)
        self.connect(self.ui.implCoroOpenButton, QtCore.SIGNAL('clicked()'), self.openImplCoro)
        self.connect(self.ui.implCoroImportButton, QtCore.SIGNAL('clicked()'), self.importImplCoro)
        self.connect(self.ui.implSagOpenButton, QtCore.SIGNAL('clicked()'), self.openImplSag)
        self.connect(self.ui.implSagImportButton, QtCore.SIGNAL('clicked()'), self.importImplSag)
        self.connect(self.ui.implPptOpenButton, QtCore.SIGNAL('clicked()'), self.openImplPpt)
        self.connect(self.ui.implPptImportButton, QtCore.SIGNAL('clicked()'), self.importImplPpt)
        self.connect(self.ui.implPdfOpenButton, QtCore.SIGNAL('clicked()'), self.openImplPdf)
        self.connect(self.ui.implPdfImportButton, QtCore.SIGNAL('clicked()'), self.importImplPdf)
        self.connect(self.ui.implPdfOpenButton2, QtCore.SIGNAL('clicked()'), self.openImplPdf2)
        self.connect(self.ui.implPdfImportButton2, QtCore.SIGNAL('clicked()'), self.importImplPdf2)
        self.connect(self.ui.implElListPdfOpenButton, QtCore.SIGNAL('clicked()'), self.openImplElListPdf)
        self.connect(self.ui.implElListPdfImportButton, QtCore.SIGNAL('clicked()'), self.importImplElListPdf)
        # Infos patient tab
        self.connect(self.ui.patProtocolCombo, QtCore.SIGNAL('currentIndexChanged(QString)'), lambda s:self.selectProtocol(s,self.ui.patSubjectCombo))
        #self._pathologyinfo = lambda :pathologypatientinfo.buildUI()
        self.connect(self.ui.patProtocolCombo, QtCore.SIGNAL('activated(QString)'),lambda :pathologypatientinfo.buildUI(self.ui.patpatGroupBox,self.currentProtocol)) #lambda :pathologypatientinfo.buildUI(self.ui.patpatGroupBox,self.currentProtocol)) self._pathologyinfo(self.ui.patpatGroupBox,self.currentProtocol)
        self.connect(self.ui.patSubjectCombo, QtCore.SIGNAL('currentIndexChanged(QString)'), self.selectPatSubject)
        self.connect(self.ui.patSubjectCombo, QtCore.SIGNAL('activated(QString)'), self.setCurrentSubject)
        self.connect(self.ui.patInfoValidate,QtCore.SIGNAL('clicked()'),self.ValidatePatInfo)
    
        # Importation SEEG tab
        self.connect(self.ui.seegProtocolCombo, QtCore.SIGNAL('currentIndexChanged(QString)'), lambda s:self.selectProtocol(s,self.ui.seegSubjectCombo))
        self.connect(self.ui.seegSubjectCombo, QtCore.SIGNAL('currentIndexChanged(QString)'), self.selectSeegSubject)
        self.connect(self.ui.seegSubjectCombo, QtCore.SIGNAL('activated(QString)'), self.setCurrentSubject)
        self.connect(self.ui.seegManipsBaseCombo, QtCore.SIGNAL('activated(QString)'), self.manipSelected)
        self.connect(self.ui.chooseSeegButton, QtCore.SIGNAL('clicked()'), self.chooseSeeg)
        self.connect(self.ui.seegImportButton, QtCore.SIGNAL('clicked()'), self.seegImport)
    
        #Preference tab
        self.connect(self.ui.prefSpmTemplateButton, QtCore.SIGNAL('clicked()'), self.setSpmTemplatePath)
        self.connect(self.ui.prefANTsButton, QtCore.SIGNAL('clicked()'), self.setANTsPath)
        self.connect(self.ui.prefSpmStandaloneButton, QtCore.SIGNAL('clicked()'), self.setSpmStandalonePath)
        self.connect(self.ui.prefFreesurferButton, QtCore.SIGNAL('clicked()'), self.setFreesurferPath)
        self.connect(self.ui.prefNoDBFileLocationButton, QtCore.SIGNAL('clicked()'),self.setNoDBFilePath)
        self.connect(self.ui.prefANTScheckbox, QtCore.SIGNAL('clicked()'), lambda: self.setPrefCoregister('ANTS'))
        self.connect(self.ui.prefSPMcheckbox, QtCore.SIGNAL('clicked()'), lambda: self.setPrefCoregister('SPM'))
        self.connect(self.ui.prefAimscheckbox, QtCore.SIGNAL('clicked()'), lambda: self.setPrefCoregister('Aims'))
        self.ui.radioButtonProject_4.toggled.connect(self.switchProjectButton)
        self.ui.radioButtonProject_3.toggled.connect(self.switchProjectButton)
        self.ui.radioButtonProject_2.toggled.connect(self.switchProjectButton)
        self.ui.radioButtonProject.toggled.connect(self.switchProjectButton)
        self.ui.radioButtonProject_4.toggled.connect(self.updatePatientCode)
        self.ui.radioButtonProject_3.toggled.connect(self.updatePatientCode)
        self.ui.radioButtonProject_2.toggled.connect(self.updatePatientCode)
        self.ui.radioButtonProject.toggled.connect(self.updatePatientCode)
        self.ui.SavePreferencespushButton.clicked.connect(self.savePreferences)
    
        #self.connect(self.ui.radioButtonProject_3, QtCore.SIGNAL('clicked()'), lambda: self.switchProjectButton)
        #self.connect(self.ui.radioButtonProject_2, QtCore.SIGNAL('clicked()'), lambda: self.switchProjectButton)
        #self.connect(self.ui.radioButtonProject, QtCore.SIGNAL('clicked()'), lambda: self.switchProjectButton)
        #self.connect(self.ui.radioButtonProject_3, QtCore.SIGNAL('clicked()'), lambda: self.updatePatientCode)
        #self.connect(self.ui.radioButtonProject_2, QtCore.SIGNAL('clicked()'), lambda: self.updatePatientCode)
        #self.connect(self.ui.radioButtonProject, QtCore.SIGNAL('clicked()'), lambda: self.updatePatientCode)
    
        # Finds the available protocols in brainvisa database and fills the comboboxes
        self.analyseBrainvisaDB()
        self.setGuide("Choose on of the tabs to import data in the BrainVisa database")
        self.enableACPCButtons(False)
        # Disable work-in-progress tabs
        self.ui.tabWidget.setTabEnabled(6, False)
        self.ui.tabWidget.setTabEnabled(7, False)
        self.ui.tabWidget.setTabEnabled(8, False)
        self.ui.tabWidget.setTabEnabled(9, False)



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


    def initialize(self):
        """Reads the preference file, checks the available patients and sequences"""
        #Localisation du fichier de preference
        prefpath = os.path.join(os.path.expanduser('~'), '.imageimport')

        # Chargement du fichier et stockage dans self.prefs
        # Basic prefs
        self.prefs = {'dicomOutputPath':os.path.join(os.path.expanduser('~'), 'dicomOut'),'niftiOutputPath':os.path.join(os.path.expanduser('~'), 'niftiOut'), 'sites':['Gre',]}
        try:
            if (os.path.exists(prefpath)):
                filein = open(prefpath, 'rb')
                #self.prefs = pickle.load(filein)
                #fout.write(json.dumps({'mars_atlas':resection_mars_atlas_info}))
                try:
                    self.prefs = json.loads(filein.read())
                except:
                    filein.close()
                    filein = open(prefpath, 'rb')
                    self.prefs = pickle.load(filein)
    
                filein.close()
        except:
            print("There is a problem with %s opening apparently"%prefpath)

        self.dicomOutputPath = self.prefs['dicomOutputPath']
        self.niftiOutputPath = self.prefs['niftiOutputPath']

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
        brainvisa_spm12_standalone_path = configuration.SPM.spm12_standalone_path
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
                elf.setSpmTemplatePath(brainvisa_spm12_path)
                self.prefs['spm']=brainvisa_spm12_path
            else:
                QtGui.QMessageBox.warning(self, u"SPM", u"SPM path is not defined (or wrong) in tab 'Préférences'\n Normalization and SPM coregistration won't work !")

        if 'spm_standalone' in self.prefs:
            if brainvisa_spm12_standalone_path != self.prefs['spm_standalone']:
                QtGui.QMessageBox.warning(self, u"SPM standalone", u"SPM standalone path different between IntrAnat and BrainVisa, strange, you should check that, by default keep the one precised in IntrAnat")
            self.setSpmStandalonePath(self.prefs['spm_standalone'])
        else:
            self.setSpmStandalonePath(brainvisa_spm12_standalone_path)

        if 'freesurfer' in self.prefs:
            if brainvisa_freesurfer_home_path != self.prefs['freesurfer']:
                QtGui.QMessageBox.warning(self, u"Freesurfer", u"Freesurfer path different between IntrAnat and BrainVisa, strange, you should check that, by default keep the one precised in IntrAnat")
            self.setFreesurferPath(self.prefs['freesurfer'])
        else:
            self.setFreesurferPath(brainvisa_freesurfer_home_path)

            if 'FileNoDBPath' in self.prefs:
                self.setNoDBFilePath(self.prefs['FileNoDBPath'])


        #self.prefs['freesurfer'] = str(self.ui.prefFreesurferEdit.text())

        if 'ants' in self.prefs:
            self.setANTsPath(self.prefs['ants'])

        if 'coregisterMethod' in self.prefs:
            self.coregMethod = self.prefs['coregisterMethod']

            if self.coregMethod == 'ANTs':
                self.ui.prefSPMcheckbox.setCheckState(False)
                self.ui.prefANTScheckbox.setCheckState(True)
                self.ui.prefAimscheckbox.setCheckState(False)
            if self.coregMethod == 'spm':
                self.ui.prefSPMcheckbox.setCheckState(True)
                self.ui.prefANTScheckbox.setCheckState(False)
                self.ui.prefAimscheckbox.setCheckState(False)
            if self.coregMethod == 'Aims':
                self.ui.prefSPMcheckbox.setCheckState(False)
                self.ui.prefANTScheckbox.setCheckState(False)
                self.ui.prefAimscheckbox.setCheckState(True)

        if 'projectSelected' in self.prefs:
            if self.prefs['projectSelected'] != []:
                projects= [self.ui.radioButtonProject,self.ui.radioButtonProject_2,self.ui.radioButtonProject_3,self.ui.radioButtonProject_4]
                projects[self.prefs['projectSelected'][0]].setChecked(True)
        else:
            self.prefs.update({'projectSelected':[2]})
            projects= [self.ui.radioButtonProject,self.ui.radioButtonProject_2,self.ui.radioButtonProject_3,self.ui.radioButtonProject_4]
            projects[self.prefs['projectSelected'][0]].setChecked(True)



        # Test if matlab is available
        if 'matlabOk' not in self.prefs or self.prefs['matlabOk'] == False:
            if matlabIsPresent():
                self.prefs['matlabOk'] = True
            else:
                QtGui.QMessageBox.warning(self, u"Matlab", u"The command 'matlab' doesn't work !\nSPM coregistration and normalization won't work")
                self.prefs['matlabOk'] = False
        # Scanner les dossiers et séquences déjà importées, en faire un dictionnaire :
        self.dicomOutContent = self.analyzeDicomOutput(self.dicomOutputPath)
        self.niftiOutContent = self.analyzeNiftiOutput(self.niftiOutputPath)

    def warningMEDIC(self):
        
        shortwarning = TimerMessageBox(5,self)
        shortwarning.exec_()


    def savePreferences(self):

        prefpath = os.path.join(os.path.expanduser('~'), '.imageimport')

        self.prefs['dicomOutputPath'] = self.dicomOutputPath
        self.prefs['niftiOutputPath'] = self.niftiOutputPath
        self.prefs['currentProtocol'] = self.currentProtocol
        self.prefs['currentSubject'] = self.currentSubject
        self.prefs['sites'] = [str(self.ui.subjectSiteCombo.itemText(item)) for item in range(self.ui.subjectSiteCombo.count())]
        if len(str(self.ui.prefSpmTemplateEdit.text()))>0:
            self.prefs['spm'] = str(self.ui.prefSpmTemplateEdit.text())
        else:
            self.prefs.pop('spm',None)
        if len(str(self.ui.prefANTsEdit.text()))>0:
            self.prefs['ants']=str(self.ui.prefANTsEdit.text())
        else:
            self.prefs.pop('ants',None)
        if len(str(self.ui.prefSpmStandaloneEdit.text()))>0:
            self.prefs['spm_standalone'] = str(self.ui.prefSpmStandaloneEdit.text())
        else:
            self.prefs.pop('spm_standalone',None)
        if len(str(self.ui.prefFreesurferEdit.text()))>0:
            self.prefs['freesurfer'] = str(self.ui.prefFreesurferEdit.text())
        else:
            self.prefs.pop('freesurfer',None)

        if len(str(self.ui.prefNoDBFileLocationEdit.text()))>0:
            self.prefs['FileNoDBPath'] = str(self.ui.prefNoDBFileLocationEdit.text())
        else:
            self.prefs.pop('FileNoDBPath',None)

        if self.ui.prefANTScheckbox.isChecked():
            self.prefs['coregisterMethod'] = 'ANTs'
        elif self.ui.prefSPMcheckbox.isChecked():
            self.prefs['coregisterMethod'] = 'spm'
        elif self.ui.prefAimscheckbox.isChecked():
            self.prefs['coregisterMethod'] = 'Aims'
        else:
            self.prefs['coregisterMethod'] = 'spm'

        self.coregMethod = self.prefs['coregisterMethod']

        projects= [self.ui.radioButtonProject.isChecked(),self.ui.radioButtonProject_2.isChecked(),self.ui.radioButtonProject_3.isChecked(),self.ui.radioButtonProject_4.isChecked()]
        project_selected = [x for x in range(len(projects)) if projects[x]==True]
    
        self.prefs['projectSelected'] = project_selected
    
        fileout = open(prefpath, 'wb')
    
        pickle.dump(self.prefs, fileout)
        #to modify to json
    
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


    def analyzeDicomOutput(self, root):
        """ Structure of the dicom output directory is supposed to be
        root/patient/sequence/dicomfiles.
        This returns the dictionary of patients and sequences already imported in there"""
        if not os.path.isdir(root):
            return {}
        return dict([(f, os.listdir(os.path.join(root, f))) for f in os.listdir(root) if os.path.isdir(os.path.join(root, f))])

    def analyzeNiftiOutput(self, root):
        """ Structure of the nifti output directory is supposed to be
        root/patient/sequence.nii.
        This returns the dictionary of patients and sequences already imported in there"""
        if not os.path.isdir(root):
            return {}
        return dict([(f, [seq.split('.nii')[0] for seq in os.listdir(os.path.join(root, f))]) for f in os.listdir(root) if os.path.isdir(os.path.join(root, f))])

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
        #wdiScan = WriteDiskItem( 'Scanner Based Referential', 'Referential', exactType=True )
        #refScan = wdiScan.findValue( image )
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
        #print "ATTRIBUTES : "+repr(att)
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
        if tab == 3 or tab is None: # Tab 3 is Import Nifti tab
            self.ui.niftiProtocolCombo.clear()
            self.ui.niftiProtocolCombo.addItems(protocols)
            self.ui.niftiProtocolCombo.setCurrentIndex(currentIndex)
            selectSubj(self.ui.niftiSubjectCombo)
        if tab == 7 or tab is None: # Tab 7 is Implantation tab
            self.ui.implProtocolCombo.clear()
            self.ui.implProtocolCombo.addItems(protocols)
            self.ui.implProtocolCombo.setCurrentIndex(currentIndex)
            selectSubj(self.ui.implSubjectCombo)
        if tab == 4 or tab is None: # Tab 4 is Registration tab
            self.ui.regProtocolCombo.clear()
            self.ui.regProtocolCombo.addItems(protocols)
            self.ui.regProtocolCombo.setCurrentIndex(currentIndex)
            selectSubj(self.ui.regSubjectCombo)
        if tab == 2 or tab is None: # Tab 2 is Patient Information tab
            #self.ui.patProtocolCombo.currentIndexChanged.disconnect(self._pathologyinfo)
            self.ui.patProtocolCombo.clear()
            self.ui.patProtocolCombo.addItems(protocols)
            #self.connect(self.ui.patProtocolCombo, QtCore.SIGNAL('currentIndexChanged(QString)'), lambda s:self.selectProtocol(s,self.ui.patSubjectCombo))
            #self._pathologyinfo = pathologypatientinfo.buildUI(self.ui.patpatGroupBox,self.currentProtocol)
            #self.connect(self.ui.patProtocolCombo, QtCore.SIGNAL('currentIndexChanged(QString)'),self._pathologyinfo) #lambda :pathologypatientinfo.buildUI(self.ui.patpatGroupBox,self.currentProtocol))
            self.ui.patProtocolCombo.setCurrentIndex(currentIndex)
            self.pathologypatientInfo = pathologypatientinfo.buildUI(self.ui.patpatGroupBox,self.currentProtocol)
            selectSubj(self.ui.patSubjectCombo)
            #self.pathologypatientInfo = pathologypatientinfo.buildUI(self.ui.patpatGroupBox,self.currentProtocol)
        if tab == 6 or tab is None: # Tab 6 is SEEG tab
            self.ui.seegProtocolCombo.clear()
            self.ui.seegProtocolCombo.addItems(protocols)
            self.ui.seegProtocolCombo.setCurrentIndex (currentIndex)
            selectSubj(self.ui.seegSubjectCombo)


    def selectProtocol(self, proto, comboToFillWithSubjects):
        """ A BrainVisa protocol was selected : query the database to get the list of subjects and put them in the provided QComboBox"""
        self.currentProtocol = str(proto)
        rdi = ReadDiskItem( 'Subject', 'Directory', requiredAttributes={'center':str(proto)} )
        subjects = list( rdi._findValues( {}, None, False ) )
        comboToFillWithSubjects.clear()#self.ui.bvSubjectCombo.clear()
        comboToFillWithSubjects.addItems(sorted([s.attributes()['subject'] for s in subjects]))

        #if self.ui.tabWidget.currentIndex() ==2:
        #   self.pathologypatientInfo = pathologypatientinfo.buildUI(self.ui.patpatGroupBox,self.currentProtocol)

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
        rdi = ReadDiskItem( 'Raw T1 MRI', 'BrainVISA volume formats', requiredAttributes={'center':str(protocol), 'subject':str(subj) } )
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
#         rdi = ReadDiskItem( 'HippoFreesurferAtlas', 'BrainVISA volume formats', requiredAttributes={'center':str(protocol), 'subject':str(subj) } )
#         images += list( rdi._findValues( {}, None, False ) )
        return images

    def selectBvSubject(self, subj):
        """ A BrainVisa subject was selected : query the database to get the list of images"""
        # Change subject
        self.clearAnatomist()
        # Display "Date : XXXX-XX-XX - Seq: T1 - Acq : T1Pre
        self.ui.bvImageList.clear()
        images = self.findAllImagesForSubject(self.ui.bvProtocolCombo.currentText(), subj)
        #name = modality + acquisition + subacquisition if exist subacquisition key else modality + acquisition
        self.ui.bvImageList.addItems(sorted([i.attributes()['modality'] + ' - '+ i.attributes()['acquisition'] + ' - ' + i.attributes()['subacquisition'] if 'subacquisition' in i.attributes().keys()\
                #i.attributes()['modality'] if 'acquisition' not in i.attributes().keys()\
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

    def patientSelectionChanged(self, t):
        """ Update the studies selected when the patient selection changes. t is the list, 'pacs' or 'dir' in the UI"""
    
        print 'patientSelectionChanged('+t+')'
        list = None
        if t == 'pacs':
            list = self.ui.pacsPatientList
            out = self.ui.pacsStudiesList
            studies = self.pacsStudies
        elif t == 'dir':
            list = self.ui.dirPatientList
            out = self.ui.dirStudiesList
            studies = self.dirStudies
    
        # Unselect all studies
        for s in out.selectedItems():
            out.setItemSelected(s, False)
    
        # Select all studies that match the selected patients
        names = [str(p.text()) for p in list.selectedItems()]
        outItems = [out.item(i) for i in xrange(out.count())]
        for s in outItems:
            if studies[self.studiesUIDbyName[str(s.text())]]['PatientsName'] in names:
                out.setItemSelected(s, True)


    def studiesSelectionChanged(self, t):
        """ Update the series selected when the studies selection changes. t is the list, 'pacs' or 'dir' in the UI"""
        list = None
        if t == 'pacs':
            list = self.ui.pacsStudiesList
            out = self.ui.pacsSeriesList
            series = self.pacsStudies
        elif t == 'dir':
            list = self.ui.dirStudiesList
            out = self.ui.dirSeriesList
            series = self.dirSeries
    
        # Unselect all studies
        for s in out.selectedItems():
            out.setItemSelected(s, False)
    
        # Select all series that match the selected studies
        names = [self.studiesUIDbyName[str(s.text())] for s in list.selectedItems()]
        outItems = [out.item(i) for i in xrange(out.count())]
        for s in outItems:
            if series[self.seriesUIDbyName[str(s.text())]]['StudyInstanceUID'] in names:
                out.setItemSelected(s, True)

    def selectOutput(self, t):
        """ Select the output directories for t='Dicom' or t='Nifti' images"""
        if t=='Dicom':
            self.dicomOutputPath = str(QtGui.QFileDialog.getExistingDirectory(self, u"Choix d'un répertoire de sortie pour les fichiers DICOM",self.dicomOutputPath ))
            self.ui.dicomOutputPathEdit.setText(self.dicomOutputPath)
        elif t == 'Nifti':
            self.niftiOutputPath = str(QtGui.QFileDialog.getExistingDirectory(self, u"Choix d'un répertoire de sortie pour les fichiers Nifti",self.niftiOutputPath ))
            self.ui.niftiOutputPathEdit.setText(self.niftiOutputPath)

    def openDicomDir(self, direc=None, updateUI = True):
        """Loads a DICOM directory (direc), organize the files and update the UI"""
        if direc is None:
            direc = QtGui.QFileDialog.getExistingDirectory(self, u"Choix d'un répertoire de fichiers DICOM")#,"/default/Path")
            direc = str(direc)
        if not os.path.exists(direc) or not os.path.isdir(direc):
            print "OpenDicomDir : %s does not exist or is not a directory"%repr(direc)
            return
        files = os.listdir(direc)
        files = [os.path.join(direc,f) for f in files if os.path.isfile(os.path.join(direc,f))]
        self.setStatus('Chargement de '+str(len(files))+" fichiers DICOM...")
        # Create executor object, connect the output function to it, then start it
        organize = lambda :organizeDicomFiles(files, True)
        thread = PythonExecutor(organize)
        self.threads.append(thread)
        self.connect(thread, QtCore.SIGNAL('dicomFilesOutput(PyQt_PyObject)'), lambda x: self.finishOpenDicomDir(updateUI, x))
        thread.start()

    def finishOpenDicomDir(self, updateUI, dicomFiles):
        self.setStatus(u'Fichiers DICOM chargés')
        print "dicomFiles ok"
        if updateUI:
            self.ui.dirSeriesList.clear()
            self.ui.dirStudiesList.clear()
            self.ui.dirPatientList.clear()
    
        series = getSeries(dicomFiles)
        print "*** Found %d series"%len(series)
    
        self.dirSeries = series
        #Afficher la liste des series
        if updateUI:
            self.seriesUIDbyName = dict((series[uid]['SeriesDate']+' '+series[uid]['SeriesTime']+' '+series[uid]['SeriesDescription']  , uid) for uid in series if 'SeriesTime' in series[uid] and 'SeriesDate' in series[uid] and 'SeriesDescription' in series[uid])
            self.ui.dirSeriesList.addItems(self.seriesUIDbyName.keys())
            self.ui.dirSeriesList.sortItems()
            # Colorer les séries déjà importées:
            try:
                for i in xrange(self.ui.dirSeriesList.count()):
                    self.ui.dirSeriesList.item(i).setBackgroundColor(QColor.fromRgb(255,0,0))
                    s = self.seriesUIDbyName[str(self.ui.dirSeriesList.item(i).text())]
                    if  patientNameDecode(s['PatientsName']) in niftiOutContent:
                        if serie2filename(s) in niftiOutContent[patientNameDecode(s['PatientsName'])]:
                            self.ui.dirSeriesList.item(i).setBackgroundColor(QColor.fromRgb(0,255,0))
            except:
                pass
    
        self.dirStudies = getStudies(dicomFiles)
        studies = self.dirStudies
        print "Studies : "+repr(studies)
        #Afficher la liste des studies
        if updateUI:
            self.studiesUIDbyName = dict((studies[uid]['name'], uid) for uid in studies.keys())
            self.ui.dirStudiesList.addItems(self.studiesUIDbyName.keys())
            self.ui.dirStudiesList.sortItems()
    
        self.dirPatients = getPatients(series)
        print "Patients : "+repr(studies)
        if updateUI:
            self.ui.dirPatientList.addItems(self.dirPatients)
            self.ui.dirPatientList.sortItems()


    ## DCMTK retrieve example :
    #	# get list of studies for patient "John Doe"
    #findscu --study -k 0008,0052=STUDY -k 0010,0010="Doe^John" -k 0020,000D server 104
    #now get list of CT series within a certain study identified by a UID returned by the last call to findscu
    #findscu --study -k 0008,0052=SERIES -k 0020,000D=1.2.276.1.2.3.4.5 -k 0008,0060=CT -k 0020,000E server 104
    #retrieve a series found by findscu. Note that we need to specify both Study and Series Instance UID.
    #movescu --study -k 0008,0052=SERIES -k 0020,000D=1.2.276.1.2.3.4.5 -k 0020,000E=1.2.276.9.8.7.6.5 server 104
    def searchPacs(self):
        """ Search the PACS for the patient Name provided"""
        # Get the list of patients'names if the provided name is llong enough (2 chars minimum)
        name = str(self.ui.patientSearchName.text())+'*'
        name.rstrip()#replace(' ', '^') --> spaces are allowed inside names, ^separates first/last/middle names
        # Remove wildcards to avoid large matches
        name.replace('?', '')
        name.replace('*', '')
        server = str(self.ui.pacsServerEdit.text())
        port = str(self.ui.pacsPortSpin.value())
        # We want study mode, with tag 0008,0052 set to study, tag 0010,0010 matching the patient, and listing all 0020,000 (studies)
        cmd = ['findscu','--study', '-k','0008,0052=STUDY', '-k' , '0010,0010="'+name+'"', '-k','0020,000D', server, port ]
    
        lines = runCmd(cmd)
        print lines
        # List of studies
        studies = []
        # In all studies get  the patient's name and add it to the list
        patientsList = []
    
        # Maximum number of patients : 10
        patientsList =patientsList[:10]
    
        # For each study that matches a patient from the list, get all series
        for s in studies:
            cmd = ['findscu','--study', '-k','0008,0052=SERIES', '-k' , '0010,0010="'+name+'"',\
                   '-k', '0020,000D='+s, '-k', '0020,000E', server, port ]
            lines = runCmd(cmd)
            # Store the results


    def importDataPacs(self):
        """ Get the selected series from the pacs to a local directory and convert to Nifti """
        # Regarder la sélection patient et series
        # importer avec findscu dans un répertoire temporaire
        # getscu ou movescu
        cmd = ['getscu', '--study', '-k', '0008,0052=SERIES', '-k', '0020,000D='+studyId, '-k', '0020,000E='+serieId, server, port]
        # puis lancer self.openDicomDir(tmpDir)


    def importDataDir(self):
        # Find selected series
        selectedSeries = dict((self.seriesUIDbyName [str(item.text())], self.dirSeries[self.seriesUIDbyName [str(item.text())]]) for item in self.ui.dirSeriesList.selectedItems())
    
        # If nothing is selected
        if not selectedSeries:
            return
    
        # Copy the dicom series to a directory
        copySeriesTo(selectedSeries, self.dicomOutputPath)
    
        # Convert to Nifti
        th = convertSeriesToNifti(selectedSeries, self.niftiOutputPath)
        self.threads.extend(th)
    
        print "Series importees : " + repr(selectedSeries.keys())
        # Regarder la sélection patient et series
        # importer avec findscu dans un répertoire temporaire
        # puis lancer copySeriesTo(series, dicomOutputPath)

    def chooseNifti(self):
        self.ui.niftiFileLabel.setText('')
        path = str(QtGui.QFileDialog.getOpenFileName(self, "Open MRI Image", "", "Images (*.nii *.nii.gz *.ima *.img *.ima.gz *.img.gz *.img.bz2 *ima.bz2 *.nii.bz2 *.mgz)"))
        if not path:
            return
        self.ui.niftiFileLabel.setText(path)
        # Display the image and remove others
        print "Loading %s"%path
        split_path = os.path.splitext(path)
        if split_path[-1] == ".mgz":
            try:
                # tmp_nii_path = "/ tmp / tmp_mgz_nii.nii"
                tmp_nii_path = getTmpFilePath('nii')
                #on reslice par rapport au t1 pre
                di = ReadDiskItem( 'Raw T1 MRI', 'BrainVISA volume formats', requiredAttributes={'center':self.currentProtocol, 'subject':self.currentSubject } )
                allT1 = list(di.findValues({},None,False))
                idxT1pre = [i for i in range(len(allT1)) if 'T1pre' in str(allT1[i])]
                if len(idxT1pre)==0:
                    print "to import freesurfer images, the T1pre has to be imported because a reslicing is necessary"
                    return
    
                #mri_convert -rl /data/brainvisa/Epilepsy/Gre_2015_BESk/t1mri/T1pre_2015-1-1/Gre_2015_BESk.nii /home/b67-belledone/freesurfer/subjects/Gre_2015_BESk/mri/aparc.a2009s+aseg.mgz /home/b67-belledone/freesurfer/subjects/Gre_2015_BESk/mri/test.nii -rt nearest -ncD
                #remplacer par launchFreesurferCommand
                self.brainvisaContext.write("mri_convert from freesurfer")
    
                launchFreesurferCommand(context, None, 'mri_convert', '-i',path,'-o',tmp_nii_path,'-rl',str(allT1[idxT1pre[0]].fullPath()),'-rt','nearest','-nc')
                #mriconvert_call ='mri_convert -i {} -o {} -rl {} -rt nearest -nc'.format(path,tmp_nii_path,str(allT1[idxT1pre[0]].fullPath()))
                #mriconvert_call ='mri_convert -i {} -o {}'.format(path,tmp_nii_path)
                #runCmd(mriconvert_call.split())
            except:
                print "conversion mgz to nii didn't work"
                return
        else:
            mri = self.a.loadObject(path)
            tmp_nii_path = path
        # Display images
        self.displayImage(tmp_nii_path, self.wins)
        # Reset other components
        self.brainCenter = None
        self.ui.brainCenterLabel.setText('')
        self.ui.niftiFileLabel.setText(tmp_nii_path)

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
    
        #filetype = {'T1':'Raw T1 MRI', 'T2':'T2 MRI', 'CT':'CT', 'TEP':'PET', 'PET':'PET', 'fMRI':'fMRI-epile', 'Statistics':'Statistic-Data', 'FLAIR':'FLAIR', 'FreesurferAtlas':'FreesurferAtlas', 'FGATIR':'FGATIR','HippoFreesurferAtlas':'HippoFreesurferAtlas'}[modality]
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
#         if filetype == 'FreesurferAtlas':
#             self.generateAmygdalaHippoMesh(proto, patient, acq, di)
            
        # self.importNiftiWorker(di, path, filetype, proto, patient)
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
        # Check that all the files exist
        for key in allFiles:
            if not os.path.isfile(allFiles[key]['file']):
                QtGui.QMessageBox.warning(self, "Error", u"FreeSurfer file not found:\n" + allFiles[key]['file'])
                return
        
        # Run copy and conversion in a separate thread
        res = ProgressDialog.call(lambda thr:self.importFSoutputWorker(FsSubjDir, subject, allFiles, diT1pre, thr), True, self, "Processing...", "Import FreeSurfer output")
        # res = self.importFSoutputWorker(FsSubjDir, subject, allFiles, diT1pre)
        

    def importFSoutputWorker(self, FsSubjDir, subject, allFiles, diT1pre, thread=None):
        # Where to copy the new files
        acq =  str(diT1pre.attributes()['acquisition']).replace('T1','FreesurferAtlas')
        write_filters = { 'center': str(self.ui.niftiProtocolCombo.currentText()), 'acquisition': str(acq), 'subject' : subject }
        # Progress bar
        if thread:
            thread.emit(QtCore.SIGNAL("PROGRESS_TEXT"), "Copying files to freesurfer database...")
        # Add all the files to the local freesurfer database
        for key in allFiles:
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

#         # Get output folder
#         wdi = WriteDiskItem('Acquisition', 'Directory')
#         diFolder = wdi.findValue(write_filters)
#         # Importing Destrieux atlas to BrainVISA database
#         wdi = WriteDiskItem('Label volume', 'NIFTI-1 image', requiredAttributes={'subject':subject, '_ontology':diFolder.attributes()["_ontology"], '_database':diFolder.attributes()["_database"]})
#         diFSDestrieux = wdi.findValue(diFolder.fullPath() + '/aparc_a2009s+aseg_' + subject + '.nii')
        
        # Create the folder if doesn't exist
        self.createItemDirs(diFSDestrieux)
        # Reslice and convert to AIMS
        launchFreesurferCommand(self.brainvisaContext, None, 'mri_convert', '-i',str(allFiles['destrieux']['file']),'-o',str(diFSDestrieux.fullPath()),'-rl',str(diT1pre.fullPath()),'-rt','nearest','-nc')
        ret = subprocess.call(['AimsFileConvert', '-i', str(diFSDestrieux.fullPath()), '-o', str(diFSDestrieux.fullPath()), '-t', 'S16'])
        # Add reference in the database (creates .minf)
        neuroHierarchy.databases.insertDiskItem(diFSDestrieux, update=True)
        # Generate amygdala and hippocampus meshes
        if thread:
            thread.emit(QtCore.SIGNAL("PROGRESS_TEXT"), "Creating hippcampus and amygdala meshes...")
        self.generateAmygdalaHippoMesh(str(self.ui.niftiProtocolCombo.currentText()), subject, acq, diFSDestrieux)
 
        
#     def importFSoutputOld(self,subject=None):
# 
#         if subject == None:
#             subject = self.currentSubject
# 
#         #find T1pre of the subject
#         rT1BV = ReadDiskItem('Raw T1 MRI', 'BrainVISA volume formats',requiredAttributes={'subject':subject})
#         allT1 = list(rT1BV.findValues({},None,False))
#         idxT1pre = [i for i in range(len(allT1)) if 'T1pre' in str(allT1[i])]
#         diT1pre = allT1[idxT1pre[0]]
# 
#         acq =  str(diT1pre.attributes()['acquisition']).replace('T1','FreesurferAtlas')
#         write_filters = { 'center': str(self.ui.niftiProtocolCombo.currentText()), 'acquisition': str(acq), 'subject' : subject }
# 
#         FreesurferNu = ReadDiskItem('Nu FreesurferAnat', 'FreesurferMGZ',requiredAttributes={'subject':subject,'_ontology':'freesurfer'})
#         Nufound = list(FreesurferNu.findValues({},None,False))
#         if not Nufound:
#             QtGui.QMessageBox.warning(self, "Error", u"Subject \"" + subject + "\" was not found in the local FreeSurfer database.")
#             return
#         Nufound = [x for x in Nufound if '.mgz' in str(x)][0]
# 
#         FreesurferRibbon = ReadDiskItem('Ribbon Freesurfer','FreesurferMGZ',requiredAttributes={'subject':subject,'_ontology':'freesurfer'})
#         Ribbonfound = list(FreesurferRibbon.findValues({},None,False))
#         Ribbonfound = [x for x in Ribbonfound if '.mgz' in str(x)][0]
# 
#         FreesurferPial = ReadDiskItem('Any Type','FreesurferPial',requiredAttributes={'subject':subject,'_ontology':'freesurfer'})
# 
#         Pialfound = list(FreesurferPial.findValues({},None,False))
#         PialLeft = [x for x in Pialfound if '/lh.' in str(x)][0]
#         PialRight = [x for x in Pialfound if '/rh.' in str(x)][0]
# 
#         PialMeshBV = WriteDiskItem('Pial','aims mesh formats')
#         diPialMeshBV = PialMeshBV.findValue(write_filters)
# 
#         FreesurferWhite = ReadDiskItem('Any Type','FreesurferWhite',requiredAttributes={'subject':subject,'_ontology':'freesurfer'})
# 
#         Whitefound = list(FreesurferWhite.findValues({},None,False))
#         WhiteLeft = [x for x in Whitefound if '/lh.' in str(x)][0]
#         WhitelRight = [x for x in Whitefound if '/rh.' in str(x)][0]
# 
#         right_grey_white_output=WriteDiskItem('Right Grey White Mask','Aims writable volume formats')
#         diright_grey_white = right_grey_white_output.findValue(write_filters)
#         left_grey_white_output=WriteDiskItem('Left Grey White Mask','Aims writable volume formats')
#         dileft_grey_white = left_grey_white_output.findValue(write_filters)
# 
#         WhiteMeshBV = WriteDiskItem('White','aims mesh formats')
#         diWhiteMeshBV = WhiteMeshBV.findValue(write_filters)
#         #diWhiteLeftMeshBV = [x for x in diWhiteMeshBV
# 
#         FreesurferThickness = ReadDiskItem('Any Type','FreesurferThickness',requiredAttributes={'subject':subject,'_ontology':'freesurfer'})
# 
#         Thicknessfound = list(FreesurferThickness.findValues({},None,False))
#         ThickLeft = [x for x in Thicknessfound if '/lh.' in str(x)][0]
#         ThickRight = [x for x in Thicknessfound if '/rh.' in str(x)][0]
# 
#         FreesurferDestrieux = ReadDiskItem('Freesurfer Cortical Parcellation using Destrieux Atlas','FreesurferMGZ',requiredAttributes={'subject':subject,'_ontology':'freesurfer'})
# 
#         Destrieuxfound = list(FreesurferDestrieux.findValues({},None,False))
# 
#         #generer tous les writediskitem
#         self.brainvisaContext.write("convert output from freesurfer")
#         launchFreesurferCommand(context, None, 'mri_convert', '-i',str(Nufound),'-o',os.path.join(os.path.dirname(str(Nufound)),'nu.nii'),'-rl',str(diT1pre.fullPath()))
#         launchFreesurferCommand(context, None, 'mri_convert', '-i',str(Ribbonfound),'-o',os.path.join(os.path.dirname(str(Ribbonfound)),'ribbon.nii'),'-rl',str(diT1pre.fullPath()),'-rt','nearest','-nc')
# 
# 
#         wdi = WriteDiskItem('FreesurferAtlas', 'NIFTI-1 image' )
#         diFSDestrieux = wdi.findValue(write_filters)
#         #create the folder if doesn't exist
#         if not os.path.isfile(os.path.dirname(str(diFSDestrieux.fullPath()))):
#             try:
#                 os.makedirs(os.sep.join(os.path.dirname(str(diFSDestrieux.fullPath())).split(os.sep)[:-1]))
#             except:
#                 #already exist probably
#                 pass
#             try:
#                 os.makedirs(os.path.dirname(str(diFSDestrieux.fullPath())))
#             except:
#                 #already exist probably
#                 pass
# 
#         launchFreesurferCommand(context, None, 'mri_convert', '-i',str(Destrieuxfound[0]),'-o',str(diFSDestrieux.fullPath()),'-rl',str(diT1pre.fullPath()),'-rt','nearest','-nc')
#         ret = subprocess.call(['AimsFileConvert', '-i', str(diFSDestrieux.fullPath()), '-o', str(diFSDestrieux.fullPath()), '-t', 'S16'])
#         #pour destrieux
#         self.generateAmygdalaHippoMesh(str(self.ui.niftiProtocolCombo.currentText()), subject, acq, diFSDestrieux)
# 
#         launchFreesurferCommand(context, None, 'mris_convert', '--to-scanner', WhiteLeft, str(diWhiteleft.fullPath()))
#         launchFreesurferCommand(context, None, 'mris_convert', '--to-scanner', WhiteRight, str(diWhiteright.fullPath()))
#         launchFreesurferCommand(context, None, 'mris_convert', '--to-scanner', PialLeft, str(diPialleft.fullPath()))
#         launchFreesurferCommand(context, None, 'mris_convert', '--to-scanner', PialRight, str(diPialright.fullPath()))
#         launchFreesurferCommand(context, None, 'mris_convert', '--to-scanner', ThickLeft, str(diPialright.fullPath()))
#         launchFreesurferCommand(context, None, 'mris_convert', '--to-scanner', ThickRight, str(diPialright.fullPath()))
# 
#         #loading the object to get its referential and assign it
#         self.a.loadObject(diT1pre)
#         refes=self.a.getReferentials()
#         for element in refes:
#             if 'T1pre' in str(element.getInfos()):
#                 ref=element
#         whiteL=self.a.loadObject(diWhiteleft)
#         whiteL.assignReferential(ref)
#         whiteR=self.a.loadObject(diWhiteright)
#         whiteR.assignReferential(ref)
#         pialL=self.a.loadObject(diPialleft)
#         pialL.assignReferential(ref)
#         pialR=self.a.loadObject(diPialright)
#         pialR.assignReferential(ref)


    def importPacs(self):
        """ Import the selected images from the pacs to BrainVisa (after the Nifti importation) TODO does not work """
    
        #self.importNifti(path, patient, modality, acq)
        selectedSeries = dict((self.seriesUIDbyName [str(item.text())], self.pacsSeries[self.seriesUIDbyName [str(item.text())]]) for item in self.ui.pacsSeriesList.selectedItems())
        self.importBV(selectedSeries)


    def importDir(self):
        """ Import the selected image from dicom directory to BrainVisa (after Nifti importation) """
        # Find selected series
        selectedSeries = dict((self.seriesUIDbyName [str(item.text())], self.dirSeries[self.seriesUIDbyName [str(item.text())]]) for item in self.ui.dirSeriesList.selectedItems())
        self.importBV(selectedSeries)

    def importBV(self, selectedSeries):
        if len(selectedSeries) != 1:
            QtGui.QMessageBox.warning(self, "Error", u"Select one (and only one ! ) sequence !")
            return
        serie = None
        for k in selectedSeries:
            serie = selectedSeries[k]
    
        if 'niftiPath' not in serie:
            QtGui.QMessageBox.warning(self, "Error", u"You have to convert the sequence in Nifti format before importing in Brainvisa")
            return
    
        acq = str(self.ui.dirSeqType.currentText() + self.ui.dirAcqType.currentText() + '_' + dateDecode(serie['SeriesDate']))
        patient = patientNameDecode(serie['PatientsName'])
        self.importNifti(serie['niftiPath'], patient, str(self.ui.dirSeqType.currentText()), acq)
    

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



    # ******************************* Implantation

    def selectImplSubject(self, subj):
        """ A BrainVisa subject was selected : query the database to get the available images"""
        rdi = ReadDiskItem( 'Electrode Implantation Coronal Image', 'BrainVISA image formats', requiredAttributes={'center':str(self.ui.implProtocolCombo.currentText()), 'subject':str(subj)} )
        coros = list( rdi._findValues( {}, None, False ) )
        coro = None
        if len(coros) == 0:
            self.ui.implCoroLabel.setText("no image")
        elif len(coros) > 1:
            QtGui.QMessageBox.warning(self, "Erreur", "Ce sujet a plusieurs implantations enregistrees ! Cette situation n'est pas geree par ce logiciel")
            print repr(coros)
        else:
            self.ui.implCoroLabel.setText(str(coros[0].fullPath()))
            coro = coros[0].fullPath()
    
        rdi = ReadDiskItem( 'Electrode Implantation Sagittal Image', 'BrainVISA image formats', requiredAttributes={'center':str(self.ui.implProtocolCombo.currentText()), 'subject':str(subj)} )
        sags = list( rdi._findValues( {}, None, False ) )
        sag = None
        if len(sags) == 0:
            self.ui.implSagLabel.setText("no image")
        elif len(sags) > 1:
            QtGui.QMessageBox.warning(self, "Erreur", "Ce sujet a plusieurs implantations enregistrees ! Cette situation n'est pas geree par ce logiciel")
            print repr(sags)
        else:
            self.ui.implSagLabel.setText(str(sags[0].fullPath()))
            sag = sags[0].fullPath()
        self.displayAnatomist(coro, sag)

    def openImplCoro(self):
        path = str(QtGui.QFileDialog.getOpenFileName(self, "Open Coronal 2D Image", "", "Images (*.gif *.jpg *.jpeg *.png)"))
        if not path:
            return
        self.ui.implCoroLabel.setText(path)
        sag = str(self.ui.implSagLabel.text())
        if sag == "":
            sag = None
        self.displayAnatomist(path, sag)

    def importImplCoro(self):
        wdi = WriteDiskItem('Electrode Implantation Coronal Image', 'JPEG image' )#'gz compressed NIFTI-1 image' )
        di = wdi.findValue( { 'center': str(self.ui.implProtocolCombo.currentText()), 'subject' : str(self.ui.implSubjectCombo.currentText()) } )
        if di is None:
            QtGui.QMessageBox.warning(self, "Erreur", "Impossible de trouver un chemin valide pour importer l'image coro dans BrainVisa")
            return
        # Copy the file : lancer l'importation standard de T1 MRI pour la conversion de format et
        destination = di.fileName()
        print "Importing file as "+destination
        # Create directories that do not exist yet
        self.createItemDirs(di)
        #
        ret = subprocess.call(['AimsFileConvert', '-i', str(self.ui.implCoroLabel.text()), '-o', str(destination)])
        if ret < 0:
            print "Erreur d'importation"
            QtGui.QMessageBox.warning(self, "Erreur", "Erreur d'importation BrainVisa / AimsFileConvert")
            return
        neuroHierarchy.databases.insertDiskItem( di, update=True )
        self.setStatus(u"Rapport d'implantation coronal importé")


    def openImplSag(self):
        path = str(QtGui.QFileDialog.getOpenFileName(self, "Open Sagittal 2D Image", "", "Images (*.gif *.jpg *.jpeg *.png)"))
        if not path:
            return
        self.ui.implSagLabel.setText(path)
        coro = str(self.ui.implCoroLabel.text())
        if coro == "":
            coro = None
        self.displayAnatomist(coro, path)


    def importImplSag(self):
        wdi = WriteDiskItem('Electrode Implantation Sagittal Image', 'JPEG image' )#'gz compressed NIFTI-1 image' )
        di = wdi.findValue( { 'center': str(self.ui.implProtocolCombo.currentText()), 'subject' : str(self.ui.implSubjectCombo.currentText()) } )
        if di is None:
            QtGui.QMessageBox.warning(self, "Erreur", "Impossible de trouver un chemin valide pour importer l'image sagittale dans BrainVisa")
            return
        # Copy the file : lancer l'importation standard de T1 MRI pour la conversion de format et
        destination = di.fileName()
        print "Importing file as "+destination
        # Create directories that do not exist yet
        self.createItemDirs(di)
        #
        ret = subprocess.call(['AimsFileConvert', '-i', str(self.ui.implCoroLabel.text()), '-o', str(destination)])
        if ret < 0:
            print "Erreur d'importation"
            QtGui.QMessageBox.warning(self, "Erreur", "Erreur d'importation BrainVisa / AimsFileConvert")
            return
        neuroHierarchy.databases.insertDiskItem( di, update=True )
        self.setStatus(u"Rapport d'implantation Sagittal importé")

    def openImplPpt(self):
        path = str(QtGui.QFileDialog.getOpenFileName(self, "Open PPT implantation schema", "", "Powerpoint file (*.ppt *.pptx)"))
        if not path:
            return
        self.ui.implPptLabel.setText(path)


    def importImplPpt(self):
        wdi = WriteDiskItem('Electrode Implantation Powerpoint report', 'Powerpoint file' )
        di = wdi.findValue( { 'center': str(self.ui.implProtocolCombo.currentText()), 'subject' : str(self.ui.implSubjectCombo.currentText()) } )
        if di is None:
            QtGui.QMessageBox.warning(self, "Erreur", "Impossible de trouver un chemin valide pour importer le rapport Powerpoint dans BrainVisa")
            return
        # Copy the file
        destination = di.fileName()
        print "Importing file as "+destination
        # Create directories that do not exist yet
        self.createItemDirs(di)
        shutil.copyfile(str(self.ui.implPptLabel.text()), str(destination))
        if ret < 0:
            print "Erreur d'importation"
            QtGui.QMessageBox.warning(self, "Erreur", "Erreur lors de la copie du fichier")
            return
        neuroHierarchy.databases.insertDiskItem( di, update=True )
        self.setStatus(u"Rapport d'implantation Powerpoint importé")

    def openImplPdf(self):
        path = str(QtGui.QFileDialog.getOpenFileName(self, "Open PDF implantation schema", "", "PDF file (*.pdf)"))
        print "Chemin du PDF %s"%path
        if not path:
            return
        self.ui.implPdfLabel.setText(path)

    def openImplPdf2(self):
        path = str(QtGui.QFileDialog.getOpenFileName(self, "Open PDF implantation schema", "", "PDF file (*.pdf)"))
        print "Chemin du PDF %s"%path
        if not path:
            return
        self.ui.implPdfLabel2.setText(path)

    def openImplElListPdf(self):
        path = str(QtGui.QFileDialog.getOpenFileName(self, "Open PDF Electrode List", "", "PDF file (*.pdf)"))
        print "Chemin du PDF %s"%path
        if not path:
            return
        self.ui.implElListPdfLabel.setText(path)


    def importImplPdf(self):
        wdi = WriteDiskItem('Electrode Implantation PDF report', 'PDF file' )
        di = wdi.findValue( { 'center': str(self.ui.implProtocolCombo.currentText()), 'subject' : str(self.ui.implSubjectCombo.currentText()), 'planning' : 'True' } )
        if di is None:
            QtGui.QMessageBox.warning(self, "Erreur", "Impossible de trouver un chemin valide pour importer le rapport PDF dans BrainVisa")
            return
        # Copy the file
        destination = di.fileName()
        print "Importing file as "+destination
        # Create directories that do not exist yet
        self.createItemDirs(di)
        try:
            shutil.copyfile(str(self.ui.implPdfLabel.text()), str(destination))
        except:
            print "Erreur d'importation"
            QtGui.QMessageBox.warning(self, "Erreur", "Erreur lors de la copie du fichier PDF")
            return
        neuroHierarchy.databases.insertDiskItem( di, update=True )
        self.setStatus(u"Rapport d'implantation PDF importé")

    def importImplPdf2(self):
        wdi = WriteDiskItem('Electrode Implantation PDF report', 'PDF file' )
        di = wdi.findValue( { 'center': str(self.ui.implProtocolCombo.currentText()), 'subject' : str(self.ui.implSubjectCombo.currentText()), 'planning' : 'False' } )
        if di is None:
            QtGui.QMessageBox.warning(self, "Erreur", "Impossible de trouver un chemin valide pour importer le rapport PDF dans BrainVisa")
            return
        # Copy the file
        destination = di.fileName()
        print "Importing file as "+destination
        # Create directories that do not exist yet
        self.createItemDirs(di)
        try:
            shutil.copyfile(str(self.ui.implPdfLabel2.text()), str(destination))
        except:
            print "Erreur d'importation"
            QtGui.QMessageBox.warning(self, "Erreur", "Erreur lors de la copie du fichier PDF")
            return
        neuroHierarchy.databases.insertDiskItem( di, update=True )
        self.setStatus(u"Rapport d'implantation PDF importé")

    def importImplElListPdf(self):
        wdi = WriteDiskItem('Electrode List PDF', 'PDF file' )
        di = wdi.findValue( { 'center': str(self.ui.implProtocolCombo.currentText()), 'subject' : str(self.ui.implSubjectCombo.currentText()), 'planning' : 'False' } )
        if di is None:
            QtGui.QMessageBox.warning(self, "Erreur", "Impossible de trouver un chemin valide pour importer la liste PDF dans BrainVisa")
            return
        # Copy the file
        destination = di.fileName()
        print "Importing file as "+destination
        # Create directories that do not exist yet
        self.createItemDirs(di)
        try:
            shutil.copyfile(str(self.ui.implElListPdfLabel.text()), str(destination))
        except:
            print "Erreur d'importation"
            QtGui.QMessageBox.warning(self, "Erreur", "Erreur lors de la copie du fichier PDF")
            return
        neuroHierarchy.databases.insertDiskItem( di, update=True )
        self.setStatus(u"Rapport d'implantation PDF importé")

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
            if acq.startswith('T1'):
                self.setStatus(u"SPM normalization %s..."%acq)
                progressThread.emit(QtCore.SIGNAL("PROGRESS_TEXT"), "SPM normalization: " + subj + "/" + image.attributes()['modality'] + "...")
                self.spmNormalize(image.fileName(), proto, patient, acq)
                self.taskfinished(u"SPM normalization done")
                # If there is a T1pre, remember the image
                if acq.find('T1pre') == 0:
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
                progressThread.emit(QtCore.SIGNAL("PROGRESS_TEXT"), "ANTS coregistration: " + subj + "/" + image.attributes()['modality'] + "...")
                temp_folder_ants = tempfile.mkdtemp('ANTs_IntrAnat') +'/'
                ants_call = 'antsRegistrationSyN.sh -d 3 -f {} -m {} -o {} -t r'.format(str(t1preImage.fullPath()),str(image.fullPath()),temp_folder_ants)
#                 thr = Executor(ants_call.split())
#                 thr.finished.connect(lambda:self.taskfinished(u"ANTs Coregister done", thr))
#                 thr.finished.connect(lambda im = image, tmp_folder=temp_folder_ants:self.setANTstrm_database(im, tmp_folder))
#                 self.threads.append(thr)
#                 thr.start()
                print("ANTs call: " + ants_call)
                # Run ANTs system command
                runCmd(ants_call.split())
                # Register transformation in the database
                self.setANTstrm_database(image, temp_folder_ants)
                
            # ===== SPM =====
            elif self.coregMethod == 'spm':
                print("Coregistration method: SPM")
                # SPM coregister
                if self.ui.regResampleCheck.isChecked():
                    progressThread.emit(QtCore.SIGNAL("PROGRESS_TEXT"), "SPM coregistration+resample: " + subj + "/" + image.attributes()['modality'] + "...")
                    # Call SPM coregister+resample
                    call = spm_coregisterReslice%("'"+str(self.prefs['spm'])+"'","{'"+str(t1preImage)+",1'}", "{'"+str(image.fileName())+",1'}")
                    spl = os.path.split(image.fileName())
                    registeredPath = os.path.join(spl[0], 'r'+spl[1])
                    matlabRun(call)
                    # Update resampled volume
                    self.setResampledToT1pre(image, registeredPath)
                    
                # SPM coregister
                else:
                    progressThread.emit(QtCore.SIGNAL("PROGRESS_TEXT"), "SPM coregistration: " + subj + "/" + image.attributes()['modality'] + "...")
                    self.setStatus('Coregistration of the image: '+ acq)
                    #self.spmCoregister(image, t1preImage)
                    
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
                    matlabRun(call)
                    # Register transformation in the database
                    self.insertTransformationToT1pre(tmpOutput, image)
                    if ('data_type' in image.attributes().keys()) and (image.attributes()['data_type'] == 'RGB'):
                        os.remove(imageFileName)
                        os.remove(imageFileName+'.minf')

    
            # ===== AIMS =====
            elif self.coregMethod == 'Aims':
                raise Exception("AIMS coregistration not supported yet.")
                
                #ret = subprocess.call(['AimsMIRegister', '-r', str(t1preImage.fullPath()), '-t', str(image.fullPath()), '--dir', tmp_trm_path, '--inv',tmp_trm_path2])
                #if ret < 0:
                #   print "coregistration error: "+ str(image.fullPath())#terminal
                #   QtGui.QMessageBox.warning(self, "Error", u"The coregistration didn't work") #utilisateur
                #   returnFreesurferAtlaspre
                #self.insertTransformationToT1pre(tmp_trm_path,image)
        self.taskfinished(u"Coregistration done")
        # Clear all the views
        self.clearAnatomist()


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
                delPath = os.path.dirname(os.path.abspath(hemis[0].fullPath()))
                removeFromDB(delPath, neuroHierarchy.databases.database(hemis[0].get("_database")))
            else:
                return
    
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
                # Get reference to current database
                db = neuroHierarchy.databases.database(hemis[idxT1pre[0]].get("_database"))
                # Delete folder default_analysis/segmentation/mesh/surface_analysis                
                delPath = os.path.dirname(os.path.abspath(hemis[idxT1pre[0]].fullPath())) + '/surface_analysis'
                removeFromDB(delPath, db)
                # Get pial meshes
                delFiles = hemis
                # Get white meshes
                rdi = ReadDiskItem('Hemisphere White Mesh', 'Anatomist mesh formats', requiredAttributes={'subject':subject, 'center':protocol})
                delFiles += list(rdi._findValues({}, None, False))
                # Get head meshe
                rdi = ReadDiskItem('Head Mesh', 'Anatomist mesh formats', requiredAttributes={'subject':subject, 'center':protocol})
                delFiles += list(rdi._findValues({}, None, False))
                # Delete all files
                for h in delFiles:
                    # Delete only the ones in T1pri
                    if 'T1pre' in h.attributes()["acquisition"]:
                        removeFromDB(h.fullPath(), db)
                        # Delete associated .minf
                        if os.path.isfile(h.fullPath() + '.minf'):
                            os.remove(h.fullPath() + '.minf')
            else:
                return

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
        
        # ===== STEP 3: MARS ATLAS =====
        # Start hip-hop
        if thread:
            thread.emit(QtCore.SIGNAL("PROGRESS_TEXT"), "Running Hip-Hop...")
            thread.emit(QtCore.SIGNAL("PROGRESS"), 50)
        self.hiphopStart(diOut.attributes()['center'], diOut.attributes()['subject'])
        

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
        #info_mat = scipy.io.loadmat(file_to_open)setANTstrm_database
        #matrix_unshaped = info_mat['AffineTransform_double_3_3']
        #matrix_shaped = matrix_unshaped.reshape(4,3)
        #matrix_shaped = numpy.roll(matrix_shaped,3)
    
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
        ProgressDialog.call(lambda thr:self.validateAcPcWorker(thr), True, self, "Processing...", "BrainVISA segmentation")
        
#         # No gado
#         if (not 'Gado' in self.mriAcPc.attributes().keys()) or (not self.mriAcPc.attributes()['Gado']):
#             self.brainvisaContext.runInteractiveProcess(lambda x='':self.hiphopStart() , 'Morphologist 2015', t1mri = self.mriAcPc, perform_normalization = False, anterior_commissure = self.AcPc['AC'],\
#                     posterior_commissure = self.AcPc['PC'], interhemispheric_point = self.AcPc['IH'], left_hemisphere_point = self.AcPc['LH'], perform_sulci_recognition = True)
#         # Gado
#         else:
#             # Morphologist step #1: Prepare subject
#             processPrepare = getProcessInstance('preparesubject')
#             processPrepare.T1mri = self.mriAcPc
#             processPrepare.Normalised = "No"
#             processPrepare.Anterior_Commissure = self.AcPc['AC']
#             processPrepare.Posterior_Commissure = self.AcPc['PC']
#             processPrepare.Interhemispheric_Point = self.AcPc['IH']
#             processPrepare.Left_Hemisphere_Point = self.AcPc['LH']
#             processPrepare.allow_flip_initial_MRI = True
#             self.brainvisaContext.runInteractiveProcess(lambda x='':self.morphologistGado1(), processPrepare)
        
        
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
            call = matlab_removeGado%("'"+str(self.prefs['spm'])+"'","'"+nobiasPre+",1'","'"+str(pathTPMseg)+",1'","'"+str(pathTPMseg)+",2'","'"+str(pathTPMseg)+",3'","'"+str(pathTPMseg)+",4'","'"+str(pathTPMseg)+",5'","'"+str(pathTPMseg)+",6'",\
                   "'"+c1Name+"'","'"+c2Name+"'","'"+c3Name+"'","'"+c4Name+"'","'"+nobiasPre+"'","'"+nogadoPre+"'")
            matlabRun(call)
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
        self.hiphopStart(self.mriAcPc.attributes()['center'], self.mriAcPc.attributes()['subject'])
        
            
            
#     def morphologistGado1(self):
#         # Morphologist step #2: Bias correction
#         processBias = getProcessInstance('t1biascorrection')
#         processBias.t1mri = self.mriAcPc
#         self.brainvisaContext.runInteractiveProcess(lambda x='':self.morphologistGado2(), processBias)

#     def morphologistGado2(self):
#         # Remove GADO with SPM
#         print self.currentSubject + ": Running segmentation to remove Gado on T1 no bias..."
#         nobiasRDI = ReadDiskItem("T1 MRI Bias Corrected", 'BrainVISA volume formats',requiredAttributes={"center":self.currentProtocol,"subject":self.currentSubject})
#         nobiasimages = list( nobiasRDI._findValues( {}, None, False ) )
#         id_pre = [x for x in range(len(nobiasimages)) if 'T1pre' in str(nobiasimages[x])]
#         nobiasPre = str(nobiasimages[id_pre[0]])
#         # Run SPM segmentation
#         pathTPMseg = os.path.join(str(self.prefs['spm']),'tpm','TPM.nii')
#         import copy
#         splittedName = nobiasPre.split('/')
#         c1Name = copy.deepcopy(splittedName)
#         c2Name = copy.deepcopy(splittedName)
#         c3Name = copy.deepcopy(splittedName)
#         c4Name = copy.deepcopy(splittedName)
#         c1Name[-1] = str("c1")+c1Name[-1]
#         c1Name = '/'.join(c1Name)
#         c2Name[-1] = str("c2")+c2Name[-1]
#         c2Name = '/'.join(c2Name)
#         c3Name[-1] = str("c3")+c3Name[-1]
#         c3Name = '/'.join(c3Name)
#         c4Name[-1] = str("c4")+c4Name[-1]
#         c4Name = '/'.join(c4Name)
#         splittedName[-1]=str("WithoutGado.nii")
#         nogadoPre = str('/'.join(splittedName))
#         call = matlab_removeGado%("'"+str(self.prefs['spm'])+"'","'"+nobiasPre+",1'","'"+str(pathTPMseg)+",1'","'"+str(pathTPMseg)+",2'","'"+str(pathTPMseg)+",3'","'"+str(pathTPMseg)+",4'","'"+str(pathTPMseg)+",5'","'"+str(pathTPMseg)+",6'",\
#                "'"+c1Name+"'","'"+c2Name+"'","'"+c3Name+"'","'"+c4Name+"'","'"+nobiasPre+"'","'"+nogadoPre+"'")
#         thr = matlabRunNB(call, lambda :self.morphologistGado3(nobiasPre, nogadoPre))
#         self.threads.append(thr)
#         thr.start()
            
#     def morphologistGado3(self, nobiasPre, nogadoPre):         
#         print self.currentSubject + ": Segmentation gado done."
#         # Replace segmented nobias image with segmented image
#         print self.currentSubject + ": Replacing nobias.nii with segmented image..."
#         nobiasBak = os.path.join(getTmpDir(),self.currentSubject + 'backup.nii')
#         cmd1 = ['mv', nobiasPre, nobiasBak]
#         cmd2 = ['cp', nogadoPre, nobiasPre]
#         line1 = runCmd(cmd1)
#         line2 = runCmd(cmd2)
#         # Force SPM volume to be saved in S16 (otherwise brainvisa4.6 crashes)
#         ret = subprocess.call(['AimsFileConvert', '-i', nobiasPre, '-o', nobiasPre, '-t', 'S16'])
#         if ret < 0:
#             QtGui.QMessageBox.warning(self, "Error", u"Error in the conversion of the segmented image to int16.")
#             return
#         
#         # Execute the rest of the Morphologist pipeline
#         morphologist = getProcessInstance('morphologist')
#         morphologist.executionNode().PrepareSubject.setSelected(False)
#         morphologist.executionNode().BiasCorrection.setSelected(False)
#         morphologist.executionNode().HistoAnalysis.setSelected(True)
#         morphologist.executionNode().BrainSegmentation.setSelected(True)
#         morphologist.executionNode().Renorm.setSelected(True)
#         morphologist.executionNode().SplitBrain.setSelected(True)
#         morphologist.executionNode().TalairachTransformation.setSelected(False)
#         morphologist.executionNode().HeadMesh.setSelected(True)
#         morphologist.executionNode().HemispheresProcessing.setSelected(True)
#         morphologist.executionNode().SulcalMorphometry.setSelected(True)
#         # Synchronous computation
#         self.brainvisaContext.runInteractiveProcess(lambda x='':self.morphologistGado4(nobiasPre, nobiasBak), morphologist, t1mri = self.mriAcPc, perform_normalization = False, anterior_commissure = self.AcPc['AC'],\
#                 posterior_commissure = self.AcPc['PC'], interhemispheric_point = self.AcPc['IH'], left_hemisphere_point = self.AcPc['LH'], perform_sulci_recognition = True)

#     def morphologistGado4(self, nobiasPre, nobiasBak):
#         # Task finishesd
#         self.taskfinished(self.currentSubject + u': BrainVISA segmentation and meshes generation')
#         # Restore initial nobias image
#         print self.currentSubject + ": Restoring original nobias.nii..."
#         cmd = ['mv', nobiasBak, nobiasPre]
#         line1 = runCmd(cmd)
#         # Compute MarsAtlas segmentation
#         self.hiphopStart()

            
    def hiphopStart(self, center, subject):
        self.taskfinished(self.currentSubject + u': Morphologist segmentation and meshes generation')
        self.setStatus(self.currentSubject + u": Starting Hip-Hop")

        Lrdi = ReadDiskItem('Labelled Cortical folds graph', 'Graph and data', requiredAttributes={ 'side': 'left', 'subject':subject, 'center':center})
        Rrdi = ReadDiskItem('Labelled Cortical folds graph', 'Graph and data', requiredAttributes={ 'side': 'right', 'subject':subject, 'center':center})
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

    #def spmRegisterPatient(self, protocol, patient, acq):
        ## Comment choisir le moment pour recaler les post avec les pre, si l'import est fait dans l'ordre post puis pre ?
        ## Solution :
        ## - a chaque importation, on regarde les acquisitions existantes pour le patient concerné
        #images = self.findAllImagesForSubject(protocol, patient)
        #nameAcq=[t.attributes()['acquisition'] for t in images]
        #na = [n.split('_')[0] for n in nameAcq] # same without the acquisition date e.g. 'T1pre' instead of 'T1pre_2001_01_01'
        #current = [t for t in images if str(t.attributes()['acquisition']) == str(acq)][0]
    
        #thrs = []
    
        ## - si on est en train d'importer la T1pre, on recale toutes les autres vers cette dernière
        #if acq.split('_')[0] == 'T1pre':
        #thrs.append(self.spmCoregister([t.fileName() for t in images if t != current], current.fileName()))
    
        ## - si on a une T1pre existante, on recale la nouvelle image vers la T1pre
        #if 'T1pre' in na and acq.split('_')[0] != 'T1pre':
        #thrs.append(self.spmCoregister([current.fileName(),], images[na.index('T1pre')].fileName()))
    
        ## - si on a une T1, on la recale vers MNI (que ce soit pre ou post) et on marque dans la base de données le recalage MNI actif
        #if acq.startswith('T1'):
        #self.setStatus("Normalisation de %s vers le referentiel MNI en cours..."%acq)
        #thrs.append(self.spmNormalize(current.fileName(), protocol, patient, acq))
    
        #return thrs


#     def spmCoregister(self, image, target):
#         """ Coregisters an image (ReadDiskItem) to a target image (filepath) - rigid registration"""
#         # Temporary txt file to store the trm transformation
#         tmpOutput = getTmpFilePath('txt')
#         #truc = subprocess.call(['AimsFileInfo', '-i', str(image.fileName())])
#         #someinfo_keys = image._minfAttributes.keys()
#         #IndexDataType = someinfo_keys.index('data_type')
#         #someinfo_values = image._minfAttributes.values()
#         imageFileName = image.fileName()
#         if  image.attributes()['data_type'] == 'RGB':
#             print "it is RGB"
#             imageFileName = getTmpFilePath('nii')
#             ret = subprocess.call(['AimsFileConvert', '-i', str(image.fileName()), '-o', str(imageFileName), '-t', 'S16'])
#             if ret < 0:
#                 print "Conversion to S16 error: "+repr(registeredPath) #terminal
#                 QtGui.QMessageBox.warning(self, "Error", u"The conversion into S16 didn't worked!") #utilisateur
#                 return
#     
#         if 'brainCenter' not in image.attributes() or 'brainCenter' not in target.attributes():
#             call = spm_coregister%("'"+str(self.prefs['spm'])+"'","'"+str(imageFileName)+",1'", "'"+str(target)+",1'", str([0, 0 ,0]), str([1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1]), str([0, 0, 0]), str([1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1]),"'"+tmpOutput+"'")
#         else:
#             call = spm_coregister%("'"+str(self.prefs['spm'])+"'","'"+str(imageFileName)+",1'", "'"+str(target)+",1'", str(image.attributes()['brainCenter']), str(image.attributes()['SB_Transform']), str(target.attributes()['brainCenter']), str(target.attributes()['SB_Transform']), "'"+tmpOutput+"'")
#     
#     
#         # TODO CHECK REGISTRATION AT THE END !
#         thr = matlabRunNB(call)#, lambda: self.taskfinished("spm coregister")) #matlabRunNB(call)
#         thr.finished.connect(lambda:self.taskfinished(u"SPM Coregister done", thr))
#         thr.finished.connect(lambda im=image, trm=tmpOutput:self.insertTransformationToT1pre(trm,im))
#         if image.attributes()['data_type'] == 'RGB':
#             thr.finished.connect(lambda im=imageFileName:os.remove(im) and os.remove(im+'.minf'))
#         self.threads.append(thr)
#         thr.start()
#         return thr


#     def spmCoregisterResample(self, image, target):
#         """ Coregisters an image (ReadDiskItem) to a target image (filepath) - rigid registration
#            !!! Replaces the original image with the resampled image in the database !!!
#            If the registration is faulty, the image must be imported again.
#         """
#         call = spm_coregisterReslice%("'"+str(self.prefs['spm'])+"'","{'"+str(target)+",1'}", "{'"+str(image.fileName())+",1'}")
#         spl = os.path.split(image.fileName())
#         registeredPath = os.path.join(spl[0], 'r'+spl[1])
#         thr = matlabRunNB(call)#, lambda: self.taskfinished("spm coregister")) #matlabRunNB(call)
#         thr.finished.connect(lambda:self.taskfinished(u"SPM Coregister-Resampling done", thr))
#         thr.finished.connect(lambda im=image, resamp = registeredPath:self.setResampledToT1pre(im, resamp))
#         self.threads.append(thr)
#         thr.start()
#         return thr

    def spm_template_t1(self):
        spm_version = checkSpmVersion(str(self.prefs['spm']))
        if spm_version == '(SPM12)':
            return "{'"+str(self.ui.prefSpmTemplateEdit.text())+os.sep+'tpm/TPM.nii'+",1'}"  #ou tpm/TPM.nii' parce que pour la normalisation dartel c'est tpm.nii /'toolbox/OldNorm/T1.nii'
        elif spm_version == '(SPM8)':
            return "{'"+str(self.ui.prefSpmTemplateEdit.text())+os.sep+'templates/T1.nii'+",1'}"

    def spmNormalize(self, image, protocol, patient, acq):
        """ Normalize one image (filepath of nifti file) to the T1 template MNI referential"""
        # Check SPM version
        spm_version = checkSpmVersion(str(self.prefs['spm']))
        if spm_version == '(SPM12)':
            print 'SPM12 used'
            call = spm12_normalise%("'"+str(self.prefs['spm'])+"'","{'"+str(image)+",1'}", "{'"+str(image)+",1'}", self.spm_template_t1())
        elif spm_version == '(SPM8)':
            print 'SPM8 used'
            call = spm8_normalise%("'"+str(self.prefs['spm'])+"'","{'"+str(image)+",1'}", "{'"+str(image)+",1'}", self.spm_template_t1())
        # Call SPM normalization
        matlabRun(call)
        # Register new files
        self.insertSPMdeformationFile(protocol, patient, acq)
        self.StatisticDataMNItoScannerBased(protocol, patient, acq)
        
#         thr = matlabRunNB(call)
#         # Try to connect the thread signals to a function here that will notice when it finishes
#         thr.finished.connect(lambda:self.taskfinished(u"SPM Normalize done", thr))
#         thr.finished.connect(lambda pr=protocol, pat=patient, a=acq:self.insertSPMdeformationFile(pr, pat, a))
#         thr.finished.connect(lambda pr=protocol, pat=patient, a=acq:self.StatisticDataMNItoScannerBased(pr, pat, a))
#         self.threads.append(thr)
#         thr.start()
#         return thr

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
        # Move the trm file to the database
        shutil.move(trmpath, transformT1.fullPath())
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
        spm_version = checkSpmVersion(str(self.prefs['spm']))
        if spm_version == '(SPM8)':
            print 'SPM8 used'
            wdi = WriteDiskItem('SPM2 normalization matrix',  'Matlab file')#'gz compressed NIFTI-1 image' )
            di = wdi.findValue( { 'center': protocol, 'subject' : patient, 'acquisition':acq } )
        elif spm_version == '(SPM12)':
            print 'SPM12 used'
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

            diT1 = ReadDiskItem( 'Raw T1 MRI', 'BrainVISA volume formats', requiredAttributes={'center':protocol, 'subject':patient } )
            allT1 = list(diT1.findValues({},None,False))
            idxT1pre = [i for i in range(len(allT1)) if 'T1pre' in str(allT1[i])]
            self.mriAcPc = allT1[idxT1pre[0]]

            spm_version = checkSpmVersion(str(self.prefs['spm']))
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
                        if spm_version == '(SPM12)':
                            print 'SPM12 used'
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
        
                        else:
                            print "doesn't work with SPM8"
                            return
    
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

        di = ReadDiskItem( 'Raw T1 MRI', 'BrainVISA volume formats', requiredAttributes={'center':protocol, 'subject':patient } )
        allT1 = list(di.findValues({},None,False))
        idxT1pre = [i for i in range(len(allT1)) if 'T1pre' in str(allT1[i])]
        T1pre = allT1[idxT1pre[0]]
        self.storeImageReferentialsAndTransforms(T1pre)
        constraints =  { 'center':protocol, 'subject':patient, 'acquisition':T1pre.attributes()['acquisition'] }

        #il faudrait mettre le referentiel à freesurferatlas la aussi ça serait fait comme ça.
        self.transfoManager.setReferentialTo(diFS, T1pre.attributes()['referential'] )

        if AmygdalaRight:
            for ii in range(len(notrightamygdalapx[0])):
                volDestrieux.setValue(0,notrightamygdalapx[3][ii],notrightamygdalapx[2][ii],notrightamygdalapx[1][ii])

            aims.write(volDestrieux, os.path.join(getTmpDir(),'rightamygdala.nii'))
            volDestrieux = aims.read(diFS.fullPath())
            wdirightamygdala =  WriteDiskItem('rightAmygdala', 'GIFTI file' )
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
            for ii in range(len(notleftamygdalapx[0])):
                volDestrieux.setValue(0,notleftamygdalapx[3][ii],notleftamygdalapx[2][ii],notleftamygdalapx[1][ii])

            aims.write(volDestrieux,os.path.join(getTmpDir(),'leftamygdala.nii'))
            volDestrieux = aims.read(diFS.fullPath())
            wdileftamygdala =  WriteDiskItem('leftAmygdala', 'GIFTI file' )
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
            for ii in range(len(notrighthippopx[0])):
                volDestrieux.setValue(0,notrighthippopx[3][ii],notrighthippopx[2][ii],notrighthippopx[1][ii])

            aims.write(volDestrieux,os.path.join(getTmpDir(),'righthippo.nii'))
            volDestrieux = aims.read(diFS.fullPath())
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

            #ret = subprocess.call(['AimsMeshCut', '-i', str(dirightHippo.fullPath()), '-o', os.path.join(getTmpDir(),'testhippocut.gii'), '-a',str(rotation[2,2]),'-b',str(rotation[2,1]),'-c',str(rotation[2,0]),'-d',str(-numpy.inner(rotation[2,:],center)),'-p',os.path.join(getTmpDir(),'testplan.nii')])
            #ret = subprocess.call(['AimsMeshCut', '-i', str(dirightHippo.fullPath()), '-o', os.path.join(getTmpDir(),'testhippocut2.gii'), '-a',str(-rotation[2,2]),'-b',str(-rotation[2,1]),'-c',str(-rotation[2,0]),'-d',str(numpy.inner(rotation[2,:],center)),'-p',os.path.join(getTmpDir(),'testplan2.nii')])

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

                hippo_vol_antero *=5301
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
            for ii in range(len(notlefthippopx[0])):
                volDestrieux.setValue(0,notlefthippopx[3][ii],notlefthippopx[2][ii],notlefthippopx[1][ii])

            aims.write(volDestrieux,os.path.join(getTmpDir(),'lefthippo.nii'))
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
                #ret = subprocess.call(['AimsMeshCut', '-i', str(dileftHippo.fullPath()), '-o', os.path.join(getTmpDir(),'testhippocut.gii'), '-a',str(rotation[2,2]),'-b',str(rotation[2,1]),'-c',str(rotation[2,0]),'-d',str(-numpy.inner(rotation[2,:],center))])
                #ret = subprocess.call(['AimsMeshCut', '-i', str(dileftHippo.fullPath()), '-o', os.path.join(getTmpDir(),'testhippocut2.gii'), '-a',str(-rotation[2,2]),'-b',str(-rotation[2,1]),'-c',str(-rotation[2,0]),'-d',str(numpy.inner(rotation[2,:],center))])

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

                    hippo_vol_antero *=1701
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


    def runFreesurferReconAll(self):

        proto = str(self.ui.regProtocolCombo.currentText())
        subj = str(self.ui.regSubjectCombo.currentText())
        images = self.findAllImagesForSubject(proto, subj)
        # Find all images, normalize the T1s, find the T1pre if there is one, read image headers, store their referentials and transformations in the DB
        t1preImage = None
        FlairpreImage = None
        T2preImage = None
        FreesurferMethod = 0
        for image in images:
            patient = image.attributes()['subject']
            acq = image.attributes()['acquisition']
            # Store Scanner-based referential and referential files in the DB
            if acq.startswith('T1'):
                # If there is a T1pre, remember the image
                if acq.find('T1pre') == 0:
                    t1preImage = image

            elif acq.startswith('FLAIR'):
                if acq.find('FLAIRpre') == 0:
                    FlairpreImage = image

            elif acq.startswith('T2'):
                if acq.find('T2pre') == 0:
                    T2preImage = image

        # No T1pre : nothing else to do
        if t1preImage is None:
            print('no T1pre image')
            return

        self.brainvisaContext.write("recon-all from freesurfer")

        if t1preImage and FlairpreImage:
            FreesurferMethod = 1
        elif t1preImage and T2preImage:
            FreesurferMethod = 2
        elif t1preImage:
            FreesurferMethod = 0

        #prevenir que freesurfer demarre dans la bare de tache
        #il faut rajouter le test si le sujet exist déjà dans freesurfer il faut d'abord le supprimer pour pouvoir le relancer.

        if FreesurferMethod == 0:
            try:
                launchFreesurferCommand(context, None, 'recon-all', '-all', '-subjid', subj ,'-i',str(t1preImage.fullPath()))
            except:
                launchFreesurferCommand(context, None, 'recon-all', '-all', '-subjid', subj)
        elif FreesurferMethod == 1:
            try:
                launchFreesurferCommand(context, None, 'recon-all', '-all', '-subjid', subj ,'-i',str(t1preImage.fullPath()),'-FLAIRpial', '-FLAIR', str(FlairpreImage.fullPath()))
            except:
                shutil.rmtree(os.path.join(configuration.freesurfer._get_subjects_dir_path(),subj))
                launchFreesurferCommand(context, None, 'recon-all', '-all', '-subjid', subj ,'-i',str(t1preImage.fullPath()),'-FLAIRpial', '-FLAIR', str(FlairpreImage.fullPath()))
        elif FreesurferMethod == 2:
            try:
                launchFreesurferCommand(context, None, 'recon-all', '-all', '-subjid', subj ,'-i',str(t1preImage.fullPath()),'-T2pial', '-T2',str(T2preImage.fullPath()))
            except:
                shutil.rmtree(os.path.join(configuration.freesurfer._get_subjects_dir_path(),subj))
                launchFreesurferCommand(context, None, 'recon-all', '-all', '-subjid', subj ,'-i',str(t1preImage.fullPath()),'-T2pial', '-T2',str(T2preImage.fullPath()))

        self.importFSoutput(subject=subj)



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


    def downloadPacs(self):
        """ Import the selected images from the pacs to the Nifti output directory """
        pass
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

    def setNoDBFilePath(self, path=None):
        if path is None:
            path = QtGui.QFileDialog.getExistingDirectory(self, u"Select path for file export not include in database")
        if path is not None:
            self.ui.prefNoDBFileLocationEdit.setText(path)

    def setSpmStandalonePath(self, path=None):
        if path is None:
            path = QtGui.QFileDialog.getExistingDirectory(self, u"Select SpmStandalone path")
        if path is not None:
            self.ui.prefSpmStandaloneEdit.setText(path)

    def setPrefCoregister(self,key):

        if key == 'ANTS':
            if self.ui.prefANTScheckbox.isChecked():
                self.ui.prefSPMcheckbox.setCheckState(False)
                self.ui.prefAimscheckbox.setCheckState(False)
    
        elif key == 'SPM':
            if self.ui.prefSPMcheckbox.isChecked():
                self.ui.prefANTScheckbox.setCheckState(False)
                self.ui.prefAimscheckbox.setCheckState(False)
    
        elif key == 'Aims':
            if self.ui.prefAimscheckbox.isChecked():
                self.ui.prefANTScheckbox.setCheckState(False)
                self.ui.prefSPMcheckbox.setCheckState(False)
    

    def switchProjectButton(self):
        projects= [self.ui.radioButtonProject.isChecked(),self.ui.radioButtonProject_2.isChecked(),self.ui.radioButtonProject_3.isChecked(),self.ui.radioButtonProject_4.isChecked()]
        project_selected = [x for x in range(len(projects)) if projects[x]==True]

        self.prefs['projectSelected'] = project_selected

    ################################## Patient Information Tab #####################################
    def selectPatSubject(self, subj):
        pass


    ################################## SEEG Import Tab #####################################
    def selectSeegSubject(self, subj):
        # Fill this from the database ReadDiskItems :  {"MANIP_20-09-2012-2":ReadDiskItem, ...}
        self.seegManips = {}
        # Update the list of manips
        self.ui.seegManipsBaseCombo.clear()
        self.ui.seegManipsBaseCombo.addItems(self.seegManips.keys())

    def manipSelected(self, manip):
        self.ui.seegAvailableProcessingCombo.clear()
        methods = seegprocessing.getProcessingMethods(seegprocessing.getManipNameFromDirName(str(manip)))
        self.ui.seegAvailableProcessingCombo.addItems([method[0] + ' - '+method[1] for method in methods])

    def chooseSeeg(self):
        """Choose a SEEG file to import"""
        path = str(QtGui.QFileDialog.getOpenFileName(self, "Open SEEG File", "", "TRC SEEG file (*.trc)"))
        print "Path of the TRC/EEG %s"%path
        if not path:
          return
        self.ui.seegFileLabel.setText(path)

    def seegImport(self, subj = None, proto = None):
        """Import a TRC file in the database, anonymize it and convert using ELAN"""
        #QtGui.QMessageBox.warning(self, u'Erreur', u"L'importation de fichier SEEG n'est pas encore implémenté !")
        path = str(self.ui.seegFileLabel.text())
        if not os.path.isfile(path):
            print path + " is not a valid file !"
            return
        date = self.ui.seegAcqDate.date()
        manip = str(self.ui.seegManipName.currentText())
        submanip = str(self.ui.seegSubmanipName.currentText())
        acq = str(self.ui.seegAcqCombo.currentText())
    
        if proto is None:
            proto = str(self.ui.seegProtocolCombo.currentText())
        if subj is None:
            subj = str(self.ui.seegSubjectCombo.currentText())
        # Find the path in the DB
        wdi = WriteDiskItem('Raw SEEG recording', 'EEG TRC format' )
        constraints =  { 'center': proto, 'experiment': str(manip), 'subject' :subj }
    
        if len(submanip) > 0:
            constraints['subexperiment'] = submanip
        if len(acq) > 0:
            constraints['expId'] = acq
        di = wdi.findValue(constraints)
        if di is None:
            print "TRC import : could not find valid path"
    
        # Copy the file
        #generate path
    
        #first I have to anonymize.
    
        if not os.path.exists(os.path.dirname(di.fullPath())):
            try:
                os.makedirs(os.path.dirname(di.fullPath()))
            except:
                print "can't generate folder"
    
        shutil.copyfile(path, di.fullPath())
    
        # Anonymize in-place   IL FAUDRAIT REDEMANDER LE NOM ET LE PRENOM DU PATIENT ?
        self.anonymizeTRC_Sys98_t4(di.fullPath(),firstname=subj, lastname=subj)
    
        #anonymizeTRC.anonymizeTRC(di.fullPath(), lastname = subj)  #firstname = str(self.ui.firstnameToRemove.text()), lastname = str(self.ui.lastnameToRemove.text())
        neuroHierarchy.databases.insertDiskItem(di, update=True )

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


    def anonymizeTRC_Sys98_t4(filepath, firstname="No", lastname="Name", nowrite=False, overwriteMontage=True, overwriteUndocumentedElectrode = True): # Anonymize Micromed's System 98 type 4 TRC files
        """ Anonymize Micromed's System 98 type 4 TRC files """
        fo=open(filepath, "r+b") # Open read-write binary mode
        # should check for a TRC header !
        fo.seek(2)
        headerTitle = fo.read(26)
        fo.seek(175)
        headerType = fo.read(1)
        trcVersion =  struct.unpack('B',headerType)[0]
        if headerTitle != "MICROMED  Brain-Quick file" or (trcVersion != 4 and trcVersion != 3):
            print "Not a MICROMED System 98 Type 3/4 TRC file -> ignoring"
            fo.close()
            return False
        fo.seek(64) # go to patient data offset in header
        # get a 22-char string, padding with spaces, convert to integers and pack as 22 unsigned chars in little endian
        if nowrite:
            fo.close()
            return True
        fo.write(struct.pack('<22B',*[ord(a) for a in lastname[:22].ljust(22,' ')]))
        # Same with 20 chars
        fo.write(struct.pack('<20B',*[ord(a) for a in firstname[:20].ljust(20,' ')]))
        # Same with date (3 unsigned chars, for example [10, 05, 72] for october 5th 1972
        fo.seek(107)
        fo.write(struct.pack('<B',1))

        if overwriteMontage:
            try:
                if trcVersion == 3:
                     sizeMontage = 3072
                elif trcVersion == 4:
                    sizeMontage = 4096
                # Read unsigned short int number of montages
                fo.seek(152)
                nbMontages = struct.unpack('H',fo.read(2))[0]
                # Reading Montage header to get offset of montage data
                fo.seek(288)
                if fo.read(8) != 'MONTAGE ':
                    raise Exception("Incorrect MONTAGE header")
                montageOffs = struct.unpack('I', fo.read(4))[0]
                # Montage name
                for i in range(nbMontages):
                    newMontageName = "Montage "+str(i)
                    fo.seek(montageOffs + sizeMontage*i+264)
                    # To print description string : desc = f.read(64); print desc[:desc.find('\x00')]
                    fo.write(struct.pack('<64B',*[ord(a) for a in newMontageName.ljust(64,'\x00')]))

                # Montage name in HISTORY : find history offset,
                fo.seek(336)
                if fo.read(8) != 'HISTORY ':
                    raise Exception ("Incorrect HISTORY header")
                historyOffs = struct.unpack('I', fo.read(4))[0]
                # There is first an area unsigned long int[MAX_SAMPLE] --> 128* 4-bytes values (the sample where themontage was changed when viewing)
                # We can skip that, then MAX_HISTORY = 128  "specific montage" structures which are identical to the montages above.
                # Description string starts at offset 264 of each montage structure and is 64 bytes long,, just like above.
                for i in range(30):
                    newMontageName = "Montage History "+str(i)
                    fo.seek(historyOffs + 128*4 + sizeMontage*i+264)
                    desc = fo.read(64)
                    # To print description string : desc = f.read(64); print desc[:desc.find('\x00')]
                    if "".join(["\x00" for j in range(64)]) != desc:
                        fo.seek(historyOffs + 128*4 + sizeMontage*i+264)
    
                # undocumented string in place of the last electrode in  name in HISTORY : find history offset,
                fo.seek(192)
                if fo.read(8) != 'LABCOD  ':
                    raise Exception ("Incorrect LABCOD header")
                elecOffs = struct.unpack('I', fo.read(4))[0]
                # There should be 640 electrode structures of 128 bytes. But the last one is not an electrode structure. It seems to be a 32 bits integer and a 64 bytes string that contains a montage name...
                if overwriteUndocumentedElectrode:
                    fo.seek(elecOffs + 639*128 + 4)
                    fo.write(struct.pack('<64B',*[ord(a) for a in "undocumented montage".ljust(64,'\x00')]))
            except:
                print "Could not overwrite Montage name"

        fo.close()
        return True



# =============================================================================
# MAIN: Main function that starts the interface
# =============================================================================
def main():
    
    # Create application
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


