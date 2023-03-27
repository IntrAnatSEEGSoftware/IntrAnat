# -*-coding:utf-8 -*

# import os, sys, pickle, urllib2, json, subprocess
import os, sys, pickle, urllib2, json, subprocess, ssl

from brainvisa import axon
from brainvisa.data.writediskitem import ReadDiskItem
from soma.qt_gui.qt_backend import uic, QtGui, QtCore

from progressbar import *
from externalprocesses import *

# Distant account
# ssh_account = 'odavid@f-tract.eu'
# ssh_account = 'odavid@gin-serv.ujf-grenoble.fr'
# ssh_account = 'davido@129.88.196.130'
ssh_account = 'davido@gin-serv.ins-amu.fr'
ssh_key = '/home/odavid/.ssh/id_rsa' # Now specified to avoid conflicts with the key in Brainvisa container (/casa/home) when executing ssh commands with pOpen

# Correction problem in the nifti file.
spm_reinitialiaze_mat = """try
    addpath(genpath(%s));
    VF = spm_vol(%s);
    YVF = spm_read_vols(VF);
    VF.private.mat0 = VF.mat;
    VF.private.mat = VF.mat;
    spm_write_vol(VF, YVF);
catch
    disp 'AN ERROR OCCURED';
end
quit;"""
# Path to SPM
spmpath = None

# Transfer files with RSYNC
def transferFileRsync(srcDir, destDir, isDelete=False):
    # cmd = ['scp', '-P', '206', '-rp', srcDir, destDir]
    if isDelete:
        cmd = ['rsync', '-avzhO', '-e ssh -p 206 -i ' + ssh_key, srcDir, destDir + '/', '--delete']
    else:
        cmd = ['rsync', '-avzhO', '-e ssh -p 206 -i ' + ssh_key, srcDir, destDir + '/']
    print('Copy: ' + ' '.join(cmd))
    subprocess.Popen(cmd, stdout=subprocess.PIPE, env = dict()).communicate()


# Transfer files with SCP
def transferFileScp(srcDir, destDir):
    cmd = ['scp', '-P', '206', '-i', ssh_key, '-rp', srcDir, destDir]
    print('Copy: ' + ' '.join(cmd))
    subprocess.Popen(cmd, stdout=subprocess.PIPE, env = dict()).communicate()


# Run SSH command remotely
def runSSH(cmdRemote):
    cmd = ['ssh', '-Y', ssh_account, '-p', '206', '-i', ssh_key, cmdRemote]
    print('Register: ' + ' '.join(cmd))
    subprocess.Popen(cmd, stdout=subprocess.PIPE, env = dict()).communicate()   


# Delete local folder
def deleteLocalFolder(localDir):
    cmd = ['rm', '-rf', localDir]
    print('Delete: ' + ' '.join(cmd))
    subprocess.Popen(cmd, stdout=subprocess.PIPE, env = dict()).communicate()


class DialogOverwrite(QtGui.QDialog):
    def __init__(self, parent, strTitle, strMessage):
        super(DialogOverwrite, self).__init__(parent)

        msgBox = QtGui.QMessageBox()
        msgBox.setText(strMessage)
        msgBox.setWindowTitle(strTitle)
        msgBox.addButton(QtGui.QPushButton('Yes'), QtGui.QMessageBox.YesRole)
        msgBox.addButton(QtGui.QPushButton('No'), QtGui.QMessageBox.NoRole)
        msgBox.addButton(QtGui.QPushButton('Only FreeSurfer'), QtGui.QMessageBox.HelpRole)
        ret = msgBox.exec_()
        if ret == 0:
            self.isOverwrite = 1
        elif ret == 1:
            self.isOverwrite = 0
        elif ret == 2:
            self.isOverwrite = 2 
        else:
            self.isOverwrite = None
            
            
# Main Class
class ftractTransfer(QtGui.QDialog):

    def __init__(self, locateData=None):
        # UI init
        QtGui.QWidget.__init__(self)
        self.ui = uic.loadUi("ftractTransfer.ui", self)
        self.setWindowTitle('ftractTransfer')

        self.ui.filterSiteCombo.currentIndexChanged.connect(self.updatePatientFilters)
        self.ui.filterYearCombo.currentIndexChanged.connect(self.updatePatientFilters)
        self.ui.patientSelectAllButton.clicked.connect(lambda:self.selectAllListItems(self.filteredPatientList, True) )
        self.ui.selectedSelectAllButton.clicked.connect(lambda:self.selectAllListItems(self.selectedPatientList, True) )
        self.ui.patientDeselectAllButton.clicked.connect(lambda:self.selectAllListItems(self.filteredPatientList, False) )
        self.ui.selectedDeselectAllButton.clicked.connect(lambda:self.selectAllListItems(self.selectedPatientList, False) )
        self.ui.radioButtonPtoL.toggled.connect(self.populateFromDB)
        self.ui.CSVcheckBox.clicked.connect(self.updatePatientFilters) #.isChecked()
        self.ui.filterAddPatientButton.clicked.connect(lambda:self.moveSelectedItemsToOtherListWidget(self.filteredPatientList, self.selectedPatientList))
        self.ui.filterRemovePatientButton.clicked.connect(lambda:self.moveSelectedItemsToOtherListWidget(self.selectedPatientList, self.filteredPatientList))
        self.ui.transferButton.clicked.connect(self.TransferPatientScp)
        # Get list of subjects
        self.populateFromDB()


    # Test if local=>playground or playground=>local
    def isLocalToPlayground(self):
        
        FTsens = [self.ui.radioButtonLtoP.isChecked(),self.ui.radioButtonPtoL.isChecked()]
        FTsens_selected = [x for x in range(len(FTsens)) if FTsens[x]==True]
        return (FTsens_selected[0] == 0)
    

    def populateFromDB(self, thread=None):
        # Local to playground
        if self.isLocalToPlayground():
            self.ui.CSVcheckBox.setEnabled(True)
            print "Update patients: Local to playground"
            rdi = ReadDiskItem( 'Subject', 'Directory',requiredAttributes={'center':'Epilepsy'})
            subjects = list( rdi._findValues( {}, None, False ) )
            self.subjects = dict([(s.attributes()['subject'], {'rdi':s, 'center':s.attributes()['center']}) for s in subjects])

            sites = ['*',] + sorted(set([s.split('_')[0] for s in self.subjects]))
            years = ['*',] + sorted(set([s.split('_')[1] for s in self.subjects if len(s.split('_')) > 1]))

            self.ui.filterSiteCombo.clear()
            self.ui.filterSiteCombo.addItems(sites)
            self.ui.filterYearCombo.clear()
            self.ui.filterYearCombo.addItems(years)
            self.ui.selectedPatientList.clear()

        # Playground to local
        else:
            print "Update patients: Playground to local"
            self.ui.CSVcheckBox.setEnabled(False)
            proxy_handler = urllib2.ProxyHandler({})
            opener = urllib2.build_opener(proxy_handler)
            urllib2.install_opener(opener)
            # Get the data
            # response = urllib2.urlopen('https://f-tract.eu:85/ftdata/brainvisaCRF/')
            # response = ProgressDialog.call(lambda x:urllib2.urlopen('https://f-tract.eu:85/ftdata/brainvisaCRF/'), True, self, "Getting list from server...", "Load patients")
            req = urllib2.Request('https://f-tract.eu:85/ftdata/brainvisaCRF/', headers={ 'X-Mashape-Key': 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' })
	    context = ssl._create_unverified_context()  # Only for gangstars
            response = ProgressDialog.call(lambda x:urllib2.urlopen(req, context=context), True, self, "Getting list from server...", "Load patients")
            if not response:
                self.ui.radioButtonLtoP.setChecked(True)
                return
            
            # Convert from json
            data = json.load(response)
            if data['status'] != u'ok':
                print "status from django : not ok "
                return

            self.subjects = []
            sites = ['*']
            years = ['*']

            for ii in data['crf'].keys():
                info = ii.split(" - ")
                date = info[1].split('/')
                self.subjects.append(ii)
                sites.append(info[0][4:])
                years.append(date[2])

            self.data = data
            self.ui.filterSiteCombo.clear()
            self.ui.filterSiteCombo.addItems(sorted(set(sites)))
            self.ui.filterYearCombo.clear()
            self.ui.filterYearCombo.addItems(sorted(set(years)))
            self.ui.selectedPatientList.clear()

        self.updatePatientFilters()


    def filterSubjectsBasic(self):

        # Look for csv
        filteredPatients = []
        subs = self.subjects.keys()
        for patname in subs:
            wdi_csv = ReadDiskItem('Final Export Dictionaries','CSV file',requiredAttributes={'subject':patname})
            di_csv = list(wdi_csv.findValues({},None,False))
            
            if len(di_csv)>0 and len(di_csv)<2:
               filteredPatients.append(patname)
            elif len(di_csv)==0:
                pass
                #print "no csv for this patient"
            elif len(di_csv)>=2:
                print "Error: More that one csv, delete one of the two before continuing."
                for f in di_csv:
                    print "       " + str(f)
                return

        if str(self.filterSiteCombo.currentText()) != '*':
            subs = [s for s in subs if s.split('_')[0] == str(self.filterSiteCombo.currentText())]
        if str(self.filterYearCombo.currentText()) != '*':
            subs = [s for s in subs if len(s.split('_')) > 1 and s.split('_')[1] == str(self.filterYearCombo.currentText())]
        subs = [s for s in subs if s not in [str(self.selectedPatientList.item(idx).text()) for idx in range(self.selectedPatientList.count())]]
        
        if self.ui.CSVcheckBox.isChecked():
            subs = [s for s in subs if s in filteredPatients]

        self.filteredPatientList.clear()
        self.filteredPatientList.addItems(sorted(subs))
        if not self.ui.CSVcheckBox.isChecked():
            self.colorFilterListWidget(self.filteredPatientList, filteredPatients)


    def colorFilterListWidget(self, listwidget, selected):
        
       for idx in range(listwidget.count()):
            it = listwidget.item(idx)
            if str(it.text()) in selected:
                it.setBackground(QtGui.QColor(150,255,150))
            else:
                it.setBackground(QtGui.QColor(255,150,150))
            

    def filterSubjectsBasicPlay(self):

        subs = self.subjects
        if str(self.filterSiteCombo.currentText()) != '*':
            subs = [s for s in subs if str(self.filterSiteCombo.currentText()) in s.split(" - ")[0]]
        if str(self.filterYearCombo.currentText()) != '*':
            subs = [s for s in subs if str(self.filterYearCombo.currentText()) in s.split(" - ")[1]]
        subs = [s for s in subs if s not in [str(self.selectedPatientList.item(idx).text()) for idx in range(self.selectedPatientList.count())]]
        self.filteredPatientList.clear()
        self.filteredPatientList.addItems(sorted(subs))

        filteredPatients = []
        for crf_idx in subs:
            if len(self.data['DestrieuxLabelling'][crf_idx]) > 0:
                filteredPatients.append(crf_idx)
            else:
                pass
        self.colorFilterListWidget(self.filteredPatientList, filteredPatients)


    def updatePatientFilters(self):

        if self.isLocalToPlayground():
            self.filterSubjectsBasic()
        else:
            self.filterSubjectsBasicPlay()


    def moveSelectedItemsToOtherListWidget(self, lwFrom, lwTo):
        """Takes the selected items of the list 'from' and adds them to the 'to' list"""
        for idx in reversed(range(lwFrom.count())):
            it = lwFrom.item(idx)
            if it.isSelected():
                lwTo.addItem(str(it.text()))
                lwFrom.takeItem(idx)
        lwTo.sortItems()


    def selectAllListItems(self, listwidget, select = True):
        """Selects or deselects all items of a QListWidget (select = False to deselect)"""
        for idx in range(listwidget.count()):
            listwidget.setItemSelected(listwidget.item(idx), select)
            

    def TransferPatientScp(self):
        # Get local directory where to save patients
        if not self.isLocalToPlayground():
            saveDir = str(QtGui.QFileDialog.getExistingDirectory(self, "Select directory"))
            if not saveDir:
                return
        else:
            saveDir = None
            
        # Get list of selected patients
        patients = []
        patientsExist = []
        for i in range(self.selectedPatientList.count()):
            patientName = str(self.selectedPatientList.item(i).text())
            patients.append(patientName)
            # Test if patient already exists locally
            if (not self.isLocalToPlayground()) and os.path.exists(saveDir + '/' + patientName.replace('/','_')):
                if (len(patientsExist) < 10):
                    patientsExist.append(patientName)
                elif (len(patientsExist) == 10):
                    patientsExist.append("...")
        
        # Overwrite existing files?
        if patientsExist:
            ret = DialogOverwrite(self, 'Overwrite',  u"The following folders already exist locally:\n" + "\n".join(patientsExist) + u"\n\nOverwrite existing files ?")
            if ret.isOverwrite == None:
                return
            else:
                isOverwrite = ret.isOverwrite
        else:
            isOverwrite = 0
        # Call transfer function in a separate thread
        ProgressDialog.call(lambda thr:self.TransferPatientScpWorker(patients, saveDir, isOverwrite, thr), True, self, "Copying patients...", "Transfer")
        # self.TransferPatientScpWorker(patients, saveDir, isOverwrite)
        
        
    def TransferPatientScpWorker(self, patients, saveDir, isOverwrite, thread=None):
        # Check number of selected patients
        nPatients = len(patients)
        if (nPatients == 0):
            return
        # Local to playground
        if self.isLocalToPlayground():
            # Destination folder
            dbDir = "/gin/data/database/03-preprocessed/Brainvisa/Epilepsy"
            # Copy subject by subject
            for i in range(nPatients):
                # Update progress bar
                if thread is not None:
                    #thread.emit(QtCore.SIGNAL("PROGRESS"), round(100*i/nPatients))
                    #thread.emit(QtCore.SIGNAL("PROGRESS_TEXT"), "Copying patient: " + patients[i] + "...")
                    # Refactored according to PyQt5                    
                    thread.progress.emit(round(100*i/nPatients))
                    thread.progress_text.emit("Copying patient: " + patients[i] + "...")
                # Get source and destination
                srcDir = str(self.subjects[patients[i]]['rdi'])
                destDir = dbDir
                # Copy subject
                transferFileRsync(srcDir, ssh_account + ':' + destDir, True)
                # Register in database
                runSSH("export PYTHONPATH=/home/davido/ft_pipeline:/home/davido/ft_database:/home/davido/ft_pipeline/scripts; python3 -c \"import tools.get_assign_CRF_from_path as s;s.scan_database_with_CRF(path='03-preprocessed/Brainvisa/Epilepsy/" + patients[i] + "', conf='brainvisa_epilepsy')\"")

        # Playground to local
        else:
            # Database folder on the server
            dbDir = "/gin/data/database/02-raw"
            # Copy subject by subject
            for i in range(nPatients):
                # Files to fix with SPM
                spmFiles = []
                # Update progress bar
                if thread is not None:
                    strWait = "Copying patient: " + patients[i] + "..."
                    #thread.emit(QtCore.SIGNAL("PROGRESS"), round(100*i/nPatients))
                    #thread.emit(QtCore.SIGNAL("PROGRESS_TEXT"), strWait)
                    thread.progress.emit(round(100*i/nPatients))
                    thread.progress_text.emit(strWait)
                # Get subject info
                infoT1 = self.data['t1paths'][patients[i]]
                infoPost = self.data['postpaths'][patients[i]]
                infoPostop = self.data['postoppaths'][patients[i]]
                if 'DestrieuxLabelling' in self.data.keys():
                    infoDestrieuxLabelling = self.data['DestrieuxLabelling'][patients[i]]
                else:
                    infoDestrieuxLabelling = []
                
                # Create local patient directory
                destDir = saveDir + '/' + patients[i].replace('/','_')
                if not os.path.exists(destDir):
                    os.mkdir(destDir)
                    
                # === PRE ===
                preDir = destDir + '/' + 'pre'
                if infoT1 and (not os.path.exists(preDir) or (isOverwrite == 1)):
                    if thread is not None:
                        #thread.emit(QtCore.SIGNAL("PROGRESS_TEXT"), strWait + ' (pre)')
                        thread.progress_text.emit(strWait + ' (pre)')
                    # Create pre folder
                    if not os.path.exists(preDir):
                        os.mkdir(preDir)
                    # Copy pre images
                    for kk in range(len(infoT1)):
                        localFile = os.path.join(preDir, os.path.split(infoT1[kk])[-1])
                        transferFileScp(ssh_account + ':/gin/data/database/' + infoT1[kk], localFile)
                        spmFiles.append(localFile)

                # === POST ===
                postDir = destDir + '/' + 'post'
                if infoPost and (not os.path.exists(postDir) or (isOverwrite == 1)):
                    if thread is not None:
                        #thread.emit(QtCore.SIGNAL("PROGRESS_TEXT"), strWait + ' (post)')
                        thread.progress_text.emit(strWait + ' (post)')
                    # Create pre folder
                    if not os.path.exists(postDir):
                        os.mkdir(postDir)
                    # Copy post images
                    for kk in range(len(infoPost)):
                        localFile = os.path.join(postDir, os.path.split(infoPost[kk])[-1])
                        transferFileScp(ssh_account + ':/gin/data/database/' + infoPost[kk], localFile)
                        spmFiles.append(localFile)

                # === POSTOP ===
                postopDir = destDir + '/' + 'postop'
                if infoPostop and (not os.path.exists(postopDir) or (isOverwrite == 1)):
                    if thread is not None:
                        #thread.emit(QtCore.SIGNAL("PROGRESS_TEXT"), strWait + ' (postop)')
                        thread.progress_text.emit(strWait + ' (postop)')
                    # Create postop folder
                    if not os.path.exists(postopDir):
                        os.mkdir(postopDir)
                    # Copy postop images
                    for kk in range(len(infoPostop)):
                        localFile = os.path.join(postopDir, os.path.split(infoPostop[kk])[-1])
                        transferFileScp(ssh_account + ':/gin/data/database/' + infoPostop[kk], localFile)
                        spmFiles.append(localFile)

                # === FREESURFER ===
                fsDirLocal = destDir + '/' + 'freesurfer'
                if infoDestrieuxLabelling and (not os.path.exists(fsDirLocal) or (isOverwrite >= 1)):
                    if thread is not None:
                        #thread.emit(QtCore.SIGNAL("PROGRESS_TEXT"), strWait + ' (freesurfer)')
                        thread.progress_text.emit(strWait + ' (freesurfer)')
                    # Delete existing local Freesurfer folder
                    deleteLocalFolder(fsDirLocal)
                    # Copy entire FreeSurfer folder
                    fsDir = os.path.dirname(os.path.dirname(infoDestrieuxLabelling[0]))
                    transferFileScp(ssh_account + ':/gin/data/database/' + fsDir, fsDirLocal)
                # Try to copy only Lausanne2008 segmentation
                elif infoDestrieuxLabelling:
                    lausanneDirLocal = fsDirLocal + '/' + 'parcellation_Lausanne2008'
                    if (not os.path.exists(lausanneDirLocal)):
                        if thread is not None:
                            #thread.emit(QtCore.SIGNAL("PROGRESS_TEXT"), strWait + ' (lausanne2008)')
                            thread.progress_text.emit(strWait + ' (lausanne2008)')
                        # Copy Lausanne subfolder
                        try:
                            lausanneDir = os.path.dirname(os.path.dirname(infoDestrieuxLabelling[0])) + '/' + 'parcellation_Lausanne2008'
                            transferFileScp(ssh_account + ':/gin/data/database/' + lausanneDir, lausanneDirLocal)
                        except:
                            pass
                    
                # Fix transformation matrix with SPM
                for kk in range(len(spmFiles)):
                    if thread is not None:
                        #thread.emit(QtCore.SIGNAL("PROGRESS_TEXT"), "Reinitialize transformations: " + str(kk+1) + "/" + str(len(spmFiles)))
                        thread.progress_text.emit("Reinitialize transformations: " + str(kk+1) + "/" + str(len(spmFiles)))
                    #check if .mat, .private.mat and .private.mat0 are simalar, if not rewrite the file.
                    print "SPM: Reinitialize .mat .private.mat and .private.mat0 in \"" + spmFiles[kk] + "\""
                    call = spm_reinitialiaze_mat%("'"+spmpath+"'", "'"+spmFiles[kk]+"'")
                    matlabRun(call)


if __name__ == "__main__":
    # Read preferences from .imageimport
    prefpath_imageimport = os.path.join(os.path.expanduser('~'), '.imageimport')
    try:
        if (os.path.exists(prefpath_imageimport)):
            filein = open(prefpath_imageimport, 'rU')
            prefs_imageimport = pickle.load(filein)
            spmpath = prefs_imageimport['spm']
            filein.close()
    except:
        pass
    if not spmpath or not os.path.exists(spmpath):
        print 'ERROR: SPM path not set.'
        print 'Open ImageImport and select the SPM path in the preferences tab.'
        sys.exit(1)
    
    # Start application
    app = QtGui.QApplication(sys.argv)
    axon.initializeProcesses()
    QtCore.pyqtRemoveInputHook()
    window = ftractTransfer()
    window.show()
    sys.exit(app.exec_())
