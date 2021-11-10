# ls /gin/data/database/03-preprocessed/Freesurfer/*/mri/aparc+aseg.vep.mgz
#
#
# 1. Identify VEP files from the DB
# 2. Check if there is already a copy in IntrAnat database
# 3. Import the missing ones in IntrAnat
# 4. Export CSV for each implantation
#
import subprocess
import glob

from brainvisa import axon
from brainvisa.data.readdiskitem import ReadDiskItem
from brainvisa.data.writediskitem import WriteDiskItem
from freesurfer.brainvisaFreesurfer import launchFreesurferCommand
from brainvisa.processes import defaultContext
from brainvisa.data import neuroHierarchy

from LocateElectrodes import LocateElectrodes

from externalprocesses import createItemDirs
from batch_csv import generateCsv

freesurferDB = "/gin/data/database/03-preprocessed/Freesurfer"
brainvisaDB = "/gin/data/database/03-preprocessed/Brainvisa"

def findT1pre(subject):
    rT1BV = ReadDiskItem('Raw T1 MRI', 'BrainVISA volume formats', requiredAttributes={'subject': subject})
    allT1 = list(rT1BV.findValues({}, None, False))
    if not allT1:
        return None
    idxT1pre = [i for i in range(len(allT1)) if 'T1pre' in str(allT1[i])]
    if not idxT1pre:
        return None
    if len(idxT1pre) > 1:
        print("WARNING, found multiple T1pre for subject ", subject, " Returning the first one\n", repr(allT1))

    diT1pre = allT1[idxT1pre[0]]
    return diT1pre

def importFsAtlas(subject, proto='Epilepsy', imgName='VEP', imgPath="default.mgz", diT1pre=None, fileVoxelType=None):
    """
        Import Freesurfer atlas volume in T1 referential
        subject name, protocol, image name (e.g. VEP), MGZ image path, T1pre ReadDiskItem, fileVoxelType (S16 or S32)
    """
    # If external call: Find T1pre of the subject
    if not diT1pre:
        diT1pre = findT1pre(subject)
        if not diT1pre:
            return [False, ['No T1pre found']]
    # Where to copy the new files
    acq = str(diT1pre.attributes()['acquisition']).replace('T1', imgName + "-")
    # Importing Destrieux atlas to BrainVISA database
    wdi = WriteDiskItem('FreesurferAtlas', 'NIFTI-1 image')
    di = wdi.findValue({'center': proto, 'acquisition': acq, 'subject': subject})
    # Create the folder if doesn't exist
    createItemDirs(di)
    # Reslice volume to match T1pre
    try:
        launchFreesurferCommand(defaultContext(), None, 'mri_convert', '-i', imgPath,
                                '-o', str(di.fullPath()), '-rl', str(diT1pre.fullPath()),
                                '-rt', 'nearest', '-nc')
    except:
        print("ERROR: ", "Could not launch Freesurfer command")
        return False, ["Could not launch Freesurfer command"]
    # Convert to AIMS
    if fileVoxelType is None:
        fileVoxelType = 'S16'  # Default value
        if 'data_type' in di.attributes():
            if di.attributes()['data_type'] in ['U32', 'S32']:  # Do not convert to S16 if type is U32 or S32
                fileVoxelType = di.attributes()['data_type']

    ret = subprocess.call(['AimsFileConvert', '-i', str(di.fullPath()),
                           '-o', str(di.fullPath()), '-t', fileVoxelType])
    # Add reference in the database (creates .minf)
    neuroHierarchy.databases.insertDiskItem(di, update=True)
    return True, []



if __name__ == '__main__':
    # Check if this is necessary to get ReadDiskItems
    app = QtGui.QApplication(sys.argv)
    axon.initializeProcesses()
    w = LocateElectrodes(app=app, loadAll=True, isGui=False)
    # Find available patients in BV database
    rdi = ReadDiskItem( 'Subject', 'Directory',requiredAttributes={'_ontology':'brainvisa-3.2.0'}) #, requiredAttributes={'center':'Epilepsy'} )
    w.allSubjects = list( rdi._findValues( {}, None, False ))
    w.currentProtocol = 'Epilepsy'
    w.subjects = [s.attributes()['subject'] for s in w.allSubjects if 'center' in s.attributes() and s.attributes()['center'] == w.currentProtocol]
    w.subjects = sorted(w.subjects)


    # 1) List subjects who have VEP files in Freesurfer database
    vepFiles = glob.glob(freesurferDB + "/*/mri/aparc+aseg.vep.mgz")
    print("Found ", len(vepFiles), " VEP files found")
    for v in vepFiles:
        veppath = v
        v = v.replace(freesurferDB + "/", "").replace("/mri/aparc+aseg.vep.mgz", "")
        # 2) Check if it is already imported
        # brainvisadb/Epilepsy/0001GRE_25112014/FreesurferAtlas/VEP-pre_2014-1-1/0001GRE_25112014-VEP-pre_2014-1-1.nii
        vepBV = glob.glob(brainvisaDB + "/Epilepsy/" + v + "/FreesurferAtlas/VEP-pre_*/" + v + "-VEP-pre_*.nii")
        if len(vepBV) > 0:
            print("Already imported ", v, " -> ", repr(vepBV))
            continue
        else:
            print("Importing ", v)
        # 3) Import it into BrainVisa
        worked, msg = importFsAtlas(v, proto='Epilepsy', imgName='VEP', imgPath=veppath, fileVoxelType='S32')
        if(worked):
            print("Successfully imported ", v)
        else:
            print("Failed to import ", v, " -> ", msg)
        print("Exporting CSV")
        isOk, errMsg = generateCsv(w, v, True, True)
        if isOk:
            print(v, ": CSV exported")
        else:
            print(v, ": CSV export failed: ", errMsg)


    app.quit()
    del app
    #
        # [
        #     "Compute MNI coordinates for all contacts",
        #     "Compute parcels for all contacts",
        #     "Compute MarsAtlas resection position",
        #     "Compute parcel metrics",
        #     "Save contact coordinates (.pts/.txt files)",
        #     "Save contact info (CSV file)",
        #     "Save contact info (BIDS .tsv)",
        #     "Save screenshots",
        #     "Save video (MP4)"],
        #     "Export", "Select options to run:",
        #     [True, True, False, False, True, True, self.bidspath is not None, False, False]
        # selOptions = [True, True, True, False, True, True, False, False, False]
        # locEl = LocateElectrodes(app=None, loadAll=False, isGui=False)
        # locEl.loadPatient(v)
        # locEl.exportAllWorker(selOptions)
