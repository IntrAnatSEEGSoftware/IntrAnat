Technical Documentation
***********************

.. contents::
   :depth: 3
..

IntrAnat

Technical documentation

Written : January 2014
Last updated : March 2021

List of software
===================

**IntrAnat** is a collection of software for Stereo-ElectroEncephaloGraphy (SEEG),
also called Intracerebral EEG and for Deep Brain Stimulation (DBS).

Main parts:

-  **ImageImport**: importing data into the IntrAnat/Brainvisa database

-  **locateElectrodes**: Display 3D imaging data and locate electrodes on
   images. Multimodal display (including atlases) and export electrode coordinates and parcels.

-  **groupDisplay**: Display electrode locations of multiple subjects in MNI space (**deprecated**)

-  **electrodeEditor**: creating a new model of electrode for *locateElectrodes*
   CombiView: visualisation de données de groupes de sujets

Software dependencies:
==========================

IntrAnat was developed and used mostly on Linux.

It should run on Windows and Mac OS X but bugs may be platform-specific.

Il uses multiple pieces of software :

-  **BrainVisa 4.6.1** (Started on 4.2) that should already be configured with a database.
   BrainVisa manages the database used by IntrAnat, provides many tools used by IntrAnat, notably
   Anatomist which displays 3D volumes and meshes.
-  **Qt4 and PyQt4** (included in BrainVisa)
-  **Python 2.6** (included in BrainVisa), and **Python 3**
-  **SPM12** (and **matlab**) fo intrasubject coregistration, and normalization to MNI referential.
-  **FreeSurfer** atlases can be imported and used to tag electrode contacts with the name of brain regions
-  **dcmtk** was used to import from DICOM servers but function was removed

| For Brainvisa 5 (released march 2021) code will be ported to Qt5 (see relevant git branch) e.g. rewriting
| *self.connect(self.electrodeList,
  QtCore.SIGNAL("currentRowChanged(int)"), self.electrodeSelect)*
| as
| *self.electrodeList.currentRowChanged[int].connect(self.electrodeSelect)*

**Matlab** is called through command line, so version changes should not be a problem
as long as Matlab binary is in the PATH.

For **SPM12**, functions called by IntrAnat were broken between versions in the past.
To update affected code, search for spm\_\* strings that contain spm calls in the code.

**Example** : spm_coregister ::

    spm_coregister = "try,VF=spm_vol(%s);VG=spm_vol(%s);x = spm_coreg(VF,VG);
    trm = spm_matrix(x(:)');trm = [trm(1:3,4)';trm(1:3,1:3)];
    dlmwrite(%s,trm, 'delimiter',' ','precision',16);
    catch, disp 'AN ERROR OCCURED'; end;quit;"

This string contains the matlab code using spm functions. It is then called
in *spmCoregister(self, image,target)*

How SPM coregister works : each Nifti image has an internal transform that goes to a
« scanner-based » referential. Sometimes this transform is not present in the file header
and SPM behavior is not always predictable (depends on versions and file header).

In principle, spm_coreg output transformation does not include scanner-based transforms.
So to get from an image to another, the transform should be::

    trm = VF.mat\spm_matrix(x(:)')*VG.mat

In IntrAnat, the matrix from one scanner-based ref to the other is directly stored and used by Anatomist
(Brainvisa's viewer)

Integration into BrainVisa
==========================

L'intégration comprend deux parties :

-  A toolbox « epilepsy » in BrainVisa **deprecated**
-  software using BrainVisa but that are called stand-alone (not from Brainvisa GUI)

Epilepsy toolbox
----------------

brainvisa-4.6.1/brainvisa/toolboxes/epilepsy/

This directory contains definitions of file types (electrode models, implantations
and so on), a description of how this data is organized in the database (inside
each subject directory), as well as a few « BrainVisa processes » that can be called from Brainvisa GUI

Definitions for the database
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Data stored for IntrAnat in BrainVisa database is defined in two steps:

-  in the **types** directory: epilepsy.py defines file formats (name
   and extension) and data types. Example for SEEG recordings using
   Micromed's TRC file format:
   Declare file format and extension::

      Format( 'EEG TRC format', 'f|*.trc' )

   Declare a data type called 'SEEG recording', which is not a subtype (so
   base type is 'Any Type')
   and which can be in multiple file formats::

      FileType( 'SEEG recording', 'Any Type', ['EEG TRC format', 'Elan EEG format'])#'ImaGIN matlab format'

   Now let's declare a type 'Raw SEEG recording', subtype of
   'SEEG recording', and which is stored in 'EEG TRC format'.::

      FileType( 'Raw SEEG recording', 'SEEG recording', 'EEG TRC format' )

   Various examples are visible in the main BrainVisa hierarchy
   (brainvisa-4.6.0/brainvisa/types/) but also in other toolboxes.
-  In *hierarchies* directory, multiple subdirectories, one for each version
   of the database structure.

   Currently, user database works with hierarchy brainvisa-3.2.0,
   in the same way as it is organized in brainvisa-4.6.0/brainvisa/hierarchies/.
   Directory *shared* corresponds to the internal database of Brainvisa (e.g. to
   store its templates, available for all studies).

   In brainvisa-4.6.0/brainvisa/toolboxes/epilepsy/hierarchies/brainvisa-3.2.0/,
   a few files are used to declare where IntrAnat-specific files will be inserted
   into the standard BrainVisa hierarchy.

   Example : in *images.py*, storing CT images in the database:
   Create a tuple *ct_content*, containing a string description representing
   the file name in the database. This name is an expression that will match
   a real directory in the database.

   Here, {acquisition} means that the directory name will be the same as the
   'acquisition' property of the same object. E.g. if a directory is called
   'postOp-2012-11-11', Brainvisa will set for each file in this directory
   a property 'acquisition' defined by the directory name.

   This property can be reused in directories and file names inside this directory.
   In this example, we will set a default value, and set this property as optional
   (e.g. we may want to store a subject's CT scan without setting an acquisition
   name if there will be only one::

      ct_content = (
      "{acquisition}", SetDefaultAttributeValue( 'acquisition',
      default_acquisition ), SetNonMandatoryKeyAttribute( 'acquisition' ),

   Then, declare the content of the directory::

      SetContent(

   A CT file (its type is declared just like 'SEEG recording') which name
   is subject's name, dash, and acquisition name. These values are properties
   that were declared previously and which values are known (e.g. declared
   through {acquisition})::

      "<subject>-<acquisition>", SetType( 'CT' ),

   A registration directory, to store referentials and geometrical transforms
   from CT imqge to other referentials::

      'registration', SetContent(
      'CT-<subject>_<acquisition>', SetType( 'Referential of CT' ),
      'CT-<subject>_<acquisition>_TO_Talairach-ACPC', SetType( 'Transform
      CT to Talairach-AC/PC-Anatomist' ),
      'CT-<subject>_<acquisition>_TO_Talairach-MNI', SetType( 'Transform CT
      to Talairach-MNI template-SPM'),
      'CT-<subject>_<acquisition>_TO_Scanner_Based', SetType(
      'Transformation to Scanner Based Referential' ),

   Here we are adding a new transform to another image of the same subject, with
   a specific modality and acquisition: those are new properties, declared with {}::

      'CT-<subject>_<acquisition>_TO_{modalityTarget}_{acquisitionTarget}',
      SetType( 'Transform CT to another image' ),
      'CT-<subject>_<acquisition>_Scanner_Based', SetType( 'Scanner Based
      Referential' ),
      ),
      )
      )

   Finally we insert this into the existing hierarchy:
   in directory '{protocol}/{subject}' we add a 'ct' directory, which has a 'modality'
   attribute with value 'ct'. We then add its content previously declared as ct_content.::

      insert( '{protocol}/{subject}',
      'ct', SetWeakAttr( 'modality', 'ct' ),
      apply( SetContent, ct_content)
      )

   Numerous other examples are available in main Brainvisa hierarchy and its toolboxes, e.g. ::

      brainvisa-4.3.0/brainvisa/hierarchies/brainvisa-3.1.0/base.py
      brainvisa-4.3.0/brainvisa/toolboxes/morphologist/hierarchies/brainvisa-3.1.0/anatomy.py

Patches
--------

Referential and transformation patches developed for IntrAnat were
integrated into BrainVisa from 4.5.0


Tools used for development
---------------------------

**Qt Designer / Qt Creator** to create Graphical User Interfaces.
It generates a .ui file that can be directly loaded from Python (e.g.
\__init_\_ function in ImageImportWindow)::

   from PyQt4 import uic
   self.ui = uic.loadUi("epilepsie-electrodes.ui", self)

Programmed using Python/PyQt with bindings included in BrainVisa.

Editors: Kate (on KDE), PyCharm.

IntrAnat Librairies
--------------------

Python files with useful functions:

-  *electrode.py* is used to manage electrode models and their display in Anatomist

-  *dicomutilities.py* (DEPRECATED) functions to read and access DICOM files and servers

-  *externalprocesses.py* allows to call external software (synchronous or asynchronous,
   with callback functions), especially to run matlab code.

-  *referentialconverter.py* defines an object storing multiple referentials and can convert
   coordinates from one referential to another.


ImageImport
--------------

*ImageImport* is used to add patients in the BrainVisa database,
import images (MRI, CT scans, PET scans...), de register these images
to normalize them with SPM12 (to MNI referential), and to run the main
segmentation process of BrainVisa to get brain meshes and sulci.
It can also import Freesurfer segmentations.

UI is defined in ImageImportWindow.ui

Main code is contained in the class *ImageImportWindow* in ImageImportWindow.py
and ImageImport.py is used to run it.

Software structure :

-  Buttons and other UI elements are connected to functions in \__init_\_ function
   of *ImageImportWindow* class::

     self.connect(self.ui.regSubjectCombo,
     QtCore.SIGNAL('currentIndexChanged(QString)'),
     self.setCurrentSubject)

   For PyQt5, newer syntax must be used::

     self.ui.regSubjectCombo.currentIndexChanged[str].connect(self.setCurrentSubject)

   In this example *regSubjectCombo* (a combo box of subjects in registration tab) emitting
   *currentIndexChanged* calls *self.setCurrentSubject* function.

-  code comments should help understand the code.

LocateElectrodes
----------------

LocateElectrodes allows to locate electrode models in post-implantation images
and export contact coordinates in multiple formats.
It can also display various data imported or computed
(e.g. cortical meshes from T1 MRI, CT scanner, T1, T2,
pre/post implantation/post-résection MRI...) and realistic 3D models of electrodes
(or even just contacts to improve visibility).

UI is defined in epilepsie-electrodes.ui

Software structure:
-  a few functions placed at the top of the file allow easier management of electrodes

- a main class defines functions linked to UI elements.


Data formats
==================

**Images :** MRI, CT, PET : Nifti (.nii) or compressed nifti (.nii.gz)

**SEEG**: TRC (Micromed) .eeg (ELAN). DEPRECATED: eeg signal is usually processed outside IntrAnat.

**Electrodes**: .elecmodel (python pickle format, could be converted to json).
   Elecmodel files are python dictionaries saved by pickle  library.
   An electrode is a list of cylinders as a dictionary::

      {'Plot1', {...}, 'Plot2':{...}, 'Element1':{...}}

   Elements are inactive parts of the electrode model, plots are electrical contacts.
   The electrode model is defined in a referential where the endpoint of the electrode
   is at 0,0,0, and the electrode is aligned to the Z axis. Z increases when we go towards
   the cable.
   Each part of the electrode is also defined by a dictionary::

     'Plot1': {'axis': 'Axe Z',
     'diameter': 0.80000000000000004,
     'length': 2.0,
     'position': [0.0, 0.0, 0.0],
     'type': 'Plot',
     'vector': [0.0, 0.0, 1.0]}

   Axis gives the main direction of the cylinder, diameter and length are in mm, position is
   the center of one circle of the cylinder, type (Plot or Element), and vector that is added to
   center to go to the other extremity of the cylinder.

   To import those files from python::

     import pickle
     f=open('Dixi-D08-15BM.elecdef')
     d=pickle.load(f)

** Electrode Implantations ** : .elecimplant (pickle or json), .pts, .txt

   Like elecmodel files, elecimplant files are python dictionaries saved by pickle
   They include deprecated '2mni' which used to store the linear transform to MNI referential,
   'ReferentialUuid' is a unique identifier of the referential used for electrode coordinates, which
   is Anatomist's native referential of the pre-operative T1 MRI.

   'electrodes' contains a list of electrodes. Each electrode is a dictionary with the following
   keys: 'entry' for coordinates of the entry point (on the skull), 'model' for electrode model name,
   'name' for the electrode name, 'target' for coordinates of the electrode endpoint::

     { '2mni': None,
       'ReferentialUuid': '2506a605-3d18-fa50-3557-a47922440c41',
       'electrodes': [{'entry': [126.82000732421875,
                                 121.16304779052734,
                                 135.00004577636719],
                       'model': 'Dixi-D08-08AM',
                       'name': 'A',
                       'target': [101.98080444335938,
                                  120.58821868896484,
                                  132.00001525878906]},
                      {'entry': …....},]
     }

Development guide
==================

Debugger
--------

Use ipython -q4thread File.py

Insert "import pdb;pdb.set_trace()" in the code where you want the debugger to start.

Use Database Browser in BrainVisa to check whether files in the database are correctly
indexed. « Update » the database from there can reset the index if it was damaged by a bug.

Adding a new data type in the database
---------------------------------------

(WARNING this may change in next major brainvisa version)

As explained (Data Formats, above) to add a new data type, its type
must be declared (for a new one), its place and name in the database too.

For example, adding a file for each patient to store the list of implanted
structures.

This file will be called structures_patientName.txt and will be in
« implantation » directory of each subject in the database. We need :

-  A **file format**, here the txt format, already declared in the database.
   If not, just add your format to
   brainvisa/toolboxes/epilepsy/types/epilepsy.py with a line like::

      Format( 'PTS format', 'f|*.pts' ) # For a .pts file format
      Format ( 'Powerpoint file', ["f|*.ppt","f|*.pptx"] ) # For a format with multiple extensions

   Already known formats may have been declared in
   brainvisa-4.6.0/brainvisa/types/\*.py (e.g. \*txt) or in other toolboxes :
   brainvisa-4.6.0/brainvisa/toolboxes/\*/types/\*.py

-  A **File type** "Implanted Structures", which is not a subtype of an existing type.
   "Right Side Implanted Structures" could be a subtype of "Implanted Structures".
   It is declared in brainvisa/toolboxes/epilepsy/types/epilepsy.py::

      FileType( '**Implanted Structures'**, 'Any Type', 'Text file' )

-  A declaration inside the **database hierarchy**: e.g. in
   brainvisa/toolboxes/epilepsy/hierarchies/brainvisa-3.2.0/electrodes.py

   In this file « implantation » is inserted into the subject directory
   Just add a line in its content, as for existing data::

     "structures_<subject>", SetType('**Implanted Structures**'),

   This means that in the implantation directory, there might be a file called
   structures_PatientName.txt which is of type 'Implanted Structures'.

Accessing data from Python
---------------------------
Use Brainvisa API::

   from brainvisa import axon
   axon.initializeProcesses()
   from brainvisa.data.readdiskitem import ReadDiskItem
   from brainvisa.data.writediskitem import WriteDiskItem

Finding all files of type 'Implanted Structures' from epileptic patients. We will use
file types and attributes from the database, e.g. '{protocol}' used to defined the directory
name in the hierarchy creates a 'protocol' attribute filled by the real name of the directory::

   rdi = ReadDiskItem( 'Implanted Structures', 'Text file' ,
   requiredAttributes={'protocol':'Epilepsy'} )
   # If you know subject name, add a constraint
   rdi2 = ReadDiskItem( 'Implanted Structures', 'Text file' ,
                         requiredAttributes={'protocol':'Epilepsy',
                         'subject' :'GRE_2021_TEST'} )
   # Getting the result as a list
   implStructures = list( rdi._findValues( {}, None, False ) )

_findValues is the one used because in 2014 there was only this internal function to access all results.
Returned objects are a list of ReadDiskItem. This object gives access to the file path and its attributes::

   implS = implStructures[0]
   print 'Subject '+implS.attributes()['subject']+'. File is here: '+ implS.fullPath()

If the file type can be loaded by Anatomist (e.g. an MRI)::

   from brainvisa import anatomist
   anatomist.loadObject(implS)

Or implS.fullPath() to read the file manually.

The logic to locate where to write the data is the same: just set the attributes andn file type you want
to write and BrainVisa will generate the path for you. To make things simpler, you can use an existing
ReadDiskItem to populate the properties and attribute. For example, we have a ReadDiskItem of type T1 MRI
and want to save our 'Implanted Structures' file for the same subject. File type and diskitem are all that
is required here::

   wdi = WriteDiskItem( 'Implanted Structures', 'Text file' )
   di = wdi.findValue({'subject':'monSujet', 'protocol':'Epilepsy'} )
   di2 = wdi.findValue(diskItemT1)
   print 'Output file: ' + di.fullPath()


Referentials and geometrical transformations
--------------------------------------------

Voxels in 3D images (MRI, CT, PET...) are located in space using a
**referential**. Most software use internally a **'native' referential**,
e.g. for Anatomist it is defined as position 0,0,0 at the center of the uppermost, righmost "deepmost" voxel.
Coordinates x,y,z are the distance in mm along the voxel matrix axes.

Unfortunately this convention varies between software.

DICOM format usually defines a transformation matrix which allow to compute the
position of image voxels in a physical referential. Converting from DICOM
to Nifti usually keeps this information in the Nifti header under the
**"scanner-based"** name. SPM uses this matrix for the coordinates it shows
when using the "display" button on an image. Loading an MRI in SPM, the matrix
can be displayed::

   a=spm_vol('GRE_2021_TEST.nii');a.mat

The BrainVisa commandline tool « AimsFileInfo GRE_2021_TEST.nii » also displays
the list of transform matrices in the Nifti file header::

   'referentials' : [ 'Scanner-based anatomical coordinates' ],
   'transformations' : [ [ -0.999992, 0, 0, 90.9604, 0, -1, 0, 134.016, 0, 0, -1, 121.85, 0, 0, 0, 1 ] ],

We call this referential **Scanner-based referential**.

BrainVisa can store in its database referentials and transforms to convert
from one to the other.
A referential has a unique identifier (UUID, as all other database objects).
Transformations are declared in the database with metadata to define from which
and to which referential they are going, which allows to load them automatically
(using TransformationManager object).

**ImageImport** thus stores the native and scanner-based referentials
and the trasnform between them for all imqges when importing.
Those files (.referential, .trm) are available in the registration directory
of all images imported into the database.

Entering the location of AC and PC will add transforms to Talairach referential.

Things get more complex while registering a post-operative MRI (or
CT, PET...) to a pre-operative MRI. SPM coregister
computes a transformation matrix from
scanner-based post-MRI referential towards pre-MRI scanner-based referential.

To convert coordinates from native T1-post referential to T1-pre native referential,
three transforms must be combined:

natif post → scanner-based post → scanner-based pre → natif pre.

**Frequent issue** : imqges processed through other software might lose
the scanner-based transformation matrix.
In that cas IntrAnat may interpret header transformations in an incorrect way.

There can be two different matrices in a Nifti header.
If None or both are called 'scanner-based', IntrAnat does not know which one
to use, and it can select a one different that another software such as SPM.
Registering this image using SPM, IntrAnat would not know which referential
was used by SPM to compute its matrix, so the results will be wrong.

**MNI referential**: SPM normalize computes both a linear transform
(a matrix, sur as .trm files) and a non-linear transform in the \_sn.mat file.

BrainVisa does not process nonlinear transforms; IntrAnat converts electrode contact
coordinates to MNI by a call to matlab and SPM. \_sn.mat transforms go from
scanner-based MRI referential to MNI referential.

BrainVisa API now contains a TransformationManager, which can look for
referentials and transforms linked to objects. It was updated from our patches
to automatically load needed transforms between the images stored by IntrAnat.

ReferentialConverter is an internal tool for IntrAnat used to declare
transformations (CA-CP, Talairach, Goetz or any linear transformation)
and to convert points coordinates from one referential to another.

Computing electrode contact (plot) coordinates:

-  electrode coordinates determined through IntrAnat's locateElectrodes
   are saved in T1-pre MRI native coordinates.
-  To export to .pts, those coordinates are transformed into scanner-based
-  Those coordinates are saved in a temporary file, matlab and SPM is used
   ton convert its nonlinear transform into a vector field (y_field.nii)
-  This vector field is used to convert coordinates to MNI in the temporary file.
-  When matlab is done the outpuut file is read again from Python and the PTS file is saved with MNI coordinates.

Adding a BrainVisa process
------------------------------

To get a process visible from BrainVisa main UI, just add a file to
brainvisa-4.6.0/brainvisa/toolboxes/epilepsy/processes.

This file must follow the standard BrainVisa model :
import some files, declare a signature (its parameters), a few variables (name...),
an init function and a run function to be called by the 'run' button in BrainVisa UI.

Here is a simple exemple. For more complex processes, use existing processes in other toolboxes::

   from neuroProcesses import \*
   import shfjGlobals
   from brainvisa import anatomist
   import glob, registration
   name = 'Anatomist Show Electrode Model' # Process name in BrainVisa UI
   userLevel = 0 # level 0 for everyone, level 1 : advanced users, level 2 experts
   roles = ('viewer',) # some process have custom roles. This one is a viewer and may be called to view this data type
   def validation(): # Validata parameter values
      anatomist.validation()
   # Here only one parameter, an electrode model file to be read
   signature = Signature(
      'model', ReadDiskItem( 'Electrode Model', 'Electrode Model format' ),
   )
   # May be used to fill some parameters with default values
   def initialization( self ):
      pass
   # This will be run when clicking on the 'run' button
   def execution( self, context ):
   a = anatomist.Anatomist()
   elec = ElectrodeEditorDialog(a)
   elec.open(self.model.fullPath()) # Using the parameter from the signature
   meshes = elecDialog.getAnatomistObjects()
   w = a.createWindow('Axial')
   a.addObjects(meshes, [w,])
   return (w, elec, meshes) # Return all objects that should not be destroyed when the function ends, here all 3D objects that are displayed

Brainvisa API
-------------

BrainVisa is mainly developed at CEA Neurospin by Denis Rivière, Yann Cointepas.
A lot of documentation is available online:

http://brainvisa.info

For python bindings, pyaims (3D coords, meshes, volumes...),
pyanatomist (load and display images and volumes in Anatomist).


Interacting with Matlab
------------------------

Function defined in externalprocesses.py exist to call matlab.

Simplest way is to write matlab code in a python string with %s for parameters
and ending by "quit;" as in ImageImport.py::

   from externalprocesses import \*
   spm_coregister = "VF=spm_vol(%s);VG=spm_vol(%s);x = spm_coreg(VF, VG);\\
   trm = spm_matrix(x(:)');trm = [trm(1:3,4)';trm(1:3,1:3)];\\
   dlmwrite(%s,trm, 'delimiter',' ','precision',16); \\
   quit;"

Set the parameters :

   call = spm_coregister%("'monFichier.img,1'", "'AutreFichier.img,1'","'fichierOutput'")

Execute it::

   matlabRun(call)

This is a **blocking function**, so the software will be locked until matlab computation is finished.

A **non-blocking call** is also available and creates a Qthread object.
Connecting to its end of execution signal allows to call a function. The object
must be stored to avoid destroyed the thread before it is finished, and the start function is used to start it::

   thr = matlabRunNB(call)
   thr.finished.connect(lambda:self.taskfinished(u"SPM Coregister ended", thr))
   self.threads.append(thr)
   thr.start()

Tu create a **temporary file** for matlab to write, one can be created::
   tempfile = getTmpFilePath('txt')

It must be **removed manually** after execution.

Parallel execution
--------------------

BrainVisa ships with soma-workflow that can manage parallel execution.
IntrAnat does not use it yet.
http://brainvisa.info/soma/soma-workflow/index.html
