IntrAnat installation for BrainVisa 5
=====================================

Requirements
------------

IntrAnat requires the installation of:

-  **BrainVISA**: Version >5.0.4
-  **Matlab**: Requires a Matlab license to run SPM12 for MNI
   normalization
-  **SPM12**: For MNI normalization and volume co-registration
-  **ANTs 2.2**: For volume co-registration (optional)
-  **FreeSurfer >6.0**: For importing the result of the FreeSurfer T1
   MRI segmentation pipeline (optional)

Supported operating systems:

-  **Ubuntu 20.04**: Reference for all the instructions below
-  **Other Linux distributions**: Should work but haven’t been tested,
   the main limiting factor for portability being Singularity (see
   below)
-  **Windows 10/WSL2**: `Linux subsystem on Windows
   10 <https://docs.microsoft.com/en-us/windows/wsl/install-win10>`__

   -  Should works with Ubuntu 20.04 installation,
   -  WSL can run only the latest versions of Matlab (>= 2018b), but
      Matlab can be installed on the Windows system and called from
      WSL/IntrAnat

Choose installation directory

Install Ubuntu packages
-----------------------

Instructions for Ubuntu 20.04.

Update your system:

.. code:: bash

   sudo apt update
   sudo apt dist-upgrade

FreeSurfer dependencies:

.. code:: bash

   sudo apt install build-essential libjpeg62 libxss1 libgomp1 tcsh bc

ANTS dependencies:

.. code:: bash

   sudo apt install cmake-curses-gui gcc g++ zlib1g-dev

BrainVISA/IntrAnat:

.. code:: bash

   sudo apt install libtinfo-dev libapt-pkg-dev git libav-tools mencoder libglm-dev 

Installation directory
----------------------

Choose and create installation directory, set environment variable:

.. code:: bash

   mkdir $HOME/IntrAnat
   export INTRANAT_INSTALL=$HOME/IntrAnat

Matlab and SPM installation
---------------------------

Install any version of Matlab. Make sure it is in the system PATH.

Install SPM12:

.. code:: bash

   cd $INTRANAT_INSTALL
   wget http://www.fil.ion.ucl.ac.uk/spm/download/restricted/eldorado/spm12.zip
   unzip spm12.zip
   rm spm12.zip

Matlab on Windows 10/WSL
~~~~~~~~~~~~~~~~~~~~~~~~

You can either install MATLAB in WSL/Ubuntu, or install the Windows
version and create a link to matlab.exe inside the Ubuntu PATH. For
example:

.. code:: bash

   sudo -s
   printf '#!/bin/bash\n/mnt/c/Program\ Files/MATLAB/R2017b/bin/matlab.exe -nodesktop -wait "$@"\nexit $?' > /usr/local/bin/matlab 
   chmod a+x /usr/local/bin/matlab

Exchanging data between the Ubuntu and Windows drives:

-  Everything must be on the Linux drive: using a network drive
   connected via SSH to access the files in Matlab
-  ``sudo vi /etc/ssh/sshd_config``
-  Change the port to 2222 and edit all the options as in `bash - How to
   SSH into WSL from Windows on the same machine - Super
   User <https://superuser.com/questions/1123552/how-to-ssh-into-wsl>`__
-  ``sudo service ssh start``
-  On Windows: Install SFTP Net Drive:
   https://www.nsoftware.com/netdrive/sftp/

   -  Connect to the drive:

      -  Server: 127.0.0.1:2222
      -  Username/password: The authentication of your Ubuntu user
      -  Drive letter: “L:”

Freesurfer installation (optional)
----------------------------------

Download FreeSurfer 7.2.0:

.. code:: bash

   cd $INTRANAT_INSTALL
   wget -O freesurfer.tgz https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/7.2.0/freesurfer-linux-ubuntu18_amd64-7.2.0.tar.gz
   tar zxvf freesurfer.tgz

Add the FreeSurfer configuration to your .bashrc:

.. code:: bash

   grep -q -F "FREESURFER_HOME=" ~/.bashrc || printf "\n# FREESURFER\nexport FREESURFER_HOME=$INTRANAT_INSTALL/freesurfer\nsource \$FREESURFER_HOME/SetUpFreeSurfer.sh\n" >> ~/.bashrc
   rm freesurfer.tgz

Get a license file for FreeSurfer
(http://surfer.nmr.mgh.harvard.edu/registration.html) and save it in
``$INTRANAT_INSTALL/freesurfer/license.txt`` Example:

::

   printf "francois.tadel@univ-grenoble-alpes.fr\n34309\n *COj3JXOXnbes\n FSCI/SECcEOfM" > $INTRANAT_INSTALL/freesurfer/license.txt

Install ANTs (optional)
-----------------------

Get and compile ANTs >= 2.2.0:

.. code:: bash

   cd $INTRANAT_INSTALL
   wget -O ANTs.tgz https://github.com/stnava/ANTs/tarball/master
   tar zxvf ANTs.tgz
   rm ANTs.tgz
   mv ANTsX-ANTs* ANTs
   cd ANTs
   mkdir build
   cd build
   ccmake ../

-  Press “c” to configure, then “c” again

-  If no errors, press “g” to generate the make files

-  Full compilation instructions: `Compiling ANTs on Linux and Mac OS ·
   ANTsX/ANTs Wiki ·
   GitHub <https://github.com/ANTsX/ANTs/wiki/Compiling-ANTs-on-Linux-and-Mac-OS>`__

-  Expect the compilation to last for a few hours

   ::

      make

Add the ANTs configuration to your .bashrc:

.. code:: bash

   grep -q -F "ANTSPATH=" ~/.bashrc || printf "\n# ANTs\nexport ANTSPATH=$INTRANAT_INSTALL/ANTs/build/bin/\nexport PATH=$INTRANAT_INSTALL/ANTs/Scripts:\$ANTSPATH:\$PATH\nexport ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=4\n" >> ~/.bashrc

BrainVisa installation
----------------------

Download BrainVisa
~~~~~~~~~~~~~~~~~~

Download the BrainVisa singularity image from its `official
website <https://brainvisa.info/web/download.html>`__.

At the time of this writing the latest image is
`brainvisa-5.0.4.sif <https://brainvisa.info/download/brainvisa-5.0.4.sif>`__

Save it into **$INTRANAT_INSTALL**

To run the image, you need to install singularity. BrainVisa provides a
package for Ubuntu and other linux distributions on the same webpage
that you should `download
here <https://brainvisa.info/download/singularity-ce_3.8.3~ubuntu-20.04_amd64.deb>`__.
If you are using Mac OS X, please follow the instructions on Brainvisa’s
web page.

Install BrainVisa
~~~~~~~~~~~~~~~~~

First, install Singularity by double-clicking on the .deb package or by
running

.. code:: bash

   sudo dpkg -i singularity-container-*.deb

You can now run the setup of BrainVisa and add it to your PATH:

.. code:: bash

   mkdir $INTRANAT_INSTALL/brainvisa-5.0.4
   singularity run -B $INTRANAT_INSTALL/brainvisa-5.0.4:/casa/setup $INTRANAT_INSTALL/brainvisa-5.0.4.sif
   echo "export INTRANAT_INSTALL=$HOME/IntrAnat" >> ~/.bashrc
   echo "export PATH=$INTRANAT_INSTALL/brainvisa-5.0.4/bin:$PATH" >> ~/.bashrc

Now that Brainvisa installation is complete, setup its database.

Create an empty directory for the database.

.. code:: bash

   mkdir $INTRANAT_INSTALL/brainvisa_db

Set the database and program paths in BrainVISA:

-  Start BrainVISA:

   .. code:: bash

      brainvisa

-  In the BrainVISA interface:

   -  If the window “Update databases” is open, click on the button
      “Update”, then close the update window.
   -  If you get the window “Welcome to BrainVISA” click on button “Open
      preferences”, otherwise select menu “BrainVISA > Preferences”

      -  In the “Databases” section, click on “Add”, then select folder
         $INTRANAT_INSTALL/brainvisa_db, then click “OK”
      -  Click again on “Add”, select Freesurfer subjects folder
         ``$INTRANAT_INSTALL/freesurfer/subjects``, click on advanced
         settings and select ontology “freesurfer”. This will add
         Freesurfer database to Brainvisa.
      -  Stil in Preferences: Set SPM path to
         ``$INTRANAT_INSTALL/spm12`` and Freesurfer path to
         ``$INTRANAT_INSTALL/freesurfer``

-  Close preferences

-  Close BrainVISA

Install IntrAnat
================

Setting container to read write

.. code:: bash

   bv

In the window, click on the button to set the container to
**read-write**. You need to do this once.

Then select Terminal to run the next commands inside the terminal. You
can also run directly bash inside the container with the following
command:

.. code:: bash

   bv bash

Define installation directory:

.. code:: bash

   export INTRANAT_INSTALL=/home/USERNAME/IntrAnat

Clone the GitHub repository:

.. code:: bash

   cd $INTRANAT_INSTALL
   git clone https://github.com/IntrAnatSEEGSoftware/IntrAnat

Checkout brainvisa_5.0.0dev branch in the git repository:

.. code:: bash

   cd $INTRANAT_INSTALL/IntrAnat/
   git checkout --track origin/brainvisa_5.0.0dev

Set up epilepsy toolbox in BrainVisa:

.. code:: bash

   ln -s  $INTRANAT_INSTALL/IntrAnat/epilepsy-toolbox/ /casa/host/install/brainvisa/toolboxes/epilepsy

Add the electrode models to the shared database:

.. code:: bash

   ln -s $INTRANAT_INSTALL/IntrAnat/electrode_models/ /casa/host/install/share/brainvisa-share-5.0
  
Start Brainvisa, set up and update the databases: run the brainvisa command, in the main window select Data management, Update databases. Make sure that the shared database and your own database are checked. Then run the process, wait until brainvisa has finished indexing the electrode models and your own data.

**Matlab license**. Inside the container the home directory is not your
user directory. To have matlab running, you need to have its license
files in the container’s home, and if matlab command is not in the path
you need to add it to the .bashrc file in the container’s home:

.. code:: bash

   ln -s /home/USERNAME/.matlab ~/
   # If matlab command is not in the PATH
   echo "export PATH=$PATH:/PATH-TO-MATLAB/" >> ~/.bashrc

==================

Start IntrAnat
==============

Create startup scripts:

.. code:: bash

   cd $INTRANAT_INSTALL
   printf "#!/bin/bash\nbv bash -c \"cd $INTRANAT_INSTALL/IntrAnat\npython ImageImport.py\"\n" > ImageImport.sh
   printf "#!/bin/bash\nsource $INTRANAT_INSTALL/brainvisa-4.6.1/bin/bv_env.sh $INTRANAT_INSTALL/brainvisa-4.6.1\ncd IntrAnat\npython locateElectrodes.py" > locateElectrodes.sh
   printf "#!/bin/bash\nsource $INTRANAT_INSTALL/brainvisa-4.6.1/bin/bv_env.sh $INTRANAT_INSTALL/brainvisa-4.6.1\ncd IntrAnat\npython groupDisplay.py" > groupDisplay.sh
   chmod a+x *.sh

Manual execution:

.. code:: bash

   cd ~/IntrAnat/IntrAnat
   bv python ImageImport.py

Or all in one line:

.. code:: bash

   cd ~/IntrAnat/IntrAnat && bv python ImageImport.py

Set program paths:

-  Open ImageImport, go to the tab “Preferences”
-  Set path to SPM12: ``$INTRANAT_INSTALL/spm12``
-  Set path to ANTs: ``$INTRANAT_INSTALL/ANTs-build``
-  Set path to FreeSurfer: ``$INTRANAT_INSTALL/freesurfer`` (should be
   set automatically if the FreeSurfer path is properly set in the
   BrainVISApreferences )
-  Set path to a BIDS database (for BIDS export)
-  Click on button “Save preferences”

Update IntrAnat from GitHub:

.. code:: bash

   cd ~/IntrAnat/IntrAnat
   git pull

Install MRIConvert(optional)
============================

MRIConvert is not needed for running IntrAnat, but is a very useful tool
for converting DICOM images into .nii files.
https://lcni.uoregon.edu/downloads/mriconvert

.. code:: bash

   cd $INTRANAT_INSTALL
   wget -O MRIConvert.tgz https://lcni.uoregon.edu/downloads/mriconvert/MRIConvert-2.1.0-x86_64-rhel.tar.gz/at_download/file
   tar zxvf MRIConvert.tgz
   rm MRIConvert.tgz
   cd MRIConvert-*
   chmod a+x install.sh
   sudo ./install.sh

Lausanne2008 parcellation(optional)
===================================

These scripts are not publicly available yet…

Install FSL:

.. code:: bash

   sudo apt-get install neurodebian
   sudo apt-get update
   sudo apt-get install fsl-complete 
   sudo pip install nypipe
   sudo pip install nibabel
   sudo pip install networkx==1.11

Edit .bashrc, add at the end:

.. code:: bash

   source /usr/share/fsl/5.0/etc/fslconf/fsl.sh
