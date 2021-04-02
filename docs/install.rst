Installation
***************

Intranat depends on multiple software. The installation process is therefore not completely
straightforward.

Requirements
============

IntrAnat requires the installation of:

- **BrainVISA**: Version > 4.6.1
- **Matlab**: Requires a Matlab license to run SPM12 for MNI normalization
- **SPM12**: For MNI normalization and volume co-registration
- **ANTs (>2.2)**: For volume co-registration (optional)
- **FreeSurfer (> 6.0)**: For importing the result of the FreeSurfer T1 MRI segmentation pipeline (optional)

Supported operating systems:

- **Ubuntu 16.04 (or newer)**: Reference for all the instructions below
- **CentOS 7**: For CentOS/Fedora/RHEL, similar packages are available for install using "sudo yum install" instead of "sudo apt install", but their name might have to be updated.
- **Other Linux distributions**: Might work but haven't been tested yet, the main limiting factor for portability being BrainVISA
- **Windows 10/WSL**: `Linux subsystem on Windows 10 / Ubuntu 16.04 <https://www.howtogeek.com/249966/how-to-install-and-use-the-linux-bash-shell-on-windows-10/>`__ Almost everything works just as with a Ubuntu 16.04, WSL can run only the latest versions of Matlab (>= 2018b), but Matlab can be installed on the Windows system and called from WSL/IntrAnat

Install Ubuntu packages
=======================

Instructions for Ubuntu 16.04 (should work for 18.04 and 20.04).

Update your system:

::

    sudo apt update
    sudo apt dist-upgrade

FreeSurfer dependencies:

::

    sudo apt install build-essential libjpeg62 libxss1 libgomp1 tcsh bc

ANTS dependencies:

::

    sudo apt install cmake-curses-gui gcc g++ zlib1g-dev

BrainVISA/IntrAnat:

::

    sudo apt install libtinfo-dev libapt-pkg-dev git libav-tools mencoder libglm-dev
    # Next instructions are not needed for ubuntu >= 18.04
    wget http://ftp.br.debian.org/debian/pool/main/a/apt/libapt-pkg4.12_1.0.9.8.4_amd64.deb
    dpkg -i libapt-pkg4.12_1.0.9.8.4_amd64.deb
    rm libapt-pkg4.12_1.0.9.8.4_amd64.deb

Install IntrAnat
================

Define installation directory:

::

    mkdir $HOME/IntrAnat
    export INTRANAT_INSTALL=$HOME/IntrAnat

Clone the GitHub repository:

::

    cd $INTRANAT_INSTALL
    git clone https://github.com/ftadel/IntrAnat

Checkout brainvisa\_4.6 branch in the git repository:

::

    cd $INTRANAT_INSTALL/IntrAnat/
    git checkout --track origin/brainvisa_4.6

Download BrainVISA
==================

The current official version of BrainVISA (4.6.1) contains bugs in
Anatomist that make IntrAnat interface unusable. We distribute
temporarily a patched version of BrainVISA 4.6.2 on an alternate
website, until the BrainVISA developers release a new version including
these bug fixes.

-  Linux Ubuntu16 (glibc 2.23):
   https://insermfrance-my.sharepoint.com/:u:/g/personal/olivier\_david\_inserm\_eu/EQnoFTBdrFhJi8mjMoy4VyUBptlS\_7yWlMrwucAfJwVhMA
-  Other Linux distributions (eg. CentOS 7):
   https://insermfrance-my.sharepoint.com/:u:/g/personal/olivier\_david\_inserm\_eu/EfGuisNIOaRGnV5EdWgPvooB2mvnNPp1CNNvEcwgb\_EWpQ
-  Windows 10/WSL/Ubuntu 16: Read the instructions on the BrainVISA
   installation page to prepare your system:
   http://brainvisa.info/web/download.html#running-on-windows-10-ubuntu-shell

A last option is to compile BrainVISA with bv\_maker on your system.

Install BrainVISA
=================

Start the installer:

::

    chmod u+x brainvisa-4.6-install
    ./brainvisa-4.6-install

-  Set the installation to: /home/USERNAME/IntrAnat/brainvisa-4.6.1.
-  Use the default options

Set up epilepsy toolbox:

::

    ln -s $INTRANAT_INSTALL/IntrAnat/epilepsy-toolbox $INTRANAT_INSTALL/brainvisa-4.6.1/brainvisa/toolboxes/epilepsy

Add the electrode models:

::

    ln -s $INTRANAT_INSTALL/IntrAnat/electrode_models $INTRANAT_INSTALL/brainvisa-4.6.1/share/brainvisa-share-4.6/electrode_models

Start Brainvisa, set up and update the databases:

::

    mkdir $INTRANAT_INSTALL/brainvisa_db
    $INTRANAT_INSTALL/brainvisa-4.6.1/BrainVISA

In the BrainVISA interface:

* In the window "Update databases", click on the button "Update", then close the figure
* If you get the window "Welcome to BrainVISA" click on button "Open preferences", otherwise select menu "BrainVISA > Preferences"
* In the "Databases" section, click on "Add", then select folder $INTRANAT\_INSTALL/brainvisa\_db, then click "OK"
* Close BrainVISA
* Add your freesurfer database: Read the help in BrainVISA's FreeSurfer toolbox

Install additional packages in BrainVISA's Python environment if
necessary:

* openpyxl:

::

    cd $INTRANAT_INSTALL/
    source $INTRANAT_INSTALL/brainvisa-4.6.1/bin/bv_env.sh $INTRANAT_INSTALL/brainvisa-4.6.1
    wget https://files.pythonhosted.org/packages/d9/dd/5952829956827de7ff36eb70877fdffd6dbfacb670fae05eb7ccba52ace7/openpyxl-2.5.5.tar.gz
    tar zxvf openpyxl-2.5.5.tar.gz
    rm openpyxl-2.5.5.tar.gz
    cd openpyxl-2.5.5/
    python setup.py install

*  jdcal
*  et\_xmlfile

Specific instructions for Windows 10/WSL:

* Delete brainvisa/lib/libxcb
* to avoid the errors "libxcb-dri3.so.0: undefined symbol: xcb\_send\_fd"
* Delete additional duplicated libraries
* ``cd $INTRANAT_INSTALL/brainvisa-4.6.1/lib``
* ``rm libxcb* libgcc_s* libpcre* libstdc++* libtinfo* libdl* libz*``

Specific instructions for Mandriva2008:
* ``rm libgcc_s* libstdc++* libdl* libz*``
* Install in brainvisa Python environment: jdcal, et\_xmlfile, openpyxl

Install FreeSurfer
==================

Download FreeSurfer 6.0 (or more recent):

::

    cd $INTRANAT_INSTALL
    wget -O freesurfer.tgz ftp://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/6.0.0/freesurfer-Linux-centos6_x86_64-stable-pub-v6.0.0.tar.gz
    tar zxvf freesurfer.tgz

Add the FreeSurfer configuration to your .bashrc:

::

    grep -q -F "FREESURFER_HOME=" ~/.bashrc || printf "\n# FREESURFER\nexport FREESURFER_HOME=$INTRANAT_INSTALL/freesurfer\nsource \$FREESURFER_HOME/SetUpFreeSurfer.sh\n" >> ~/.bashrc
    rm freesurfer.tgz

Get a license file for FreeSurfer
(http://surfer.nmr.mgh.harvard.edu/registration.html) and save it in
``$INTRANAT_INSTALL/freesurfer/license.txt`` Example:

::

    printf "francois.tadel@univ-grenoble-alpes.fr\n34309\n *COj3JXOXnbes\n FSCI/SECcEOfM" > $INTRANAT_INSTALL/freesurfer/license.txt

Install ANTs
============

Get and compile ANTs >= 2.2.0:

::

    cd $INTRANAT_INSTALL
    wget -O ANTs.tgz https://github.com/stnava/ANTs/tarball/master
    tar zxvf ANTs.tgz
    rm ANTs.tgz
    mv ANTsX-ANTs* ANTs
    cd ANTs
    mkdir build
    cd build
    ccmake ../

-  Press "c" to configure, then "c" again
-  If no errors, press "g" to generate the make files
-  Full compilation instructions:
   https://github.com/ANTsX/ANTs/wiki/Compiling-ANTs-on-Linux-and-Mac-OS
-  Expect the compilation to last for a few hours

   ::

       make

Add the ANTs configuration to your .bashrc:

::

    grep -q -F "ANTSPATH=" ~/.bashrc || printf "\n# ANTs\nexport ANTSPATH=$INTRANAT_INSTALL/ANTs/build/bin/\nexport PATH=$INTRANAT_INSTALL/ANTs/Scripts:\$ANTSPATH:\$PATH\nexport ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=4\n" >> ~/.bashrc

Install MATLAB
==============

Install any version of Matlab. Make sure it is in the system PATH.

Install SPM12:

::

    cd $INTRANAT_INSTALL
    wget http://www.fil.ion.ucl.ac.uk/spm/download/restricted/eldorado/spm12.zip
    unzip spm12.zip
    rm spm12.zip

Set the program paths in BrainVISA:

* Start BrainVISA: ``$INTRANAT_INSTALL/brainvisa-4.6.1/BrainVISA``
* Open menu BrainVISA > Preferences: Set SPM path to ``$INTRANAT_INSTALL/spm12``
* Close BrainVISA

Matlab on Windows 10/WSL
------------------------

You can either install MATLAB in WSL/Ubuntu, or install the Windows
version and create a link to matlab.exe inside the Ubuntu PATH. For
example:

::

    sudo -s
    printf '#!/bin/bash\n/mnt/c/Program\ Files/MATLAB/R2017b/bin/matlab.exe -nodesktop -wait "$@"\nexit $?' > /usr/local/bin/matlab
    chmod a+x /usr/local/bin/matlab

Exchanging data between the Ubuntu and Windows drives:

* Everything must be on the Linux drive: using a network drive connected via SSH to access the files in Matlab
* ``sudo vi /etc/ssh/sshd_config``
* Change the port to 2222 and edit all the options as in https://superuser.com/questions/1123552/how-to-ssh-into-wsl
* ``sudo service ssh start``
* On Windows: Install SFTP Net Drive: https://www.nsoftware.com/netdrive/sftp/
* Connect to the drive:
   * Server: 127.0.0.1:2222
   * Username/password: The authentication of your Ubuntu user
   * Drive letter: "L:"

Install MRIConvert
==================

MRIConvert is not needed to run IntrAnat, but is a very useful tool for converting DICOM images into .nii files used by IntrAnat.
https://lcni.uoregon.edu/downloads/mriconvert

::

    cd $INTRANAT_INSTALL
    wget -O MRIConvert.tgz https://lcni.uoregon.edu/downloads/mriconvert/MRIConvert-2.1.0-x86_64-rhel.tar.gz/at_download/file
    tar zxvf MRIConvert.tgz
    rm MRIConvert.tgz
    cd MRIConvert-*
    chmod a+x install.sh
    sudo ./install.sh

Lausanne2008 parcellation
=========================

These scripts are not publicly available yet...

Install FSL:

::

    sudo apt-get install neurodebian
    sudo apt-get update
    sudo apt-get install fsl-complete
    sudo pip install nypipe
    sudo pip install nibabel
    sudo pip install networkx==1.11

Edit .bashrc, add at the end:

::

    source /usr/share/fsl/5.0/etc/fslconf/fsl.sh




Running IntrAnat
================

Create startup scripts:

::

    cd $INTRANAT_INSTALL
    printf "#!/bin/bash\nsource $INTRANAT_INSTALL/brainvisa-4.6.1/bin/bv_env.sh $INTRANAT_INSTALL/brainvisa-4.6.1\ncd IntrAnat\npython ImageImport.py" > ImageImport.sh
    printf "#!/bin/bash\nsource $INTRANAT_INSTALL/brainvisa-4.6.1/bin/bv_env.sh $INTRANAT_INSTALL/brainvisa-4.6.1\ncd IntrAnat\npython locateElectrodes.py" > locateElectrodes.sh
    printf "#!/bin/bash\nsource $INTRANAT_INSTALL/brainvisa-4.6.1/bin/bv_env.sh $INTRANAT_INSTALL/brainvisa-4.6.1\ncd IntrAnat\npython groupDisplay.py" > groupDisplay.sh
    chmod a+x *.sh

Manual execution:

::

    cd ~/IntrAnat
    source brainvisa-4.6.1/bin/bv_env.sh
    cd IntrAnat
    python ImageImport.py

Or all in one line:

::

    cd ~/IntrAnat && source brainvisa-4.6.1/bin/bv_env.sh && cd IntrAnat && python ImageImport.py

Set program paths:

* Open ImageImport, go to the tab "Preferences"
* Set path to SPM12: ``$INTRANAT_INSTALL/spm12``
* Set path to ANTs: ``$INTRANAT_INSTALL/ANTs-build``
* Set path to FreeSurfer: ``$INTRANAT_INSTALL/freesurfer`` (should be set automatically if the FreeSurfer path is properly set in the BrainVISApreferences )
* Click on button "Save preferences"

Update IntrAnat from GitHub:

::

    cd ~/IntrAnat/IntrAnat
    git pull


