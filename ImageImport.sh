#!/bin/bash
cd /Users/admin/Softwares/BrainVISA-4.5/bin/
. /Users/admin/Softwares/BrainVISA-4.5/bin/bv_env.sh /Users/admin/Softwares/BrainVISA-4.5
cd /Users/admin/Documents/GIT/IntrAnatElectrodes/epilepsie/
echo "#############################################################################" >> intranat-imageImport-`whoami`.log  
date >> intranat-imageImport-`whoami`.log
python ImageImport.py  >> intranat-imageImport-`whoami`.log

