#!/bin/bash
cd /Users/admin/Softwares/BrainVISA-4.5/bin/
. /Users/admin/Softwares/BrainVISA-4.5/bin/bv_env.sh /Users/admin/Softwares/BrainVISA-4.5
cd /Users/admin/Documents/GIT/IntrAnatElectrodes/epilepsie/
echo "#############################################################################" >> intranat-locateElectrodes-`whoami`.log  
date >> intranat-locateElectrodes-`whoami`.log
python locateElectrodes.py >> intranat-locateElectrodes-`whoami`.log 

