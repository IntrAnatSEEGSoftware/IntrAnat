#!/bin/bash
cd /home/lnppeeg/software/brainvisa-4.3.0/bin/
. /home/lnppeeg/software/brainvisa-4.3.0/bin/bv_env.sh
cd /home/lnppeeg/prog/electrophysiology/epilepsie/
echo "#############################################################################" >> intranat-locateElectrodes-`whoami`.log  
date >> intranat-locateElectrodes-`whoami`.log
python groupDisplay.py >> intranat-locateElectrodes-`whoami`.log 

