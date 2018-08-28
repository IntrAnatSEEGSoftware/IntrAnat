# -*- coding: utf-8 -*-
#  This software and supporting documentation are distributed by
#      Institut des Neurosciences de Grenoble - INSERM U836
#      France
#
# This software is governed by the CeCILL license version 2 under
# French law and abiding by the rules of distribution of free software.
# You can  use, modify and/or redistribute the software under the 
# terms of the CeCILL license version 2 as circulated by CEA, CNRS
# and INRIA at the following URL "http://www.cecill.info". 
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license version 2 and that you accept its terms.

include( 'base' )
include( 'anatomy' )


# SEEG DATA
seeg_content = (

  "<subject>_{experiment}_{expNumber} ", SetType( 'SEEG Experiment' ),
     SetContent(
       "<subject>_<experiment>_{subexperiment}__{expId}", SetType( 'Raw SEEG recording' ), 
       "<subject>_<experiment>__{expId}", SetType( 'Raw SEEG recording' ), 
       "<subject>_<experiment>_{subexperiment}", SetType( 'Raw SEEG recording' ), 
       "<subject>_<experiment>", SetType( 'Raw SEEG recording' ), 
       "<seegAcq>", SetType( 'Raw SEEG recording' ),
       "{seegProcessing}", SetType( 'SEEG processing directory' ),
     ),
  "<subject>_{experiment}", SetType( 'SEEG Experiment' ),
     SetContent(
       "<subject>_<experiment>_{subexperiment}__{expId}", SetType( 'Raw SEEG recording' ), 
       "<subject>_<experiment>__{expId}", SetType( 'Raw SEEG recording' ), 
       "<subject>_<experiment>_{subexperiment}", SetType( 'Raw SEEG recording' ), 
       "<subject>_<experiment>", SetType( 'Raw SEEG recording' ), 
       "<seegAcq>", SetType( 'Raw SEEG recording' ),
       "{seegProcessing}", SetType( 'SEEG processing directory' ),
       #"{seegAcq}", SetType( 'Elan EEG' ), 
       #"{seegAcq}", SetType( 'Trigger'),
       #"{localizer}", SetType ('JP Localizer'),
       #"*_trials", SetType ('SEEG Trials Directory' ),
         #SetContent(
           #"*", SetType ('JP image'),
   #),
       #"*", SetType ('BrainTV Film'),
       #"*", SetType ('JP image'),
     ),

  #"indexepi_{date}_{location}_{experiment}_{expNumber}", SetDefaultAttributeValue( 'expNumber', 1 ), SetNonMandatoryKeyAttribute( 'expNumber' ),
     #SetType( 'Epilepsy Index SEEG Experiment' ), SetWeakAttr( 'experimentType', 'epilepsy index' ),
     #SetContent(
       #"EEG_{seegAcq}", SetType( 'Raw SEEG recording' ),
       #"EEG_{seegAcq}", SetType( 'ImaGIN matlab files' ),
       #"baselineEEG_{seegAcq}", SetType( 'ImaGIN matlab files' ),
       #"EI_EEG_{seegAcq}__{freqMin}_{freqMax}_{window}_{onset}", SetType( 'ImaGIN matlab files' ),
       #"EI_Group__{freqMin}_{freqMax}_{window}_{onset}", SetType( 'ImaGIN matlab files' ),
       #"SPM_EI_EEG_{seegAcq}__{freqMin}_{freqMax}_{window}_{onset}", SetType( 'EpiIndex Directory' ),
         #SetContent(
           #"*", SetType( '4D Volume' ),
           #"*", SetType( 'ImaGIN matlab files' ),
	 #),
       #"Results", SetType( 'Results Directory' ), SetContent( "*", SetType( 'PDF Report' ), )
     #),

  #"stim_{date}_{location}_{experiment}_{expNumber}", SetDefaultAttributeValue( 'expNumber', 1 ), SetNonMandatoryKeyAttribute( 'expNumber' ),
     #SetType( 'Stimulation SEEG Experiment' ), SetWeakAttr( 'experimentType', 'stimulation' ),
     #SetContent(
       #"EEG_{seegAcq}", SetType( 'Raw SEEG recording' ),
       #"EEG_{seegAcq}", SetType( 'ImaGIN matlab files' ),
     #),
)

# Microelectrode recordings data
mer_content = (
  "<subject>_{experiment}_{expNumber}", SetType('MER Experiment'), SetContent(
       "<subject>_<experiment>_{subexperiment}__{expId}", SetType('Raw MER recording'),
       "<subject>_<experiment>__{expId}", SetType('Raw MER recording'),
       "<subject>_<experiment>_{subexperiment}", SetType('Raw MER recording'),
       "<subject>_<experiment>", SetType('Raw MER recording'),
  ),
  "<subject>_{experiment}", SetType('MER Experiment'), SetContent(
       "<subject>_<experiment>_{subexperiment}__{expId}", SetType('Raw MER recording'),
       "<subject>_<experiment>__{expId}", SetType('Raw MER recording'),
       "<subject>_<experiment>_{subexperiment}", SetType('Raw MER recording'),
       "<subject>_<experiment>", SetType('Raw MER recording'),
  ),
)

# Eye movements
eye_content = (
  "<subject>_{experiment}_{expNumber}", SetType('Eye Tracking Experiment'), SetContent(
     "<subject>_<experiment>_{subexperiment}__{expId}", SetType('Raw Eye Tracking'),
     "<subject>_<experiment>__{expId}", SetType('Raw Eye Tracking'),
     "<subject>_<experiment>_{subexperiment}", SetType('Raw Eye Tracking'),
     "<subject>_<experiment>", SetType('Raw Eye Tracking'),
  ),
  "<subject>_{experiment}", SetType('Eye Tracking Experiment'), SetContent(
     "<subject>_<experiment>_{subexperiment}__{expId}", SetType('Raw Eye Tracking'),
     "<subject>_<experiment>__{expId}", SetType('Raw Eye Tracking'),
     "<subject>_<experiment>_{subexperiment}", SetType('Raw Eye Tracking'),
     "<subject>_<experiment>", SetType('Raw Eye Tracking'),
  ),
)

# Including in the main database hierarchy
insert( '{center}/{subject}',
  'seeg', SetWeakAttr( 'modality', 'seeg' ),
    apply( SetContent, seeg_content)
)

insert( '{center}/{subject}',
  'mer', SetWeakAttr( 'modality', 'mer' ),
    apply( SetContent, mer_content)
)

insert( '{center}/{subject}',
  'eye', SetWeakAttr( 'modality', 'eye' ),
    apply( SetContent, eye_content)
)

insert( '{center}/{subject}', 'seeg', SetWeakAttr( 'modality', 'seeg' ), SetContent("<subject>", SetType('Electrodes SEEG Labels'),SetWeakAttr( 'sEEG Labels', 'default' )))


