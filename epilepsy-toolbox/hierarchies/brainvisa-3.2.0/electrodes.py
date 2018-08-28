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

insert( '{center}/{subject}','implantation', SetType('Electrode Implantation Directory'), SetContent(
  "<subject>_{implantation_session}", SetType('Electrode Implantation'), # In case there is more than one
  "<subject>", SetType('Electrode Implantation'),  SetWeakAttr( 'implantation_session', 'default' ),
  "<subject>_{implantation_session}", SetType('Electrodes Labels'), # In case there is more than one
  "<subject>", SetType('Electrodes Labels'),  SetWeakAttr( 'implantation_session', 'default' ),
  "<subject>", SetType('Electrode Implantation PTS'), SetWeakAttr( 'no_ref_name', 'True', 'implantation_session', 'default' ),
  "<subject>_{ref_name}", SetType('Electrode Implantation PTS'),  SetWeakAttr( 'implantation_session', 'default' ),
  "<subject>_{ref_name}_Pos", SetType('Electrode Implantation Position TXT'),  SetWeakAttr( 'implantation_session', 'default' ),
  "<subject>_{ref_name}_Name", SetType('Electrode Implantation Name TXT'), SetWeakAttr( 'implantation_session', 'default' ),
  "<subject>_Pos", SetType('Electrode Implantation Position TXT'), SetWeakAttr( 'no_ref_name', 'True', 'implantation_session', 'default' ),
  "<subject>_Name", SetType('Electrode Implantation Name TXT'), SetWeakAttr( 'no_ref_name', 'True', 'implantation_session', 'default' ),
  "<subject>_{implantation_session}", SetType('Electrode Implantation PTS'), SetWeakAttr( 'no_ref_name', 'True' ),
  "<subject>_{implantation_session}_{ref_name}", SetType('Electrode Implantation PTS'), 
  "<subject>_{implantation_session}_{ref_name}_Pos", SetType('Electrode Implantation Position TXT'), 
  "<subject>_{implantation_session}_{ref_name}_Name", SetType('Electrode Implantation Name TXT'),
  "<subject>_{implantation_session}_Pos", SetType('Electrode Implantation Position TXT'), SetWeakAttr( 'no_ref_name', 'True' ),
  "<subject>_{implantation_session}_Name", SetType('Electrode Implantation Name TXT'), SetWeakAttr( 'no_ref_name', 'True' ),
  "<subject>_sag", SetType('Electrode Implantation Sagittal Image'), SetWeakAttr( 'implantation_session', 'default' ),
  "<subject>_coro", SetType('Electrode Implantation Coronal Image'), SetWeakAttr( 'implantation_session', 'default' ),
  "<subject>", SetType('Electrode Implantation Powerpoint report'), SetWeakAttr( 'implantation_session', 'default' ),
  "<subject>", SetType('Electrode Implantation PDF report'), SetWeakAttr( 'implantation_session', 'default', 'planning', 'True' ),
  "<subject>_SI", SetType('Electrode Implantation PDF report'), SetWeakAttr( 'implantation_session', 'default', 'planning', 'False' ),
  "<subject>_elec", SetType('Electrode List PDF'), SetWeakAttr( 'implantation_session', 'default' ),
  "<subject>_sag_{implantation_session}", SetType('Electrode Implantation Sagittal Image'), # In case there is more than one
  "<subject>_coro_{implantation_session}", SetType('Electrode Implantation Coronal Image'), # In case there is more than one
  "<subject>_{implantation_session}", SetType('Electrode Implantation Powerpoint report'), # In case there is more than one
  "<subject>_{implantation_session}", SetType('Electrode Implantation PDF report'), SetWeakAttr('planning', 'True' ), 
  "<subject>_SI_{implantation_session}", SetType('Electrode Implantation PDF report'), SetWeakAttr( 'implantation_session', 'default', 'planning', 'False' ),
  "<subject>_elec_{implantation_session}", SetType('Electrode List PDF'), SetWeakAttr( 'implantation_session', 'default' ),
  
  )
)

# Protocol-specific electrode models
insertFirst( '{center}', 'electrode_models', SetType('Electrode Model Directory'), SetContent(
  '{model_name}', SetType('Electrode Model')  
  )
)

