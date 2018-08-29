# -*- coding: utf-8 -*-
#  This software and supporting documentation are distributed by
#      INSERM U836 - Institut des Neurosciences de Grenoble
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

include( 'builtin' )
include( 'registration' )
include( 'anatomy' )

####################### File Formats ##########################
Format( 'Electrode Implantation format', 'f|*.elecimplant' )
Format( 'Electrode Label Format', 'f|*.eleclabel')
Format( 'Electrode sEEG Label Format', 'f|*.seegeleclabel')
Format( 'Resection json', 'f|*.resection')
Format( 'Electrode Model format', 'f|*.elecdef' )
Format( 'Subject Information format', 'f|*.subjectinfo')
Format( 'Patient Template format', 'f|*.patienttemplate')
Format( 'PTS format', 'f|*.pts' )
Format( 'EEG TRC format', 'f|*.TRC' )
Format( 'Elan EEG format', ['f|*.eeg', 'f|*.eeg.ent', 'f|*.notes.txt', 'f|*.elec.out'] )
Format ('Elan trigger format', 'f|*.pos')
Format ('JP Localizer format', 'f|*.par')
Format ('Blackrock MER format', ['f|*.ccf', 'f|*.nev', 'f|*.ns5'] )
#Format ('ImaGIN matlab format', ['f|*.mat', 'f|*.dat', 'f|*.hdr', 'f|*.txt'] )
Format ( 'Eyelink format', 'f|*.edf' )
Format ( 'PDF file', "f|*.pdf" )
Format ( 'Powerpoint file', ["f|*.ppt","f|*.pptx"] )
Format( 'Atlas metrics', 'f|*.atlasmetrics')

####################### File types ##########################
# Images
# FileType( 'T2 MRI', '3D Volume' ) -> already in builtin.py

FileType( 'CT', '3D Volume' )
#FileType( 'PET', '3D Volume' ) # Or should it be 4D ? # Already defined in /usr/local/brainvisa/toolboxes/nuclearimaging/types/OLDnuclear_imaging.py
FileType( 'fMRI-epile', '3D Volume') #'4D Volume'
FileType( 'Statistic-Data', '3D Volume')
FileType('T1 SPM resampled in MNI', '3D Volume')
FileType('FLAIR', '3D Volume')
FileType('FGATIR','3D Volume')
FileType('Resection', 'Label Volume')
FileType('ROI IntrAnat', 'ROI Graph')
FileType('Final Export Dictionaries','CSV file')
FileType('FreesurferAtlas','3D Volume')
FileType('HippoFreesurferAtlas','3D Volume')
FileType('DTIVolume','3D Volume')

FileType('leftAmygdala','Mesh')
FileType('rightAmygdala','Mesh')
FileType('leftHippo','Mesh')
FileType('rightHippo','Mesh')
FileType('leftanteroHippocampus','Mesh')
FileType('leftposteroHippocampus','Mesh')
FileType('rightanteroHippocampus','Mesh')
FileType('rightposteroHippocampus','Mesh')

FileType('leftHippocampusNII','3D Volume')
#FileType('leftposteroHippocampusNII','3D Volume')
FileType('rightHippocampusNII','3D Volume')
#FileType('rightposteroHippocampusNII','3D Volume')

FileType('Screenshot of Mars Atlas' ,'Any Type','PNG image')
#FileType('GIF of Mars Atlas','Any Type', 'GIF image')
#FileType('GIF of Electrodes','Any Type', 'GIF image')
FileType('Parcels Infos','Any Type', 'Atlas metrics')
FileType('Outliers','Any Type', 'Atlas metrics')
FileType('MP4 of Mars Atlas', 'Any Type','MP4 film')
FileType('MP4 of Electrodes','Any Type', 'MP4 film')

#Patient
FileType('SubjectInfo','Any Type', 'Subject Information format')
FileType('PatientInfoTemplate','Any Type', 'Patient Template format')

# Electrodes
FileType( 'Electrode Model', 'Any Type', 'Electrode Model format' )
FileType( 'Electrode Implantation', 'Any Type', 'Electrode Implantation format' )
FileType( 'Electrode Implantation PTS', 'Any Type', 'PTS format' )
FileType( 'Electrode Implantation TXT', 'Any Type', 'Text file' )
FileType( 'Electrode Implantation Position TXT', 'Electrode Implantation TXT' )
FileType( 'Electrode Implantation Name TXT', 'Electrode Implantation TXT' )
FileType( 'Electrodes Labels', 'Any Type', 'Electrode Label Format' )
FileType( 'Electrodes SEEG Labels', 'Any Type', 'Electrode sEEG Label Format' )

# Scanned images of the implantation (surgeon drawing)
FileType( 'Electrode Implantation Image', '2D Image' )
FileType( 'Electrode Implantation Sagittal Image', 'Electrode Implantation Image')
FileType( 'Electrode Implantation Coronal Image', 'Electrode Implantation Image')
FileType( 'Electrode Implantation report', 'Any Type', ['PDF file', 'Powerpoint file'])
FileType( 'Electrode Implantation PDF report', 'Electrode Implantation report', 'PDF file')
FileType( 'Electrode Implantation Powerpoint report', 'Electrode Implantation report', 'Powerpoint file')
FileType( 'Electrode List PDF', 'Any Type', 'PDF file')

FileType( 'Electrode Implantation Directory', 'Directory', 'Directory' )
FileType( 'Electrode Model Directory', 'Directory', 'Directory' )

# Referentials and transformations
FileType( 'Referential of Electrode', 'Referential' )
FileType( 'Electrode Implantation transformation', 'Transformation Matrix', 'Transformation Matrix' )

FileType( 'Transform Raw T1 MRI to another image', 'Transformation Matrix')
FileType( 'SPM normalization inverse deformation field', 'Any Type','BrainVISA volume formats' ) #Most probably Nifti, 5-dim
FileType( 'SPM normalization deformation field', 'Any Type','BrainVISA volume formats' ) #Most probably Nifti, 5-dimx,y,z,t,vectorDirection (three values)

FileType( 'Referential of T2 MRI', 'Referential' )
FileType( 'Transform T2 MRI to Talairach-AC/PC-Anatomist', 'Transformation Matrix')
FileType( 'Transform T2 MRI to Talairach-MNI template-SPM', 'SPM normalization matrix')
FileType( 'Transform T2 MRI to another image', 'Transformation Matrix')

FileType( 'Referential of CT', 'Referential' )
FileType( 'Transform CT to Talairach-AC/PC-Anatomist', 'Transformation Matrix')
FileType( 'Transform CT to Talairach-MNI template-SPM', 'SPM normalization matrix')
FileType( 'Transform CT to another image', 'Transformation Matrix')

FileType( 'Referential of PET', 'Referential' )
FileType( 'Transform PET to Talairach-AC/PC-Anatomist', 'Transformation Matrix')
FileType( 'Transform PET to Talairach-MNI template-SPM', 'SPM normalization matrix')
FileType( 'Transform PET to another image', 'Transformation Matrix')

FileType( 'Referential of fMRI-epile', 'Referential' )
FileType( 'Transform fMRI-epile to Talairach-AC/PC-Anatomist', 'Transformation Matrix')
FileType( 'Transform fMRI-epile to Talairach-MNI template-SPM', 'SPM normalization matrix')
FileType( 'Transform fMRI-epile to another image', 'Transformation Matrix')

FileType( 'Referential of Statistic-Data', 'Referential' )
FileType( 'Transform Statistic-Data to Talairach-AC/PC-Anatomist', 'Transformation Matrix')
FileType( 'Transform Statistic-Data to Talairach-MNI template-SPM', 'SPM normalization matrix')
FileType( 'Transform Statistic-Data to another image', 'Transformation Matrix')

FileType( 'Referential of FreesurferAtlas', 'Referential' )
FileType( 'Transform FreesurferAtlas to Talairach-AC/PC-Anatomist', 'Transformation Matrix')
FileType( 'Transform FreesurferAtlas to Talairach-MNI template-SPM', 'SPM normalization matrix')
FileType( 'Transform FreesurferAtlas to another image', 'Transformation Matrix')

FileType( 'Referential of HippoFreesurferAtlas', 'Referential' )
FileType( 'Transform HippoFreesurferAtlas to Talairach-AC/PC-Anatomist', 'Transformation Matrix')
FileType( 'Transform HippoFreesurferAtlas to Talairach-MNI template-SPM', 'SPM normalization matrix')
FileType( 'Transform HippoFreesurferAtlas to another image', 'Transformation Matrix')

FileType( 'Referential of FLAIR', 'Referential' )
FileType( 'Transform FLAIR to Talairach-AC/PC-Anatomist', 'Transformation Matrix')
FileType( 'Transform FLAIR to Talairach-MNI template-SPM', 'SPM normalization matrix')
FileType( 'Transform FLAIR to another image', 'Transformation Matrix')

FileType( 'Referential of FGATIR', 'Referential' )
FileType( 'Transform FGATIR to Talairach-AC/PC-Anatomist', 'Transformation Matrix')
FileType( 'Transform FGATIR to Talairach-MNI template-SPM', 'SPM normalization matrix')
FileType( 'Transform FGATIR to another image', 'Transformation Matrix')

FileType( 'Referential of DTIVolume', 'Referential' )
FileType( 'Transform DTIVolume to Talairach-AC/PC-Anatomist', 'Transformation Matrix')
FileType( 'Transform DTIVolume to Talairach-MNI template-SPM', 'SPM normalization matrix')
FileType( 'Transform DTIVolume to another image', 'Transformation Matrix')


FileType( 'Referential of Resection', 'Referential' )
FileType( 'Transform Resection to Talairach-AC/PC-Anatomist', 'Transformation Matrix')
FileType( 'Transform Resection to Talairach-MNI template-SPM', 'SPM normalization matrix')
FileType( 'Transform Resection to another image', 'Transformation Matrix')

FileType( 'Referential of ROI IntrAnat', 'Referential' )
FileType( 'Transform ROI IntrAnat to Talairach-AC/PC-Anatomist', 'Transformation Matrix')
FileType( 'Transform ROI IntrAnat to Talairach-MNI template-SPM', 'SPM normalization matrix')
FileType( 'Transform ROI IntrAnat to another image', 'Transformation Matrix')

# SEEG
FileType( 'SEEG Experiment', 'Directory', 'Directory' )
FileType( 'Cognitive SEEG Experiment', 'SEEG Experiment' )
FileType( 'Epilepsy Index SEEG Experiment', 'SEEG Experiment' )
FileType( 'Stimulation SEEG Experiment', 'SEEG Experiment' )

FileType( 'SEEG processing directory', 'Directory', 'Directory' )

FileType( 'SEEG Trials Directory', 'Directory', 'Directory' )

FileType( 'SEEG recording', 'Any Type', ['EEG TRC format', 'Elan EEG format'])#'ImaGIN matlab format'
FileType( 'Raw SEEG recording', 'SEEG recording', 'EEG TRC format' )

FileType( 'EpiIndex Directory', 'Directory', 'Directory' )
FileType( 'Results Directory', 'Directory', 'Directory' )
FileType( 'ImaGIN matlab files', 'Any Type', ['Matlab file', 'Text Data Table'] )
FileType( 'PDF Report', 'Any Type', 'PDF file' )

FileType( 'Elan EEG', 'SEEG recording', 'Elan EEG format' )
FileType( 'Trigger', 'Any Type', 'Elan trigger format' )

FileType( 'JP Localizer', 'Any Type', 'JP Localizer format' )
FileType( 'JP image', 'Any Type', 'JPEG image' )
FileType( 'BrainTV Film', 'Any Type', 'AVI film' )

# Micro Electrode Recordings*
FileType( 'MER Experiment', 'Directory', 'Directory' )
FileType( 'Raw MER recording', 'Any Type', 'Blackrock MER format' )

# Eye movements
FileType( 'Eye Tracking Experiment', 'Directory', 'Directory' )
FileType( 'Raw Eye Tracking', 'Any Type', 'Eyelink format' )

# Resection
FileType( 'Resection Description' ,'Any Type','Resection json')
