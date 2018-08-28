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



t2mri_content = (
  "{acquisition}", SetDefaultAttributeValue( 'acquisition', default_acquisition ), SetNonMandatoryKeyAttribute( 'acquisition' ),
     SetContent(
       "<subject>-<acquisition>", SetType( 'T2 MRI' ),
       'registration', SetContent(
        'T2-<subject>_<acquisition>', SetType( 'Referential of T2 MRI' ),
        'T2-<subject>_<acquisition>_TO_Talairach-ACPC', SetType( 'Transform T2 MRI to Talairach-AC/PC-Anatomist' ),
        'T2-<subject>_<acquisition>_TO_Talairach-MNI', SetType( 'Transform T2 MRI to Talairach-MNI template-SPM'),
        'T2-<subject>_<acquisition>_TO_Scanner_Based', SetType( 'Transformation to Scanner Based Referential' ),
        'T2-<subject>_<acquisition>_TO_{modalityTarget}_{acquisitionTarget}', SetType( 'Transform T2 MRI to another image' ),
        'T2-<subject>_<acquisition>_Scanner_Based', SetType( 'Scanner Based Referential' ),
       ),
     )
)

flair_content = (
  "{acquisition}", SetDefaultAttributeValue( 'acquisition', default_acquisition ), SetNonMandatoryKeyAttribute( 'acquisition' ),
     SetContent(
       "<subject>-<acquisition>", SetType( 'FLAIR' ),
       'registration', SetContent(
        'FLAIR-<subject>_<acquisition>', SetType( 'Referential of FLAIR' ),
        'FLAIR-<subject>_<acquisition>_TO_Talairach-ACPC', SetType( 'Transform FLAIR to Talairach-AC/PC-Anatomist' ),
        'FLAIR-<subject>_<acquisition>_TO_Talairach-MNI', SetType( 'Transform FLAIR to Talairach-MNI template-SPM'),
        'FLAIR-<subject>_<acquisition>_TO_Scanner_Based', SetType( 'Transformation to Scanner Based Referential' ),
        'FLAIR-<subject>_<acquisition>_TO_{modalityTarget}_{acquisitionTarget}', SetType( 'Transform FLAIR to another image' ),
        'FLAIR-<subject>_<acquisition>_Scanner_Based', SetType( 'Scanner Based Referential' ),
       ),
     )
)

fgatir_content = (
  "{acquisition}", SetDefaultAttributeValue( 'acquisition', default_acquisition ), SetNonMandatoryKeyAttribute( 'acquisition' ),
     SetContent(
       "<subject>-<acquisition>", SetType( 'FGATIR' ),
       'registration', SetContent(
        'FGATIR-<subject>_<acquisition>', SetType( 'Referential of FGATIR' ),
        'FGATIR-<subject>_<acquisition>_TO_Talairach-ACPC', SetType( 'Transform FGATIR to Talairach-AC/PC-Anatomist' ),
        'FGATIR-<subject>_<acquisition>_TO_Talairach-MNI', SetType( 'Transform FGATIR to Talairach-MNI template-SPM'),
        'FGATIR-<subject>_<acquisition>_TO_Scanner_Based', SetType( 'Transformation to Scanner Based Referential' ),
        'FGATIR-<subject>_<acquisition>_TO_{modalityTarget}_{acquisitionTarget}', SetType( 'Transform FGATIR to another image' ),
        'FGATIR-<subject>_<acquisition>_Scanner_Based', SetType( 'Scanner Based Referential' ),
       ),
     )
)


ct_content = (
  "{acquisition}", SetDefaultAttributeValue( 'acquisition', default_acquisition ), SetNonMandatoryKeyAttribute( 'acquisition' ),
     SetContent(
       "<subject>-<acquisition>", SetType( 'CT' ),
       'registration', SetContent(
       'CT-<subject>_<acquisition>', SetType( 'Referential of CT' ),
       'CT-<subject>_<acquisition>_TO_Talairach-ACPC', SetType( 'Transform CT to Talairach-AC/PC-Anatomist' ),
       'CT-<subject>_<acquisition>_TO_Talairach-MNI', SetType( 'Transform CT to Talairach-MNI template-SPM'),
       'CT-<subject>_<acquisition>_TO_Scanner_Based', SetType( 'Transformation to Scanner Based Referential' ),
       'CT-<subject>_<acquisition>_TO_{modalityTarget}_{acquisitionTarget}', SetType( 'Transform CT to another image' ),
       'CT-<subject>_<acquisition>_Scanner_Based', SetType( 'Scanner Based Referential' ),
       ),
     )
)

pet_content = (
  "{acquisition}", SetDefaultAttributeValue( 'acquisition', default_acquisition ), SetNonMandatoryKeyAttribute( 'acquisition' ),
     SetContent(
       "<subject>-<acquisition>", SetType( 'PET' ),
       'registration', SetContent(
       'PET-<subject>_<acquisition>', SetType( 'Referential of PET' ),
       'PET-<subject>_<acquisition>_TO_Talairach-ACPC', SetType( 'Transform PET to Talairach-AC/PC-Anatomist' ),
       'PET-<subject>_<acquisition>_TO_Talairach-MNI', SetType( 'Transform PET to Talairach-MNI template-SPM'),
       'PET-<subject>_<acquisition>_TO_Scanner_Based', SetType( 'Transformation to Scanner Based Referential' ),
       'PET-<subject>_<acquisition>_TO_{modalityTarget}_{acquisitionTarget}', SetType( 'Transform PET to another image' ),
       'PET-<subject>_<acquisition>_Scanner_Based', SetType( 'Scanner Based Referential' ),
       ),
     )
)

fMRI_content  = (
  "{acquisition}", SetDefaultAttributeValue( 'acquisition', default_acquisition ), SetNonMandatoryKeyAttribute( 'acquisition' ),
     SetContent(
       "<subject>-<acquisition>_{subacquisition}", SetType( 'fMRI-epile' ),
       'registration', SetContent(
       'fMRI-<subject>_<acquisition>_{subacquisition}', SetType( 'Referential of fMRI-epile' ),
       'fMRI-<subject>_<acquisition>_{subacquisition}_TO_Talairach-ACPC', SetType( 'Transform fMRI-epile to Talairach-AC/PC-Anatomist' ),
       'fMRI-<subject>_<acquisition>_{subacquisition}_TO_Talairach-MNI', SetType( 'Transform fMRI-epile to Talairach-MNI template-SPM'),
       'fMRI-<subject>_<acquisition>_{subacquisition}_TO_Scanner_Based', SetType( 'Transformation to Scanner Based Referential' ),
       'fMRI-<subject>_<acquisition>_{subacquisition}_TO_{modalityTarget}_{acquisitionTarget}', SetType( 'Transform fMRI-epile to another image' ),
       'fMRI-<subject>_<acquisition>_{subacquisition}_Scanner_Based', SetType( 'Scanner Based Referential' ),
       ),
     )
)

dtivolume_content = (
  "{acquisition}", SetDefaultAttributeValue( 'acquisition', default_acquisition ), SetNonMandatoryKeyAttribute( 'acquisition' ),
     SetContent(
       "<subject>-<acquisition>", SetType( 'DTIVolume' ),
       'registration', SetContent(
        'DTIVolume-<subject>_<acquisition>', SetType( 'Referential of DTIVolume' ),
        'DTIVolume-<subject>_<acquisition>_TO_Talairach-ACPC', SetType( 'Transform DTIVolume to Talairach-AC/PC-Anatomist' ),
        'DTIVolume-<subject>_<acquisition>_TO_Talairach-MNI', SetType( 'Transform DTIVolume to Talairach-MNI template-SPM'),
        'DTIVolume-<subject>_<acquisition>_TO_Scanner_Based', SetType( 'Transformation to Scanner Based Referential' ),
        'DTIVolume-<subject>_<acquisition>_TO_{modalityTarget}_{acquisitionTarget}', SetType( 'Transform DTIVolume to another image' ),
        'DTIVolume-<subject>_<acquisition>_Scanner_Based', SetType( 'Scanner Based Referential' ),
       ),
     )
)

Statistics_content  = (
  "{acquisition}", SetDefaultAttributeValue( 'acquisition', default_acquisition ), SetNonMandatoryKeyAttribute( 'acquisition' ),
     SetContent(
       "<subject>-<acquisition>_{subacquisition}", SetType( 'Statistic-Data' ),
       'registration', SetContent(
       'Statistics-<subject>_<acquisition>_{subacquisition}', SetType( 'Referential of Statistic-Data' ),
       'Statistics-<subject>_<acquisition>_{subacquisition}_TO_Talairach-ACPC', SetType( 'Transform Statistic-Data to Talairach-AC/PC-Anatomist' ),
       'Statistics-<subject>_<acquisition>_{subacquisition}_TO_Talairach-MNI', SetType( 'Transform Statistic-Data to Talairach-MNI template-SPM'),
       'Statistics-<subject>_<acquisition>_{subacquisition}_TO_Scanner_Based', SetType( 'Transformation to Scanner Based Referential' ),
       'Statistics-<subject>_<acquisition>_{subacquisition}_TO_{modalityTarget}_{acquisitionTarget}', SetType( 'Transform Statistic-Data to another image' ),
       'Statistics-<subject>_<acquisition>_{subacquisition}_Scanner_Based', SetType( 'Scanner Based Referential' ),
       ),
     )
)

FreeSurferAtlas_content  = (
  "{acquisition}", SetDefaultAttributeValue( 'acquisition', default_acquisition ), SetNonMandatoryKeyAttribute( 'acquisition' ),
     SetContent(
       "<subject>-<acquisition>", SetType( 'FreesurferAtlas' ),
       'registration', SetContent(
       'FreeSurferAtlas-<subject>_<acquisition>', SetType( 'Referential of FreesurferAtlas' ),
       'FreeSurferAtlas-<subject>_<acquisition>_TO_Talairach-ACPC', SetType( 'Transform FreesurferAtlas to Talairach-AC/PC-Anatomist' ),
       'FreeSurferAtlas-<subject>_<acquisition>_TO_Talairach-MNI', SetType( 'Transform FreesurferAtlas to Talairach-MNI template-SPM'),
       'FreeSurferAtlas-<subject>_<acquisition>_TO_Scanner_Based', SetType( 'Transformation to Scanner Based Referential' ),
       'FreeSurferAtlas-<subject>_<acquisition>_TO_{modalityTarget}_{acquisitionTarget}', SetType( 'Transform FreesurferAtlas to another image' ),
       'FreeSurferAtlas-<subject>_<acquisition>_Scanner_Based', SetType( 'Scanner Based Referential' ),
       ),
     )
)

HippoFreesurferAtlas_content  = (
  "{acquisition}", SetDefaultAttributeValue( 'acquisition', default_acquisition ), SetNonMandatoryKeyAttribute( 'acquisition' ),
     SetContent(
       "<subject>-<acquisition>", SetType( 'HippoFreesurferAtlas' ),
       'registration', SetContent(
       'HippoFreesurferAtlas-<subject>_<acquisition>', SetType( 'Referential of HippoFreesurferAtlas' ),
       'HippoFreesurferAtlas-<subject>_<acquisition>_TO_Talairach-ACPC', SetType( 'Transform HippoFreesurferAtlas to Talairach-AC/PC-Anatomist' ),
       'HippoFreesurferAtlas-<subject>_<acquisition>_TO_Talairach-MNI', SetType( 'Transform HippoFreesurferAtlas to Talairach-MNI template-SPM'),
       'HippoFreesurferAtlas-<subject>_<acquisition>_TO_Scanner_Based', SetType( 'Transformation to Scanner Based Referential' ),
       'HippoFreesurferAtlas-<subject>_<acquisition>_TO_{modalityTarget}_{acquisitionTarget}', SetType( 'Transform HippoFreesurferAtlas to another image' ),
       'HippoFreesurferAtlas-<subject>_<acquisition>_Scanner_Based', SetType( 'Scanner Based Referential' ),
       ),
     )
)

resection_content = (
  "{acquisition}", SetDefaultAttributeValue( 'acquisition', default_acquisition ), SetNonMandatoryKeyAttribute( 'acquisition' ),
     SetContent(
       "<subject>-<acquisition>", SetType( 'Resection' ),
       "ROI-<subject>-<acquisition>", SetType( 'ROI IntrAnat' ),
       "Info-<subject>-<acquisition>", SetType('Resection Description'),
       ),
)

#roi_intranat_content = (
#  "{acquisition}", SetDefaultAttributeValue( 'acquisition', default_acquisition ), SetNonMandatoryKeyAttribute( 'acquisition' ),
#     SetContent(
#
#       'registration', SetContent(
#       'ROI IntrAnat-<subject>_<acquisition>', SetType( 'Referential of ROI IntrAnat' ),
#       'ROI IntrAnat-<subject>_<acquisition>_TO_Talairach-ACPC', SetType( 'Transform ROI IntrAnat to Talairach-AC/PC-Anatomist' ),
#       'ROI IntrAnat-<subject>_<acquisition>_TO_Talairach-MNI', SetType( 'Transform ROI IntrAnat to Talairach-MNI template-SPM'),
#       'ROI IntrAnat-<subject>_<acquisition>_TO_Scanner_Based', SetType( 'Transformation to Scanner Based Referential' ),
#       'ROI IntrAnat-<subject>_<acquisition>_TO_{modalityTarget}_{acquisitionTarget}', SetType( 'Transform Resection to another image' ),
#       'ROI IntrAnat-<subject>_<acquisition>_Scanner_Based', SetType( 'Scanner Based Referential' ),
#       ),
#     )
#)

# SEEG DATA
#seeg_content = (
#
#  "<subject>_{experiment}_{expNumber} ", SetType( 'SEEG Experiment' ),
#     SetContent(
#       "<subject>_<experiment>_{subexperiment}__{expId}", SetType( 'Raw SEEG recording' ),
#       "<subject>_<experiment>__{expId}", SetType( 'Raw SEEG recording' ),
#       "<subject>_<experiment>_{subexperiment}", SetType( 'Raw SEEG recording' ),
#       "<subject>_<experiment>", SetType( 'Raw SEEG recording' ),
#       "<seegAcq>", SetType( 'Raw SEEG recording' ),
#       "{seegProcessing}", SetType( 'SEEG processing directory' ),
#     ),
#  "<subject>_{experiment}", SetType( 'SEEG Experiment' ),
#     SetContent(
#       "<subject>_<experiment>_{subexperiment}__{expId}", SetType( 'Raw SEEG recording' ),
#       "<subject>_<experiment>__{expId}", SetType( 'Raw SEEG recording' ),
#       "<subject>_<experiment>_{subexperiment}", SetType( 'Raw SEEG recording' ),
#       "<subject>_<experiment>", SetType( 'Raw SEEG recording' ),
#       "<seegAcq>", SetType( 'Raw SEEG recording' ),
#       "{seegProcessing}", SetType( 'SEEG processing directory' ),
#     ),
#)



insert( '{center}/{subject}',
  't2mri', SetWeakAttr( 'modality', 't2mri' ),
    apply( SetContent, t2mri_content)
)

insert( '{center}/{subject}',
  'flair', SetWeakAttr( 'modality', 'flair' ),
    apply( SetContent, flair_content)
)

insert( '{center}/{subject}',
  'fgatir', SetWeakAttr( 'modality', 'fgatir' ),
    apply( SetContent, fgatir_content)
)

insert( '{center}/{subject}',
  'ct', SetWeakAttr( 'modality', 'ct' ),
    apply( SetContent, ct_content)
)

insert( '{center}/{subject}',
  'pet', SetWeakAttr( 'modality', 'pet' ),
    apply( SetContent, pet_content)
)

insert( '{center}/{subject}',
  'fmri_epile', SetWeakAttr( 'modality', 'fmri_epile' ),
    apply( SetContent, fMRI_content)
)

insert( '{center}/{subject}',
  'DTI', SetWeakAttr( 'modality', 'dtivolume' ),
    apply( SetContent, dtivolume_content)
)

insert( '{center}/{subject}',
  'Statistic-Data', SetWeakAttr( 'modality', 'statistic_data' ),
    apply( SetContent, Statistics_content)
)

insert( '{center}/{subject}',
  'FreesurferAtlas', SetWeakAttr( 'modality', 'freesurfer_atlas' ),
    apply( SetContent, FreeSurferAtlas_content)
)

insert( '{center}/{subject}',
  'HippoFreesurferAtlas', SetWeakAttr( 'modality', 'hippofreesurfer_atlas' ),
    apply( SetContent, HippoFreesurferAtlas_content)
)

insert( '{center}/{subject}',
  'Resection', SetWeakAttr( 'modality', 'resection' ),
    apply( SetContent,resection_content)
)

#Antonella
insert( '{center}/{subject}/QualityControl', '<subject>_ScreenshotMarsAtlas',SetType('Screenshot of Mars Atlas'))
insert( '{center}/{subject}/QualityControl', '<subject>_MP4MarsAtlas',SetType('MP4 of Mars Atlas'))
insert( '{center}/{subject}/QualityControl', '<subject>_MP4Electrodes',SetType('MP4 of Electrodes'))
insert( '{center}/{subject}/QualityControl', '<subject>_InfosParcels',SetType('Parcels Infos'))
#insert( '{center}/{subject}/QualityControl', '<subject>_MP4MarsAtlas',SetType('MP4 of Mars Atlas'))

#####

#insert( '{center}/{subject}',
#  'ROI IntrAnat', SetWeakAttr( 'modality', 'roi_intranat' ),
#    apply( SetContent,roi_intranat_content)
#)

insertFirst( '{center}/{subject}/t1mri/{acquisition}/registration', 'T1-<subject>_<acquisition>_TO_{modalityTarget}_{acquisitionTarget}', SetType( 'Transform Raw T1 MRI to another image' )
)

insert( '{center}/{subject}/t1mri/{acquisition}', 'y_<subject>_inverse', SetType( 'SPM normalization inverse deformation field' )
)
insert( '{center}/{subject}/t1mri/{acquisition}', 'y_<subject>', SetType( 'SPM normalization deformation field' )
)

insert( '{center}/{subject}/t1mri/{acquisition}','w<subject>',SetType('T1 SPM resampled in MNI'))

insert( '{center}/{subject}/t1mri/{acquisition}/default_analysis/segmentation/mesh','leftAmygdala<subject>',SetType('leftAmygdala'))
insert( '{center}/{subject}/t1mri/{acquisition}/default_analysis/segmentation/mesh','rightAmygdala<subject>',SetType('rightAmygdala'))
insert( '{center}/{subject}/t1mri/{acquisition}/default_analysis/segmentation/mesh','leftanteroHippocampus<subject>',SetType('leftanteroHippocampus'))
insert( '{center}/{subject}/t1mri/{acquisition}/default_analysis/segmentation/mesh','leftposteroHippocampus<subject>',SetType('leftposteroHippocampus'))
insert( '{center}/{subject}/t1mri/{acquisition}/default_analysis/segmentation/mesh','rightanteroHippocampus<subject>',SetType('rightanteroHippocampus'))
insert( '{center}/{subject}/t1mri/{acquisition}/default_analysis/segmentation/mesh','rightposteroHippocampus<subject>',SetType('rightposteroHippocampus'))
insert( '{center}/{subject}/t1mri/{acquisition}/default_analysis/segmentation/mesh','rightHippocampus<subject>',SetType('rightHippo'))
insert( '{center}/{subject}/t1mri/{acquisition}/default_analysis/segmentation/mesh','leftHippocampus<subject>',SetType('leftHippo'))

insert( '{center}/{subject}/t1mri/{acquisition}/default_analysis/segmentation','leftHippocampus<subject>',SetType('leftHippocampusNII'))
#insert( '{center}/{subject}/t1mri/{acquisition}/default_analysis/segmentation/mesh','leftposteroHippocampus<subject>',SetType('leftposteroHippocampusNII'))
insert( '{center}/{subject}/t1mri/{acquisition}/default_analysis/segmentation','rightHippocampus<subject>',SetType('rightHippocampusNII'))
#insert( '{center}/{subject}/t1mri/{acquisition}/default_analysis/segmentation/mesh','rightposteroHippocampus<subject>',SetType('rightposteroHippocampusNII'))

