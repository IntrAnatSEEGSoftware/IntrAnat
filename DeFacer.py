#!/usr/bin/env python

import argparse, getopt, sys

from soma import aims, aimsalgo
from brainvisa import axon
#from brainvisa import anatomist
from brainvisa.data import neuroHierarchy

from brainvisa.data.readdiskitem import ReadDiskItem
from brainvisa.data.writediskitem import WriteDiskItem
from PyQt4 import QtGui, QtCore, uic, Qt
import numpy
import os, subprocess

import pdb

def main(PatientName):
    axon.initializeProcesses()
    print("MAIN")
    rdi = ReadDiskItem('Raw T1 MRI', 'aims readable volume formats', requiredAttributes={'subject':PatientName} )
    images = list( rdi._findValues( {}, None, False ) )
    
    if len(images)<1:
        print("Patient not find in the DataBase")
        return
    
    #on cherche la T1pre dans images avant de continuer a incrementer
    impre = [images[i] for i in range(len(images)) if "T1pre" in images[i].attributes()['acquisition']] 
    impost = [images[i] for i in range(len(images)) if "T1post" in images[i].attributes()['acquisition'] if not "T1postOP" in images[i].attributes()['acquisition']]
    impostOp = [images[i] for i in range(len(images)) if "T1postOp" in images[i].attributes()['acquisition']]
    
    rdi = ReadDiskItem( 'T2 MRI', 'BrainVISA volume formats', requiredAttributes={'subject':PatientName} )
    images += list( rdi._findValues( {}, None, False ) )
    rdi = ReadDiskItem( 'CT', 'BrainVISA volume formats', requiredAttributes={'subject':PatientName} )
    images += list( rdi._findValues( {}, None, False ) )
    CTim = list( rdi._findValues( {}, None, False ) )
    rdi = ReadDiskItem( 'PET', 'BrainVISA volume formats', requiredAttributes={'subject':PatientName} )
    images += list( rdi._findValues( {}, None, False ) )
    PETim = list( rdi._findValues( {}, None, False ) )
    rdi = ReadDiskItem( 'fMRI-epile', 'BrainVISA volume formats', requiredAttributes={'subject':PatientName} )
    images += list( rdi._findValues( {}, None, False ) )
    fMRIim = list( rdi._findValues( {}, None, False ) )
    rdi = ReadDiskItem( 'FLAIR', 'BrainVISA volume formats', requiredAttributes={'subject':PatientName} )
    images += list( rdi._findValues( {}, None, False ) )
    FLAIRim = list( rdi._findValues( {}, None, False ) )
    rdi = ReadDiskItem( 'FGATIR', 'BrainVISA volume formats', requiredAttributes={'subject':PatientName} )
    images += list( rdi._findValues( {}, None, False ) )
    FGATIRim= list( rdi._findValues( {}, None, False ) )

    rdi3 = ReadDiskItem( 'Head Mesh', 'Anatomist mesh formats', requiredAttributes={'subject':PatientName})
    head = list(rdi3._findValues( {}, None, False ) )
    
    voxel_size_T1 = [impre[0].attributes()['voxel_size'][0], impre[0].attributes()['voxel_size'][1], impre[0].attributes()['voxel_size'][2], 1.0]
    headVol=aims.Volume(*impre[0].attributes()['volume_dimension']+[10],dtype='S16')
    headVol.header()['voxel_size']= voxel_size_T1
    
    
    if len(head)>1:
      headPre=[head[i] for i in range(len(head)) if "T1pre" in head[i].attributes()['acquisition']]
    else:
      headPre = head
    headMesh = aims.read(headPre[0].fullName())
    vert = headMesh.vertex()
    varr = numpy.array(vert)
    norm = numpy.array(headMesh.normal())
    varr += norm * 3 # push vertices 2mm away along normal
    vert.assign(varr)
    headMesh.updateNormals()
    impreVol = aims.read(impre[0].fullName())
    aims.write(headMesh,'/tmp/headMeshInflate.gii')
    aims.SurfaceManip.rasterizeMesh(headMesh,headVol,1)
    
    #posCenterBrain = impre[0].attributes()['brainCenter']
    aims.write(headVol,'/tmp/testheadVolnoDil.nii')
    headVolnoDil = aims.read('/tmp/testheadVolnoDil.nii',20)
    dilateVol = aimsalgo.AimsMorphoDilation(headVolnoDil,12)
    aims.write(dilateVol,'/tmp/testheadVolDil.nii')
    
    #ret = subprocess.call(['cartoLinearComb','-i','/tmp/testheadVolDil.nii','-j',str(impre[0].fullName()),'-f','I1 * I2','-o','/tmp/testCarto.nii'])

    
    tofloodFill = aims.read('/tmp/testheadVolDil.nii')
    aims.floodFill(tofloodFill,[10,10,10],32767)
    #aims.floodFill( headVol, posCenterBrain, 1)
    tofloodFillArray =  numpy.asarray(tofloodFill)
    sup0=numpy.where(tofloodFillArray>0)
    pos0=numpy.where(tofloodFillArray==0)
    tofloodFillArray[sup0]=0
    tofloodFillArray[pos0]=1
    aims.write(tofloodFill,'/tmp/floodFill.nii')
    
    matimPre = numpy.asarray(impreVol)
    mask = numpy.asarray(tofloodFill)
    matimPre[numpy.where(mask==0)]=0
    
    aims.write(impreVol,'/tmp/testDeface.nii')
    
    wdiTransform = ReadDiskItem('Transform Raw T1 MRI to another image', 'Transformation matrix', exactType=True, requiredAttributes = {'subject':PatientName})
    diTransform = list(wdiTransform.findValues({}, None, False ))

    wdiTransform = ReadDiskItem('Transform PET to another image', 'Transformation matrix', exactType=True, requiredAttributes = {'subject':PatientName})
    diTransform += list(wdiTransform.findValues({}, None, False ))
    
    wdiTransform = ReadDiskItem('Transform FLAIR to another image', 'Transformation matrix', exactType=True, requiredAttributes = {'subject':PatientName})
    diTransform += list(wdiTransform.findValues({}, None, False ))
    
    wdiTransform = ReadDiskItem('Transform FGATIR to another image', 'Transformation matrix', exactType=True, requiredAttributes = {'subject':PatientName})
    diTransform += list(wdiTransform.findValues({}, None, False ))
    
    wdiTransform = ReadDiskItem('Transform fMRI-epile to another image', 'Transformation matrix', exactType=True, requiredAttributes = {'subject':PatientName})
    diTransform += list(wdiTransform.findValues({}, None, False ))
    
    wdiTransform = ReadDiskItem('Transform CT to another image', 'Transformation matrix', exactType=True, requiredAttributes = {'subject':PatientName})
    diTransform += list(wdiTransform.findValues({}, None, False ))  
    
    wdiTransform = ReadDiskItem('Transform T2 MRI to another image', 'Transformation matrix', exactType=True, requiredAttributes = {'subject':PatientName})
    diTransform += list(wdiTransform.findValues({}, None, False ))
 

    for t in diTransform:
      print(t.attributes()['modality'])
      if t.attributes()['modality'] == 't1mri' and 'postOp' in t.attributes()['acquisition']:
        trmpostop_to_pre = t
      elif t.attributes()['modality'] == 't1mri' and 'post' in t.attributes()['acquisition'] and 'postOp' not in t.attributes()['acquisition']:
          trmpost_to_pre = t
      elif t.attributes()['modality'] == 'PET':
          trmPET_to_pre = t
      elif t.attributes()['modality'] == 'ct':
          trmCT_to_pre = t
      elif t.attributes()['modality'] == 'flair':
          trmFLAIR_to_pre = t
      else:
        print("this modality will not be defaced 1: %s"%str(t.attributes()['modality']))
        print(str(t))
        pdb.set_trace()  
        
    
    
    
    wdiTransform2 = ReadDiskItem('Transformation to Scanner Based Referential', 'Transformation matrix', exactType=True, requiredAttributes = {'subject':PatientName})
    diTransform2 = list(wdiTransform2.findValues({}, None, False ))

    for t in diTransform2:
      if t.attributes()['modality'] == 't1mri' and 'postOp' in t.attributes()['acquisition']:
         trmpostop_to_SB = t
         #transfo_postop_to_SB = aims.read(trmpostop_to_SB.fullPath()).toMatrix()
      elif t.attributes()['modality'] == 't1mri' and 'pre' in t.attributes()['acquisition']:
         trmpre_to_SB = t
      elif t.attributes()['modality'] == 't1mri' and 'post' in t.attributes()['acquisition'] and 'postOp' not in t.attributes()['acquisition']:
         trmpost_to_SB = t   
      elif t.attributes()['modality'] == 'PET':
         trmPET_to_SB = t
      elif t.attributes()['modality'] == 'ct':
         trmCT_to_SB = t
      elif t.attributes()['modality'] == 'flair':
         trmFLAIR_to_SB = t
      else:
         print("this modality will not be defaced 2: %s"%str(t.attributes()['modality']))
         print(str(t))
         pdb.set_trace()
      #transfo_pre_to_SB = aims.read(trmpre_to_SB.fullPath()).toMatrix()

    import copy
    
    if 'trmpostop_to_SB' in locals():
       print("t1mri post op")
       pdb.set_trace()
       #faut foutre toutes les transfo dans le tmp sinon ca fout la merde lors d'update database
       trmpre_to_SBinvpath = trmpre_to_SB.fullPath().split('/')
       trmpre_to_SBinvpath[-1] = 'inv'+trmpre_to_SBinvpath[-1]
       trmpostop_to_pre_path = copy.deepcopy(trmpre_to_SBinvpath)
       trmpostop_to_pre_path[-1] = 'postop_to_pre.trm'
       trmpre_to_postop_path = copy.deepcopy(trmpre_to_SBinvpath)
       trmpre_to_postop_path[-1] = 'pre_to_postop.trm'
       trmpre_to_SBinvpath = os.path.join('/tmp',trmpre_to_SBinvpath[-1])
       trmpostop_to_pre_path = os.path.join('/tmp',trmpostop_to_pre_path[-1])
       trmpre_to_postop_path = os.path.join('/tmp',trmpre_to_postop_path[-1])
       
       ret = subprocess.call(['AimsInvertTransformation','-i',trmpre_to_SB.fullPath(),'-o', trmpre_to_SBinvpath])
       ret = subprocess.call(['AimsComposeTransformation', '-o',trmpostop_to_pre_path, trmpre_to_SBinvpath, trmpostop_to_pre.fullPath(), trmpostop_to_SB.fullPath()])
       ret = subprocess.call(['AimsInvertTransformation','-i',trmpostop_to_pre_path,'-o',trmpre_to_postop_path])
       
       ret = subprocess.call(['AimsResample', '-i', '/tmp/floodFill.nii', '-m', trmpre_to_postop_path, '-o', '/tmp/maskT1postop.nii', '-t', 'n', '-r',str(impostOp[0].fullPath())])
       
       aimsPostop = aims.read(impostOp[0].fullPath())
       aimsmask = aims.read('/tmp/maskT1postop.nii')
       Postopim = numpy.asarray(aimsPostop)
       Mask = numpy.asarray(aimsmask)
       
       Postopim[numpy.where(Mask==0)]=0
       
       aims.write(aimsPostop,'/tmp/T1postopDefaced.nii')
           
    if 'trmpost_to_pre' in locals():
       print("t1mri post")
       #faut foutre toutes les transfo dans le tmp sinon ca fout la merde lors d'update database
       trmpre_to_SBinvpath = trmpre_to_SB.fullPath().split('/')
       trmpre_to_SBinvpath[-1] = 'inv'+trmpre_to_SBinvpath[-1]
       trmpost_to_pre_path = copy.deepcopy(trmpre_to_SBinvpath)
       trmpost_to_pre_path[-1] = 'postop_to_pre.trm'
       trmpre_to_post_path = copy.deepcopy(trmpre_to_SBinvpath)
       trmpre_to_post_path[-1] = 'pre_to_postop.trm'
       trmpre_to_SBinvpath = os.path.join('/tmp',trmpre_to_SBinvpath[-1])
       trmpost_to_pre_path = os.path.join('/tmp',trmpost_to_pre_path[-1])
       trmpre_to_post_path = os.path.join('/tmp',trmpre_to_post_path[-1])
       
       ret = subprocess.call(['AimsInvertTransformation','-i',trmpre_to_SB.fullPath(),'-o', trmpre_to_SBinvpath])
       ret = subprocess.call(['AimsComposeTransformation', '-o',trmpost_to_pre_path, trmpre_to_SBinvpath, trmpost_to_pre.fullPath(), trmpost_to_SB.fullPath()])
       ret = subprocess.call(['AimsInvertTransformation','-i',trmpost_to_pre_path,'-o',trmpre_to_post_path])
       ret = subprocess.call(['AimsResample', '-i', '/tmp/floodFill.nii', '-m', trmpre_to_post_path, '-o', '/tmp/maskT1post.nii', '-t', 'n', '-r',str(impost[0].fullPath())])
       
       aimsPost = aims.read(impost[0].fullPath())
       aimsmask = aims.read('/tmp/maskT1post.nii')
       Postim = numpy.asarray(aimsPost)
       Mask = numpy.asarray(aimsmask)
       
       Postim[numpy.where(Mask==0)]=0
       
       aims.write(aimsPost,'/tmp/T1postDefaced.nii')
    
    if 'trmFLAIR_to_pre' in locals():
       print("FLAIR")
       trmpre_to_SBinvpath = trmpre_to_SB.fullPath().split('/')
       trmpre_to_SBinvpath[-1] = 'inv'+trmpre_to_SBinvpath[-1]
       trmFLAIR_to_pre_path = copy.deepcopy(trmpre_to_SBinvpath)
       trmFLAIR_to_pre_path[-1] = 'FLAIR_to_pre.trm'
       trmpre_to_FLAIR_path = copy.deepcopy(trmpre_to_SBinvpath)
       trmpre_to_FLAIR_path[-1] = 'pre_to_FLAIR.trm'
       trmpre_to_SBinvpath = os.path.join('/tmp',trmpre_to_SBinvpath[-1])
       trmFLAIR_to_pre_path = os.path.join('/tmp',trmFLAIR_to_pre_path[-1])
       trmpre_to_FLAIR_path = os.path.join('/tmp',trmpre_to_FLAIR_path[-1])
       
       ret = subprocess.call(['AimsInvertTransformation','-i',trmpre_to_SB.fullPath(),'-o', trmpre_to_SBinvpath])
       ret = subprocess.call(['AimsComposeTransformation', '-o',trmFLAIR_to_pre_path, trmpre_to_SBinvpath, trmFLAIR_to_pre.fullPath(), trmFLAIR_to_SB.fullPath()])
       ret = subprocess.call(['AimsInvertTransformation','-i',trmFLAIR_to_pre_path,'-o',trmpre_to_FLAIR_path])
       ret = subprocess.call(['AimsResample', '-i', '/tmp/floodFill.nii', '-m', trmpre_to_FLAIR_path, '-o', '/tmp/maskT1FLAIR.nii', '-t', 'n', '-r',str(FLAIRim[0].fullPath())])
       
       aimsFLAIR= aims.read(FLAIRim[0].fullPath())
       aimsmask = aims.read('/tmp/maskT1FLAIR.nii')
       FLAIRmat = numpy.asarray(aimsFLAIR)
       Mask = numpy.asarray(aimsmask)
       
       FLAIRmat[numpy.where(Mask==0)]=0
       
       aims.write(aimsFLAIR,'/tmp/FLAIRDefaced.nii')       
       
       
    if 'trmPet_to_pre' in locals():
       print("PET")
       trmpre_to_SBinvpath = trmpre_to_SB.fullPath().split('/')
       trmpre_to_SBinvpath[-1] = 'inv'+trmpre_to_SBinvpath[-1]
       trmPET_to_pre_path = copy.deepcopy(trmpre_to_SBinvpath)
       trmPET_to_pre_path[-1] = 'PET_to_pre.trm'
       trmpre_to_PET_path = copy.deepcopy(trmpre_to_SBinvpath)
       trmpre_to_PET_path[-1] = 'pre_to_PET.trm'
       trmpre_to_SBinvpath = os.path.join('/tmp',trmpre_to_SBinvpath[-1])
       trmPET_to_pre_path = os.path.join('/tmp',trmPET_to_pre_path[-1])
       trmpre_to_PET_path = os.path.join('/tmp',trmpre_to_PET_path[-1])
       
       ret = subprocess.call(['AimsInvertTransformation','-i',trmpre_to_SB.fullPath(),'-o', trmpre_to_SBinvpath])
       ret = subprocess.call(['AimsComposeTransformation', '-o',trmPET_to_pre_path, trmpre_to_SBinvpath, trmPET_to_pre.fullPath(), trmPET_to_SB.fullPath()])
       ret = subprocess.call(['AimsInvertTransformation','-i',trmPET_to_pre_path,'-o',trmpre_to_PET_path])
       ret = subprocess.call(['AimsResample', '-i', '/tmp/floodFill.nii', '-m', trmpre_to_PET_path, '-o', '/tmp/maskT1PET.nii', '-t', 'n', '-r',str(PETim[0].fullPath())])
       
       aimsPET= aims.read(PETim[0].fullPath())
       aimsmask = aims.read('/tmp/maskT1PET.nii')
       PETmat = numpy.asarray(aimsPET)
       Mask = numpy.asarray(aimsmask)
       
       PETmat[numpy.where(Mask==0)]=0
       
       aims.write(aimsPET,'/tmp/PETDefaced.nii')          
       
    if 'trmCT_to_pre' in locals():
       print("CT")
       trmpre_to_SBinvpath = trmpre_to_SB.fullPath().split('/')
       trmpre_to_SBinvpath[-1] = 'inv'+trmpre_to_SBinvpath[-1]
       trmCT_to_pre_path = copy.deepcopy(trmpre_to_SBinvpath)
       trmCT_to_pre_path[-1] = 'CT_to_pre.trm'
       trmpre_to_CT_path = copy.deepcopy(trmpre_to_SBinvpath)
       trmpre_to_CT_path[-1] = 'pre_to_CT.trm'
       trmpre_to_SBinvpath = os.path.join('/tmp',trmpre_to_SBinvpath[-1])
       trmCT_to_pre_path = os.path.join('/tmp',trmCT_to_pre_path[-1])
       trmpre_to_CT_path = os.path.join('/tmp',trmpre_to_CT_path[-1])
       
       ret = subprocess.call(['AimsInvertTransformation','-i',trmpre_to_SB.fullPath(),'-o', trmpre_to_SBinvpath])
       ret = subprocess.call(['AimsComposeTransformation', '-o',trmCT_to_pre_path, trmpre_to_SBinvpath, trmCT_to_pre.fullPath(), trmCT_to_SB.fullPath()])
       ret = subprocess.call(['AimsInvertTransformation','-i',trmCT_to_pre_path,'-o',trmpre_to_CT_path])
       ret = subprocess.call(['AimsResample', '-i', '/tmp/floodFill.nii', '-m', trmpre_to_CT_path, '-o', '/tmp/maskT1CT.nii', '-t', 'n', '-r',str(CTim[0].fullPath())])
       
       aimsCT = aims.read(CTim[0].fullPath())
       aimsmask = aims.read('/tmp/maskT1CT.nii')
       CTmat = numpy.asarray(aimsCT)
       Mask = numpy.asarray(aimsmask)
       
       CTmat[numpy.where(Mask==0)]=0
       
       aims.write(aimsCT,'/tmp/CTDefaced.nii')          
    
    
    pdb.set_trace()


app = QtGui.QApplication(sys.argv)  

parser = argparse.ArgumentParser(description='Deface all images found for a patient - need the head mask made by morphologist and coregistration matrix')
parser.add_argument('-i','--input_patient',help='Input Patient Name')
#parser.add_argument('-m','--mask_dir',help='Mask directory')
#parser.add_argument('-t','--mask_type',help='Mask type: could be "", SSMMI_CompositeVersor2DBSplineNCC, SSMMI_VersorOnlyNCC or MAN')
#parser.add_argument('-o','--output_dir',help='Output directory')
#parser.add_argument('-S','--suffix',help='Suffix of image filename, such as filename ended by *_Final_${suffix}_iteration_${iteration}.nii.gz, where ${iteration} is given by input flag -I or --iteration. suffix="uni_bcorr_reo" (not denoised images) or suffix="nlm_uni_bcorr_reo" (denoised images)')
#parser.add_argument('-I','--iteration',help='Reconstruction iteration, needed to load the corresponding input images being histogram equalized')

args = parser.parse_args()


if (args.input_patient!=None):
    print 'Input directory: '+str(args.input_patient)
    main(str(args.input_patient))
else:
    print("Usage: %s -i input_patient_name")%sys.argv[0]
    #print("Usage: %s -i input_directory -m mask_directory -t mask_type -o output_directory -S image_file_suffix -I iteration" % sys.argv[0])
    sys.exit(2)
