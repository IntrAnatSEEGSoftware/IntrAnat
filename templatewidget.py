#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Widget to choose or create a common template for group studies
# Examples include MNI referential with Colin27 of MNI152 templates, custom SPM/Dartel templates...
#
# (c) Inserm U836 2012-2014 - Manik Bhattacharjee
#
# License GNU GPL v3
#
#
from soma.qt_gui.qt_backend import QtGui, QtCore, uic

import sys, os


from brainvisa import axon
from brainvisa.data.readdiskitem import ReadDiskItem
from brainvisa.data.writediskitem import WriteDiskItem
import brainvisa.registration as registration
import numpy

from externalprocesses import *

########## SPM calls
# Convert SPM normalization _sn.mat to vector field
spm_SnToField8 = """try, spm('defaults', 'FMRI');spm_jobman('initcfg');
clear matlabbatch;
FileNameSN = '%s';
matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.sn2def.matname{1}=FileNameSN;
matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.sn2def.vox=[NaN NaN NaN];
matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.sn2def.bb=NaN*ones(2,3);
matlabbatch{1}.spm.util.defs.comp{1}.inv.space{1}=['%s' ',1'];
matlabbatch{1}.spm.util.defs.ofname='%s';
matlabbatch{1}.spm.util.defs.fnames='';
matlabbatch{1}.spm.util.defs.savedir.saveusr{1}=spm_str_manip(FileNameSN,'h');
matlabbatch{1}.spm.util.defs.interp=1;
spm_jobman('run',matlabbatch);catch, disp 'AN ERROR OCCURED'; end;quit;""" # %(FileNameSN, FileSource, ofname) -> '_sn.mat' file and source image FileSource (normalized with the _sn). For the Database, we want y_<subject>_inverse.nii, so we need ofname = '<subject>_inverse' --> Maybe should provide also the output dir ? Right now, same as _sn.mat

# API changed in SPM12...
spm_SnToField12 = """try, spm('defaults', 'FMRI');spm_jobman('initcfg');
clear matlabbatch;
FileNameSN = '%s';
matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.sn2def.matname{1}=FileNameSN;
matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.sn2def.vox=[NaN NaN NaN];
matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.sn2def.bb=NaN*ones(2,3);
matlabbatch{1}.spm.util.defs.comp{1}.inv.space = {'%s'};
matlabbatch{1}.spm.util.defs.out{1}.savedef.ofname = '%s';
matlabbatch{1}.spm.util.defs.out{1}.savedef.savedir.saveusr{1}=spm_str_manip(FileNameSN,'h');
spm_jobman('run',matlabbatch);catch, disp 'AN ERROR OCCURED'; end;quit;"""

spm_SnToField = spm_SnToField12

# Using a vector field, transform points MRI Scanner-based coordinates to MNI coordinates
spm_ConvertPointsToMNI = """try,
P='%s';
P1=spm_vol([P ',1,1']);
P2=spm_vol([P ',1,2']);
P3=spm_vol([P ',1,3']);
[V1,XYZ]=spm_read_vols(P1);
V2=spm_read_vols(P2);
V3=spm_read_vols(P3);

%% Apply tranformation to electrodes
wPosElectrode=%s;
for i1=1:size(PosElectrode,1)
    D=(XYZ(1,:)-PosElectrode(i1,1)).^2+(XYZ(2,:)-PosElectrode(i1,2)).^2+(XYZ(3,:)-PosElectrode(i1,3)).^2;
    [tmp,order]=sort(D);
    tmp=tmp(1:18);      %% cubic neighborhood
    order=order(1:18);
    W=1./tmp;           %% weight inverse to distance
    if sum(isinf(W))>0
        W=[1 zeros(1,length(W)-1)];
    end
    wPosElectrode(i1,:)=[sum(V1(order).*W)./sum(W) sum(V2(order).*W)./sum(W) sum(V3(order).*W)./sum(W)];
end
dlmwrite(%s,wPosElectrode, 'delimiter',' ','precision',16);
;catch, disp 'AN ERROR OCCURED'; end;quit;
""" # %(field_inverse.nii,matrixOfCoords(plots, [1,2,3]), fileOut)

# Read deformation field y_<subject>_inverse.nii, apply the vector field to scanner-based coordinates of electrodes
spm_normalizePoints = """
try, P='%s';
P1=spm_vol([P ',1,1']);
P2=spm_vol([P ',1,2']);
P3=spm_vol([P ',1,3']);
[V1,XYZ]=spm_read_vols(P1);
V2=spm_read_vols(P2);
V3=spm_read_vols(P3);

%% Apply tranformation to electrodes
PosElectrode = dlmread('%s');
wPosElectrode=PosElectrode;
for i1=1:size(PosElectrode,1)
D=(XYZ(1,:)-PosElectrode(i1,1)).^2+(XYZ(2,:)-PosElectrode(i1,2)).^2+(XYZ(3,:)-PosElectrode(i1,3)).^2;
[tmp,order]=sort(D);
tmp=tmp(1:18);      %%  cubic neighborhood
order=order(1:18);
W=1./tmp;           %%  weight inverse to distance
if sum(isinf(W))>0
W=[1 zeros(1,length(W)-1)];
end
wPosElectrode(i1,:)=[sum(V1(order).*W)./sum(W) sum(V2(order).*W)./sum(W) sum(V3(order).*W)./sum(W)];
end
dlmwrite('%s',wPosElectrode,'precision',18);
catch, disp 'AN ERROR OCCURED'; end;quit;
"""


class TemplateMRI:
  def __init__(self, name="", refId=-1, volRefId=-1, anatomist=None):
    self.name = name
    self.refId = refId
    if refId != -1 and volRefId == -1:
      self.volRefId = refId
    else:
      self.volRefId = volRefId
    self.referentialDiskItem = None
    self.referentialAnatomist = None
    if self.refId != -1:
      try:
        self.referentialDiskItem = registration.getTransformationManager().referential(self.refId)
        if anatomist is not None:
          self.referentialAnatomist = anatomist.createReferential(self.referentialDiskItem)
      except:
        pass
    self.volumes = None
    self.findVolumesInDB(volRefId)

  def findVolumesInDB(self, volRefId=None):
    if self.volumes:
      return self.volumes
    if volRefId is None:
      volRefId = self.volRefId
    rdi = ReadDiskItem( 'anatomical Template', 'aims readable volume formats', requiredAttributes={'referential':str(volRefId)} )
    self.volumes = list( rdi._findValues( {}, None, False ) )
    return self.volumes

  def normalizeCoordinates(self, coords, referential):
    """ Fonction that normalizes a list of (x,y,z) to the template referential if sufficient data is available for the provided referential (MUST BE REIMPLEMENTED BY EACH KIND OF TEMPLATE) """
    print "TemplateMRI : normalization not implemented for generic class !"
    return []
    pass
  def denormalizeCoordinates(self, coords, referential):
    """ Fonction that denormalizes a list of (x,y,z) from the template referential if sufficient data is available to the provided referential (MUST BE REIMPLEMENTED BY EACH KIND OF TEMPLATE) """
    print "TemplateMRI : denormalization not implemented for generic class !"
    return []
    pass

  def resampleVolume(referential, readPath, writePath):
    """ Resamples the provided volume (Nifti...) in the template referential """
    print "TemplateMRI : resampling not implemented for generic class !"
    return []
    pass


class TemplateMNI(TemplateMRI):

  def __init__(self, anatomist = None):
    TemplateMRI.__init__(self, name="MNI", refId=registration.talairachMNIReferentialId, volRefId='19bfee8e-51b1-4d9e-8721-990b9f88b12f', anatomist=anatomist) #registration.talairachMNIReferentialId

  def normalizeCoordinates(self, coords, refId):
    """Normalize the coordinates in MNI referential. refId is the uuid of a referential of the T1pre used to store the coords"""
    return self.convertT1ScannerBasedToMni(coords, refId)

  def getT1preMniTransform(self, refId):
    """Returns the path of the transformation to MNI (vector field) and compute it if necessary (from _sn.mat)"""
    # Get referential file
    transfoManager = registration.getTransformationManager()
    refDiskitem = transfoManager.referential( refId )
    # Find _sn.mat
    rdi = ReadDiskItem( 'SPM2 normalization matrix', 'Matlab file' )
    rdiT1 = ReadDiskItem( 'Raw T1 MRI', 'aims readable volume formats')
    diT1 = rdiT1.findValue(refDiskitem)
    di = rdi.findValue(diT1)
    if di is None:
      print "SPM deformation _sn.mat not found in database"
      return None
    # Convert to field
    wdi = WriteDiskItem( 'SPM normalization inverse deformation field', 'NIFTI-1 image' )
    diField = wdi.findValue(di)
    if diField is None:
      print "Cannot find path to save MNI vector field in the DB"
      return None
    #For a file /database/y_SubjectName_inverse.nii, get SubjectName_inverse
    ofname = os.path.basename(diField.fullPath()).lstrip('y_').rsplit('.',1)[0]
    if not os.path.exists(diField.fullPath()):
      matlabRun(spm_SnToField%(str(di.fullPath()), str(diT1.fullPath()),  ofname) )
    else:
      print "Deformation field already present : not recomputed in %s"%diField.fullPath()
    if os.path.exists(diField.fullPath()):
      return diField.fullPath()
    else:
      print "Matlab did not convert the MNI transform to vector field !"
      return None

  def convertT1ScannerBasedToMni(self, points, refId):
    """Converts an array of points [x,y,z] in scanner-based coords to MNI coords if deformation field is available"""
    field = self.getT1preMniTransform(refId)
    if field is None:
      print "MNI deformation field not found"
      return None
    tmpOutput = getTmpFilePath('csv')
    arr = numpy.asarray(points)#([ [1,2,3], [4,5,6], [7,8,9] ])
    numpy.savetxt(tmpOutput, arr, delimiter=",")
    print "Launching SPM NORMALIZE POINTS with %s"%tmpOutput
    matlabRun(spm_normalizePoints % (field, tmpOutput, tmpOutput))
    out = numpy.loadtxt(tmpOutput, delimiter=",")
    os.remove(tmpOutput)
    if numpy.array_equal(out, arr):
      print "Points to MNI : Error, result read is identical to input"
      return None
    if out.shape != arr.shape:
      print "Points to MNI : Error, result (%s) has not the same number of elements as input (%s)"%(repr(out),repr(arr))
      return None
    return out.tolist()



class TemplateWidget(QtGui.QWidget):
  def __init__(self, app=None):
    QtGui.QWidget.__init__(self)




if __name__ == "__main__":
  app = QtGui.QApplication(sys.argv)
  axon.initializeProcesses()
  from brainvisa.data.readdiskitem import ReadDiskItem
  from brainvisa.data.writediskitem import WriteDiskItem
  window = TemplateWidget()
  window.show()
  sys.exit(app.exec_())
