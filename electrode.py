#! /usr/bin/env python
# -*- coding: utf-8 -*-
# (c) Inserm U836 2013 - Manik Bhattacharjee
#
# License GNU GPL v3

import os, pickle
from soma import aims
from soma.qt_gui.qt_backend import QtGui
import pdb

# TODO : should use JSON, not pickle, but Brainvisa-4.3.0 pack does not contain json module for python

# Each electrode is a group of cylinders.
#
# This can load electrode models and create the necessary objects in Anatomist to display them
# In that use case, you should :
# el = ElectrodeModel(theCurrentAnatomist)
# el.open('fileOdModel.elecdef')
# el.setDisplayReferential(electrodeReferential)
# cylinders = el.getDisplayed()
# meshes = el.getAnatomistObjects()
class ElectrodeModel:
  
  def __init__(self, anato = None, typeColors = {}, modelPath=None, dispMode = None, dispParams=None):
    # Init of variables -> typeColors is a dict {'Plot':QColor(r,g,b), 'Other':QColor(R,G,B),...}
    self.cylinders = {}
    self.dispMode = 'real'
    self.dispParams = {}
    self.displayed = {}
    self.typeColors = typeColors
    self.a = anato
    if modelPath is not None:
      self.open(modelPath, dispMode, dispParams)
   
  # Load an Electrode Definition file
  def open(self, path=None, dispMode = None, dispParams=None,bipole=None):
    if path is None or not os.path.isfile(path):
      print 'No electrode model at '+repr(path)
      return False

    filein = open(path, 'rb')
    self.cylinders = pickle.load(filein)
    filein.close()
    
    # Reset colors
    self.currentColorHue = 0
    self.typeColors = {}
    # If shape is not defined, it's a cylinder
    for name in self.cylinders:
      if 'shape' not in self.cylinders[name]:
        self.cylinders[name]['shape']='cylinder'
      # If center of elements was not saved, compute it
      if 'center' not in self.cylinders[name]:
        pos = self.cylinders[name]['position']
        v = self.cylinders[name]['length']
        vec = self.cylinders[name]['vector']
        self.cylinders[name]['center'] = [pos[0] + v*vec[0]/2.0, pos[1]+v*vec[1]/2.0, pos[2]+v*vec[2]/2.0]
    # Create the meshes
    if dispMode:
      self.dispMode = dispMode
      if dispParams:
        self.dispParams = dispParams
    self.updateDisplay()
    return True

  # Save an Electrode Definition file
  def save(self, path = None):
    if path is None or not os.path.isfile(path):
      print "Cannot write electrode definition : invalid file path  " + repr(path)
      return False
    fileout = open(path, 'wb')
    pickle.dump(self.cylinders, fileout)
    fileout.close()
    return True

    
  # Display functions
  
  def setDisplayMode(self, mode='real', parameters={}):
	""" Mode display : 'real', 'sphere', 'off'... parameters is {'diameter: 2mm} for the sphere mode"""
	self.dispMode = mode
	self.dispParams = parameters
	
	
  def displayCylinder(self, name):
    # If it is already there, remove it
    self.undisplayCylinder(name)
    self.displayed[name] = None
    t = self.cylinders[name]['type']
    p = self.cylinders[name]['position']
    v = self.cylinders[name]['vector']
    r = self.cylinders[name]['diameter']/2.0
    newCyl = None

    if self.dispMode == 'off':
      self.displayed[name] = {'mesh':None, 'type':t}
      return
    if self.dispMode not in ['real', 'sphere', 'bipole']:
      print "Unknown dispMode %s for electrode !  Using 'real'..."%repr(self.dispMode)
      self.dispMode = 'real'
      
    if self.dispMode == 'real':
      pEnd = p[0] + v[0]*self.cylinders[name]['length'], p[1] + v[1]*self.cylinders[name]['length'], p[2] + v[2]*self.cylinders[name]['length']
      newCyl = self.a.toAObject(aims.SurfaceGenerator.cylinder(aims.Point3df(p[0], p[1], p[2]), aims.Point3df(pEnd), r, r, 24, True, True))
      self.a.releaseObject(newCyl)
    elif self.dispMode == 'sphere':
      diam = 2.0
      if self.dispParams.has_key('diameter'):
        diam = float(self.dispParams['diameter'])
      if t == 'Plot': # Ignore the other parts
        pCenter = (p[0] + v[0]*self.cylinders[name]['length']/2.0, p[1] + v[1]*self.cylinders[name]['length']/2.0, p[2] + v[2]*self.cylinders[name]['length']/2.0)
        newCyl = self.a.toAObject(aims.SurfaceGenerator.sphere(aims.Point3df(pCenter[0], pCenter[1], pCenter[2]), diam, 32))
        self.a.releaseObject(newCyl)
        #newCyl = self.a.toAObject(aims.SurfaceGenerator.cube(aims.Point3df(pCenter[0], pCenter[1], pCenter[2]), 2.0))
    elif self.dispMode == 'bipole':
        #pCenter = (p[0] + v[0]*self.cylinders[name]['length']/2.0, p[1] + v[1]*self.cylinders[name]['length']/2.0, p[2] + v[2]*self.cylinders[name]['length']/2.0)
        #to change one day
        newCyl = self.a.toAObject(aims.SurfaceGenerator.ellipse(aims.Point3df(p[0],p[1],p[2]),2.5,1.5,50))
        self.a.releaseObject(newCyl)
        #aims.SurfaceGenerator.ellipse()
        
    # Automatic color for an unknown type
    if t not in self.typeColors:
      self.currentColorHue = (self.currentColorHue + 40) % 256
      self.typeColors[t] = QtGui.QColor.fromHsv(self.currentColorHue, 245, 220, 255);
    color = self.typeColors[t] 
    if newCyl is not None:
      self.a.setMaterial(newCyl, diffuse=[color.redF(), color.greenF(), color.blueF(), color.alphaF()]) 
    self.displayed[name] = {'mesh':newCyl, 'type':t}
    #print "Adding %s mesh for %s : %s"%(t,name,repr(newCyl))

  def undisplayCylinder(self, name):
    if name in self.displayed:
      print "electrode : undisplay cylinder "+name
      if self.displayed[name]['mesh'] is not None:
        self.a.deleteObjects(self.displayed[name]['mesh'])
        print "UNDISPLAY cylinder : DELETED %s"%name
        self.displayed[name]['mesh'] = None # CURRENT

  def updateDisplay(self):
    self.clearDisplay()
    # Sort objects by name so that Elements/Plot always have the same order, hence same colors
    for name in sorted(self.cylinders, reverse=True):
      self.displayCylinder(name)

  def updateDisplayCylinder(self, name):
    self.undisplayCylinder(name)
    self.displayCylinder(name)

  def clearDisplay(self):
    # Destroy all anatomist objects and reset "displayed" list
    if len(self.displayed) != 0:
      meshes = [self.displayed[name]['mesh'] for name in self.displayed if self.displayed[name]['mesh'] is not None]
      print "electrode : Removing all meshes"
      #traceback.print_stack(limit=4)
      try:
        self.a.deleteObjects(meshes) # Does not work from locateElectrodes.py... CURRENT
      except:
        pass
    for n in self.displayed:
      self.displayed[n]['mesh']=None
    self.displayed={}
    
  def getDisplayed(self):
    return self.displayed
    
  def getCylinder(self, name):
    return self.cylinders[name]
    
  def getCylinders(self):
    return self.cylinders
    
  def getPlots(self):
    return dict([(k,p) for k,p in self.cylinders.iteritems() if p['type'] == 'Plot'])
    
  def countPlots(self):
    return len(self.getPlots())
  
  def countCylinders(self):
    return len(self.cylinders)

  def getAnatomistObjects(self):
	return [self.displayed[n]['mesh'] for n in self.displayed if self.displayed[n]['mesh'] is not None]
    
  def plotMeshes(self):
	return [self.displayed[n]['mesh'] for n in self.displayed if self.displayed[n]['type'] == 'Plot' and self.displayed[n]['mesh'] is not None]
    
  def setDisplayReferential(self, referential):
    for name in self.displayed:
      if self.displayed[name]['mesh'] is not None:
        self.a.assignReferential(referential, self.displayed[name]['mesh'])
      
  def setTypeColors(self, typeColors):
    self.typeColors = typeColors
    # Automatic color for an unknown type
    for name in self.cylinders:
      t = self.cylinders[name]['type']
      if t not in self.typeColors:
        self.currentColorHue = (self.currentColorHue + 40) % 256
        self.typeColors[t] = QtGui.QColor.fromHsv(self.currentColorHue, 245, 220, 255);
      color = self.typeColors[t]
      pdb.set_trace()
      self.a.setMaterial(newCyl, diffuse=[color.redF(), color.greenF(), color.blueF(), color.alphaF()])
    
  def displaySelect(self, name):
    g = self.a.getDefaultWindowsGroup()
    g.setSelection(self.displayed[name]['mesh'])


