#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# (c) Inserm U836 2012-2014 - Manik Bhattacharjee
#
# License GNU GPL v3
#


import sys, os, pickle
from soma.qt_gui.qt_backend import QtGui, QtCore, uic

import anatomist.direct.api as anatomist
from soma import aims

import pdb


# Small dialog to create a 3D electrode template for use in BrainVisa/Anatomist
#
# Each electrode is a group of cylinders.
# This software should provide an easy way to create new models of electrodes with realtime visualization
#
# It can also be used from another program that needs electrode models : it can load electrod models and
# create the necessary objects in Anatomist to display them along other objects.
# In that use case, you should :
# dialog = ElectrodeEditorDialog(theCurrentAnatomist)
# dialog.open()
# dialog.setDisplayReferential(electrodeReferential)
# cylinders = dialog.getDisplayed()
# meshes = dialog.getAnatomistObjects()
#
class ElectrodeEditorDialog(QtGui.QWidget):
  
  def __init__(self, anato = None, app=None):
    
    # UI init
    QtGui.QWidget.__init__(self)
    self.ui = uic.loadUi("electrodeEditor.ui", self)
    self.setWindowTitle('Electrode editor')
    
    # Init of variables
    self.app = app
    self.cylinders = {}
    self.displayed = {}
    #self.elementIndex = 1
    self.currentColorHue = 0
    self.typeColors = {}
    
    # Linking UI elements to functions
    self.connect(self.exitButton, QtCore.SIGNAL('clicked()'), self.quit)
    self.connect(self.openButton, QtCore.SIGNAL('clicked()'), self.open)
    self.connect(self.saveButton, QtCore.SIGNAL('clicked()'), self.save)
    self.connect(self.addButton, QtCore.SIGNAL('clicked()'), self.addCylinder)
    self.connect(self.updateButton, QtCore.SIGNAL('clicked()'), self.updateCylinder)
    self.connect(self.deleteButton, QtCore.SIGNAL('clicked()'), self.deleteCylinder)
    self.connect(self.axisCombo,  QtCore.SIGNAL('currentIndexChanged(QString)'), self.axisChanged)
    #self.connect(self.typeCombo, QtCore.SIGNAL('currentIndexChanged(QString)'), self.updateName)
    
    self.cylinderList.sortItems(QtCore.Qt.AscendingOrder)
    self.connect(self.cylinderList, QtCore.SIGNAL("currentItemChanged(QListWidgetItem*,QListWidgetItem*)"), self.cylinderListClick) #itemClicked
    
    # Anatomist windows
    if anato is None:
      self.a = anatomist.Anatomist('-b' )
    else:
      self.a = anato  
    layout = QtGui.QHBoxLayout( self.anatomistBox )
    self.axWindow = self.a.createWindow( 'Axial', no_decoration=True )
    self.sagWindow = self.a.createWindow( 'Sagittal', no_decoration=True )
    self.corWindow = self.a.createWindow( 'Coronal', no_decoration=True )
    self.wins = [self.axWindow, self.sagWindow, self.corWindow]
    [w.setHasCursor(0) for w in self.wins]
    self.sagWindow.setParent( self.anatomistBox )
    self.axWindow.setParent( self.anatomistBox )
    self.corWindow.setParent( self.anatomistBox )
    layout.addWidget( self.sagWindow.getInternalRep() )
    layout.addWidget( self.axWindow.getInternalRep() )
    layout.addWidget( self.corWindow.getInternalRep() )
    
    
  # Enable or disable axis values widgets if the axisCombo changes  
  def axisChanged(self, value):
    value = str(value)
    if value == "Custom...":
      self.xValue.setEnabled(True)
      self.yValue.setEnabled(True)
      self.zValue.setEnabled(True)
    else:
      self.xValue.setEnabled(False)
      self.yValue.setEnabled(False)
      self.zValue.setEnabled(False)
      
      
  # Load an Electrode Definition file
  def open(self, model=None):
    if model is not None:
      if not os.path.isfile(model):
        path = '/home/manik/prog/electrophysiology/epilepsie/'+str(model)
      else:
        path = model
    else:
      path = QtGui.QFileDialog.getOpenFileName(self, "Open Electrode Definition File", "", "Electrode Definition files (*.elecdef)")
    if not os.path.exists(path):
      print 'No electrode model at '+str(path)
      return

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
        self.cylinders[name]['center'] = [pos[0] + v*vec[0]/2.0, self.posYValue.value()+v*vec[1]/2.0, self.posZValue.value()+v*vec[2]/2.0]
    # Display the cylinders
    self.updateCylinderList()
    self.updateDisplay()

  # Save an Electrode Definition file
  def save(self):
    path = QtGui.QFileDialog.getSaveFileName(self, "Save Electrode Definition File", "", "Electrode Definition files (*.elecdef)")
    if path is None or str(path) == "":
      return
    fileout = open(path, 'wb')
    pickle.dump(self.cylinders, fileout)
    fileout.close()

  # Get the real orientation vector from the UI
  def getAxisVector(self):
    # Valeur réelle de l'axe
    ax = str(self.axisCombo.currentText())
    if  ax == "Custom...":
      vec = [self.xValue.value(),self.yValue.value(), self.zValue.value()]
    elif ax == "Axe X":
      vec = [1.0,0.0,0.0]
    elif ax == "Axe Y":
      vec = [0.0,1.0,0.0]
    elif ax == "Axe Z":
      vec = [0.0,0.0,1.0]
    else:
      QtGui.QMessageBox.warning(self, "Erreur", "L'axe ("+ax+") est inconnu !")
      return None
    return ax, vec


  # Add a new cylinder to the current electrode
  def addCylinder(self):
    # Safety check to avoid overwriting a cylinder
    name = str(self.nameEdit.text())
    if name in self.cylinders:
      QtGui.QMessageBox.warning(self, "Erreur", "Un élément avec ce nom existe déjà ! Changez le nom !")
      return
      
    # Store the cylinder
    (ax, vec) = self.getAxisVector()
    v = self.lengthValue.value()
    self.cylinders[name] = {'axis':ax,
                            'vector':vec,
                            'position':[self.posXValue.value(), self.posYValue.value(), self.posZValue.value()], 
                            'diameter': self.diameterValue.value(),
                            'length': v, 
                            'center': [self.posXValue.value()+v*vec[0]/2.0, self.posYValue.value()+v*vec[1]/2.0, self.posZValue.value()+v*vec[2]/2.0],
                            'type': str(self.typeCombo.currentText()), 
			    'shape':'cylinder'}
                            
    # add an item to the cylinderList
    item = QtGui.QListWidgetItem(name,self.cylinderList)
    
    # Display the cylinder
    self.displayCylinder(name)
    
    # Prepare for the next cylinder
    self.posXValue.setValue(self.posXValue.value() + vec[0]*self.lengthValue.value())
    self.posYValue.setValue(self.posYValue.value() + vec[1]*self.lengthValue.value())
    self.posZValue.setValue(self.posZValue.value() + vec[2]*self.lengthValue.value())
    
    ## Set an element name that is not used yet
    #while ("Element "+str(self.elementIndex)) in self.cylinders:
      #self.elementIndex += 1
    #self.nameEdit.setText("Element "+str(self.elementIndex))
    
    if self.typeCombo.currentIndex()==0:
        self.typeCombo.setCurrentIndex(1)
    elif self.typeCombo.currentIndex()==1:
        self.typeCombo.setCurrentIndex(0)

    self.updateName()
  
  # Update the currently selected cylinder
  def updateCylinder(self):
    # TODO if the name changes ?
    name = str(self.cylinderList.currentItem().text())
    (ax, vec) = self.getAxisVector()
    v = self.lengthValue.value()
    self.cylinders[name] = {'axis':ax,
                            'vector':vec,
                            'position':[self.posXValue.value(), self.posYValue.value(), self.posZValue.value()], 
                            'diameter': self.diameterValue.value(),
                            'length': v,
                            'center': [self.posXValue.value()+v*vec[0]/2.0, self.posYValue.value()+v*vec[1]/2.0, self.posZValue.value()+v*vec[2]/2.0],
                            'type': str(self.typeCombo.currentText()), 
			    'shape':'cylinder'}
    self.displayCylinder(name)
  
  # Display a cylinder selected from the list
  def selectCylinder(self, name):
    cyl = self.cylinders[name]
    self.axisCombo.setCurrentIndex(self.axisCombo.findText(cyl['axis']))
    v = cyl['vector']
    self.xValue.setValue(v[0])
    self.yValue.setValue(v[1])
    self.zValue.setValue(v[2])
    v = cyl['position']
    self.posXValue.setValue(v[0])
    self.posYValue.setValue(v[1])
    self.posZValue.setValue(v[2])
    
    self.diameterValue.setValue(cyl['diameter'])
    self.lengthValue.setValue(cyl['length'])
    self.typeCombo.setCurrentIndex(self.typeCombo.findText(cyl['type']))

    self.displaySelect(name)

  # Delete the selected 
  def deleteCylinder(self):
    item = self.cylinderList.takeItem(self.cylinderList.currentRow())
    self.undisplayCylinder(str(item.text()))
    del self.cylinders[str(item.text())]
    item = None
    
  # Refills the cylinderlist from the stored data
  def updateCylinderList(self):
    self.cylinderList.clear()
    item = None
    for name in sorted(self.cylinders):
      item = QtGui.QListWidgetItem(name,self.cylinderList)
    self.cylinderList.setCurrentItem(item)
    
    
  # Click on an item and column in the list
  def cylinderListClick(self, item,prevItem=None):
    name = str(item.text())
    self.selectCylinder(name)

  def closeEvent(self, event):
    self.quit()     
    
  def quit(self):
        reply = QtGui.QMessageBox.question(self, 'Message',
            "Are you sure to quit?", QtGui.QMessageBox.Yes | 
            QtGui.QMessageBox.No, QtGui.QMessageBox.No)

        if reply == QtGui.QMessageBox.Yes:
            self.app.quit()
        else:
            pass

  # Display functions
  def displayCylinder(self, name):
    # If it is already there, remove it
    self.undisplayCylinder(name)
    self.displayed[name] = None
    t = self.cylinders[name]['type']
    p = self.cylinders[name]['position']
    v = self.cylinders[name]['vector']
    r = self.cylinders[name]['diameter']/2.0
    pEnd = p[0] + v[0]*self.cylinders[name]['length'], p[1] + v[1]*self.cylinders[name]['length'], p[2] + v[2]*self.cylinders[name]['length']
    #print "New cylinder mesh at %.2f, %.2f, %.2f,  radius=%.2f" % (p[0], p[1], p[2], r)
    newCyl = self.a.toAObject(aims.SurfaceGenerator.cylinder(aims.Point3df(p[0], p[1], p[2]), aims.Point3df(pEnd), r, r, 24, True, True))
    
    # Couleur automatique par catégorie
    if t not in self.typeColors:
      self.currentColorHue = (self.currentColorHue + 40) % 256
      self.typeColors[t] = QtGui.QColor.fromHsv(self.currentColorHue, 245, 220, 255);
    color = self.typeColors[t] 
    self.a.setMaterial(newCyl, diffuse=[color.redF(), color.greenF(), color.blueF(), color.alphaF()])
    self.a.addObjects(newCyl, self.wins)
    self.displayed[name] = {'mesh':newCyl, 'type':t}

  def undisplayCylinder(self, name):
    if name in self.displayed:
      self.a.removeObjects(self.displayed[name]['mesh'], self.wins)
      self.a.deleteObjects(self.displayed[name]['mesh'])

  def updateDisplay(self):
    self.clearDisplay()
    # Add the new ones
    for name in self.cylinders:
      self.displayCylinder(name)

  def updateDisplayCylinder(self, name):
    # Modify
    self.a.deleteObjects([self.displayed[name]['mesh'],])
    self.displayCylinder(name)

  def clearDisplay(self):
    # Destroy all anatomist objects and reset "displayed" list
    if len(self.displayed) != 0:
      meshes = [self.displayed[name]['mesh'] for name in self.displayed]
      self.a.removeObjects(meshes, self.wins)
      self.a.deleteObjects(meshes)
    self.displayed={}
    
  def getDisplayed(self):
    return self.displayed
    
  def getCylinders(self):
    return self.cylinders
    
  def getAnatomistObjects(self):
    return [self.displayed[n]['mesh'] for n in self.displayed]
    
  def plotMeshes(self):
    return [self.displayed[n]['mesh'] for n in self.displayed if self.displayed[n]['type'] == 'Plot']
    
  def setDisplayReferential(self, referential):
    for name in self.displayed:
      self.a.assignReferential(referential, self.displayed[name]['mesh'])

  def displaySelect(self, name):
    g = self.a.getDefaultWindowsGroup()
    if name in self.displayed:
      g.setSelection(self.displayed[name]['mesh'])
    else:
      print "Cannot find %s in self.displayed -> cannot light up selection"%name
      
  def updateName(self):
      
    element_number = 1
    last_length = 0
    if self.typeCombo.currentIndex()==0:
        for i_index in range(len(self.cylinders.keys())):
            if self.cylinders[self.cylinders.keys()[i_index]]['type'] == 'Plot':
                element_number += 1
        element_name = 'Plot'+str(element_number)
        if 'Plot'+str(element_number-1) in self.cylinders.keys():
            last_length = self.cylinders['Plot'+str(element_number-1)]['length']
    elif self.typeCombo.currentIndex()==1:
        for i_index in range(len(self.cylinders.keys())):
            if self.cylinders[self.cylinders.keys()[i_index]]['type'] == 'Tube':
                element_number += 1
        element_name = 'Element '+str(element_number)
        if 'Element '+str(element_number-1) in self.cylinders.keys():
            last_length = self.cylinders['Element '+str(element_number-1)]['length']
    #while ("Element "+str(self.elementIndex)) in self.cylinders:
      #self.elementIndex += 1
    self.nameEdit.setText(element_name)
    if last_length !=0:
      self.lengthValue.setValue(last_length)
      

# Fonction principale qui lance l'interface
def main(noapp=0):
    app = None
    if noapp == 0:
      app = QtGui.QApplication(sys.argv)
    window = ElectrodeEditorDialog(app = app)
    window.show()
    if noapp == 0:
      sys.exit(app.exec_())

if __name__ == "__main__":
	main()

