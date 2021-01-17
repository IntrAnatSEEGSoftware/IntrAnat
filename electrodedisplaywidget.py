#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Widget to select and display electrode plots on a common referential/template
#
# (c) Inserm U836 2012-2014 - Manik Bhattacharjee
#
# License GNU GPL v3
#
#
from soma.qt_gui.qt_backend import QtGui, QtCore, uic

import sys, pickle, shutil, traceback, os, json, re, numpy, csv


from brainvisa import axon
from brainvisa.configuration import neuroConfig
neuroConfig.gui = True
from brainvisa import anatomist
from soma import aims
import brainvisa.registration as registration

from locateElectrodes import createElectrode, getPlotsCenters, getPlots, getPlotsNames, createBipole
from referentialconverter import ReferentialConverter
from templatewidget import TemplateMRI, TemplateMNI
from brainvisa.data.readdiskitem import ReadDiskItem
from brainvisa.data.writediskitem import WriteDiskItem
from brainvisa.data import neuroHierarchy
from readSulcusLabelTranslationFile import *
from readFreesurferLabelFile import *
# from externalprocesses import PythonExecutor
from readFunctionalTractography import *
from scipy import spatial as sc_sp
from collections import OrderedDict
from locateElectrodes import natural_keys

from bipoleSEEGColors import bipoleSEEGColors
from control_ftract2 import *

import pdb

def loadElectrodeModels():
    """Load electrode models from the database"""
    models = {}
    rdiEM = ReadDiskItem('Electrode Model', 'Electrode Model format')
    result = list (rdiEM._findValues( {}, None, False ) )
    for e in result:
      #WARNING if a model is available from multiple protocols, it will use only one
      models[str(e.attributes()['model_name'])] = e
    return models


class ElectrodeDisplayWidget(QtGui.QWidget):
  def __init__(self, app=None, ana = None,dataSubjects = None):
    QtGui.QWidget.__init__(self)
    uic.loadUi("groupPlots.ui", self)
    self.subjects = []
    self.subjItems = []
    self.implantations = {}
    self.plotsData = {}
    self.testDataSubjects = dataSubjects
    self.taskCounter = 0
    self.tasks = []
    self.meshes = {}
    self.bipolesmeshes = {}
    self.dispMode = 'off'#'sphere'
    self.dispParams = {'diameter':2.0}
    self.transfoManager = registration.getTransformationManager()
    # Get ReferentialConverter (for Talairach, AC-PC...)
    self.refConv = ReferentialConverter()

    self.electrodeModels = loadElectrodeModels()
    self.addSelectionButton.clicked.connect(self.addSelection)
    self.removeSelectionButton.clicked.connect(self.removeSelection)
    self.removePlotsNotTRC.clicked.connect(self.removeNotTRC)
    self.removePlotsLeftSide.clicked.connect(lambda :self.removePlotsLeftRight('Left'))
    self.removePlotsRightSide.clicked.connect(lambda :self.removePlotsLeftRight('Right'))
    self.addAroundButton.clicked.connect(self.selectAround)
    #self.normalizeButton.clicked.connect(self.normalizeCoords)
    #self.saveNormalizedButton.clicked.connect(self.saveNormalizedCoords)
    self.selectionList.itemDoubleClicked.connect(self.updatePlotSelected)
    self.generateStatsButton.clicked.connect(self.generateStatisticsContacts)
    
    #fill the combo possibility.
    loca = ['*']
    parcels_namesMA= readSulcusLabelTranslationFile('parcels_label_name.txt')
    loca.extend(parcels_namesMA.values())
    loca.sort()
    self.AddMAparcels2SelectioncomboBox.clear()
    self.AddMAparcels2SelectioncomboBox.addItems(loca)
    self.AddMAparcels2SelectioncomboBox.currentIndexChanged.connect(self.AddMAparcels2selection)
    
    self.radioButtonbothHemi.toggled.connect(self.changeBothRightDisplay)
    #self.radioButtonAllRight.toggled.connect(self.changeBothRightDisplay)
    
    self.radioButtonContactDisplay.toggled.connect(self.contactSEEGDisplay)
    #self.radioButtonsEEGResults.toggled.connect(self.contactSEEGDisplay)
    
    pix = QtGui.QPixmap('/home/b67-belledone/Desktop/epilepsie-manik/Logo-F-TRACT.xpm' )
    anatomist.anatomist.cpp.IconDictionary.instance().addIcon('ftract_control', pix)
    ad = anatomist.anatomist.cpp.ActionDictionary.instance()
    #control = ensemble d'action
    ad.addAction( 'fTract_Action', StimulateResults )
    cd = anatomist.anatomist.cpp.ControlDictionary.instance()
    cd.addControl( 'ftract_control', ControlFtract, 25 )
    cm = anatomist.anatomist.cpp.ControlManager.instance()   
    cm.addControl('QAGLWidget3D','','ftract_control')


    
    # Anatomist windows and objects
    if ana == None:
       self.a = anatomist.Anatomist('-b' )
    else:
      self.a = ana

    layout = QtGui.QHBoxLayout( self.viewWidget )
    self.axWindow = self.a.createWindow( 'Axial' )#, no_decoration=True )
    self.axWindow.setParent(self.viewWidget)
    layout.addWidget( self.axWindow.getInternalRep() )

    self.sagWindow = self.a.createWindow( 'Sagittal' )#, no_decoration=True )
    self.sagWindow.setParent(self.viewWidget)
    layout.addWidget( self.sagWindow.getInternalRep() )


    self.axWindow.internalRep.otherwindow = self.sagWindow
    self.windows = [self.axWindow, self.sagWindow]
    #pdb.set_trace()
    #self.axWindow.connect()



    self.templates = {'MNI':TemplateMNI(self.a)}
    #self.templateCombo.clear()
    #self.templateCombo.addItems(sorted(self.templates.keys()))
    self.setTemplate(self.templates['MNI'])
    self.templReferential = None
    #self.templateCombo.currentIndexChanged.connect(self.templateChanged)

    self.subjectList.itemSelectionChanged.connect(self.subjectSelectionChanged)
    self.electrodeList.itemSelectionChanged.connect(self.electrodeSelectionChanged)
    self.selectionList.itemSelectionChanged.connect(self.selectedSelectionChanged)
    self.addDisplayButton.clicked.connect(self.displayImage)
    self.addMNIImageToDisplayList.clicked.connect(self.addMNIImagetoList)
    self.addMNIMeshTextToDisplayList.clicked.connect(self.addMNIMeshTexttoList)


  def setStatus(self, text):
    self.statusLabel.setText(str(text))

  def incTaskCounter(self):
    self.taskCounter = self.taskCounter + 1
    self.setStatus(u"Tasks in progress : "+str(self.taskCounter))
    return self.taskCounter

  def decTaskCounter(self):
    self.taskCounter = self.taskCounter - 1
    self.setStatus(u"Tasks in progress : "+str(self.taskCounter))
    return self.taskCounter

#   def startTask(self, taskFunction):
#     pe = PythonExecutor(taskFunction)
#     self.tasks.append(pe)
#     # Remove it from the list of threads when finished
#     pe.finished.connect(lambda th=pe:self.taskFinished(th))
#     self.incTaskCounter()
#     pe.start()

#   def taskFinished(self, thread):
#     self.tasks.remove(thread)
#     self.decTaskCounter()

  def setTemplate(self, templ):
    """Set the template used as a common referential"""
    # Un nom, des données (IRM ?) un identifiant de référentiel pour le refconv ?
    self.template = templ
    if self.template.referentialAnatomist:
      self.a.assignReferential(self.template.referentialAnatomist, self.windows)
    if self.template.volumes:
      self.displayCombo.clear()
      #self.displayCombo.addItems(["Image "+str(i) for i in range(len(self.template.volumes))])
      self.displayCombo.addItems([os.path.split(im.fullPath())[1] for im in self.template.volumes])

  #def templateChanged(self, tpl):
    #"""The combo box to select the template was changed"""
    ## Read combo, setTemplate, remove selection or reset display (no spheres, no images)...
    #pass

  def displayImage(self):
    try:
      #try:
        #self.a.removeObjects([self.currentImage,],self.windows) #self.axWindow.removeObjects(self.currentImage)
        ##self.sagWindow.removeObjects(self.currentImage)
      #except:
          #pass
      self.currentImage = self.a.loadObject(self.template.volumes[self.displayCombo.currentIndex()])
      self.a.addObjects([self.currentImage], self.windows)
    except:
      print "Could not add selected image"
      pdb.set_trace()

  def subjectSelectionChanged(self):
    """Subject selection changed, update electrode list selection"""
    selected = [str(item.text()) for item in self.subjectList.selectedItems()]
    for i in xrange(self.electrodeList.count()):
      selec = False
      for s in selected:
        if str(self.electrodeList.item(i).text()).startswith(s):
          selec = True
          break
      self.electrodeList.item(i).setSelected(selec)

  def electrodeSelectionChanged(self):
    """Electrode selection changed, update plot list selection"""
    selected = [str(item.text()) for item in self.electrodeList.selectedItems()]
    for i in xrange(self.plotList.count()):
      selec = False
      for s in selected:
        if str(self.plotList.item(i).text()).startswith(s):
          selec = True
          break
      self.plotList.item(i).setSelected(selec)

  #def plotSelectionChanged(self):
  #  pass

  def selectedSelectionChanged(self):
    """In the selected plots list, the selected items changed -> update the view"""
    # Update Anatomist selection -> select all meshes for the selected plots
    g = self.a.getDefaultWindowsGroup()
    g.setSelection([self.meshes[str(item.text())] for item in self.selectionList.selectedItems()])

  def plotDataFromFullName(self, name):
    """Get the data from a plot using its full name (e.g. Gre_2014_DUPj : A 2)"""
    (sub, elec, plot) = self.plotNameFromFullPlotName(name)
    return self.plotsData[sub][elec][plot]

  def fullPlotName(self, subj, elec, plot):
    """Compute fully qualified plot name from subject, electrode, plot names"""
    return subj + ' : ' + elec + ' ' + plot

  def plotNameFromFullPlotName(self, name):
    """Get subject, electrode, plot names from the displayed name Subject : Electrode Plot (e.g. Gre_2014_DUPj : A 2)"""
    (sub, elecplot) = name.split(' : ')
    (elec, plot) = elecplot.split()
    return (sub, elec, plot)

  #def normalizeCoords(self):
    #"""Get plot coordinates in the selected template referential and store these coordinates in self.plotData"""

    #def miniFunc(myself, coordsSB, s):
      #"""Internal mini function to launch in a thread"""
      #refId = myself.implantations[s]['ReferentialUuid']
      #normCoords = myself.template.normalizeCoordinates([coordsSB[s][el][p] for el in coordsSB[s] for p in coordsSB[s][el]], refId)
      ##pdb.set_trace()
      #idx = 0
      #if normCoords is None:
        #print "normCoords is None for subject %s : could not convert coordinates to template referential"%s
        #return
      ##import pdb; pdb.set_trace()
      #for el in coordsSB[s]:
        #for p in coordsSB[s][el]:
          #myself.plotsData[s][el][p][myself.template.name] = normCoords[idx]
          #idx = idx + 1
    ##fin de la fonction miniFunc

    #print "go back to locateElectrode for now, we don't manage call to normalisation from locateElectrode for now"
    #return
    #coordsSB = {}
    #for s in self.subjects:
      #coordsSB[s] = {}
      #for elec, plots in self.plotsData[s].iteritems():
         #sb = dict([(p,plots[p]['Scanner-based']) for p in plots if self.template.name not in plots[p].keys()])
         #if len(sb) > 0:
           #coordsSB[s][elec] = sb
      ## Compute template coords for this subject
      #self.startTask(lambda myself=self, cSB=coordsSB, suj=s:miniFunc(myself, cSB, suj))


  #def saveNormalizedCoords(self):
    #"""Should save all computed coordinates of plots in elecimplant file for all subject
       #with the timestamp of the original data, to avoid recomputation"""
    #for s,rdi in zip(self.subjects, self.subjItems):
      #self.saveImplantation(s, rdi)
    #return

  def selectAround(self):
    """Find in the list of selected plots the ones near the linked cursor"""
    # Get linked cursor coords (in template referential)
    pTempl = self.a.linkCursorLastClickedPosition(self.template.referentialAnatomist).items()[:3]
    # Get accepted radius
    r2 = self.radiusSpin.value()**2
    meshes = []
    # Compute distance to all selected plots
    pdb.set_trace() #need to check if need to take absolute value in x in mni
    for i in xrange(self.selectionList.count()):
      fullname = str(self.selectionList.item(i).text())
      (sub, elec, plot) = self.plotNameFromFullPlotName(fullname)
      coords = self.plotsData[sub][elec][plot][self.template.name][:3]
      dist2 = (coords[0]-pTempl[0])**2 + (coords[1]-pTempl[1])**2 + (coords[2]-pTempl[2])**2
      # Select item if in range, deselect if not
      self.selectionList.item(i).setSelected(dist2 <= r2)
      if fullname in self.meshes and dist2 <= r2:
        meshes.append(self.meshes[fullname])
    # Select them in Anatomist
    g = self.a.getDefaultWindowsGroup()
    g.setSelection(meshes)

  def setSubjects(self, names, diskitems):
    """Sets the list of subjects (and corresponding readdiskitems in the database) for the widget"""
    self.subjects = sorted(names)
    self.subjItems = diskitems
    self.loadImplantations()
    self.plotsData = dict([(s, self.getPlotDataFromImplantation(s)) for s in self.subjects])
    #remove NonType from plotsData
    pat_to_remove = []
    for jj,kk in self.plotsData.iteritems():
      if kk is None:
         #self.plotsData.pop(jj,None)
         pat_to_remove.append(jj)

    for ii in range(len(pat_to_remove)):
      self.plotsData.pop(pat_to_remove[ii],None)

    if len(self.plotsData.keys()) == 0:
      print "No data to show"
      return
    self.updateUIplots()

  def addSelection(self):
    """Add the selected plots/subjects/electrodes to the selection """
    current = [str(self.selectionList.item(i).text()) for i in xrange(self.selectionList.count())] # FIXME pas juste les selected ! Tous les items
    new = [str(s.text()) for s in self.plotList.selectedItems() if str(s.text()) not in current]

    # Display the new ones
    meshes = []
    invalid = set()
    for n in new:
      if self.template.name in self.plotDataFromFullName(n):
        mesh =  self.displaySphereAt(self.plotDataFromFullName(n)[self.template.name], self.plotDiameter(), self.template.referentialAnatomist, color=(0.0,0.9,0.1,1.0),name = n)
        self.meshes[n] = mesh
        meshes.append(mesh)
      else:
        invalid.add(n) # If there, coordinates are not available in the right template referential
    if len(invalid) > 0:
      print "Some plots were not added to the selection, because normalized coordinates were not available for them"
      new = list(set(new) - invalid)
    self.selectionList.addItems(new)
    self.a.addObjects(meshes, self.windows)

  def plotDiameter(self):
    """Returns the diameter of the spheres used to display plots"""
    return 2.0

  def displaySphereAt(self, center, diameter, referential, color, name = None):
    """Returns a spherical mesh (anatomist object) with color = [1.0,0.0,0.0,1.0] for a red, not transparent sphere"""
    mesh = self.a.toAObject(aims.SurfaceGenerator.sphere(aims.Point3df(center[0], center[1], center[2]), diameter, 54))
    if name is not None:
        mesh.setName(name)
    self.a.setMaterial(mesh, diffuse=color)#[color.redF(), color.greenF(), color.blueF(), color.alphaF()] #sortir le setMaterial(diffuse=color) et le mettre à la fin de la boucle for des fonctions qui l'appelle ?
    self.a.assignReferential(referential, mesh)
    return mesh

  def removeSelection(self):
    """Remove plots from the selected list"""
    removable = [str(s.text()) for s in self.selectionList.selectedItems()]
    meshes = [self.meshes[r] for r in removable if r in self.meshes]
    for r in removable:
      if r in self.meshes:
        del self.meshes[r]
    self.a.removeObjects(meshes, self.windows)
    self.a.deleteObjects(meshes)
    # Reverse loop to remove from the bottom (avoids messing up the index)
    for idx in reversed(range(self.selectionList.count())):
      if str(self.selectionList.item(idx).text()) in removable:
        self.selectionList.takeItem(idx)

  def removeNotTRC(self):
    """Remove plots which are not registered in the TRC"""
    all_items=[str(self.selectionList.item(i).text()) for i in range(self.selectionList.count())]
    #on remplace les ' par des p dans all_items
    all_items = [all_items[x].replace("'","p") for x in range(len(all_items))]
    full_list_trc=[]
    for subj in self.subjects:
       #check if exist TRC in DB for this subject
       rdi = ReadDiskItem('Raw SEEG recording', 'EEG TRC format' )
       di = rdi.findValue({'subject':subj})
       if di is not None:
         data_micromed = neo.MicromedIO(filename = str(di)).read_segment() #all data
         taille=len(data_micromed.analogsignals)
         ##Normalisation name (between analogsignals' name and plots' name)
         number=['01', '02', '03', '04', '05', '06', '07', '08', '09']
         noms=[]
         for i in range(taille):
            name=data_micromed.analogsignals[i].name
            name=name.upper()
            if name[len(name)-2:] in number:
              name=name[:len(name)-2]+name[len(name)-1:]
            setattr(data_micromed.analogsignals[i],'name',name)
            noms+=[data_micromed.analogsignals[i].name]

         #re.sub(noms[16],re.findall('\d+',noms[16])[0],subj + " : " + "Plot"+re.findall('\d+',noms[16])[0])
         #noms_remade= [subj + " : " +  re.findall('\S+(?<![\d_])',noms[x])[0] + " Plot"+re.findall('\d+',noms[x])[0] for x in taille]
         noms_remade= [subj + " : " +  re.findall('\S+(?<![\d_])',noms[x])[0] + " Plot"+re.findall('\d+',noms[x])[0] for x in range(len(noms)) if len(re.findall('\d+',noms[x])) > 0]
         #on remplace les ' par des p dans noms_remade
         [full_list_trc.append(noms_remade[x].replace("'","p")) for x in range(len(noms_remade))]


    to_keep=[all_items[x] for x in range(len(all_items)) for y in range(len(full_list_trc)) if all_items[x]==full_list_trc[y]]
    to_remove=list(set(all_items)-set(to_keep))

    meshes = [self.meshes[r] for r in to_remove if r in self.meshes]
    for r in to_remove:
      if r in self.meshes:
        del self.meshes[r]
    self.a.removeObjects(meshes, self.windows)
    self.a.deleteObjects(meshes)
    # Reverse loop to remove from the bottom (avoids messing up the index)
    for idx in reversed(range(self.selectionList.count())):
       if str(self.selectionList.item(idx).text()) in to_remove:
         self.selectionList.takeItem(idx)



  def removePlotsLeftRight(self, side):


     all_items=[str(self.selectionList.item(i).text()) for i in range(self.selectionList.count())]

     to_remove = []
     for ii in all_items:
       MNI_pos = self.plotDataFromFullName(ii)['MNI']

       if side == 'Left':
          if MNI_pos[0] >= 0:
              to_remove.append(ii)
       elif side == 'Right':
           if MNI_pos[0] <=0:
              to_remove.append(ii)

     meshes = [self.meshes[r] for r in to_remove if r in self.meshes]
     for r in to_remove:
       if r in self.meshes:
         del self.meshes[r]

     self.a.removeObjects(meshes,self.windows)
     self.a.deleteObjects(meshes)

     for idx in reversed(range(self.selectionList.count())):
       if str(self.selectionList.item(idx).text()) in to_remove:
         self.selectionList.takeItem(idx)



  def updateUIplots(self):
    self.subjectList.clear()
    self.subjectList.addItems(self.subjects)
    allElecs = sorted([s + ' : ' + el for s in self.plotsData.keys() for el in self.plotsData[s].keys()])
    self.electrodeList.clear()
    self.electrodeList.addItems(allElecs)
    allPlots = sorted([self.fullPlotName(s, el, pl) for s in self.plotsData.keys() for el in self.plotsData[s].keys() for pl in self.plotsData[s][el].keys() ])
    self.plotList.clear()
    self.plotList.addItems(allPlots)

  def t1pre2ScannerBased(self, subject):
    """ Returns a triplet of Anatomist objects (native T1pre referential, scanner-base T1pre referential, Transformation from T1pre referential to T1pre Scanner-Based referential) """
    rdi = ReadDiskItem('Transformation to Scanner Based Referential', 'Transformation matrix', exactType=True,\
      requiredAttributes={'modality':'t1mri', 'subject':subject})
    allTransf = list (rdi._findValues( {}, None, False ) )
    for trsf in allTransf:
      if trsf.attributes()['acquisition'].startswith(u'T1pre'):
        print repr(trsf.attributes())
        srcrDiskItem = self.transfoManager.referential( trsf.attributes()['source_referential'] )
        srcr = self.a.createReferential(srcrDiskItem)
        dstrDiskItem = self.transfoManager.referential(trsf.attributes()['destination_referential'])
        self.t1pre2ScannerBasedId = trsf.attributes()['destination_referential']
        dstr = self.a.createReferential(dstrDiskItem)
        return  (srcr, dstr, self.a.loadTransformation(trsf.fullPath(), srcr, dstr))
    return None

  def loadImplantations(self):
    self.implantations = dict([(s,self.loadImplantation(rdi)) for s,rdi in zip(self.subjects, self.subjItems)])

  def loadImplantation(self, rdiSuj):
    rdi = ReadDiskItem( 'Electrode implantation', 'Electrode Implantation format')
    impl = rdi.findValue(rdiSuj)

    if not impl:
      print "Cannot find implantation for %s"%rdiSuj.attributes()['subject']
      return {}
    if (os.path.exists(str(impl))):
	    filein = open(str(impl), 'rb')
	    try:
	       dic = json.loads(filein.read())
	    except:
	       filein.close()
	       filein = open(str(impl), 'rb')
	       dic = pickle.load(filein)

	    filein.close()
	    
    #we load eleclabel now if exist
    rdi_eleclabel = ReadDiskItem('Electrodes Labels','Electrode Label Format')
    impl_label = rdi_eleclabel.findValue(rdiSuj)
    
    if not impl_label:
        print("Cannot find implantation label for %s"%rdiSuj.attributes()['subject'])
        pass
    else:
        if (os.path.exists(str(impl_label))):
           filein = open(str(impl_label),"rb")
           try:
             dic2 = json.loads(filein.read())
           except:
             filein.close()
             filein.open(str(impl_label),"rb")
             dic2 = pickle.load(filein)
           
           filein.close()
           dic.update({'label':dic2['plots_label']})       
    
    return dic

    #print "Exception while reading implantation file for %s"%rdiSuj.attributes()['subject']
    #return {}

  def saveImplantation(self, subj, rdiSubj):
    wdi = WriteDiskItem( 'Electrode implantation', 'Electrode Implantation format')
    impl = wdi.findValue(rdiSubj)
    if impl is None:
      print "Could not find electrode implantation file to save to (%s) !"%subj
      return
    try:
      #import pdb; pdb.set_trace()
      fileout = open(impl.fullPath()+'.temporary', 'wb')
      content = self.implantations[subj]
      content['plotsData-timestamp'] = content['timestamp']
      content['plotsData'] = self.plotsData[subj]
      fileout.write(json.dumps(content))
      #pickle.dump(content, fileout)
      fileout.close()
      #to modify to json
      shutil.move(impl.fullPath()+'.temporary', impl.fullPath())
      neuroHierarchy.databases.insertDiskItem( impl, update=True )
    except:
      print "Exception while writing implantation file for %s"%subj
      traceback.print_exc(file=sys.stdout)
      return



  def getPlotDataFromImplantation(self, subj):
    els = self.subjectElectrodes(subj)
    #self.addElectrode(e['name'], e['model'], e['target'], e['entry'], refId)
    if 'plotsData' in self.implantations[subj]:
      if self.implantations[subj]['plotsData-timestamp'] == self.implantations[subj]['timestamp']:
        print "Using pre-recorded plots coordinates for %s"%subj
        return self.implantations[subj]['plotsData']
      else:
        print "PlotsData timestamp was invalid for %s"%subj
    res = {}
    for e in els:
      res[e['name']] = self.getPlotsFromElectrode(e, subj)

    print(subj)
    if "plotsMNI" in self.implantations[subj].keys():
      info_plotsMNI = dict(self.implantations[subj]['plotsMNI'])
      for kk,vv in res.iteritems():
        for ll,ww in res[kk].iteritems():
            ww.update({"MNI":info_plotsMNI[kk+"%02d"%int(ll[4:])]})
            if "label" in self.implantations[subj].keys():
                try:
                 ww.update({"label":self.implantations[subj]["label"][kk+"%02d"%int(ll[4:])]})
                except:
                 pdb.set_trace()
    else:
   		print "Error MNI Coordinates"
   		QtGui.QMessageBox.warning(self, "Error", "MNI coordinates haven't been generated for electrode contacts of Subject: {}\nThey have to be generated using locateElectrodes".format(subj))
   		return
    
    return res

  def getPlotsFromElectrode(self, el, subj):
    print "Creating electrode model for %s"%subj
    traceback.print_exc(file=sys.stdout)
    (nativeRef, sbRef, t1pre2ScannerBased) = self.t1pre2ScannerBased(subj)
    (newRef, transf, elecModel) = createElectrode(el['target'], el['entry'], nativeRef, ana=self.a, model = self.electrodeModels[str(el['model'])].fullPath(), dispMode = self.dispMode, dispParams = self.dispParams)
    plots = getPlots(elecModel)
    pNames = getPlotsNames(elecModel)
    return dict([(n, {'internal':plots[n]['center'], 'native':list(transf.transform(plots[n]['center'])), 'Scanner-based':list(t1pre2ScannerBased.transform(transf.transform(plots[n]['center'])))}) for n in pNames])

  def getSubjectImplantation(self, subj):
    """Returns electrode implantation data for the subject (dictionary from the elecimplant file)"""
    return self.implantations[subj]

  def subjectElectrodes(self, subj):
    """Returns the list available electrodes for the subject"""
    impl = self.getSubjectImplantation(subj)
    if 'electrodes' in impl:
      return impl['electrodes']
    return []

  def subjectElectrodesNames(self, subj):
    """Returns the list of names of available electrodes for the subject"""
    return [el['name'] for el in self.subjectElectrodes()]

  def subjectPlot(self, subj, electrode):
    """Returns the list of plots of the chosen electrode for the subject"""
    impl = self.getSubjectImplantation(subj)
    if 'electrodes' in impl:
      els = [e for e in impl['electrodes'] if e['name'] == electrode]
      if len(els) > 0:
        els[0] # (e['name'], e['model'], e['target'], e['entry'], refId)
    return None

  def getPlotsCoordinates(self, subj, referential = None, electrode = None, plot = None):
    """Returns the coordinates of the centers of all plots of the subject, or only for the given electrodes/plots
       If referential is None, coordinates are returned in the native referential (T1pre of the subject)
    """
    return None

  def updatePlotSelected(self, item = None):
    try:
      if item is not None:
        xyz = self.plotDataFromFullName(str(item.text()))[self.template.name]
        # setLinkedCursor uses window referential : must apply transform before setting the position
        self.windows[0].moveLinkedCursor(xyz)
      else:
          print "Error moving the cursor to the contact2"
    except Exception as e:
      print "Error moving the cursor to the contact"
      #pdb.set_trace()

  def addMNIImagetoList(self,path_fichier = None):

      if path_fichier is None or path_fichier is False:
        fichier = QtGui.QFileDialog.getOpenFileName(self, "Opening file: ", "", "(*.nii *.gii *.img *.nii.gz)")
      elif not os.path.isfile(path_fichier):
        fichier = QtGui.QFileDialog.getOpenFileName(self, "Opening file: ", "", "(*.nii *.gii *.img *.nii.gz)")
      else:
        fichier = path_fichier
          
      image_mni_ref = self.a.loadObject(self.template.volumes[0])
      #image_mni_ref.loadReferentialFromHeader()
      #try:
          #self.a.removeObjects([self.currentImage,],self.windows)
      #except:
          #pass
      self.currentImage = self.a.loadObject(str(fichier))
      #self.currentImage.loadReferentialFromHeader()
      self.a.execute('LoadReferentialFromHeader', objects=[image_mni_ref,self.currentImage])
      all_trans = self.a.getTransformations()
      trans_from_vols = []
      tm=registration.getTransformationManager()
      tm.referential(registration.talairachMNIReferentialId)
      for vol in (self.currentImage, image_mni_ref):
          trans_from_vol = [t for t in all_trans if t.source() == vol.referential and not t.isGenerated()]
          # hope trans_from_vol1 contains just one transform
          # but if there are several, try to select the one going to
          # scanner-based
          if len(trans_from_vol) > 1:
            trans_from_vol_filt = [t for t in trans_from_vol if t.destination().header()['name'].startswith('Scanner-based anatomical coordinates')]
            if len(trans_from_vol_filt) == 1:
               trans_from_vol = trans_from_vol_filt
          if len(trans_from_vol) == 0:
            raise RuntimeError('could not find a non-ambiguous transform')
          elif len(trans_from_vol) > 1:
            print "There is more than one available transformation ... we take the first one and pray"
            trans_from_vol[0] = trans_from_vol_filt[0]
          trans_from_vols.append(trans_from_vol)

      #pdb.set_traceI()
      trans_from_vol1, trans_from_vol2 = trans_from_vols
      self.template.volumes.append(str(fichier))
      self.a.execute('LoadTransformation',origin=trans_from_vol1[0].destination(),destination=trans_from_vol2[0].destination(),matrix=[0, 0, 0, 1, 0, 0,  0, 1, 0,  0, 0, 1])
      self.a.addObjects([self.currentImage,], self.windows)
      self.displayCombo.addItems([os.path.split(str(fichier))[1]])

  def addMNIMeshTexttoList(self):
      
      #fichierMesh = QtGui.QFileDialog.getOpenFileName(self, "Opening mesh (surface corresponding to the texture): ", "", "(*.gii)")
      
      ##check if the file exist
      #if not os.path.isfile(fichierMesh):
          #print("the file doesn't exist")
          #return
      
      #self.addMNIImagetoList(fichierMesh)
      
      #ask for a texture gii or a functionalTractography file
      texture_info = QtGui.QMessageBox(self)
      texture_info.setText("Choose the type of texture format (gii or csv to generate the gii)")
      texture_info.setWindowTitle("texture format")
      gii_button = texture_info.addButton(QtGui.QPushButton('.gii'),QtGui.QMessageBox.AcceptRole)
      csvfuncTract_button =texture_info.addButton(QtGui.QPushButton('.csv functionalTractography'),QtGui.QMessageBox.AcceptRole)
      #center_seg.setWindowModality(QtCore.Qt.NonModal)
      texture_info.show()
      texture_info.exec_()
      #reply = texture_info.buttonRole(texture_info.clickedButton())
      if str(texture_info.clickedButton().text())=='.gii':
          print("texture already gii generated")
          fichierTexture =  QtGui.QFileDialog.getOpenFileName(self, "Opening texture (corresponding to the mesh): ", "", "(*.gii)")
          pdb.set_trace()
      elif str(texture_info.clickedButton().text())=='.csv functionalTractography':
          print("have to generate the gii texture from the csv data, functionalTractography csv model")
          fichierCSV =  QtGui.QFileDialog.getOpenFileName(self, "Opening functional tractography data: ", "", "(*.csv)")
          if not os.path.isfile(fichierCSV):
              print("the file doesn't exist")
          full_data = readFunctionalTractography(fichierCSV)

          #ask where to save the data
          path_to_save = QtGui.QFileDialog.getExistingDirectory(self,'Directory to save the mesh')
          
          BrodmannParcels = aims.read('MNI_Atlases/rbrodmann.nii')
          BrodmannParcelsArrayData = BrodmannParcels.arraydata()
          left_white = aims.read('MNI_Brainvisa/t1mri/T1pre_1900-1-3/default_analysis/segmentation/mesh/Gre_2016_MNI1_Lwhite.gii')
          right_white = aims.read('MNI_Brainvisa/t1mri/T1pre_1900-1-3/default_analysis/segmentation/mesh/Gre_2016_MNI1_Rwhite.gii')         
          
          list_remove = set(['Patient', 'Atlas'])
          list_condi_max = set(['ValueNb','PeakDelayMed','PeakDelaySTD','Probability'])
          list_condi_present = set(full_data.keys())
          
          #condi_intersect = list(list_condi_max & list_condi_present)
          condi_intersect = sorted(list(list_condi_present-list_remove), key=lambda s: s.lower())
          nb_time = len(condi_intersect)

          orderTexture = dict([(i,condi_intersect[i]) for i in range(len(condi_intersect))]) #{0:'ValueNb',1:'PeakDelayMed',2:'PeakDelaySTD',3:'Probability'}

          if full_data['Atlas'] == 'MarsAtlas':
            #read the marsAtlas parcellation     
            left_MA = aims.read('MNI_Brainvisa/t1mri/T1pre_1900-1-3/default_analysis/segmentation/mesh/surface_analysis/Gre_2016_MNI1_Lwhite_parcels_marsAtlas.gii')
            right_MA = aims.read('MNI_Brainvisa/t1mri/T1pre_1900-1-3/default_analysis/segmentation/mesh/surface_analysis/Gre_2016_MNI1_Rwhite_parcels_marsAtlas.gii')

            #for i in range(len(condi_intersect)):
            #    left_white.vertex(i).assign(left_white.vertex(0))
            #    left_white.normal(i).assign(left_white.normal(0))
            #    left_white.polygon(i).assign(left_white.polygon(0))
            #    right_white.vertex(i).assign(right_white.vertex(0))
            #    right_white.normal(i).assign(right_white.normal(0))
            #    right_white.polygon(i).assign(right_white.polygon(0))
              
            aims.write(left_white,str(path_to_save) + os.path.sep + 'left_white.gii')
            aims.write(right_white,str(path_to_save) + os.path.sep + 'right_white.gii')

            try:
              os.mkdir(str(path_to_save)+os.path.sep+'Texture')
            except:
              pass          

            #faire un test si all au lieu des noms de parcels.
            for i_parcels_stimulated in full_data[orderTexture[0]].keys():
              
              if len( full_data[orderTexture[0]][i_parcels_stimulated]) > 0:
                new_TimeSurfTextLeft = aims.TimeTexture('FLOAT')
                new_TimeSurfTextRight = aims.TimeTexture('FLOAT')
            
                for i in range(nb_time):   #assign the value
                  textnowLeft = new_TimeSurfTextLeft[i]
                  textnowRight = new_TimeSurfTextRight[i]
                    
                  textnowLeft.reserve(len(left_white.vertex(0))) #left_white.vertex(0)))
                  textnowRight.reserve(len(right_white.vertex(0)))
                  marsatlas_label = readSulcusLabelTranslationFile('parcels_label_name.txt')
              
                  #gauche #control lateral lorsqu'étude contro/ipsi
                  for iter_vert in range(len(left_white.vertex(0))):
                    #marsatlas_label[left_MA[0].arraydata()[iter_vert]]
                    #if isinstance(full_data[orderTexture[i]]['L_VCcm'][marsatlas_label[left_MA[0].arraydata()[iter_vert]]], (str, unicode)):
                    
                    #join right and left or not 
                    #for now we assume that we are using marsatlas
                    actual_marsatlas_parcels = []

                    if i_parcels_stimulated.startswith('L_') or i_parcels_stimulated.startswith('R_'):
                        try:
                            actual_marsatlas_parcels = [marsatlas_label[left_MA[0].arraydata()[iter_vert]]]
                        except:
                            pass #faudrait mieux mettre si c'est == 0 alors c'est un vertex qui n'a pas de correspondance marsatlas.
                    elif i_parcels_stimulated == 'All':
                        try:
                            actual_marsatlas_parcels = [marsatlas_label[left_MA[0].arraydata()[iter_vert]]]
                        except:
                            pass #faudrait mieux mettre si c'est == 0 alors c'est un vertex qui n'a pas de correspondance marsatlas.                        
                    else:
                        try:
                          actual_marsatlas_parcels = [[marsatlas_label[left_MA[0].arraydata()[iter_vert]]][0][2:]]
                        except:
                          pass
                    
                    try:
                        if left_MA[0].arraydata()[iter_vert] == 0:
                          textnowLeft.append(-4)
                        else:
                           if full_data[orderTexture[i]][i_parcels_stimulated].keys()[0].startswith('i_') or full_data[orderTexture[i]][i_parcels_stimulated].keys()[0].startswith('c_'):
                             if full_data[orderTexture[i]][i_parcels_stimulated]['c_'+actual_marsatlas_parcels[0]] == 'NaN':
                               textnowLeft.append(-4)
                             else:
                               textnowLeft.append(float(full_data[orderTexture[i]][i_parcels_stimulated]['c_'+actual_marsatlas_parcels[0]]))
                           else:
                             if full_data[orderTexture[i]][i_parcels_stimulated][actual_marsatlas_parcels[0]] == 'NaN':
                               textnowLeft.append(-4)
                             else:
                               textnowLeft.append(float(full_data[orderTexture[i]][i_parcels_stimulated][actual_marsatlas_parcels[0]]))  
                    except:
                        textnowLeft.append(-4)
                        #pdb.set_trace()
                        
                  #puis droite #ipsi lateral lorsqu'étude contro_ipsi
                  for iter_vert in range(len(right_white.vertex(0))):
                    
                    actual_marsatlas_parcels = []

                    if i_parcels_stimulated.startswith('L_') or i_parcels_stimulated.startswith('R_'):
                        try:
                          actual_marsatlas_parcels = [marsatlas_label[right_MA[0].arraydata()[iter_vert]]]
                        except:
                          pass
                    elif i_parcels_stimulated == 'All':
                        try:
                            actual_marsatlas_parcels = [marsatlas_label[right_MA[0].arraydata()[iter_vert]]]
                        except:
                            pass #faudrait mieux mettre si c'est == 0 alors c'est un vertex qui n'a pas de correspondance marsatlas.                      
                    else:
                        try:
                          actual_marsatlas_parcels = [[marsatlas_label[right_MA[0].arraydata()[iter_vert]]][0][2:]]
                        except:
                          pass
                    
                    try:
                        if right_MA[0].arraydata()[iter_vert] == 0:
                          textnowRight.append(-4)
                        else:
                          if full_data[orderTexture[i]][i_parcels_stimulated].keys()[0].startswith('i_') or full_data[orderTexture[i]][i_parcels_stimulated].keys()[0].startswith('c_'):
                            if full_data[orderTexture[i]][i_parcels_stimulated]['i_'+actual_marsatlas_parcels[0]] == 'NaN':
                              textnowRight.append(-4)
                            else:
                              textnowRight.append(float(full_data[orderTexture[i]][i_parcels_stimulated]['i_'+actual_marsatlas_parcels[0]]))
                          else:
                            if full_data[orderTexture[i]][i_parcels_stimulated][actual_marsatlas_parcels[0]] == 'NaN':
                              textnowRight.append(-4)
                            else:
                              textnowRight.append(float(full_data[orderTexture[i]][i_parcels_stimulated][actual_marsatlas_parcels[0]]))
                              
                    except:
                        textnowRight.append(-4)
                        #pdb.set_trace()                    
          
            
                aims.write(new_TimeSurfTextLeft,str(path_to_save) + os.path.sep + 'Texture' + os.path.sep + ('%s_left.gii')%i_parcels_stimulated)
                aims.write(new_TimeSurfTextRight,str(path_to_save) + os.path.sep + 'Texture' + os.path.sep + ('%s_right.gii')%i_parcels_stimulated)
              
              else:
                print(('No Data for %s')%i_parcels_stimulated)
                
            obj1 = self.a.loadObject('MNI_Brainvisa/t1mri/T1pre_1900-1-3/default_analysis/segmentation/mesh/Gre_2016_MNI1_Lwhite.gii')
            obj2 = self.a.loadObject('MNI_Brainvisa/t1mri/T1pre_1900-1-3/default_analysis/segmentation/mesh/surface_analysis/Gre_2016_MNI1_Lwhite_parcels_marsAtlas.gii')
            obj1.loadReferentialFromHeader()
            obj2.setPalette(palette = 'marsatlas')
            MarsAtlas_fusion_obj = self.a.fusionObjects([obj1, obj2], method='FusionTexSurfMethod')
            self.a.addObjects(MarsAtlas_fusion_obj,self.axWindow)
            self.currentImage = [obj1,obj2]

          
          if full_data['Atlas'] == 'Brodmann':
              voxel_size_T1 = [BrodmannParcels.getVoxelSize()[0], BrodmannParcels.getVoxelSize()[1], BrodmannParcels.getVoxelSize()[2], 1.0]
              sizeOutputnii = BrodmannParcels.getSize().list()
              sizeOutputnii[-1] = len(orderTexture.keys())

              for i_parcels_stimulated in full_data[orderTexture[0]].keys():
                  #di.setMinf('ColorPalette','Blue-Red-fusion')
                  volToGenerate = aims.Volume(*sizeOutputnii,dtype = 'float')
                  volToGenerate.header()['voxel_size']=voxel_size_T1
                  volToGenerate.fill(-4)
                  
                  for i_texture in orderTexture.keys():

                     for i_parcels_result in full_data[orderTexture[i_texture]][i_parcels_stimulated].keys():
                        #fait chier le droite gauche
                        if full_data[orderTexture[0]][i_parcels_stimulated][i_parcels_result]=='NaN':
                          volToGenerate.arraydata()[numpy.where(BrodmannParcelsArrayData==float(i_parcels_result))]=-4
                        else:  
                          volToGenerate.arraydata()[numpy.where(BrodmannParcelsArrayData==float(i_parcels_result))]=full_data[orderTexture[i_texture]][i_parcels_stimulated][i_parcels_result]
                        
                     aims.write(volToGenerate,str(path_to_save)+os.path.sep+'%s.nii'%(str(int(float(i_parcels_stimulated)))))

            
      print("done")
      pdb.set_trace()
      #try:
      #    self.a.removeObjects([self.currentImage,],self.windows)
      #except:
      #    pass

         
      ##self.displayCombo.addItems([os.path.split(str(fichier))[1]])
      #obj3 = self.a.loadObject(str(path_to_save) + os.path.sep + 'left_white_multipletime.gii')
      #obj4 = self.a.loadObject(str(path_to_save) + os.path.sep + 'Texture' + os.path.sep + ('%s_left.gii')%i_parcels_stimulated)
      #obj3.loadReferentialFromHeader()
      #obj4.setPalette(palette = 'Blue-Red-fusion')
      #FunctioTracto_fusion_obj = self.a.fusionObjects([obj3,obj4],method='FusionTexSurfMethod')
      #self.a.addObjects(FunctioTracto_fusion_obj,self.sagWindow)
          
          

          
      #pdb.set_trace()
      #self.displayCombo.addItems([os.path.split(str(fichier))[1]])
          
      #textureContacts = aims.TimeTexture()
      
  def doubleClickedFunctionalTractography(self):
      
      pdb.set_trace()
      
  
  def generateStatisticsContacts(self):
      
      #get selected contacts
      current = [str(self.selectionList.item(i).text()) for i in xrange(self.selectionList.count())]
      
      #il me faut un dictionnaire avec toutes les parcels mars atlas et un dictionnaire avec toutes les parcels freesurfer possible.
      dict_marsatlas = {}
      dict_freesurfer = {}
      dict_dispersion_MA = {}
      dict_dispersion_FS = {}
      parcels_names = readSulcusLabelTranslationFile('parcels_label_name.txt')
      freesurfer_parcel_names = readFreesurferLabelFile('freesurfer_label.txt')
      missing_marsatlas = []
      missing_freesurfer = []
      all_patients = []
      dict_MNI_PatientName = {}
      
      for ii in parcels_names.values():
          dict_marsatlas.update({ii:[]})
          dict_dispersion_MA.update({ii:{}})
 
      
      for ii in freesurfer_parcel_names.values():
          dict_freesurfer.update({ii[0]:[]})
          dict_dispersion_FS.update({ii[0]:{}})
         
      #je parcours toutes les "current", je regarde leur parcels et j'ajoute la position mni à la list de cette parcels.          
      for ii in current:
          (sub, elec, plot) = self.plotNameFromFullPlotName(ii)
          dataplot = self.plotDataFromFullName(ii)
          if 'MarsAtlas' in dataplot['label'].keys():
              if dataplot['label']['MarsAtlas'][1] != u'not in a mars atlas parcel':
                 dict_marsatlas[dataplot['label']['MarsAtlas'][1]].append(dataplot['MNI'])
          else:
              #signaler que certains patients n'ont pas marsatlas de généré et que ça va "fausser" les résultats
              #print("plot %s without marsAtlas parcellation estimated"%(ii))
              if sub not in missing_marsatlas:
                  missing_marsatlas.append(sub)
          if 'Freesurfer' in dataplot['label'].keys():
              if dataplot['label']['Freesurfer'][1] != u'not in a freesurfer parcel':
                 try: 
                   dict_freesurfer[dataplot['label']['Freesurfer'][1]].append(dataplot['MNI'])
                 except:
                   print sub
                   print "probleme avec ce patient"
                   pass
          else:
              #signaler que certains patient n'ont pas freesurfer de généré et que ça va "fausser" les résultats
              #print("plot %s without FreeSurfer parcellation estimated"%(ii))
              if sub not in missing_freesurfer:
                  missing_freesurfer.append(sub)
          if sub not in all_patients:
              all_patients.append(sub)   
          dict_MNI_PatientName.update({str(dataplot['MNI']):sub})  
              
      #now I calculate the dispersion per parcels: (and I'ld like to normalized it but don't know how)
      for iter_MA in dict_marsatlas.keys():
        points_array = numpy.array(dict_marsatlas[iter_MA])
        #array_median = repmat(numpy.median(points_array,axis=0),len(points_array),1)
        #diff = points_array - array_median
        #list_dist_median =   sc_sp.distance.cdist([numpy.median(points_array,axis=0)],points_array)
        if len(points_array)>1:
            dict_dispersion_MA[iter_MA].update({'nb contact':len(dict_marsatlas[iter_MA]),'average point':numpy.mean(points_array,axis=0),'median point':numpy.median(points_array,axis=0)})
            found_outlier = self.is_outlier(points_array,thresh = 3)
            pos_outlier = numpy.where(found_outlier==True)
            if len(pos_outlier[0]) == 0:
               dict_dispersion_MA[iter_MA].update({'outlier position':None})
            else:
               #[dict_MNI_PatientName[str(points_array[pos_outlier[0]][i].tolist())] for i in range(len(points_array[pos_outlier[0]]))]
               dict_dispersion_MA[iter_MA].update({'outlier position':points_array[pos_outlier[0]]})
               dict_dispersion_MA[iter_MA].update({'outlier name':[dict_MNI_PatientName[str(points_array[pos_outlier[0]][i].tolist())] for i in range(len(points_array[pos_outlier[0]]))]})
        elif len(points_array)==1:
            dict_dispersion_MA[iter_MA].update({'nb contact': 1,'average point':numpy.mean(points_array,axis=0),'median point':numpy.median(points_array,axis=0)})
        else:
            dict_dispersion_MA[iter_MA].update({'nb contact': 0})

 
      for iter_FS in dict_freesurfer.keys():
        points_array = numpy.array(dict_freesurfer[iter_FS])
        #array_median = repmat(numpy.median(points_array,axis=0),len(points_array),1)
        #diff = points_array - array_median
        #list_dist_median =   sc_sp.distance.cdist([numpy.median(points_array,axis=0)],points_array)
        if len(points_array)>1:
            dict_dispersion_FS[iter_FS].update({'nb contact':len(dict_freesurfer[iter_FS]),'average point':numpy.mean(points_array,axis=0),'median point':numpy.median(points_array,axis=0)})
            found_outlier = self.is_outlier(points_array,thresh = 3)
            pos_outlier = numpy.where(found_outlier==True)
            if len(pos_outlier[0]) == 0:
               dict_dispersion_FS[iter_FS].update({'outlier position':None})
            else:
               dict_dispersion_FS[iter_FS].update({'outlier position':points_array[pos_outlier[0]]})
               dict_dispersion_FS[iter_FS].update({'outlier name':[dict_MNI_PatientName[str(points_array[pos_outlier[0]][i].tolist())] for i in range(len(points_array[pos_outlier[0]]))]})
        elif len(points_array)==1:
            dict_dispersion_FS[iter_FS].update({'nb contact':1,'average point':numpy.mean(points_array,axis=0),'median point':numpy.median(points_array,axis=0)})
        else:
            dict_dispersion_FS[iter_FS].update({'nb contact':0})

      #ecrire le tout dans un csv
      fileName = QtGui.QFileDialog.getSaveFileName(self, 'Dialog Title', '/', '*.csv') #str(QtGui.QFileDialog.getExistingDirectory(self, "Select Directory"))
      fileName=str(fileName)
      testcsv = fileName.split('.')
      if len(testcsv)>0:
          if testcsv[1] != 'csv':
              print("error, the extension should be .csv")
              return
      else:
          fileName = fileName + '.csv'
          
      with open(fileName, 'w') as csvfile:
         writer = csv.writer(csvfile, delimiter='\t')
         writer.writerow([u'Group Analysis'])
         listwrite_allpatients = [u'Patients']
         for ii in all_patients:
             listwrite_allpatients.append(ii)
         writer.writerow(listwrite_allpatients) #writer.writerow([u'Patients',all_patients])
         listwrite_missingMA = [u'missing marsAtlas info for the following patients (analysis done without their data):']
         for ii in missing_marsatlas:
             listwrite_missingMA.append(ii)
         writer.writerow(listwrite_missingMA)#writer.writerow([u'missing marsAtlas info for the following patients (analysis done without their data):',missing_marsatlas])
         listwrite_missingFS = [u'missing FreeSurfer info for the following patients (analysis done without their data):']
         for ii in missing_freesurfer:
             listwrite_missingFS.append(ii)
         writer.writerow(listwrite_missingFS)
         #writer.writerow([u'missing FreeSurfer info for the following patients (analysis done without their data):',missing_freesurfer])
         writer.writerow([u'MarsAtlas parcellation analysis'])
         writer.writerow([u'parcel name',u'nb contact', u'average point', u'median point', u'outliers positions',u'patient names of the outliers'])
        
         dictMA_sorted_tmp = OrderedDict(sorted(dict_dispersion_MA.items()))
         for kk,vv in dictMA_sorted_tmp.iteritems():
            if vv['nb contact']>0:
              listwrite = [kk]
              listwrite.append(vv['nb contact'])
              if 'average point' in vv.keys():
                listwrite.append([float(format(vv['average point'][i],'.3f')) for i in range(3)])
                listwrite.append([float(format(vv['median point'][i],'.3f')) for i in range(3)])
              if 'outlier position' in vv.keys():
                if vv['outlier position'] is not None:
                  listwrite.append(vv['outlier position'])
                  listwrite.append(vv['outlier name'])
              writer.writerow(listwrite)
         
         #writer.writerow(dict_dispersion_MA)
         writer.writerow([])
         writer.writerow([u'Freesurfer parcellation analysis'])
         writer.writerow([u'parcel name',u'nb contact', u'average point', u'median point', u'outliers positions',u'patient names of the outliers'])
         dictFS_sorted_tmp = OrderedDict(sorted(dict_dispersion_FS.items()))
         for kk,vv in dictFS_sorted_tmp.iteritems():
              if vv['nb contact']>0:
                listwrite = [kk]
                listwrite.append(vv['nb contact'])
                if 'average point' in vv.keys():
                  listwrite.append([float(format(vv['average point'][i],'.3f')) for i in range(3)])
                  listwrite.append([float(format(vv['median point'][i],'.3f')) for i in range(3)])
                if 'outlier position' in vv.keys():
                  if vv['outlier position'] is not None:
                    listwrite.append(vv['outlier position'])
                    listwrite.append(vv['outlier name'])
                writer.writerow(listwrite)
         
      print("csv done")    

  def is_outlier(self,list_points, thresh=3):
      
      #to be changed to the MAD
      if len(list_points.shape) == 1:
         list_points = list_points[:,None]
      median = numpy.median(list_points, axis=0)
      diff = numpy.sum((list_points - median)**2, axis=-1)
      diff = numpy.sqrt(diff)
      med_abs_deviation = numpy.median(diff)

      modified_z_score = 0.6745 * diff / med_abs_deviation

      return modified_z_score > thresh

  def AddMAparcels2selection(self):
      
      print "select contacts according to marsatlas parcels"
      fullPlot_List = [str(self.plotList.item(idx).text()) for idx in range(self.plotList.count())]
      
      parcels_names = readSulcusLabelTranslationFile('parcels_label_name.txt')
      
      dict_plotMA = {}
      for ii in parcels_names.values():
          dict_plotMA.update({ii:[]})
          
      for ii in fullPlot_List:
          (sub, elec, plot) = self.plotNameFromFullPlotName(ii)
          dataplot = self.plotDataFromFullName(ii)
          if 'MarsAtlas' in dataplot['label'].keys():
              if dataplot['label']['MarsAtlas'][1] != u'not in a mars atlas parcel':         
                 dict_plotMA[dataplot['label']['MarsAtlas'][1]].append(ii)
      
      if str(self.AddMAparcels2SelectioncomboBox.currentText()) =='*':
          print('unselect all')
          for i in xrange(self.plotList.count()):
            selec = False
            self.plotList.item(i).setSelected(selec)
      else:
          list_plot2select= dict_plotMA[str(self.AddMAparcels2SelectioncomboBox.currentText())]
          for i in xrange(self.plotList.count()):
            selec = False
            if str(self.plotList.item(i).text()) in list_plot2select:
                selec = True
            self.plotList.item(i).setSelected(selec)
            
            
  def changeBothRightDisplay(self):
      
      all_items=[str(self.selectionList.item(i).text()) for i in range(self.selectionList.count())]
      refBothHemi = self.template.referentialAnatomist
      newRef = self.a.createReferential()
      
      transf = self.a.createTransformation([0,0,0,-1,0,0,0,1,0,0,0,1], origin = newRef, destination = refBothHemi)
      meshesLeft = []
      for ii in all_items:
          MNI_pos = self.plotDataFromFullName(ii)['MNI']
          
          if MNI_pos[0] >=0:
              #on fait classique
              pass
          elif MNI_pos[0] < 0:
              
            if self.radioButtonbothHemi.isChecked():
              meshesLeft.append(ii)
              self.a.assignReferential(refBothHemi,self.meshes[ii])
          
          
            elif self.radioButtonAllRight.isChecked():
              meshesLeft.append(ii)
              pdb.set_trace()
              self.a.assignReferential(newRef,self.meshes[ii])
              
      #self.a.assignReferential(newRef,meshesLeft) 


     #to_remove = []
     #for ii in all_items:
       #MNI_pos = self.plotDataFromFullName(ii)['MNI']

       #if side == 'Left':
          #if MNI_pos[0] >= 0:
              #to_remove.append(ii)
       #elif side == 'Right':
           #if MNI_pos[0] <=0:
              #to_remove.append(ii)

     #meshes = [self.meshes[r] for r in to_remove if r in self.meshes]
     #for r in to_remove:
       #if r in self.meshes:
         #del self.meshes[r]

     #self.a.removeObjects(meshes,self.windows)
     #self.a.deleteObjects(meshes)
      
      
  def contactSEEGDisplay(self):
      
      if self.radioButtonContactDisplay.isChecked():
        pdb.set_trace()
        try:
            self.bipoleSEEGColors.close()
        except:
            pass
        self.a.removeObjects([self.bipolesmeshes[x] for x in self.bipolesmeshes.keys()],self.windows)
        #current = [str(self.selectionList.item(i).text()) for i in xrange(self.selectionList.count())] # FIXME pas juste les selected ! Tous les items
        #new = [str(s.text()) for s in self.plotList.selectedItems() if str(s.text()) not in current]

        ## Display the new ones
        #meshes = []
        #invalid = set()
        #for n in new:
          #if self.template.name in self.plotDataFromFullName(n):
            #mesh =  self.displaySphereAt(self.plotDataFromFullName(n)[self.template.name], self.plotDiameter(), self.template.referentialAnatomist, color=(0.0,0.9,0.1,1.0),name = n)
            #self.meshes[n] = mesh
            #meshes.append(mesh)
          #else:
            #invalid.add(n) # If there, coordinates are not available in the right template referential
        #if len(invalid) > 0:
          #print "Some plots were not added to the selection, because normalized coordinates were not available for them"
          #new = list(set(new) - invalid)
        #self.selectionList.addItems(new)
        
        self.a.addObjects([self.meshes[x] for x in self.meshes.keys()], self.windows)        
          
      
      elif self.radioButtonsEEGResults.isChecked():
   
         meshes = [self.meshes[x] for x in self.meshes.keys()]
          
          #for ind_mesh in self.meshes.keys():
              #del self.meshes[ind_mesh]
         self.a.removeObjects(meshes,self.windows)
          #self.a.deleteObjects(meshes)
         
         
         #il faut générer les bipoles
         current = [str(self.selectionList.item(i).text()) for i in xrange(self.selectionList.count())]
         info_contact={}
         for pindex in range(0,len(current)):
          (sub, elec, plot) = self.plotNameFromFullPlotName(current[pindex])
          try:
              if sub not in info_contact.keys():
                info_contact.update({sub:[]})
                info_contact[sub].append(elec+"%02d"%int(plot.split('Plot')[1]))
              else:
                info_contact[sub].append(elec+"%02d"%int(plot.split('Plot')[1]))
          except:
              pdb.set_trace()
         info_bipole = {}
         for subj in info_contact.keys():
             
             contacts_sorted = sorted(info_contact[subj],key=natural_keys)
             if subj not in info_bipole.keys():
                 info_bipole.update({subj:{}})
             for contact_index in range(1,len(contacts_sorted)):
                 previous_contact = "".join([i for i in contacts_sorted[contact_index-1] if not i.isdigit()])
                 previous_number = "".join([i for i in contacts_sorted[contact_index-1] if i.isdigit()])
                 current_contact = "".join([i for i in contacts_sorted[contact_index] if not i.isdigit()])
                 current_number = "".join([i for i in contacts_sorted[contact_index] if i.isdigit()])
                 
                 if previous_contact == current_contact:
                     if int(current_number) - int(previous_number) == 1:
                       bipole = contacts_sorted[contact_index] + ' - ' + contacts_sorted[contact_index-1]
                       mni_bipole = ((numpy.array(self.plotsData[subj][current_contact]['Plot'+str(int(current_number))]['MNI']) + numpy.array(self.plotsData[subj][current_contact]['Plot'+str(int(previous_number))]['MNI']))/2).tolist() 
                       info_bipole[subj].update({bipole:mni_bipole})
                     else:
                       #find previous contact
                       previous_number = int(current_number)-1
                       next_number = int(current_number)+1
                       bipole_with_previous = contacts_sorted[contact_index] + ' - ' + "%s%02d"%(current_contact,int(previous_number))
                       bipole_wiht_next = "%s%02d"%(current_contact,int(next_number)) + ' - ' + contacts_sorted[contact_index]
                       if 'Plot'+str(previous_number) in self.plotsData[subj][current_contact].keys():
                         if bipole_with_previous not in info_bipole[subj].keys():
                            mni_bipole = ((numpy.array(self.plotsData[subj][current_contact]['Plot'+str(int(current_number))]['MNI']) + numpy.array(self.plotsData[subj][current_contact]['Plot'+str(int(previous_number))]['MNI']))/2).tolist()
                            info_bipole[subj].update({bipole_with_previous:mni_bipole})
                             
                       if 'Plot'+str(next_number) in self.plotsData[subj][current_contact].keys():
                         if bipole_wiht_next not in self.plotsData[subj][current_contact].keys():
                              mni_bipole = ((numpy.array(self.plotsData[subj][current_contact]['Plot'+str(int(current_number))]['MNI']) + numpy.array(self.plotsData[subj][current_contact]['Plot'+str(int(next_number))]['MNI']))/2).tolist()
                              info_bipole[subj].update({bipole_wiht_next:mni_bipole})

                       
         list_to_show = []
         for subj in info_bipole.keys():
           if 'seeg_label_all' in self.testDataSubjects[subj].keys():
               
             for index_bip in info_bipole[subj].keys():

               rdiEM = ReadDiskItem('Electrode Model', 'Electrode Model format')
               listEM = list(rdiEM.findValues({},None,False))
               matches = filter((lambda x: u"bipole" in str(x)), listEM)
               if subj + ' : ' + index_bip not in self.bipolesmeshes.keys():
                 #la faudrait que je vérifie si le bipole existe vraiment parce que sinon je crée énormément de mesh pour rien ...
                 mesh = self.displaySphereAt(info_bipole[subj][index_bip],self.plotDiameter(), self.template.referentialAnatomist, color=(0.0,0.0,0.0,0.0),name = subj + ' : ' + index_bip) #je devrais sortir le color de la ?
                 mesh.setMaterial(front_face='counterclockwise')
                 self.bipolesmeshes.update({subj + ' : ' + index_bip:mesh})
               list_to_show.append(subj + ' : ' + index_bip)
           
           #else:
             #print "No SEEG stim results"
             #QtGui.QMessageBox.warning(self, "Error", "Stim report has not been generated for the subject: {}".format(subj))
             
         
         self.a.addObjects([self.bipolesmeshes[x] for x in list_to_show], self.windows)
         
         self.bipoleSEEGColors=bipoleSEEGColors(self,indv_pat = False, group_subsample = list_to_show)
         self.bipoleSEEGColors.show()
      #removable = [str(s.text()) for s in self.selectionList.selectedItems()]
      #meshes = [self.meshes[r] for r in removable if r in self.meshes]
      #for r in removable:
        #if r in self.meshes:
        #del self.meshes[r]
      #self.a.removeObjects(meshes, self.windows)
      #self.a.deleteObjects(meshes)
      

if __name__ == "__main__":
  app = QtGui.QApplication(sys.argv)
  axon.initializeProcesses()
  from brainvisa.data.readdiskitem import ReadDiskItem
  from brainvisa.data.writediskitem import WriteDiskItem
  window = ElectrodeDisplayWidget()
  window.show()
  sys.exit(app.exec_())
