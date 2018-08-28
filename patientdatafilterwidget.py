#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Widget to select patients from the DB using available data as filters : needs an access
# to brainvisa database that is already set up ! (brainvisa axon.initializeProcesses()
#
# (c) Inserm U836 2012-2014 - Manik Bhattacharjee
#
# License GNU GPL v3
#
#
from soma.qt_gui.qt_backend import QtGui, QtCore, uic

import sys, json, pickle, numpy, csv
import math
from collections import OrderedDict
import pdb

from soma import aims
from brainvisa import axon
from brainvisa.data.readdiskitem import ReadDiskItem
from brainvisa.data.writediskitem import WriteDiskItem
from readElecLocalCSVFile import readElecLocalCSVFile
from readSulcusLabelTranslationFile import *

class PatientDataFilterWidget (QtGui.QWidget):
  def __init__(self, app=None):
    QtGui.QWidget.__init__(self)
    self.ui = uic.loadUi("patientDataFilter.ui", self)
    self.pushedWidgets = []
    self.subjects = {}
    #self.selectedSubjects=[]
    self.filters=[]
    self.filtersAtlases = {}
    self.pushedWidget = []
    self.greenDotIcon = QtGui.QIcon('greendot.png')
    self.filterProtocolCombo.currentIndexChanged.connect(self.updatePatientFilters)
    self.filterSiteCombo.currentIndexChanged.connect(self.updatePatientFilters)
    self.filterYearCombo.currentIndexChanged.connect(self.updatePatientFilters)
    self.filterLocaMarsAtlas.currentIndexChanged.connect(self.updatePatientFilters)
    
    self.filterCognitionLoca.currentIndexChanged.connect(self.updatePatientFilters)
    self.filterCognitionType.currentIndexChanged.connect(self.updatePatientFilters)
    
    self.filterResection.currentIndexChanged.connect(self.updatePatientFilters)

    self.filterStimulationLoca.currentIndexChanged.connect(self.updatePatientFilters)
    self.filterStimulationFreq.currentIndexChanged.connect(self.updatePatientFilters)    
    self.filterStimulationType.currentIndexChanged.connect(self.updatePatientFilters)
    
    self.FilterAtlascomboBox.currentIndexChanged.connect(self.changeAtlasfiltration)

    self.filterAddPatientButton.clicked.connect(lambda:self.moveSelectedItemsToOtherListWidget(self.filteredPatientList, self.selectedPatientList))
    self.filterRemovePatientButton.clicked.connect(lambda:self.moveSelectedItemsToOtherListWidget(self.selectedPatientList, self.filteredPatientList))
    self.filterAddGreenPatientButton.clicked.connect(lambda:self.moveFilteredItemsToOtherListWidget(self.filteredPatientList, self.selectedPatientList, self.filters, reverse=False))
    self.filterRemoveRedPatientButton.clicked.connect(lambda:self.moveFilteredItemsToOtherListWidget(self.selectedPatientList, self.filteredPatientList, self.filters, reverse=True))

    self.patientSelectAllButton.clicked.connect(lambda:self.selectAllListItems(self.filteredPatientList, True) )
    self.patientDeselectAllButton.clicked.connect(lambda:self.selectAllListItems(self.filteredPatientList, False))
    self.selectedSelectAllButton.clicked.connect(lambda:self.selectAllListItems(self.selectedPatientList, True))
    self.selectedDeselectAllButton.clicked.connect(lambda:self.selectAllListItems(self.selectedPatientList, False))

    self.t1preButton.clicked.connect(lambda: self.toggleFilterWidget('T1pre', self.t1preButton))
    self.t1postButton.clicked.connect(lambda: self.toggleFilterWidget('T1post', self.t1postButton))
    self.t2preButton.clicked.connect(lambda: self.toggleFilterWidget('T2pre', self.t2preButton))
    self.t2postButton.clicked.connect(lambda: self.toggleFilterWidget('T2post', self.t2postButton))
    self.petpreButton.clicked.connect(lambda: self.toggleFilterWidget('PETpre', self.petpreButton))
    self.ctpostButton.clicked.connect(lambda: self.toggleFilterWidget('CTpost', self.ctpostButton))

    self.implantationButton.clicked.connect(lambda: self.toggleFilterWidget('implantation', self.implantationButton))
    self.mniButton.clicked.connect(lambda: self.toggleFilterWidget('MNI', self.mniButton))
    self.mniParcelsButton.clicked.connect(self.mniParcelsGeneration)
    self.seegButton.clicked.connect(self.toggleSeeg)
    self.seegManipCombo.currentIndexChanged.connect(self.seegManipChanged)

    self.t1preButton.setStyleSheet("background-color: rgb(90, 90, 90);")
    self.t1postButton.setStyleSheet("background-color: rgb(90, 90, 90);")
    self.t2preButton.setStyleSheet("background-color: rgb(90, 90, 90);")
    self.t2postButton.setStyleSheet("background-color: rgb(90, 90, 90);")
    self.petpreButton.setStyleSheet("background-color: rgb(90, 90, 90);")
    self.ctpostButton.setStyleSheet("background-color: rgb(90, 90, 90);")
    self.implantationButton.setStyleSheet("background-color: rgb(90, 90, 90);")
    self.mniButton.setStyleSheet("background-color: rgb(90, 90, 90);")
    self.seegButton.setStyleSheet("background-color: rgb(90, 90, 90);")


    self.populateFromDB()


  def getSelectedPatientsNames(self):
    """Returns a list of names (string) of all selected patients"""
    return sorted([str(self.selectedPatientList.item(i).text()) for i in xrange(self.selectedPatientList.count())])

  def getSelectedPatients(self):
    """ Returns a list of ReadDiskItems representing the selected subjects"""
    return [self.subjects[p]['rdi'] for p in self.getSelectedPatientsNames()]

  def findExistingVolumeAcquisitions(self, volType, subjectName, subjectData):
    """For a subjectName, find all acquisitions in the DB of the type volType ('Raw T1 MRI') and store them in subjectData[acquisition]"""
    rdi = ReadDiskItem( volType, 'aims readable volume formats', requiredAttributes={'subject':subjectName, 'center':subjectData['center']} )
    volumes = (list( rdi._findValues( {}, None, False ) ))
    for v in volumes:
      try:
        acq = v.attributes()['acquisition'].split('_')[0]
        subjectData[str(acq)] = v
      except:
        print "Acquisition not defined for "+repr(v.fullPath())


  def populateFromDB(self):
    """Fills the list of patients and necessary information about patients to be able to filter"""
    rdi = ReadDiskItem( 'Subject', 'Directory',requiredAttributes={'_ontology':'brainvisa-3.2.0'} )
    subjects = list( rdi._findValues( {}, None, False ) )
    self.subjects = dict([(s.attributes()['subject'], {'rdi':s, 'center':s.attributes()['center']}) for s in subjects])

    protocols = list(set([s.attributes()['center'] for s in subjects]))

    # Fill the combos
    self.filterProtocolCombo.clear()
    self.filterProtocolCombo.addItems(['*',]+sorted(protocols))

    # Update the basic filters
    sites = ['*',] + sorted(set([s.split('_')[0] for s in self.subjects]))
    years = ['*',] + sorted(set([s.split('_')[1] for s in self.subjects if len(s.split('_')) > 1]))

    loca = ['*']
    loca_FS = ['*']
    loca_BM = ['*']
    cogniType = ['*']
    resec = ['*']
    resec_FS = ['*']
    resec_BM = ['*']
    stimFreq = ['*']
    stimType = ['*']

    self.filterSiteCombo.clear()
    self.filterSiteCombo.addItems(sites)
    self.filterYearCombo.clear()
    self.filterYearCombo.addItems(years)

    #self.filterLocaMarsAtlas.addItems(loca)
    #self.filterCognition.addItems(cogni)

    # Get info for all filters

    seegManips = set()
    for s,subj in self.subjects.iteritems():
      print(s)
      # Find the T1s/T2s/CT/PET
      for modality in ['T2 MRI', 'CT', 'PET', 'Raw T1 MRI']:
        self.findExistingVolumeAcquisitions(modality, s, subj)

      # Find MNI transform if T1pre is there
      if 'T1pre' in subj:
        rdi = ReadDiskItem( 'SPM2 normalization matrix', 'Matlab file' )
        di = rdi.findValue(subj['T1pre'])
        rdi2 = ReadDiskItem('SPM normalization deformation field', 'NIFTI-1 image')
        di2 = rdi2.findValue(subj['T1pre'])
        if di is not None:
          subj['MNI'] = di
        if di2 is not None:
          subj['MNI'] = di2

      # Find Electrode Implantation
      rdi = ReadDiskItem( 'Electrode implantation', 'Electrode Implantation format', requiredAttributes={'subject':s, 'center':subj['center']} )
      elecs = list(rdi._findValues( {}, None, False ) )
      if len(elecs) > 0:
        subj['implantation'] = elecs[0]

      else:
        wditxtmnipos = ReadDiskItem('Electrode Implantation Position TXT', 'Text file', requiredAttributes={'ref_name':'MNI'})
        dimnipos=wditxtmnipos.findValue(subj['rdi'])
        wditxtmniname = ReadDiskItem('Electrode Implantation Name TXT', 'Text file', requiredAttributes={'ref_name':'MNI'})
        dimniname = wditxtmniname.findValue(subj['rdi'])
        
        if dimnipos is not None:
          subj['implantationMNI'] = (dimnipos,dimniname)
        else:                   
          wdi_csvnew =  ReadDiskItem('Final Export Dictionaries','CSV file')
          rdi_csvnew= wdi_csvnew.findValue(subj['rdi'])
          subj['implantationCSV'] = rdi_csvnew
      
      # Find Electrode Label
      rdi_elec = ReadDiskItem('Electrodes Labels','Electrode Label Format',requiredAttributes={'subject':s, 'center':subj['center']})
      list_di_elec = list(rdi_elec.findValues( {}, None, False ) )
      #loca_patient = []
      loca_patient = {}
      loca_patient_FS = {}
      loca_patient_BM = {}
      loca_bipole_patient = {}
      loca_bipole_patient_FS = {}
      loca_bipole_patient_BM = {}
      
      info_label_elec = None
      
      if len(list_di_elec) > 0:
          di_elec = rdi_elec.findValue(subj['T1pre'])
          #subj['elec_label'] = di_elec
          #read the json file
          fin = open(di_elec.fullPath(),'r')
          info_label_elec = json.loads(fin.read())
          fin.close()
          if 'plots_by_label' in info_label_elec.keys():
            for kk,vv in info_label_elec['plots_by_label'].iteritems():
                if len(vv):
                    loca.append(kk)
                    #loca_patient.append(kk)
                    loca_patient.update({kk:vv})
                    
          if 'plots_by_label_FS' in info_label_elec.keys():
             for kk,vv in info_label_elec['plots_by_label_FS'].iteritems():
                if len(vv):
                    loca_FS.append(kk)
                    #loca_patient_FS.append(kk)
                    loca_patient_FS.update({kk:vv})
          if 'plots_by_label_BM' in info_label_elec.keys():
             for kk,vv in info_label_elec['plots_by_label_BM'].iteritems():
                if len(vv):
                    loca_BM.append(kk)
                    #loca_patient_BM.append(kk)
                    loca_patient_BM.update({kk:vv})
                    
          #on fait la meme chose avec les bipoles pour les stims (puis la cognition)
          if 'plots_bipolar_by_label' in info_label_elec.keys():
            for kk,vv in info_label_elec['plots_bipolar_by_label'].iteritems():
                if len(vv):
                    #loca_patient.append(kk)
                    loca_bipole_patient.update({kk:vv})              
          if 'plots_bipolar_by_label_FS' in info_label_elec.keys():
            for kk,vv in info_label_elec['plots_bipolar_by_label_FS'].iteritems():
                if len(vv):
                    #loca_patient.append(kk)
                    loca_bipole_patient_FS.update({kk:vv})     
          if 'plots_bipolar_by_label_BM' in info_label_elec.keys():
            for kk,vv in info_label_elec['plots_bipolar_by_label_BM'].iteritems():
                if len(vv):
                    #loca_patient.append(kk)
                    loca_bipole_patient_BM.update({kk:vv})               
          
          #add freesurfer and broadmann        
          subj['rdi_elec_label'] = list_di_elec[0]
      else:
        subj['rdi_elec_label'] = None      
      
      subj['elec_label'] = loca_patient
      subj['elec_label_FS'] = loca_patient_FS
      subj['elec_label_BM'] = loca_patient_BM
      subj['bipole_label'] = loca_bipole_patient
      subj['bipole_label_FS'] = loca_bipole_patient_FS
      subj['bipole_label_BM'] = loca_bipole_patient_BM


      rdi_resec_label = ReadDiskItem('Resection Description','Resection json',requiredAttributes={'subject':s, 'center':subj['center']})
      list_di_resec_label = list(rdi_resec_label.findValues( {}, None, False ) )
      resec_patient = []
      resec_patient_FS = []
      if len(list_di_resec_label) > 0:
        di_resec_label = rdi_resec_label.findValue(subj['T1pre'])
        fin = open(di_resec_label.fullPath(),'r')
        info_label_resec = json.loads(fin.read())
        fin.close()

        if 'mars_atlas' in info_label_resec.keys():
          for kk,vv in info_label_resec['mars_atlas'].iteritems():
              if len(vv):
                  resec.append(kk)
                  resec_patient.append(kk)
                
        #add freesurfer
        if 'Freesurfer' in info_label_resec.keys():
            if 'Freesurfer not calculated' in info_label_resec['Freesurfer']:
                pass
            else:
                for kk,vv in info_label_resec['Freesurfer'].iteritems():
                    if len(vv):
                        resec_FS.append(kk)
                        resec_patient_FS.append(kk)
        
        subj['rdi_resec_label'] = list_di_resec_label[0]
      else:
        subj['rdi_resec_label'] = None
        
      subj['resec_label'] = resec_patient
      subj['resec_label_FS'] = resec_patient_FS

      # Find SEEG manips
      rdi = ReadDiskItem( 'Raw SEEG recording', ['EEG TRC format', 'Elan EEG format'], requiredAttributes={'subject':s, 'center':subj['center']} )
      records = list(rdi._findValues( {}, None, False ) )
      if len(records) > 0:
        for rec in records:
          subj['seeg_'+str(rec.attributes()['experiment'])] = rec
        manips = [str(rec.attributes()['experiment']) for rec in records]
        seegManips.update(manips)

      rdi_stimresult_label = ReadDiskItem('Electrodes SEEG Labels','Electrode sEEG Label Format',requiredAttributes={'subject':s, 'center':subj['center']})
      list_di_stimresult_label = list(rdi_stimresult_label.findValues( {}, None, False ) )
      seeg_label_loca = []
      seeg_label_type = []
      seeg_label_freq = []
      seeg_label_all = []
      
      
      if len(list_di_stimresult_label)>0:
        
        by_default_list = [u'Motor',u'Sensitive',u'Sensory',u'Vegetative',u'Emotional',u'Experiencial','Superior functions']
        by_default_list2 = [u'no response', u'pathological', u'seizure sensation']
        full_list = by_default_list + by_default_list2
        
        di_stimresult_label = rdi_stimresult_label.findValue(subj['T1pre'])
        fin = open(di_stimresult_label.fullPath(),'r')
        info_seeg_label = json.loads(fin.read())
        fin.close()
        
        #la il faut que je retourne les infos dans tous les sens c'est la merde ...
        for kk,vv in info_seeg_label['contacts'].iteritems():
            #seeg_label_loca = 

            for freqVal in info_seeg_label['contacts'][kk]['cell'].keys():
                if freqVal not in stimFreq:
                    stimFreq.append(freqVal)
                    #seeg_label_freq.append(freqVal)
                list_inter = [x for x in by_default_list if info_seeg_label['contacts'][kk]['cell'][freqVal][x]['value'] != 0]
                list_inter2 = info_seeg_label['contacts'][kk]['cell'][freqVal][u'Type of response']['fontcolor']
                max_index = list_inter2[1:].index(max(list_inter2[1:]))
                try:
                    if list_inter2[max_index+1] == 0:
                      resp1 = by_default_list2[0]
                    elif max_index == 0:
                      resp1 = by_default_list2[2]
                    elif max_index == 2:
                      resp1 = by_default_list2[1] 
                except:
                    pdb.set_trace()
                seeg_label_all.append((kk.title(),freqVal,resp1,list_inter))
                try:
                  if 'MarsAtlas' in info_label_elec['plots_label_bipolar'][kk.title()].keys():
                    seeg_label_all[-1] = seeg_label_all[-1] + (info_label_elec['plots_label_bipolar'][kk.title()]['MarsAtlas'][1],)
                  else:
                    seeg_label_all[-1] = seeg_label_all[-1] + (u'not estimated',)
                except:
                    pdb.set_trace()

                if 'Freesurfer' in info_label_elec['plots_label_bipolar'][kk.title()].keys():
                  seeg_label_all[-1] = seeg_label_all[-1] + (info_label_elec['plots_label_bipolar'][kk.title()]['Freesurfer'][1],)
                else:
                  seeg_label_all[-1] = seeg_label_all[-1] + (u'not estimated',)
                  
                if 'Broadmann' in info_label_elec['plots_label_bipolar'][kk.title()].keys():
                  seeg_label_all[-1] = seeg_label_all[-1] + (info_label_elec['plots_label_bipolar'][kk.title()]['Broadmann'][1],)
                else:
                  seeg_label_all[-1] = seeg_label_all[-1] + (u'not estimated',)  
                    #[True if info_seeg_label['contacts'][kk]['cell'][freqVal][x]['value'] != 0 else False for x in by_default_list]
        #pdb.set_trace()
          
          
        subj['seeg_label_all'] = seeg_label_all
        
    if 'full_list' in locals():
        stimInter = set(stimType)
        full_listInter = set(full_list)
        newstimtype = stimInter.union(full_listInter)
        stimType = sorted(list(newstimtype))
    #stimType = stimType + full_list
    
    #remove duplicate #la on récupère les clés ? ou après ? où est-ce que je gère les infos des contacts ?
    locas = list(set(loca))
    locas.sort()
    locas_FS = list(set(loca_FS))
    locas_FS.sort()
    locas_BM = list(set(loca_BM))
    locas_BM.sort()
    
    
    
    self.filterCognitionLoca.clear()
    self.filterCognitionType.clear()

    
    self.filterCognitionLoca.addItems(locas)
    self.filterCognitionType.addItems(cogniType)
    
    self.filterStimulationLoca.clear()
    self.filterStimulationType.clear()
    self.filterStimulationFreq.clear()
    
    self.filterStimulationLoca.addItems(locas)
    self.filterStimulationType.addItems(stimType)
    self.filterStimulationFreq.addItems(stimFreq)
    
    resecs = list(set(resec))
    resecs.sort()
    resecs_FS = list(set(resec_FS))
    resecs_FS.sort()
    resecs_BM = list(set(resec_BM))
    resecs_BM.sort()

    self.filtersAtlases.update({'locaMA':locas,'locaFS':locas_FS,'locaBM':locas_BM,'resecMA':resecs,'resecFS':resecs_FS,'resecBM':resecs_BM}) 
    
    self.filterLocaMarsAtlas.clear()
    self.filterLocaMarsAtlas.addItems(locas)

    #self.filterCognition.clear()
    #self.filterCognition.addItems(cognis)

    self.filterResection.clear()
    self.filterResection.addItems(resecs)

    if len(seegManips) > 0:
      self.seegManipCombo.clear()
      self.seegManipCombo.addItems(sorted(list(seegManips)))

    self.updatePatientFilters()

  def toggleSeeg(self):
    """If SEEG is clicked, select the manip as a filter"""
    manip = str(self.seegManipCombo.currentText())
    if manip is None or manip == "":
      return
    self.toggleFilterWidget('seeg_'+manip, self.seegButton)
    if 'seeg_'+manip in self.filters:
      self.seegManipCombo.setItemIcon(self.seegManipCombo.currentIndex(), self.greenDotIcon)
    else:
      self.seegManipCombo.setItemIcon(self.seegManipCombo.currentIndex(), None)

  def seegManipChanged(self, idx):
    self.toggleWidget(self.seegButton, 'seeg_' + str(self.seegManipCombo.currentText()) in self.filters)


  def toggleFilterWidget(self, filterName, widg):
    """Adds/Removes a filter and push/unpush the widget"""
    onOrOff = self.toggleFilter(filterName)
    self.toggleWidget(widg, onOrOff)

  def toggleWidget(self, widg, on = True):
    """Sets the background color of the widget to green and add it to pushedWidgets, or reverse if it is already pushed"""
    if not on:
      widg.setStyleSheet("background-color: rgb(90, 90, 90);")
      self.pushedWidget.remove(widg)
    else:
      widg.setStyleSheet("background-color: rgb(120, 255, 120);")
      #widg.setStyleSheet("border: 2px solid #88ff88;")
      self.pushedWidget.append(widg)

  def toggleFilter(self, filterName, updateUi = True):
    """Sets or unsets a filter, then update the ui"""
    if filterName in self.filters:
      self.filters.remove(filterName)
    else:
      self.filters.append(filterName)
    if updateUi:
      self.updatePatientFilters()
    return filterName in self.filters

  def getFilters(self):
    return self.filters

  def matchFilters(self, subj, filters):
    return all([subj.has_key(f) for f in filters])

  def filterSubjectsBasic(self):
    """Filtering subject list using basic factors (protocol, year, site)"""
    subs = self.subjects.keys()
    if str(self.filterProtocolCombo.currentText()) != '*':
      subs = [s for s in subs if self.subjects[s]['center'] == str(self.filterProtocolCombo.currentText())]
    if str(self.filterSiteCombo.currentText()) != '*':
      subs = [s for s in subs if s.split('_')[0] == str(self.filterSiteCombo.currentText())]
    if str(self.filterYearCombo.currentText()) != '*':
      subs = [s for s in subs if len(s.split('_')) > 1 and s.split('_')[1] == str(self.filterYearCombo.currentText())]
    if self.FilterAtlascomboBox.currentIndex() == 0:
      if str(self.filterLocaMarsAtlas.currentText()) != '*':
        subs = [s for s in subs if 'elec_label' in self.subjects[s] and str(self.filterLocaMarsAtlas.currentText()) in self.subjects[s]['elec_label'].keys()]
      if str(self.filterResection.currentText()) != '*':
        subs = [s for s in subs if 'resec_label' in self.subjects[s] and  str(self.filterResection.currentText()) in self.subjects[s]['resec_label']]
    elif self.FilterAtlascomboBox.currentIndex() == 1:
      if str(self.filterLocaMarsAtlas.currentText()) != '*':
        subs = [s for s in subs if 'elec_label_FS' in self.subjects[s] and str(self.filterLocaMarsAtlas.currentText()) in self.subjects[s]['elec_label_FS'].keys()]
      if str(self.filterResection.currentText()) != '*':
        subs = [s for s in subs if 'resec_label_FS' in self.subjects[s] and  str(self.filterResection.currentText()) in self.subjects[s]['resec_label_FS']]       
    elif self.FilterAtlascomboBox.currentIndex() == 2:
      if str(self.filterLocaMarsAtlas.currentText()) != '*':
        subs = [s for s in subs if 'elec_label_BM' in self.subjects[s] and str(self.filterLocaMarsAtlas.currentText()) in self.subjects[s]['elec_label_BM'].keys()]
      if str(self.filterResection.currentText()) != '*':
        subs = [s for s in subs if 'resec_label_BM' in self.subjects[s] and  str(self.filterResection.currentText()) in self.subjects[s]['resec_label_BM']]
    if str(self.filterStimulationLoca.currentText())!= '*' or str(self.filterStimulationFreq.currentText()) != '*' or str(self.filterStimulationType.currentText()) != '*' :
       subs = [s for s in subs if 'seeg_label_all' in self.subjects[s]]
       CondiLoca = str(self.filterStimulationLoca.currentText()) == '*'
       CondiFreq = str(self.filterStimulationFreq.currentText()) == '*'
       CondiType = str(self.filterStimulationType.currentText()) == '*'
       subs_toremove = []
       inter_stim = {}
       for ind_sub in subs:
           isthereloc = [x for x in self.subjects[ind_sub]['seeg_label_all'] if (str(self.filterStimulationLoca.currentText())==x[4+self.FilterAtlascomboBox.currentIndex()] or CondiLoca) and (str(self.filterStimulationFreq.currentText())==x[1] or CondiFreq) and (str(self.filterStimulationType.currentText())==x[2] or str(self.filterStimulationType.currentText()) in x[3] or CondiType)] 
           if len(isthereloc)==0:
               subs_toremove.append(ind_sub)
           else:
               inter_stim.update({ind_sub:isthereloc})
           #for bipole_info in self.sub   
       [subs.remove(x) for x in subs_toremove]
    
        
        
        
    if str(self.filterCognitionLoca.currentText()) != '*':
       print "no cognition load in the data base for now" 
        
    subs = [s for s in subs if s not in [str(self.selectedPatientList.item(idx).text()) for idx in range(self.selectedPatientList.count())]]
    self.filteredPatientList.clear()
    self.filteredPatientList.addItems(sorted(subs))

  def colorFilterListWidget(self, listwidget, selected):
    for idx in range(listwidget.count()):
      it = listwidget.item(idx)
      if str(it.text()) in selected:
        it.setBackground(QtGui.QColor(150,255,150))
      else:
        it.setBackground(QtGui.QColor(255,150,150))

  def updatePatientFilters(self):
    """Colors the 'patient list' and 'selected patient list' items in green if they match the filters and updates the content of patient's list"""
    # Apply basic filters to the main patient list first and update the UI
    self.filterSubjectsBasic()
    # Filter the list of subjects and filter them
    filteredPatients = [s for s,subj in self.subjects.iteritems() if self.matchFilters(subj, self.filters)]
    # Color the patients according to the filters
    self.colorFilterListWidget(self.filteredPatientList, filteredPatients)
    self.colorFilterListWidget(self.selectedPatientList, filteredPatients)

  def moveSelectedItemsToOtherListWidget(self, lwFrom, lwTo):
    """Takes the selected items of the list 'from' and adds them to the 'to' list"""
    for idx in reversed(range(lwFrom.count())):
      it = lwFrom.item(idx)
      if lwFrom.isItemSelected(it):
        lwTo.addItem(str(it.text()))
        lwFrom.takeItem(idx)
    lwTo.sortItems()

  def moveFilteredItemsToOtherListWidget(self, lwFrom, lwTo, filters, reverse=False):
    """Takes the items from lwFrom (QListWidget) and move those that match the filters to lwTo. If reverse is true, the items that don't match will move"""
    filtered = [s for s,subj in self.subjects.iteritems() if self.matchFilters(subj, filters)]
    for idx in reversed(range(lwFrom.count())):
      it = lwFrom.item(idx)
      if (not reverse and str(it.text()) in filtered) or (reverse and str(it.text()) not in filtered):
        lwTo.addItem(str(it.text()))
        lwFrom.takeItem(idx)

    lwTo.sortItems()

  def selectAllListItems(self, listwidget, select = True):
    """Selects or deselects all items of a QListWidget (select = False to deselect)"""
    for idx in range(listwidget.count()):
      listwidget.setItemSelected(listwidget.item(idx), select)

  def mniParcelsGeneration(self):
      
    #patient selected
    PatsSelected = self.getSelectedPatientsNames()
    subject_to_perform = []
    subjects_ok = []
    
    for ind_pat in PatsSelected:
        
        #first we look if it's not already done
        if 'rdi_elec_label' in self.subjects[ind_pat].keys():
          if self.subjects[ind_pat]['rdi_elec_label'] is not None:
            fin = open(str(self.subjects[ind_pat]['rdi_elec_label']),'r')
            info_label_elec = json.loads(fin.read())
            fin.close()
            if 'plots_label' in info_label_elec.keys():
                if "BroadmannDilate" in info_label_elec['plots_label'][info_label_elec['plots_label'].keys()[0]].keys():               
                  subjects_ok.append(ind_pat)
                  print("BrodmannDilate already estimated for this patient %s"%ind_pat)
                else:
                  subject_to_perform.append(ind_pat)    
            else:
                subject_to_perform.append(ind_pat)
          else:
            subject_to_perform.append(ind_pat)
        else:
          subject_to_perform.append(ind_pat)
          
    #if not we look if MNI position have been estimated already
    subject_possible_to_perform = {}
    subject_not_performable = []
        
    for ind_pat in subject_to_perform:
        
        if 'implantation' in self.subjects[ind_pat].keys():

            fin = open(str(self.subjects[ind_pat]['implantation']),'r')
            try:    
                info_implant = json.loads(fin.read())
            except:
                fin.close()
                fin = open(str(self.subjects[ind_pat]['implantation']),'r')
                info_implant = pickle.load(fin)
            fin.close()
            
            if 'plotsMNI' in info_implant.keys():
                subject_possible_to_perform.update({ind_pat:info_implant['plotsMNI']})
            else:
                print("MNI estimation has to be done in locateElectrodes for %s, to sensitive to duplicate these functions at differente place"%ind_pat)
                subject_not_performable.append(ind_pat)
        
        elif 'implantationMNI' in self.subjects[ind_pat].keys():
            #pdb.set_trace()
            MNI_pos = numpy.loadtxt(str(self.subjects[ind_pat]['implantationMNI'][0]))
            a=open(str(self.subjects[ind_pat]['implantationMNI'][1]),"rU")
            MNI_name=a.read().splitlines()
            a.close()
            info = [[MNI_name[x],MNI_pos[x].tolist()] for x in range(len(MNI_name))]
            subject_possible_to_perform.update({ind_pat:info})
            #MNI_name = 
            #subject_possible_to_perform.update({ind_pat:info_implant['plotsMNI']})
        
        elif 'implantationCSV' in self.subjects[ind_pat].keys():
        #    pdb.set_trace()
            print("refaire %s"%str(ind_pat))
            #subject_possible_to_perform.update({ind_pat:info_implant['plotsMNI']})
        else:
            subject_not_performable.append(ind_pat)

    
    matrix_MNI_Nativ = numpy.matrix([[  -1.,    0.,    0.,   90.],[0.,   -1.,    0.,   91.],[0.,    0.,   -1.,  109.],[0.,    0.,    0.,    1.]])    
    #la blague maintenant c'est d'ouvrir un csv et de rajouter des infos mni sans perdre les infos "patients specifique"
    for ind_pat in subject_possible_to_perform.keys():
        
        print(ind_pat)
        plotMNI_sorted = sorted(subject_possible_to_perform[ind_pat], key=lambda plot_number: plot_number[0])
        plotMNI_Dict = dict(plotMNI_sorted)
        #make bipolar montage


        plot_dict_MNI_Native = {}
        for vv,kk in plotMNI_Dict.iteritems():
          inter_pos = [kk[0], kk[1], kk[2], 1]
          inter_pos = numpy.matrix(inter_pos).reshape([4,1])
          result_pos = numpy.dot(matrix_MNI_Nativ,inter_pos)
          plot_dict_MNI_Native.update({vv:[result_pos.tolist()[0][0],result_pos.tolist()[1][0],result_pos.tolist()[2][0]]})
        
        info_plotMNI_bipolaire= []
        for pindex in range(1,len(plotMNI_sorted)):
          previous_contact = "".join([i for i in plotMNI_sorted[pindex-1][0] if not i.isdigit()])
          current_contact = "".join([i for i in plotMNI_sorted[pindex][0] if not i.isdigit()])
          if previous_contact == current_contact:
               info_plotMNI_bipolaire.append((plotMNI_sorted[pindex][0]+' - '+ plotMNI_sorted[pindex-1][0],(numpy.array(plot_dict_MNI_Native[plotMNI_sorted[pindex][0]])+numpy.array(plot_dict_MNI_Native[plotMNI_sorted[pindex-1][0]]))/2 ))
               
        info_plotMNI_bipolaireSB= {}
        for pindex in range(1,len(plotMNI_sorted)):
          previous_contact = "".join([i for i in plotMNI_sorted[pindex-1][0] if not i.isdigit()])
          current_contact = "".join([i for i in plotMNI_sorted[pindex][0] if not i.isdigit()])
          if previous_contact == current_contact:
               info_plotMNI_bipolaireSB.update({plotMNI_sorted[pindex][0]+' - '+ plotMNI_sorted[pindex-1][0]:(numpy.array(plotMNI_Dict[plotMNI_sorted[pindex][0]])+numpy.array(plotMNI_Dict[plotMNI_sorted[pindex-1][0]]))/2 })       
               
        info_plotMNI_bipolaire_Dict =dict(info_plotMNI_bipolaire)
        #open MNI filtersAtlases
        #pdb.set_trace()
        
        Hammers_parcels_names = readSulcusLabelTranslationFile('parcels_label_name_Hammers.txt')
        AAL_parcels_names = readSulcusLabelTranslationFile('parcels_label_name_AAL.txt')
        AALDilate_parcels_names = readSulcusLabelTranslationFile('parcels_label_name_AALDilate.txt')
        vol_AAL = aims.read('MNI_Atlases/rAALSEEG12.nii')
        vol_AALDilate = aims.read('MNI_Atlases/rAALSEEG12Dilate.nii')
        vol_BroadmannDilate = aims.read('MNI_Atlases/rBrodmannSEEG3spm12.nii')
        vol_Broadmann = aims.read('MNI_Atlases/rbrodmann.nii')
        vol_Hammers = aims.read('MNI_Atlases/rHammersSEEG12.nii')
        
        sphere_size = 3 #en mm
        nb_voxel_sphere_MNI = [sphere_size, sphere_size, sphere_size] 
        #c'est là que ça va pas ...
        #il faut que je recharge celui du patient
        plots_label = {}
        if 'rdi_elec_label' in self.subjects[ind_pat].keys():
          if self.subjects[ind_pat]['rdi_elec_label'] is not None:
            fin = open(str(self.subjects[ind_pat]['rdi_elec_label']),'r')
            info_label_elec = json.loads(fin.read())
            fin.close()
            plots_label = info_label_elec['plots_label']
          else:
            plots_label = {}  
        else:
          plots_label = {}    
            
        for pindex in range(len(plotMNI_sorted)):

         plot_pos_pix_MNI = [round(plot_dict_MNI_Native[plotMNI_sorted[pindex][0]][i]) for i in range(3)]
         
         ##MNI Atlases
         ##AAL
         voxel_within_sphere_AAL = [round(vol_AAL.value(plot_pos_pix_MNI[0]+vox_i,plot_pos_pix_MNI[1]+vox_j,plot_pos_pix_MNI[2]+vox_k)) for vox_k in range(-nb_voxel_sphere_MNI[2],nb_voxel_sphere_MNI[2]+1) for vox_j in range(-nb_voxel_sphere_MNI[1],nb_voxel_sphere_MNI[1]+1) for vox_i in range(-nb_voxel_sphere_MNI[0],nb_voxel_sphere_MNI[0]+1) if math.sqrt(vox_i**2+vox_j**2+vox_k**2) < sphere_size]
         voxel_to_keepAAL = [x for x in voxel_within_sphere_AAL if x != 0 and not math.isnan(x)]
         
         ##AALDilate
         voxel_within_sphere_AALdilate = [round(vol_AALDilate.value(plot_pos_pix_MNI[0]+vox_i,plot_pos_pix_MNI[1]+vox_j,plot_pos_pix_MNI[2]+vox_k)) for vox_k in range(-nb_voxel_sphere_MNI[2],nb_voxel_sphere_MNI[2]+1) for vox_j in range(-nb_voxel_sphere_MNI[1],nb_voxel_sphere_MNI[1]+1) for vox_i in range(-nb_voxel_sphere_MNI[0],nb_voxel_sphere_MNI[0]+1) if math.sqrt(vox_i**2+vox_j**2+vox_k**2) < sphere_size]
         voxel_to_keepAALDilate = [x for x in voxel_within_sphere_AALdilate if x != 0 and not math.isnan(x)]
         
         ##Broadmann
         voxel_within_sphere_Broadmann = [round(vol_Broadmann.value(plot_pos_pix_MNI[0]+vox_i,plot_pos_pix_MNI[1]+vox_j,plot_pos_pix_MNI[2]+vox_k)) for vox_k in range(-nb_voxel_sphere_MNI[2],nb_voxel_sphere_MNI[2]+1) for vox_j in range(-nb_voxel_sphere_MNI[1],nb_voxel_sphere_MNI[1]+1) for vox_i in range(-nb_voxel_sphere_MNI[0],nb_voxel_sphere_MNI[0]+1) if math.sqrt(vox_i**2+vox_j**2+vox_k**2) < sphere_size]
         voxel_to_keepBroadmann = [x for x in voxel_within_sphere_Broadmann if x != 0 and not math.isnan(x)]
         
         #Brodmann dilate
         voxel_within_sphere_Broadmanndilate = [round(vol_BroadmannDilate.value(plot_pos_pix_MNI[0]+vox_i,plot_pos_pix_MNI[1]+vox_j,plot_pos_pix_MNI[2]+vox_k)) for vox_k in range(-nb_voxel_sphere_MNI[2],nb_voxel_sphere_MNI[2]+1) for vox_j in range(-nb_voxel_sphere_MNI[1],nb_voxel_sphere_MNI[1]+1) for vox_i in range(-nb_voxel_sphere_MNI[0],nb_voxel_sphere_MNI[0]+1) if math.sqrt(vox_i**2+vox_j**2+vox_k**2) < sphere_size]
         voxel_to_keepBroadmannDilate = [x for x in voxel_within_sphere_Broadmanndilate if x != 0 and not math.isnan(x)]
         
         ##Hammers
         voxel_within_sphere_Hammers = [round(vol_Hammers.value(plot_pos_pix_MNI[0]+vox_i,plot_pos_pix_MNI[1]+vox_j,plot_pos_pix_MNI[2]+vox_k)) for vox_k in range(-nb_voxel_sphere_MNI[2],nb_voxel_sphere_MNI[2]+1) for vox_j in range(-nb_voxel_sphere_MNI[1],nb_voxel_sphere_MNI[1]+1) for vox_i in range(-nb_voxel_sphere_MNI[0],nb_voxel_sphere_MNI[0]+1) if math.sqrt(vox_i**2+vox_j**2+vox_k**2) < sphere_size]
         voxel_to_keepHammers = [x for x in voxel_within_sphere_Hammers if x != 0 and not math.isnan(x)]

         ##prendre le label qui revient le plus en dehors de zero (au cas où il y en ait plusieurs)
         from collections import Counter

         if not voxel_to_keepAAL:
            label_AAL_name = "not in a AAL parcel" 
            label_AAL = round(vol_AAL.value(plot_pos_pix_MNI[0],plot_pos_pix_MNI[1],plot_pos_pix_MNI[2]))
         else:
            most_common,num_most_common = Counter(voxel_to_keepAAL).most_common(1)[0]
            label_AAL = most_common
            label_AAL_name = AAL_parcels_names[label_AAL]
             
         if not voxel_to_keepAALDilate:
            label_AALDilate_name = "not in a AALDilate parcel" 
            label_AALDilate = round(vol_AALDilate.value(plot_pos_pix_MNI[0],plot_pos_pix_MNI[1],plot_pos_pix_MNI[2]))
         else:
            most_common,num_most_common = Counter(voxel_to_keepAALDilate).most_common(1)[0]
            label_AALDilate = most_common
            label_AALDilate_name = AALDilate_parcels_names[label_AALDilate]
         
         if not voxel_to_keepBroadmann:
            label_Broadmann_name = "not in a Broadmann parcel" 
            label_Broadmann = round(vol_Broadmann.value(plot_pos_pix_MNI[0],plot_pos_pix_MNI[1],plot_pos_pix_MNI[2]))
         else:
            most_common,num_most_common = Counter(voxel_to_keepBroadmann).most_common(1)[0]
            label_Broadmann = most_common
            #label_Broadmann_name = unicode(label_Broadmann)
 
            if plot_pos_pix_MNI[0]>90:
                label_Broadmann_name = unicode(label_Broadmann+100)
            else:
                label_Broadmann_name = unicode(label_Broadmann)
                
         if not voxel_to_keepBroadmannDilate:
            label_BroadmannDilate_name = "not in a Broadmann parcel" 
            label_BroadmannDilate = round(vol_BroadmannDilate.value(plot_pos_pix_MNI[0],plot_pos_pix_MNI[1],plot_pos_pix_MNI[2]))
         else:
            most_common,num_most_common = Counter(voxel_to_keepBroadmannDilate).most_common(1)[0]
            label_BroadmannDilate = most_common
            #label_Broadmann_name = unicode(label_Broadmann)
            if plot_pos_pix_MNI[0]>90:
                label_BroadmannDilate_name = unicode(label_BroadmannDilate+100)
            else:
                label_BroadmannDilate_name = unicode(label_BroadmannDilate-48)          
             
         if not voxel_to_keepHammers:
            label_Hammers_name = "not in a Hammers parcel" 
            label_Hammers = round(vol_Hammers.value(plot_pos_pix_MNI[0],plot_pos_pix_MNI[1],plot_pos_pix_MNI[2]))
         else:    
            most_common,num_most_common = Counter(voxel_to_keepHammers).most_common(1)[0]
            label_Hammers = most_common
            label_Hammers_name = Hammers_parcels_names[label_Hammers]
                        

         #la j'ecrase les infos marsatlas etc ... sans raison
         #plots_label[plotMNI_sorted[pindex][0]]={'MarsAtlas':(0,'not calculated'),'Freesurfer':(0,'not calculated'),'Hippocampal Subfield':(0,'not calculated'),'GreyWhite':(0,'not calculated'),'AAL':(label_AAL,label_AAL_name),'AALDilate':(label_AALDilate,label_AALDilate_name),'Broadmann':(label_Broadmann,label_Broadmann_name),'Hammers':(label_Hammers,label_Hammers_name),'Resection':(0,'not calculated')}
         try:
             plots_label[plotMNI_sorted[pindex][0]].update({'AAL':(label_AAL,label_AAL_name),'AALDilate':(label_AALDilate,label_AALDilate_name),'Broadmann':(label_Broadmann,label_Broadmann_name),'BroadmannDilate':(label_BroadmannDilate,label_BroadmannDilate_name),'Hammers':(label_Hammers,label_Hammers_name)})
         except:
             plots_label[plotMNI_sorted[pindex][0]]={'MarsAtlas':(0,'not calculated'),'Freesurfer':(0,'not calculated'),'Hippocampal Subfield':(0,'not calculated'),'GreyWhite':(0,'not calculated'),'AAL':(label_AAL,label_AAL_name),'AALDilate':(label_AALDilate,label_AALDilate_name),'Broadmann':(label_Broadmann,label_Broadmann_name),'BroadmannDilate':(label_BroadmannDilate,label_BroadmannDilate_name),'Hammers':(label_Hammers,label_Hammers_name),'Resection':(0,'not calculated')}
             #je remets l'ancienne ligne ?

        plot_name = [x[0] for x in plotMNI_sorted]
        plots_by_label_BM = dict([(Lab,[p for p in plot_name if plots_label[p]['Broadmann'][1]==Lab]) for Lab in [unicode("%1.1f"%x) for x in range(0,100)]])
        plots_by_label_HM = dict([(Lab,[p for p in plot_name if plots_label[p]['Hammers'][1]==Lab]) for Lab in Hammers_parcels_names.values()])
        plots_by_label_AAL = dict([(Lab,[p for p in plot_name if plots_label[p]['AAL'][1]==Lab]) for Lab in AAL_parcels_names.values()])
        plots_by_label_AALDilate = dict([(Lab,[p for p in plot_name if plots_label[p]['AALDilate'][1]==Lab]) for Lab in AALDilate_parcels_names.values()])

        sphere_size_bipole = 5
        nb_voxel_sphere_MNI = [sphere_size, sphere_size, sphere_size]
        
        plots_label_bipolar = {}
        if 'rdi_elec_label' in self.subjects[ind_pat].keys():
          if self.subjects[ind_pat]['rdi_elec_label'] is not None:
            fin = open(str(self.subjects[ind_pat]['rdi_elec_label']),'r')
            info_label_elec = json.loads(fin.read())
            fin.close()
            plots_label_bipolar = info_label_elec['plots_label_bipolar']
          else:
            pass #pdb.set_trace()  
        else:
          pass #pdb.set_trace() 
          
        for pindex in range(len(info_plotMNI_bipolaire)): 
         plot_pos_pix_MNI = [round(info_plotMNI_bipolaire[pindex][1][i]) for i in range(3)]
         
         ##MNI Atlases
         ##AAL
         voxel_within_sphere_AAL = [round(vol_AAL.value(plot_pos_pix_MNI[0]+vox_i,plot_pos_pix_MNI[1]+vox_j,plot_pos_pix_MNI[2]+vox_k)) for vox_k in range(-nb_voxel_sphere_MNI[2],nb_voxel_sphere_MNI[2]+1) for vox_j in range(-nb_voxel_sphere_MNI[1],nb_voxel_sphere_MNI[1]+1) for vox_i in range(-nb_voxel_sphere_MNI[0],nb_voxel_sphere_MNI[0]+1) if math.sqrt(vox_i**2+vox_j**2+vox_k**2) < sphere_size]
         voxel_to_keepAAL = [x for x in voxel_within_sphere_AAL if x != 0 and not math.isnan(x)]
         
         ##AALDilate
         voxel_within_sphere_AALdilate = [round(vol_AALDilate.value(plot_pos_pix_MNI[0]+vox_i,plot_pos_pix_MNI[1]+vox_j,plot_pos_pix_MNI[2]+vox_k)) for vox_k in range(-nb_voxel_sphere_MNI[2],nb_voxel_sphere_MNI[2]+1) for vox_j in range(-nb_voxel_sphere_MNI[1],nb_voxel_sphere_MNI[1]+1) for vox_i in range(-nb_voxel_sphere_MNI[0],nb_voxel_sphere_MNI[0]+1) if math.sqrt(vox_i**2+vox_j**2+vox_k**2) < sphere_size]
         voxel_to_keepAALDilate = [x for x in voxel_within_sphere_AALdilate if x != 0 and not math.isnan(x)]
         
         ##Broadmann
         voxel_within_sphere_Broadmann = [round(vol_Broadmann.value(plot_pos_pix_MNI[0]+vox_i,plot_pos_pix_MNI[1]+vox_j,plot_pos_pix_MNI[2]+vox_k)) for vox_k in range(-nb_voxel_sphere_MNI[2],nb_voxel_sphere_MNI[2]+1) for vox_j in range(-nb_voxel_sphere_MNI[1],nb_voxel_sphere_MNI[1]+1) for vox_i in range(-nb_voxel_sphere_MNI[0],nb_voxel_sphere_MNI[0]+1) if math.sqrt(vox_i**2+vox_j**2+vox_k**2) < sphere_size]
         voxel_to_keepBroadmann = [x for x in voxel_within_sphere_Broadmann if x != 0 and not math.isnan(x)]
         
         #Brodmann dilate
         voxel_within_sphere_Broadmanndilate = [round(vol_BroadmannDilate.value(plot_pos_pix_MNI[0]+vox_i,plot_pos_pix_MNI[1]+vox_j,plot_pos_pix_MNI[2]+vox_k)) for vox_k in range(-nb_voxel_sphere_MNI[2],nb_voxel_sphere_MNI[2]+1) for vox_j in range(-nb_voxel_sphere_MNI[1],nb_voxel_sphere_MNI[1]+1) for vox_i in range(-nb_voxel_sphere_MNI[0],nb_voxel_sphere_MNI[0]+1) if math.sqrt(vox_i**2+vox_j**2+vox_k**2) < sphere_size]
         voxel_to_keepBroadmannDilate = [x for x in voxel_within_sphere_Broadmanndilate if x != 0 and not math.isnan(x)]
         
         ##Hammers
         voxel_within_sphere_Hammers = [round(vol_Hammers.value(plot_pos_pix_MNI[0]+vox_i,plot_pos_pix_MNI[1]+vox_j,plot_pos_pix_MNI[2]+vox_k)) for vox_k in range(-nb_voxel_sphere_MNI[2],nb_voxel_sphere_MNI[2]+1) for vox_j in range(-nb_voxel_sphere_MNI[1],nb_voxel_sphere_MNI[1]+1) for vox_i in range(-nb_voxel_sphere_MNI[0],nb_voxel_sphere_MNI[0]+1) if math.sqrt(vox_i**2+vox_j**2+vox_k**2) < sphere_size]
         voxel_to_keepHammers = [x for x in voxel_within_sphere_Hammers if x != 0 and not math.isnan(x)]

         ##prendre le label qui revient le plus en dehors de zero (au cas où il y en ait plusieurs)
         from collections import Counter

         if not voxel_to_keepAAL:
            label_AAL_name = "not in a AAL parcel" 
            label_AAL = round(vol_AAL.value(plot_pos_pix_MNI[0],plot_pos_pix_MNI[1],plot_pos_pix_MNI[2]))
         else:
            most_common,num_most_common = Counter(voxel_to_keepAAL).most_common(1)[0]
            label_AAL = most_common
            label_AAL_name = AAL_parcels_names[label_AAL]
             
         if not voxel_to_keepAALDilate:
            label_AALDilate_name = "not in a AALDilate parcel" 
            label_AALDilate = round(vol_AALDilate.value(plot_pos_pix_MNI[0],plot_pos_pix_MNI[1],plot_pos_pix_MNI[2]))
         else:
            most_common,num_most_common = Counter(voxel_to_keepAALDilate).most_common(1)[0]
            label_AALDilate = most_common
            label_AALDilate_name = AALDilate_parcels_names[label_AALDilate]
         
         if not voxel_to_keepBroadmann:
            label_Broadmann_name = "not in a Broadmann parcel" 
            label_Broadmann = round(vol_Broadmann.value(plot_pos_pix_MNI[0],plot_pos_pix_MNI[1],plot_pos_pix_MNI[2]))
         else:
            most_common,num_most_common = Counter(voxel_to_keepBroadmann).most_common(1)[0]
            label_Broadmann = most_common
            #label_Broadmann_name = unicode(label_Broadmann)
            if plot_pos_pix_MNI[0]>90:
                label_Broadmann_name = unicode(label_Broadmann+100)
            else:
                label_Broadmann_name = unicode(label_Broadmann)
             
         if not voxel_to_keepBroadmannDilate:
            label_BroadmannDilate_name = "not in a Broadmann parcel" 
            label_BroadmannDilate = round(vol_BroadmannDilate.value(plot_pos_pix_MNI[0],plot_pos_pix_MNI[1],plot_pos_pix_MNI[2]))
         else:
            most_common,num_most_common = Counter(voxel_to_keepBroadmannDilate).most_common(1)[0]
            label_BroadmannDilate = most_common
            #label_Broadmann_name = unicode(label_Broadmann)
            if plot_pos_pix_MNI[0]>90:
                label_BroadmannDilate_name = unicode(label_BroadmannDilate+100)
            else:
                label_BroadmannDilate_name = unicode(label_BroadmannDilate-48) 
         
         
         if not voxel_to_keepHammers:
            label_Hammers_name = "not in a Hammers parcel" 
            label_Hammers = round(vol_Hammers.value(plot_pos_pix_MNI[0],plot_pos_pix_MNI[1],plot_pos_pix_MNI[2]))
         else:    
            most_common,num_most_common = Counter(voxel_to_keepHammers).most_common(1)[0]
            label_Hammers = most_common
            label_Hammers_name = Hammers_parcels_names[label_Hammers]
                        

         #plots_label_bipolar[info_plotMNI_bipolaire[pindex][0]]={'MarsAtlas':(0,'not calculated'),'Freesurfer':(0,'not calculated'),'Hippocampal Subfield':(0,'not calculated'),'GreyWhite':(0,'not calculated'),'AAL':(label_AAL,label_AAL_name),'AALDilate':(label_AALDilate,label_AALDilate_name),'Broadmann':(label_Broadmann,label_Broadmann_name),'Hammers':(label_Hammers,label_Hammers_name),'Resection':(0,'not calculated')}
         try:
           plots_label_bipolar[info_plotMNI_bipolaire[pindex][0]].update({'AAL':(label_AAL,label_AAL_name),'AALDilate':(label_AALDilate,label_AALDilate_name),'Broadmann':(label_Broadmann,label_Broadmann_name),'BroadmannDilate':(label_BroadmannDilate,label_BroadmannDilate_name),'Hammers':(label_Hammers,label_Hammers_name)})
         except:
           plots_label_bipolar[info_plotMNI_bipolaire[pindex][0]]={'MarsAtlas':(0,'not calculated'),'Freesurfer':(0,'not calculated'),'Hippocampal Subfield':(0,'not calculated'),'GreyWhite':(0,'not calculated'),'AAL':(label_AAL,label_AAL_name),'AALDilate':(label_AALDilate,label_AALDilate_name),'Broadmann':(label_Broadmann,label_Broadmann_name),'BroadmannDilate':(label_BroadmannDilate,label_BroadmannDilate_name),'Hammers':(label_Hammers,label_Hammers_name),'Resection':(0,'not calculated')}
           
        plot_name_bip = [x[0] for x in info_plotMNI_bipolaire]
        plots_bipolar_by_label_BM = dict([(Lab,[p for p in plot_name_bip if plots_label_bipolar[p]['Broadmann'][1]==Lab]) for Lab in [unicode("%1.1f"%x) for x in range(0,100)]])
        plots_bipolar_by_label_HM = dict([(Lab,[p for p in plot_name_bip if plots_label_bipolar[p]['Hammers'][1]==Lab]) for Lab in Hammers_parcels_names.values()])
        plots_bipolar_by_label_AAL = dict([(Lab,[p for p in plot_name_bip if plots_label_bipolar[p]['AAL'][1]==Lab]) for Lab in AAL_parcels_names.values()])
        plots_bipolar_by_label_AALDilate = dict([(Lab,[p for p in plot_name_bip if plots_label_bipolar[p]['AALDilate'][1]==Lab]) for Lab in AALDilate_parcels_names.values()])
        
        wdi_csv =  ReadDiskItem('Final Export Dictionaries','CSV file',requiredAttributes={'subject':ind_pat})
        rdi_csv = list(wdi_csv.findValues({},None,False))

        if len(rdi_csv)>0:
            #print c'est la galère
            info_previous_csv = readElecLocalCSVFile(infile=str(rdi_csv[0]))
            
        else:
            
            info_previous_csv = None
            
        wdi_csvnew =  WriteDiskItem('Final Export Dictionaries','CSV file')
        rdi_csvnew= wdi_csvnew.findValue(self.subjects[ind_pat]['rdi'])
        
        with open(str(rdi_csvnew), 'w') as csvfile:
          #fieldnames=['MarsAtlas','GreyWhite','Resection']
          writer = csv.writer(csvfile, delimiter='\t')
          writer.writerow([u'Contacts Positions'])
          if info_previous_csv is not None:
              info_mniMarsAtlas = info_previous_csv[0]['MarsAtlas']
              info_mniFreesurfer = info_previous_csv[0]['Freesurfer']
              info_mniHippoFS =  info_previous_csv[0]['HippoSubfieldFreesurfer']
          else:
              info_mniMarsAtlas = 'not performed'
              info_mniFreesurfer = 'not performed'
              info_mniHippoFS = 'not performed'
              
          writer.writerow([u'Use of MNI Template','MarsAtlas',info_mniMarsAtlas,'Freesurfer',info_mniFreesurfer,'HippoSubfieldFreesurfer',info_mniHippoFS])

          #add a row with "MNI or Patient for MarsAtlas and Freesurfer
          try:
            list_to_write = set(info_previous_csv[1]['monopolar'][info_previous_csv[1]['monopolar'].keys()[0]].keys())
          except:
            list_to_write = set([])
            
          list_by_default = set([u'contact','MarsAtlas', 'MarsAtlasFull', 'Freesurfer', 'Hippocampal Subfield','GreyWhite','AAL', 'AALDilate', 'Broadmann', 'BroadmannDilate', 'Hammers', 'Resection', 'MNI','T1pre Scanner Based'])
          diff_list = list(list_to_write.difference(list_by_default))
          full_list = [u'contact','MarsAtlas', 'MarsAtlasFull', 'Freesurfer', 'Hippocampal Subfield','GreyWhite', 'AAL', 'AALDilate', 'Broadmann', 'BroadmannDilate', 'Hammers', 'Resection', 'MNI','T1pre Scanner Based']
          full_list.extend(diff_list)
          writer.writerow(full_list)

          try:
              dict_sorted_tmp = OrderedDict(sorted(info_previous_csv[1]['monopolar'].items()))
          except:
              dict_sorted_tmp = OrderedDict(sorted(plots_label.items()))
              new_dict_sorted_tmp = {}
              for kk,vv in dict_sorted_tmp.items():
                  new_dict_sorted_tmp.update({kk:{}})
                  [new_dict_sorted_tmp[kk].update({jj:ll[1]}) for jj,ll in vv.items()]
              dict_sorted_tmp = OrderedDict(sorted(new_dict_sorted_tmp.items()))
          
                
          try:
            for kk,vv in dict_sorted_tmp.iteritems():
              listwrite = [kk]
              if 'MarsAtlas' in vv.keys():
                listwrite.append(vv['MarsAtlas'])
              else:
                listwrite.append('not performed')
              if 'MarsAtlasFull' in vv.keys():
                listwrite.append(vv['MarsAtlasFull'])
              else:
                listwrite.append('not performed')
              if 'Freesurfer' in vv.keys():  
                listwrite.append(vv['Freesurfer'])
              else:
                listwrite.append('not performed')
              if 'Hippocampal Subfield' in vv.keys():             
                listwrite.append(vv['Hippocampal Subfield'])
              else:
                listwrite.append('not performed')
              if 'GreyWhite' in vv.keys():
                listwrite.append(vv['GreyWhite'])
              else:
                listwrite.append('not performed')
              if 'AAL' in vv.keys():  
                listwrite.append(vv['AAL'])
              else:
                listwrite.append(plots_label[kk]['AAL'][1])
              if 'AALDilate' in vv.keys():
                listwrite.append(vv['AALDilate'])
              else:
                listwrite.append(plots_label[kk]['AALDilate'][1])
              #if 'Broadmann' in vv.keys():
              #  listwrite.append(vv['Broadmann'])
              #else:
              listwrite.append(plots_label[kk]['Broadmann'][1])
              #if 'BroadmannDilate' in vv.keys():
              #listwrite.append(vv['BroadmannDilate'])
              #else:
              listwrite.append(plots_label[kk]['BroadmannDilate'][1])    
              if 'Hammers' in vv.keys():
                listwrite.append(vv['Hammers'])
              else:
                listwrite.append(plots_label[kk]['Hammers'][1])
              if 'Resection' in vv.keys():
                listwrite.append(vv['Resection'])
              else:
                listwrite.append('not performed')
              #[listwrite.append(x[1]) for x in vv.values()]
              listwrite.append([float(format(plotMNI_Dict[kk][i],'.3f')) for i in range(3)])
              listwrite.append('not performed, csv generated from groupDisplay')
              if len(full_list)>12:
                for i_supp in range(len(full_list)-14):
                    listwrite.append(vv[full_list[14+i_supp]])
              writer.writerow(listwrite)
          except:
            pdb.set_trace()            
        
          writer.writerow([])
          writer.writerow([])

          try:
              dict_sorted_tmp = OrderedDict(sorted(info_previous_csv[1]['bipolar'].items()))
          except:
              dict_sorted_tmp = OrderedDict(sorted(plots_label_bipolar.items()))
              new_dict_sorted_tmp = {}
              for kk,vv in dict_sorted_tmp.items():
                  new_dict_sorted_tmp.update({kk:{}})
                  [new_dict_sorted_tmp[kk].update({jj:ll[1]}) for jj,ll in vv.items()]
              dict_sorted_tmp = OrderedDict(sorted(new_dict_sorted_tmp.items()))
          
          #dict_sorted_tmp = OrderedDict(sorted(info_previous_csv[1]['bipolar'].items()))

          for kk,vv in dict_sorted_tmp.iteritems():
            #pdb.set_trace()
            listwrite = [kk]
            if 'MarsAtlas' in vv.keys():
              listwrite.append(vv['MarsAtlas'])
            else:
              listwrite.append('not performed')
            if 'MarsAtlasFull' in vv.keys():
              listwrite.append(vv['MarsAtlasFull'])
            else:
              listwrite.append('not performed')  
            if 'Freesurfer' in vv.keys():  
              listwrite.append(vv['Freesurfer'])
            else:
              listwrite.append('not performed')
            if 'Hippocampal Subfield' in vv.keys():             
              listwrite.append(vv['Hippocampal Subfield'])
            else:
              listwrite.append('not performed')
            if 'GreyWhite' in vv.keys():
              listwrite.append(vv['GreyWhite'])
            else:
              listwrite.append('not performed')
            if 'AAL' in vv.keys():  
              listwrite.append(vv['AAL'])
            else:
              listwrite.append(plots_label_bipolar[kk]['AAL'][1])
            if 'AALDilate' in vv.keys():
              listwrite.append(vv['AALDilate'])
            else:
              listwrite.append(plots_label_bipolar[kk]['AALDilate'][1])
            #if 'Broadmann' in vv.keys():
            #  listwrite.append(vv['Broadmann'])
            #else:
            listwrite.append(plots_label_bipolar[kk]['Broadmann'][1])
            #if 'BroadmannDilate' in vv.keys():
            #  listwrite.append(vv['BroadmannDilate'])
            #else:
            listwrite.append(plots_label_bipolar[kk]['BroadmannDilate'][1])  
                
            if 'Hammers' in vv.keys():
              listwrite.append(vv['Hammers'])
            else:
              listwrite.append(plots_label_bipolar[kk]['Hammers'][1])
            if 'Resection' in vv.keys():
              listwrite.append(vv['Resection'])
            else:
              listwrite.append('not performed')
            
            #[listwrite.append(x[1]) for x in vv.values()]
            listwrite.append([float(format(info_plotMNI_bipolaireSB[kk][i],'.3f')) for i in range(3)])
            listwrite.append('not performed, csv generated from groupDisplay')
            if len(full_list)>12:
                for i_supp in range(len(full_list)-14):
                    listwrite.append(vv[full_list[14+i_supp]])
            writer.writerow(listwrite)

          writer.writerow([])
          writer.writerow([])

        if info_previous_csv is not None:
          if len(info_previous_csv[2])>0:
            with open(str(rdi_csvnew), 'a') as csvfile:
              writer = csv.writer(csvfile, delimiter='\t')
              writer.writerow([u'Resection Information'])
          
              for kk,vv in info_previous_csv[2].iteritems():
                #writer.writerow([kk])
                if type(vv) == type(float()):
                  listwrite = [kk,vv]
                  writer.writerow(listwrite)
                else:
                 writer.writerow([kk])
                 for ll,bb in vv.iteritems():
                   listwrite = [ll, format(float(bb),'.1f')]
                   writer.writerow(listwrite)

          else:
            pass
        else:
           pass 
        #neuroHierarchy.databases.insertDiskItem(di, update=True )
        print("generation of csv done for %s"%ind_pat)
    print("everything done")    

  def changeAtlasfiltration(self):
      
    
    #self.filtersAtlases.update({'locaMA':locas,'locaFS':locas_FS,'locaBM':locas_BM,'cogni':cognis,'resecMA':resecs,'resecFS':resecs_FS,'resecBM':resecs_BM})
    if self.FilterAtlascomboBox.currentIndex() == 0:
      self.filterLocaMarsAtlas.clear()
      self.filterLocaMarsAtlas.addItems(self.filtersAtlases['locaMA'])

      #self.filterCognition.clear()
      #self.filterCognition.addItems(self.filtersAtlases['cogni'])

      self.filterResection.clear()
      self.filterResection.addItems(self.filtersAtlases['resecMA'])
      
      self.filterCognitionLoca.clear()
      self.filterCognitionLoca.addItems(self.filtersAtlases['locaMA'])
    
      self.filterStimulationLoca.clear()
      self.filterStimulationLoca.addItems(self.filtersAtlases['locaMA'])

      
    elif self.FilterAtlascomboBox.currentIndex() == 1:
      
      self.filterLocaMarsAtlas.clear()
      self.filterLocaMarsAtlas.addItems(self.filtersAtlases['locaFS'])

      self.filterResection.clear()
      self.filterResection.addItems(self.filtersAtlases['resecFS'])       
      #change font color depending on patient who has the information about the filtersAtlases        

      self.filterCognitionLoca.clear()
      self.filterCognitionLoca.addItems(self.filtersAtlases['locaFS'])
    
      self.filterStimulationLoca.clear()
      self.filterStimulationLoca.addItems(self.filtersAtlases['locaFS'])
      
    elif self.FilterAtlascomboBox.currentIndex() == 2:  
        
      self.filterLocaMarsAtlas.clear()
      self.filterLocaMarsAtlas.addItems(self.filtersAtlases['locaBM'])
 
      self.filterResection.clear()
      self.filterResection.addItems(self.filtersAtlases['resecBM'])  
      
      self.filterCognitionLoca.clear()
      self.filterCognitionLoca.addItems(self.filtersAtlases['locaBM'])
    
      self.filterStimulationLoca.clear()
      self.filterStimulationLoca.addItems(self.filtersAtlases['locaBM'])      
      #change font color depending on patient who has the information about the filtersAtlases
      
      
if __name__ == "__main__":
  app = QtGui.QApplication(sys.argv)
  axon.initializeProcesses()
  from brainvisa.data.readdiskitem import ReadDiskItem
  from brainvisa.data.writediskitem import WriteDiskItem
  QtCore.pyqtRemoveInputHook()
  window = PatientDataFilterWidget()
  window.show()
  sys.exit(app.exec_())
