# -*-coding:utf-8 -*

import sys, os, pickle, glob, numpy, re, string, time, subprocess, json, copy, csv

from soma.qt_gui.qt_backend import QtGui, QtCore, uic, Qt

from numpy import *
from math import sqrt
from collections import OrderedDict

from soma import aims
from brainvisa import axon

#from soma.aims.spmnormalizationreader import readSpmNormalization
from brainvisa import anatomist
from brainvisa.data import neuroHierarchy
import brainvisa.registration as registration

from externalprocesses import *
from referentialconverter import ReferentialConverter



import pdb


#Main Class
class bipoleSEEGColors(QtGui.QDialog):
    
    def __init__(self,locateData = None,indv_pat=True,group_subsample=None):
        
        QtGui.QWidget.__init__(self)
        self.ui = uic.loadUi("bipoleLabelsDisplay.ui", self)
        self.setWindowTitle('Stimulation Observation - NOT FOR MEDICAL USE')    
        self.locaData =locateData #Data from locatesElectrode.py
        self.a=self.locaData.a #Anatomist Object   
        
        self.ui.lowHz_radiobutton.toggled.connect(self.updateBipoleDisplay)
        self.ui.highHz_radiobutton.toggled.connect(self.updateBipoleDisplay)
        self.ui.bothHz_radiobutton.toggled.connect(self.updateBipoleDisplay)
        
        self.ui.fontColor_radiobutton.toggled.connect(self.updateBipoleDisplay)
        self.ui.backColor_radiobutton.toggled.connect(self.updateBipoleDisplay)

        by_default_list = [u'Motor',u'Sensitive',u'Sensory',u'Vegetative',u'Emotional',u'Experiencial','Superior functions']
        self.final_list_color = by_default_list #detect differencies ? but how because there is other column that should not be taken
        model= QtGui.QStandardItemModel(len(self.final_list_color), 1)
        for index,i_list in enumerate(self.final_list_color):
          item = QtGui.QStandardItem(i_list)
          item.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled)
          item.setData(QtCore.Qt.Unchecked, QtCore.Qt.CheckStateRole)
          model.setItem(index+1, 0, item)
        self.ui.list_backColor_Condition.setModel(model) 
        
        self.ui.list_backColor_Condition.setEnabled(False)
        self.ui.list_backColor_Condition.clicked.connect(self.updateBipoleDisplay)
        
        self.mode_indv = indv_pat
        self.subsample_group = group_subsample
        self.updateBipoleDisplay()
        
    
    def updateBipoleDisplay(self):
        
        
        infoHztoshow = [self.ui.lowHz_radiobutton.isChecked(),self.ui.highHz_radiobutton.isChecked(),self.ui.bothHz_radiobutton.isChecked()]
        infoHzselected = [x for x in range(len(infoHztoshow)) if infoHztoshow[x] == True]
        colortoshow = [self.ui.fontColor_radiobutton.isChecked(),self.ui.backColor_radiobutton.isChecked()]
        colorselected = [x for x in range(len(colortoshow)) if colortoshow[x] == True]
        
        if self.mode_indv:
          bipoleObjects= dict([(self.locaData.bipoles[x]['name'],self.locaData.bipoles[x]['elecModel'].getAnatomistObjects()[0]) for x in range(len(self.locaData.bipoles))])
          dataSEEG = self.locaData.bipoleLabels['contacts']
        else:
          bipoleObjects = self.locaData.bipolesmeshes
          dataSEEG = {}
          datatomodify = {}
          for index_subsample in self.subsample_group:
              (subject,bipole) = index_subsample.split(' : ')
              if subject not in datatomodify.keys():
                  datatomodify.update({subject:[]})
              datatomodify[subject].append(bipole)   
          for subj in datatomodify.keys():
              #if subj not in dataSEEG.keys():
                  #dataSEEG.update({subj:{}})
              for index_bipole in datatomodify[subj]:
                 matches_bipoles = [x for x in range(len(self.locaData.testDataSubjects[subj]['seeg_label_all'])) if self.locaData.testDataSubjects[subj]['seeg_label_all'][x][0] == index_bipole]
                 #new_name =  subj + ' : ' + self.locaData.testDataSubjects[subj]['seeg_label_all'][index_cont][0]
                 if len(matches_bipoles)>0:
                     for index_matches_bipoles in matches_bipoles:
                         new_name =  subj + ' : ' + self.locaData.testDataSubjects[subj]['seeg_label_all'][index_matches_bipoles][0]
                         if new_name not in dataSEEG.keys():
                           dataSEEG.update({new_name:{'cell':{}}})
                         dataSEEG[new_name]['cell'].update({self.locaData.testDataSubjects[subj]['seeg_label_all'][index_matches_bipoles][1]:{}})
                         for x in self.final_list_color:
                           if x in self.locaData.testDataSubjects[subj]['seeg_label_all'][index_matches_bipoles][3]:
                              dataSEEG[new_name]['cell'][self.locaData.testDataSubjects[subj]['seeg_label_all'][index_matches_bipoles][1]].update({x:{'value':'something'}})
                           else:
                              dataSEEG[new_name]['cell'][self.locaData.testDataSubjects[subj]['seeg_label_all'][index_matches_bipoles][1]].update({x:{'value':0}})
                                 
                         if self.locaData.testDataSubjects[subj]['seeg_label_all'][index_matches_bipoles][2] == u'no response':
                            dataSEEG[new_name]['cell'][self.locaData.testDataSubjects[subj]['seeg_label_all'][index_matches_bipoles][1]].update({u'Type of response':{'fontcolor':[255,0,0,0]}})
                         elif self.locaData.testDataSubjects[subj]['seeg_label_all'][index_matches_bipoles][2] == u'seizure sensation':
                            dataSEEG[new_name]['cell'][self.locaData.testDataSubjects[subj]['seeg_label_all'][index_matches_bipoles][1]].update({u'Type of response':{'fontcolor':[255,255,0,0]}})
                         elif self.locaData.testDataSubjects[subj]['seeg_label_all'][index_matches_bipoles][2] == u'pathological':
                            dataSEEG[new_name]['cell'][self.locaData.testDataSubjects[subj]['seeg_label_all'][index_matches_bipoles][1]].update({u'Type of response':{'fontcolor':[255,0,0,255]}})    
                               
        
        if colorselected[0] == 1:
            self.ui.list_backColor_Condition.setEnabled(True)
            default_backColor = {'Sensory': [255, 0, 176, 80], 'Sensitive': [255, 0, 112, 192], 'Superior functions': [255, 148, 82, 0], 'Motor': [255, 255, 0, 0], 'Emotional': [255, 112, 48, 160], 'Vegetative': [255, 255, 192, 0], 'Experiencial': [255, 255, 133, 255]}
            backColor = {}
            Model = self.ui.list_backColor_Condition.model()
            for index in range(Model.rowCount()-1):
                Item = Model.item(index+1)
                if Item.checkState():
                    if self.mode_indv:
                      backColor.update({str(Item.text()):self.locaData.bipoleLabels['title'][str(Item.text())]['title_backcolor'][:]})
                    else:
                      backColor.update({str(Item.text()):default_backColor[str(Item.text())]})  
                else:
                    backColor.update({str(Item.text()):[255,150,150,150]})
                    
            for index_contacts in range(len(dataSEEG.keys())):
                
                #if self.locaData.bipoleLabels['contacts'].keys()[index_contacts] == u'R02 - R01':
                
                if  (u'1 Hz' in dataSEEG[dataSEEG.keys()[index_contacts]]['cell'].keys()) and (infoHzselected[0] == 0):#(self.locaData.bipoleLabels['contacts'][self.locaData.bipoleLabels['contacts'].keys()[index_contacts]]['cell'].keys()[0] == u'1 Hz' ) and (infoHzselected[0] == 0 or infoHzselected[0] == 2):
                   
                   #enComp = [self.locaData.bipoleLabels['contacts'][self.locaData.bipoleLabels['contacts'].keys()[index_contacts]]['cell']['1 Hz'][x]['value'] for x in self.final_list_color]
                   enComp = [True if dataSEEG[dataSEEG.keys()[index_contacts]]['cell']['1 Hz'][x]['value'] != 0 else False for x in self.final_list_color]
                   try:
                      indexTrue = enComp.index(True)
                      indexesTrue = [i for i,j in enumerate(enComp) if j == True]
                      selected_color = [kk for kk,vv in backColor.items() if vv != [255,150,150,150]]
                      valueselectedcolor = [self.final_list_color.index(selected_color[i]) for i in range(len(selected_color))]
                      commonSelect = list(set(indexesTrue) & set(valueselectedcolor))
                      if len(commonSelect)>1:
                        colorn = QtGui.QColor.fromRgb(255,255,255,255)
                        bipoleObjects[dataSEEG.keys()[index_contacts]].setMaterial(diffuse=[colorn.redF(), colorn.greenF(), colorn.blueF(), colorn.alphaF()])
                      else:    
                        colorn = QtGui.QColor.fromRgb(backColor[self.final_list_color[commonSelect[0]]][1],backColor[self.final_list_color[commonSelect[0]]][2],backColor[self.final_list_color[commonSelect[0]]][3],backColor[self.final_list_color[commonSelect[0]]][0])
                        bipoleObjects[dataSEEG.keys()[index_contacts]].setMaterial(diffuse=[colorn.redF(), colorn.greenF(), colorn.blueF(), colorn.alphaF()])
                   except:
                      #no response after stimulation
                      colorn = QtGui.QColor.fromRgb(0,0,0,255)
                      bipoleObjects[dataSEEG.keys()[index_contacts]].setMaterial(diffuse=[colorn.redF(), colorn.greenF(), colorn.blueF(), colorn.alphaF()]) 
                
                elif (u'1 Hz' not in dataSEEG[dataSEEG.keys()[index_contacts]]['cell'].keys()) and (infoHzselected[0] == 0):
                    colorn = QtGui.QColor.fromRgb(150,150,150,255)
                    bipoleObjects[dataSEEG.keys()[index_contacts]].setMaterial(diffuse=[colorn.redF(), colorn.greenF(), colorn.blueF(), colorn.alphaF()])
                
                elif (u'50 Hz' in dataSEEG[dataSEEG.keys()[index_contacts]]['cell'].keys()) and (infoHzselected[0] == 1):
                   enComp = [True if dataSEEG[dataSEEG.keys()[index_contacts]]['cell']['50 Hz'][x]['value'] != 0 else False for x in self.final_list_color]
                   try:
                      indexTrue = enComp.index(True)
                      indexesTrue = [i for i,j in enumerate(enComp) if j == True]
                      selected_color = [kk for kk,vv in backColor.items() if vv != [255,150,150,150]]
                      valueselectedcolor = [self.final_list_color.index(selected_color[i]) for i in range(len(selected_color))]
                      commonSelect = list(set(indexesTrue) & set(valueselectedcolor))
                      if len(commonSelect)>1:
                        colorn = QtGui.QColor.fromRgb(255,255,255,255)
                        bipoleObjects[dataSEEG.keys()[index_contacts]].setMaterial(diffuse=[colorn.redF(), colorn.greenF(), colorn.blueF(), colorn.alphaF()])                          
                      else:    
                        colorn = QtGui.QColor.fromRgb(backColor[self.final_list_color[commonSelect[0]]][1],backColor[self.final_list_color[commonSelect[0]]][2],backColor[self.final_list_color[commonSelect[0]]][3],backColor[self.final_list_color[commonSelect[0]]][0])
                        bipoleObjects[dataSEEG.keys()[index_contacts]].setMaterial(diffuse=[colorn.redF(), colorn.greenF(), colorn.blueF(), colorn.alphaF()])
                   except:
                      #no response after stimulation
                      colorn = QtGui.QColor.fromRgb(0,0,0,255)
                      bipoleObjects[dataSEEG.keys()[index_contacts]].setMaterial(diffuse=[colorn.redF(), colorn.greenF(), colorn.blueF(), colorn.alphaF()])
                
                elif (u'50 Hz' not in dataSEEG[dataSEEG.keys()[index_contacts]]['cell'].keys()) and (infoHzselected[0] == 1):
                     colorn = QtGui.QColor.fromRgb(150,150,150,255)
                     bipoleObjects[dataSEEG.keys()[index_contacts]].setMaterial(diffuse=[colorn.redF(), colorn.greenF(), colorn.blueF(), colorn.alphaF()])
                
                elif infoHzselected[0] == 2:
                    if u'1 Hz' in dataSEEG[dataSEEG.keys()[index_contacts]]['cell'].keys():
                        enComp1Hz = [True if dataSEEG[dataSEEG.keys()[index_contacts]]['cell']['1 Hz'][x]['value'] != 0 else False for x in self.final_list_color]
                    else:
                        enComp1Hz = [False for i in range(len(self.final_list_color))]
                        
                    if u'50 Hz' in dataSEEG[dataSEEG.keys()[index_contacts]]['cell'].keys():
                        enComp50Hz = [True if dataSEEG[dataSEEG.keys()[index_contacts]]['cell']['50 Hz'][x]['value'] != 0 else False for x in self.final_list_color] 
                    else:
                        enComp50Hz = [False for i in range(len(self.final_list_color))]
                        
                    enCompBoth = numpy.logical_or(enComp1Hz,enComp50Hz).tolist()
                    try:
                      indexTrue = enCompBoth.index(True)
                      indexesTrue = [i for i,j in enumerate(enCompBoth) if j == True]
                      selected_color = [kk for kk,vv in backColor.items() if vv != [255,150,150,150]]
                      valueselectedcolor = [self.final_list_color.index(selected_color[i]) for i in range(len(selected_color))]
                      commonSelect = list(set(indexesTrue) & set(valueselectedcolor))
                      if len(commonSelect)>1:
                        colorn = QtGui.QColor.fromRgb(255,255,255,255)
                        bipoleObjects[dataSEEG.keys()[index_contacts]].setMaterial(diffuse=[colorn.redF(), colorn.greenF(), colorn.blueF(), colorn.alphaF()])                          
                      else:    
                        colorn = QtGui.QColor.fromRgb(backColor[self.final_list_color[commonSelect[0]]][1],backColor[self.final_list_color[commonSelect[0]]][2],backColor[self.final_list_color[commonSelect[0]]][3],backColor[self.final_list_color[commonSelect[0]]][0])
                        bipoleObjects[dataSEEG.keys()[index_contacts]].setMaterial(diffuse=[colorn.redF(), colorn.greenF(), colorn.blueF(), colorn.alphaF()])
                    except:
                      #no response after stimulation
                      colorn = QtGui.QColor.fromRgb(0,0,0,255)
                      bipoleObjects[dataSEEG.keys()[index_contacts]].setMaterial(diffuse=[colorn.redF(), colorn.greenF(), colorn.blueF(), colorn.alphaF()])
                
                #pdb.set_trace()   #color = QtGui.QColor.fromRgb(250,0,0,100) r g b alpha     
                #self.a.setMaterial(newCyl, diffuse=[color.redF(), color.greenF(), color.blueF(), color.alphaF()]) newCyl <anatomist.cpp.ASurface_3 object at 0xa0ff0e0>

        else:
            self.ui.list_backColor_Condition.setEnabled(False)
            for index_contacts in range(len(dataSEEG.keys())):
                #self.locaData.bipoleLabels['contacts'][self.locaData.bipoleLabels['contacts'].keys()[1]]['line']['fontcolor']
                if (u'1 Hz' in dataSEEG[dataSEEG.keys()[index_contacts]]['cell'].keys()) and (infoHzselected[0] == 0):
                  fontColorLine = dataSEEG[dataSEEG.keys()[index_contacts]]['cell']['1 Hz'][u'Type of response']['fontcolor']
                  if fontColorLine[1] == 0 and fontColorLine[2] == 0 and fontColorLine[3]==0:
                    fontColorLine = [255,0,0,0]
                  colorn = QtGui.QColor.fromRgb(fontColorLine[1],fontColorLine[2],fontColorLine[3],fontColorLine[0])
                  bipoleObjects[dataSEEG.keys()[index_contacts]].setMaterial(diffuse=[colorn.redF(), colorn.greenF(), colorn.blueF(), colorn.alphaF()])
                elif (u'1 Hz' not in dataSEEG[dataSEEG.keys()[index_contacts]]['cell'].keys()) and (infoHzselected[0] == 0):
                  colorn = QtGui.QColor.fromRgb(150,150,150,255)
                  bipoleObjects[dataSEEG.keys()[index_contacts]].setMaterial(diffuse=[colorn.redF(), colorn.greenF(), colorn.blueF(), colorn.alphaF()])
                
                elif (u'50 Hz' in dataSEEG[dataSEEG.keys()[index_contacts]]['cell'].keys()) and (infoHzselected[0] == 1):
                  fontColorLine = dataSEEG[dataSEEG.keys()[index_contacts]]['cell']['50 Hz'][u'Type of response']['fontcolor']
                  if fontColorLine[1] == 0 and fontColorLine[2] == 0 and fontColorLine[3]==0:
                    fontColorLine = [255,0,0,0]
                  colorn = QtGui.QColor.fromRgb(fontColorLine[1],fontColorLine[2],fontColorLine[3],fontColorLine[0])
                  bipoleObjects[dataSEEG.keys()[index_contacts]].setMaterial(diffuse=[colorn.redF(), colorn.greenF(), colorn.blueF(), colorn.alphaF()])
                elif (u'50 Hz' not in dataSEEG[dataSEEG.keys()[index_contacts]]['cell'].keys()) and (infoHzselected[0] == 1):
                  colorn = QtGui.QColor.fromRgb(150,150,150,255)
                  bipoleObjects[dataSEEG.keys()[index_contacts]].setMaterial(diffuse=[colorn.redF(), colorn.greenF(), colorn.blueF(), colorn.alphaF()])                  
                
                if infoHzselected[0]==2:
                    if (u'1 Hz' in dataSEEG[dataSEEG.keys()[index_contacts]]['cell'].keys()):
                        fontColorLine1 = dataSEEG[dataSEEG.keys()[index_contacts]]['cell']['1 Hz'][u'Type of response']['fontcolor']
                        if fontColorLine1[1] == 0 and fontColorLine1[2] == 0 and fontColorLine1[3]==0:
                          fontColorLine1 = [255,0,0,0]
                    else:
                        fontColorLine1 = [0,0,0,0]
                        
                    if (u'50 Hz' in dataSEEG[dataSEEG.keys()[index_contacts]]['cell'].keys()):
                        fontColorLine50 = dataSEEG[dataSEEG.keys()[index_contacts]]['cell']['50 Hz'][u'Type of response']['fontcolor']
                        if fontColorLine50[1] == 0 and fontColorLine50[2] == 0 and fontColorLine50[3]==0:
                          fontColorLine50 = [255,0,0,0]                        
                    else:
                       fontColorLine50 = [0,0,0,0]
                       
                    fontColorBoth =numpy.minimum((numpy.array(fontColorLine50) + numpy.array(fontColorLine1)),255).tolist()
                    colorn = QtGui.QColor.fromRgb(fontColorBoth[1],fontColorBoth[2],fontColorBoth[3],fontColorBoth[0])
                    bipoleObjects[dataSEEG.keys()[index_contacts]].setMaterial(diffuse=[colorn.redF(), colorn.greenF(), colorn.blueF(), colorn.alphaF()])  
    
    
    def closeEvent(self, event):
        self.quit(event)        

    def quit(self, event=None):
        reply = QtGui.QMessageBox.question(self, 'Quit',"Quit bipole colors editor ?", QtGui.QMessageBox.Yes |QtGui.QMessageBox.No, QtGui.QMessageBox.No)
        if reply == QtGui.QMessageBox.Yes:
            if event is None:
                self.a.quit()
            else:
                event.accept()
        else:
            event.ignore()
        
        
if __name__=="__main__":
    a=QtGui.QApplication(sys.argv)        
