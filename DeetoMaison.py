# -*-coding:utf-8 -*

#Import of Usefull Libraries
import struct, time, re, os, json, subprocess
from soma.qt_gui.qt_backend.Qt import QtGui, QtCore, uic, Qt
import pdb
from numpy import *
from math import sqrt

from soma import aims

from brainvisa import anatomist
from brainvisa.data.writediskitem import ReadDiskItem, WriteDiskItem

from referentialconverter import ReferentialConverter
from collections import Counter


import matplotlib.pyplot as plt
from scipy import signal, stats
from locateContacts import *
import copy

#Main Class  
class DeetoMaison(QtGui.QDialog):
    
    def __init__(self, locateData=None):
        # UI init
        QtGui.QWidget.__init__(self)
        self.ui = uic.loadUi("DeetoMaison.ui", self)
        self.setWindowTitle('Automatic electrodes localization')
        self.locaData =locateData  #Data from locatesElectrode.py
        self.a=self.locaData.a #Anatomist Object
        self.dicMeshes={} 
        self.CT=True

        model= QtGui.QStandardItemModel(len(self.locaData.electrodes), 1)
        for i,area in enumerate([self.locaData.electrodes[i]['name'] for i in range(len(self.locaData.electrodes))]):
          item = QtGui.QStandardItem(area)
          item.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled)
          item.setData(QtCore.Qt.Unchecked, QtCore.Qt.CheckStateRole)
          model.setItem(i+1, 0, item)
        self.listViewElectrodes.setModel(model)
        
        self.connect(self.pushButtonOk, QtCore.SIGNAL('clicked()'), self.verifSelec)
        self.connect(self.pushButtonCheckElec, QtCore.SIGNAL('clicked()'), self.check)
        self.connect(self.pushButtonUncheckElec, QtCore.SIGNAL('clicked()'), self.uncheck)
        self.connect(self.pushButtonElecAppear, QtCore.SIGNAL('clicked()'), self.printElec)
        self.connect(self.pushButtonRemoveElec, QtCore.SIGNAL('clicked()'), self.removeElec)
        
       
  
    def verifSelec(self):
        serpentin=False
        if self.checkBoxSerpentins.checkState()==QtCore.Qt.Checked:
            serpentin=True
        
        #dans une sous fonction appelé par un connect
        Model=self.listViewElectrodes.model()
        #boucle for sur les index Item =Model.item(index)
        index=1
        elecSelec=[]
        while index<Model.rowCount():
            Item =Model.item(index)
            if Item.checkState()==QtCore.Qt.Checked:
                elecSelec.append(self.locaData.electrodes[index-1]['name'])
            index+=1 
            
        T1post=None   
        diCT = ReadDiskItem( 'CT', 'BrainVISA volume formats', requiredAttributes={'center':self.locaData.brainvisaPatientAttributes['center'], 'subject':self.locaData.brainvisaPatientAttributes['subject'] } )
        path = list(diCT.findValues({}, None, False ))
        #check if it is a CT post
        idCTpost = [x for x in range(len(path)) if 'CTpost' in str(path[x])]
        try:
            path=path[0].fullPath()
            volCT = aims.read(path)
            npCT = volCT.arraydata()
        except: 
            diMRI = ReadDiskItem( 'Raw T1 MRI', 'BrainVISA volume formats', requiredAttributes={'center':self.locaData.brainvisaPatientAttributes['center'], 'subject':self.locaData.brainvisaPatientAttributes['subject'] } )
            pathMRI = list(diMRI.findValues({}, None, False ))         
            id_post = [x for x in range(len(pathMRI)) if 'post' in str(pathMRI[x]) and not 'postOp' in str(pathMRI[x])]
            pathMRI=pathMRI[0].fullPath()
            volCT = aims.read(pathMRI)
            npCT = volCT.arraydata()
            if id_post!=[]:    
                T1post = pathMRI[id_post[0]]
                 
        
        #si il y a ni ct post ni mri post return
        #si ct et mri post prendre ct
        if T1post is None and idCTpost==[]:
            warning = QtGui.QMessageBox(self)
            warning.setText("No CT post or MRI post found")
            warning.setWindowTitle("Warning")
            return
            
        #Recupere taille des voxels (de l'image post) CT if CT
        if idCTpost!=[]:
            volume = aims.read(str(self.locaData.diskItems['CTpost']))
            size=volume.getVoxelSize()
        ##if MRI
        if T1post is not None:
            self.CT=False
            volume = aims.read(str(self.locaData.diskItems['T1post']))
            size=volume.getVoxelSize()
        #Recupere taille des voxels T1
        volume = aims.read(str(self.locaData.diskItems['T1pre']))
        sizeT1=volume.getVoxelSize()
        
        #recupere tous les plots   
        contacts=self.locaData.getAllPlotsCentersT1preRef()
        triTempo=sorted(contacts.keys())
        tri=[]
        for el in elecSelec:
            trit=[x for x in triTempo if x[:len(el)]==el and x[len(el):len(el)+1]=="-"]
            for el in trit:
                tri.append(el)
        electrodes={}
        i=0
        j=0
        nbContacts={}
        listeCle=[]
        
        if self.CT==True:
            #transfo natif to SB t1
            rdi = ReadDiskItem('Transformation to Scanner Based Referential', 'Transformation matrix', exactType=True,  requiredAttributes={'modality':'t1mri','subject':self.locaData.brainvisaPatientAttributes['subject'], 'center':self.locaData.brainvisaPatientAttributes['center']})
            diNat2SBT1pre = rdi.findValue(self.locaData.diskItems['T1pre'])

            #transfo SB CT vers SB T1
            wdiTransform = ReadDiskItem('Transform CT to another image', 'Transformation matrix', exactType=True, requiredAttributes = {'subject':self.locaData.brainvisaPatientAttributes['subject'], 'center':self.locaData.currentProtocol }) 
            diSBCT2SBT1 = wdiTransform.findValue(self.locaData.diskItems['CTpost'] )
        
            #transfo natif to SB CT
            rdi = ReadDiskItem('Transformation to Scanner Based Referential', 'Transformation matrix', exactType=True, requiredAttributes={'modality':'ct','subject':self.locaData.brainvisaPatientAttributes['subject'], 'center':self.locaData.brainvisaPatientAttributes['center']})
            diNat2SBCT = rdi.findValue(self.locaData.diskItems['CTpost'])

            #path temporaire pour calculer les matrices

            trmpre_to_SBinvpath = diNat2SBT1pre.fullPath().split('/')
            trmpre_to_SBinvpath[-1] = 'inv'+trmpre_to_SBinvpath[-1]
            trmCTpost_to_T1pre_path = copy.deepcopy(trmpre_to_SBinvpath)
            trmCTpost_to_T1pre_path[-1] = 'CTpost_to_T1pre.trm'
            trmT1pre_to_CTpost_path = copy.deepcopy(trmpre_to_SBinvpath)
            trmT1pre_to_CTpost_path[-1] = 'T1pre_to_CTpost.trm'
            trmpre_to_SBinvpath = '/'.join(trmpre_to_SBinvpath)
            trmCTpost_to_T1pre_path = '/'.join(trmCTpost_to_T1pre_path)
            trmT1pre_to_CTpost_path = '/'.join(trmT1pre_to_CTpost_path)


            ret = subprocess.call(['AimsInvertTransformation','-i',diNat2SBT1pre.fullPath(),'-o', trmpre_to_SBinvpath])
            ret = subprocess.call(['AimsComposeTransformation', '-o',trmCTpost_to_T1pre_path, trmpre_to_SBinvpath, diSBCT2SBT1.fullPath(), diNat2SBCT.fullPath()])
            ret = subprocess.call(['AimsInvertTransformation','-i',trmCTpost_to_T1pre_path,'-o',trmT1pre_to_CTpost_path])
          
            #transfo final natif pre vers natif postOp
            transfo_pre_to_postop = aims.read(trmT1pre_to_CTpost_path).toMatrix()
            transfo_pre_to_postopInv = aims.read(trmCTpost_to_T1pre_path).toMatrix()
        
        else:
            #transfo natif to SB t1
            rdi = ReadDiskItem('Transformation to Scanner Based Referential', 'Transformation matrix', exactType=True,  requiredAttributes={'modality':'t1mri','subject':self.locaData.brainvisaPatientAttributes['subject'], 'center':self.locaData.brainvisaPatientAttributes['center']})
            diNat2SBT1pre = rdi.findValue(self.locaData.diskItems['T1pre'])
            
            #transfo SB T1 post vers SB T1 pre
            wdi = ReadDiskItem('Transform Raw T1 MRI to another image', 'Transformation matrix', exactType=True, requiredAttributes = {'subject':self.locaData.brainvisaPatientAttributes['subject'], 'center':self.locaData.currentProtocol })
            diSBCT2SBT1 = wdi.findValue(self.locaData.diskItems['T1post'])
            
            #transfo natif to SB T1post
            rdi = ReadDiskItem('Transformation to Scanner Based Referential', 'Transformation matrix', exactType=True, requiredAttributes={'modality':'t1mri','subject':self.locaData.brainvisaPatientAttributes['subject'], 'center':self.locaData.brainvisaPatientAttributes['center']})
            diNat2SBCT = rdi.findValue(self.locaData.diskItems['T1post'])
            
            
            trmpre_to_SBinvpath = diNat2SBT1pre.fullPath().split('/')
            trmpre_to_SBinvpath[-1] = 'inv'+trmpre_to_SBinvpath[-1]
            trmCTpost_to_T1pre_path = copy.deepcopy(trmpre_to_SBinvpath)
            trmCTpost_to_T1pre_path[-1] = 'CTpost_to_T1pre.trm'
            trmT1pre_to_CTpost_path = copy.deepcopy(trmpre_to_SBinvpath)
            trmT1pre_to_CTpost_path[-1] = 'T1pre_to_CTpost.trm'
            trmpre_to_SBinvpath = '/'.join(trmpre_to_SBinvpath)
            trmCTpost_to_T1pre_path = '/'.join(trmCTpost_to_T1pre_path)
            trmT1pre_to_CTpost_path = '/'.join(trmT1pre_to_CTpost_path)


            ret = subprocess.call(['AimsInvertTransformation','-i',diNat2SBT1pre.fullPath(),'-o', trmpre_to_SBinvpath])
            ret = subprocess.call(['AimsComposeTransformation', '-o',trmCTpost_to_T1pre_path, trmpre_to_SBinvpath, diSBCT2SBT1.fullPath(), diNat2SBCT.fullPath()])
            ret = subprocess.call(['AimsInvertTransformation','-i',trmCTpost_to_T1pre_path,'-o',trmT1pre_to_CTpost_path])
          
            #transfo final natif pre vers natif postOp
            transfo_pre_to_postop = aims.read(trmT1pre_to_CTpost_path).toMatrix()
            transfo_pre_to_postopInv = aims.read(trmCTpost_to_T1pre_path).toMatrix()

            
            
        d={}
        dicPointsFinal={}
        while i<len(tri):
          k=1  
          t=100
          ent=0
          while k<t:
            if tri[i][k]=="-":
                cle=tri[i][:k]
                triElec=[x for x in tri if x[:k]==cle and x[len(cle):len(cle)+1]=="-"]
                t=k
            k+=1
          listeCle.append(cle)
          tar=contacts[triElec[0]]
          tar=tuple(tar)
          for el in triElec:
              plot=el[-2:]
              if plot[-2]!="t":
                  ent=contacts[el]   
          if ent==0:
              ent=contacts[triElec[-1]]                  
          ent=tuple(ent)
          
          tar=list(tar)
          ent=list(ent)
          tar.append(1)
          ent.append(1)
          ent = array([ent])
          tar = array([tar])
        
          #calcul des vecteurs dans le natif postOp
          ent = transfo_pre_to_postop.dot(ent.T)
          tar = transfo_pre_to_postop.dot(tar.T)
          newtar=[x[0] for x in tar]
          del newtar[-1]
          newent=[x[0] for x in ent]
          del newent[-1]
          tar=tuple(newtar)
          ent=tuple(newent)
          i=tri.index(triElec[-1])+1
          dic={"entry":ent,"target":tar}
          nbContact={cle:len(triElec)}
          nbContacts.update(nbContact)
          electrode={cle:dic}
          electrodes.update(electrode)
          p12=vecteur(ent,tar)
          p12Norm=norme3D(p12)
          tritemp=[1,1]
          jo=1
          dist={0:p12Norm}
          while len(tritemp)==2:
            tritemp=[x for x in triElec if (int(x[-1])==jo or int(x[-1])==jo+1)and x[-2]=="t"]
            if len(tritemp)==2:
                p12=vecteur(contacts[tritemp[0]],contacts[tritemp[1]])
                p12Norm=norme3D(p12)
                nor={jo:p12Norm}
                dist.update(nor)
                jo+=1
          to=0 
          tritemp=[x for x in triElec if (int(x[-1])==9 and x[-2]=="t")or (int(x[-1])==0 and x[-2]!="t")]

          if len(tritemp)==2:
              p12=vecteur(contacts[tritemp[0]],contacts[tritemp[1]])
              p12Norm=norme3D(p12)
              nor={jo:p12Norm}
              dist.update(nor)
              jo+=1   
          while len(tritemp)==2:
            tritemp=[x for x in triElec if (int(x[-1])==to or int(x[-1])==to+1)and x[-2]!="t"]
            if len(tritemp)==2:
                p12=vecteur(contacts[tritemp[0]],contacts[tritemp[1]])
                p12Norm=norme3D(p12)
                nor={jo:p12Norm}
                dist.update(nor)
            jo+=1
            to+=1
          points={}  
          for el in triElec:  
            if re.findall("([0-9])",el[-1])!=None and el[-2]=="t":
                point={int(el[-1]):contacts[el]} 
            else:
                point={int(el[-2:]):contacts[el]}
            for key, vect in point.items():
                vect=list(vect)
                vect.append(1)
                vect = array([vect])
        
                #calcul des vecteurs dans le natif postOp
                vect = transfo_pre_to_postop.dot(vect.T)
                
                vect=[x[0] for x in vect]
                del vect[-1]
                vect=tuple(vect)
                point[key]=vect
            points.update(point)      
          dicPoints={cle:points}  
          dicPointsFinal.update(dicPoints)
          norme={cle:dist}
          d.update(norme)
        approxElec={}   
      
        try:
            brainMask = ReadDiskItem('Brain Mask', 'aims readable volume formats',requiredAttributes={'subject':self.locaData.brainvisaPatientAttributes['subject'], 'center':self.locaData.currentProtocol })
            diBrain = list(brainMask.findValues({}, None, False ))
            diBrain0=diBrain[0].fullPath()
            volBrainMask = aims.read(diBrain0)
            brainMaskArray = volBrainMask.arraydata() 
        except:
            brainMaskArray=None
              
        
        while j<len(electrodes):
          foo=locateContacts(electrodes[listeCle[j]]["target"],electrodes[listeCle[j]]["entry"],npCT,volCT,nbContacts[listeCle[j]],size[0],size[1],size[2],d[listeCle[j]],transfo_pre_to_postopInv,brainMaskArray,sizeT1,dicPointsFinal[listeCle[j]],serpentin,transfo_pre_to_postop,self.CT)
          dic={listeCle[j]:foo}
          approxElec.update(dic)
          print "electrode: ",listeCle[j],"number: " ,j+1, "of", len(electrodes)  
          j+=1
        self.approxElec=approxElec   
    
    def printElec(self):  
        diameter=1
        listPlots=[]
        for eleckeys,elec in self.approxElec.items():
            name=eleckeys
            for plot in elec.keys():
                #tu parcours approxElec pour générer les maillage      
                mesh = self.a.toAObject(aims.SurfaceGenerator.sphere(aims.Point3df(elec[plot][0], elec[plot][1], elec[plot][2]), diameter, 64))
                mesh.setName(name+str(plot+1))
                listPlots.append(mesh)
                self.a.setMaterial(mesh, diffuse=(0.0,0.9,0.1,1.0))#[color.redF(), color.greenF(), color.blueF(), color.alphaF()]
                self.a.assignReferential(self.locaData.preReferential(), mesh)
                self.a.addObjects(mesh, self.locaData.wins[0])
                #self.locaData.DeetoEstimated = approxElec
            meshes={name:listPlots}    
            self.dicMeshes.update(meshes)
        
    def removeElec(self):
        #dans une sous fonction appelé par un connect
        Model=self.listViewElectrodes.model()
        #boucle for sur les index Item =Model.item(index)
        index=1
        elecSelec=[]
        listEl=[]
        while index<Model.rowCount():
            Item =Model.item(index)
            if Item.checkState()==QtCore.Qt.Checked:
                for el in self.dicMeshes[str(Item.text())]:
                    self.a.removeObjects([el,], self.locaData.wins) # Remove from windows
                    listEl.append(el)
            index+=1 
        self.a.deleteObjects(listEl) # CURRENT    
        
        
    def check(self):
        model = self.listViewElectrodes.model()
        for index in range(model.rowCount()-1):
           item = model.item(index+1)
           if item.isCheckable():
                  item.setCheckState(QtCore.Qt.Checked)
        
    def uncheck(self):
        model = self.listViewElectrodes.model()
        for index in range(model.rowCount()-1):
           item = model.item(index+1)
           if item.isCheckable():
                  item.setCheckState(QtCore.Qt.Unchecked)

    def closeEvent(self, event):
        self.quit(event)

    def quit(self, event=None):
        reply = QtGui.QMessageBox.question(self, 'Message',"Annuler la segmentation?", QtGui.QMessageBox.Yes |QtGui.QMessageBox.No, QtGui.QMessageBox.No)
        if reply == QtGui.QMessageBox.Yes:
            if event is None:
                self.a.quit()
            else:
                event.accept()
        else:
            event.ignore()

if __name__=="__main__":
    a=QtGui.QApplication(sys.argv)