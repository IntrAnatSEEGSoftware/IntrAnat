#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os, pickle,numpy, re, string, json, csv, pdb


from soma.qt_gui.qt_backend import QtGui, QtCore, uic, Qt



# self.connect(self.ImportTheoriticalImplantation,QtCore.SIGNAL('clicked()'),self.importRosaImplantation)


def importRosaImplantation(locate_env):
   fichierRosa =  QtGui.QFileDialog.getOpenFileName(locate_env, "Select a file containing electrode coordinates: ", "", "(*.ros )")
   electrodes = ecritureCSVetJSON(fichierRosa)
   
   return electrodes


def ecritureCSVetJSON(fichier):
    target=[]
    nomElec=[]
    entry=[]
    electrodes={}
    kv=""
    i=0
    j=0
    
    a=open(fichier,"rU")
    fjson=open("electrodes.json","w")
    
    line=a.readlines()
    pdb.set_trace()
    line0=[nb for nb in line if len(nb)>50]
    line1=[nb for nb in line0 if (re.findall("([a-z])",nb[0].lower())!=None) and (nb[2]==" " or nb[1]==" ")]

    
    while j<len(line):
        line3=line[j]
        if line3[:12]=="[TRAJECTORY]":
            line3=line[j+1]
            break
        j+=1
        
    if j==len(line):
        raise ValueError("Pas de \"[TRAJECTORY]\" trouve dans le fichier")
        
    for el in line1:
        nm=el.split(" ")
        nomElec.append(nm[0])
        tar=nm[8],nm[9],nm[10]
        target.append(tar)
        ent=nm[4],nm[5],nm[6]
        entry.append(ent)   
    
    while i<len(line1):
        dic={"entry":entry[i],"target":target[i]}
        electrode={nomElec[i]:dic}
        electrodes.update(electrode)
        i+=1
        

    if int(line3)!=i:
        raise ValueError("Electrodes comptees et nombre d'electrodes inscrit different")
        
    with open("electrodes.csv","w") as csvfile:        
       fjson.write(str(electrodes))    
        
    with open("electrodes.csv","w") as csvfile:
       fcsv=csv.writer(csvfile, delimiter='\t')
       fcsv.writerow(["Electrodes", "Entry", "Target"])
       fcsv.writerow([i])
       
       for key, value in electrodes.items():
           kv=key+" "+str(value.values()).replace('\'','')
           kv=kv.replace('[','')
           kv=kv.replace(']','')
           kv=kv.replace(")","")
           kv=kv.split(" (")
           kv[1]=kv[1][:-1]
           fcsv.writerow(kv)
           
    a.close()
    fjson.close()
    
    return electrodes