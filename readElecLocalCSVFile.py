#
# From Electrodes localisation, read CSV file coming from INTRANAT
#

import csv, sys, os, re
from datetime import date
import pdb

def readElecLocalCSVFile(infile = None):    

    data = csv.reader(open(infile), delimiter = '\t')
    contents = []
    for row in data:
        contents.append(row)

    saveSpaceL = 0
    #find the string 'contact' to identify the starting line
    for i1 in range(len(contents[1])):
        for i2 in range(10):
            try:
                if contents[i2][i1] == 'contact':
                    savec = i1
                    savel = i2
                if contents[i2][i1] == 'Use of MNI Template':
                    saveSpaceL = i2
            except:
                continue

    props = {}

    if saveSpaceL == 0:
        props['MarsAtlas'] = 'specific'
        props['Freesurfer'] = 'empty'
        props['HippoSubfieldFreesurfer'] = 'empty'
    else:
        if contents[saveSpaceL][contents[saveSpaceL].index('MarsAtlas')+1] == 'True':
            props['MarsAtlas'] = 'MNI'
        else:
            props['MarsAtlas'] = 'specific'

        if contents[saveSpaceL][contents[saveSpaceL].index('Freesurfer')+1] =='True':
            props['Freesurfer'] = 'MNI'
        else:
            props['Freesurfer'] = 'specific'

        if contents[saveSpaceL][contents[saveSpaceL].index('HippoSubfieldFreesurfer')+1] == 'True':
            props['HippoSubfieldFreesurfer'] = 'MNI'
        else:
            props['HippoSubfieldFreesurfer'] = 'specific'

    startl = savel + 1
    #extract monopolar and bipolar contacts, assuming each list is separated by a two-lines blank
    stimpolarity = ['monopolar', 'bipolar']
    ContactsInfos = {'monopolar':{},'bipolar':{}}
    i = -2
    for pol in stimpolarity:
        i = i + 2
        while contents[startl+i]:
            # standard name : 'Ap12' (not A'12, Ap21, a'1a'2...)
            nametmp = contents[startl+i][savec]
            
            testpol = re.search(" - ",nametmp)
            if testpol is not None:
               polarity = 'bipolar'
            else:
               polarity = 'monopolar'

            resection = []
            if contents[startl+i][contents[savel].index('Resection')] == 'not in resection':
                resection = False
            elif not contents[startl+i][contents[savel].index('Resection')] == 'resection not calculated':
                resection = True 

            #if not Contact.objects.filter(CRF = crf, name = name, polarity = polarity):
            if nametmp not in ContactsInfos.keys():
              ContactsInfos[polarity].update({nametmp:{}})
            #newcontact = Contact(CRF = crf, name = name, polarity = polarity, resection = resection, segmentation = greywhite)

            for index_column in contents[savel]:
                if index_column != 'contact':
                    
                    ContactsInfos[polarity][nametmp].update({index_column:contents[startl+i][contents[savel].index(index_column)]})

            i += 1
    
    ResectionInfos={}
    i=i+2
    try:
       while contents[startl+i]:
          if contents[startl+i]==['Resection Information']:
              actualAtlas = ''
          elif contents[startl+i]==['mars_atlas']:
              actualAtlas = 'mars_atlas'
              ResectionInfos.update({'mars_atlas':{}})
          elif contents[startl+i]==['Freesurfer']:
              actualAtlas = 'Freesurfer'
              ResectionInfos.update({'Freesurfer':{}})
          elif contents[startl+i][0] == 'Volume resection (mm3):' :
              ResectionInfos.update({'Volume resection (mm3):':contents[startl+i][1]})
          else:
              ResectionInfos[actualAtlas].update({contents[startl+i][0]:contents[startl+i][1]})
          i += 1
       
    except:
       pass
   
    return (props,ContactsInfos,ResectionInfos)
    
            #try:
              #ContactsInfos[nametmp].update({'MarsAtlas':contents[startl+i][contents[savel].index('MarsAtlas')]})
            #except:
              #pass
          
            #try:
              #ContactsInfos[nametmp].update({'Freesurfer':contents[startl+i][contents[savel].index('Freesurfer')]})
            #except:
              #pass          

            #try:
              #ContactsInfos[nametmp].update({'HippoSubfieldFreeSurfer':contents[startl+i][contents[savel].index('HippoSubfieldFreeSurfer')]})
            #except:
              #pass 

            #try:
              #ContactsInfos[nametmp].update({'MNI':contents[startl+i][contents[savel].index('MNI')]})
            #except:
              #pass 
            
            #try:
              #ContactsInfos[nametmp].update({'AAL':contents[startl+i][contents[savel].index('AAL')]})
            #except:
              #pass
          
            #try:
              #ContactsInfos[nametmp].update({'AALDilate':contents[startl+i][contents[savel].index('AALDilate')]})
            #except:
              #pass          

            #try:
              #ContactsInfos[nametmp].update({'Broadmann':contents[startl+i][contents[savel].index('Broadmann')]})
            #except:
              #pass 

            #try:
              #ContactsInfos[nametmp].update({'Hammers':contents[startl+i][contents[savel].index('Hammers')]})
            #except:
              #pass 