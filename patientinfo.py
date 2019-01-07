# -*- coding: utf-8 -*-
#
# Functions to create UI, read and write patient data
# such as type of epilepsy, disease history...
#, and dictionary of available processing methods for each manip type
#
# (c) Inserm U836 2012-2014 - Manik Bhattacharjee
#
# License GNU GPL v3

from soma.qt_gui.qt_backend import QtGui, QtCore
import pdb, os, json
from brainvisa.data.writediskitem import ReadDiskItem, WriteDiskItem

# Types : Text, Line (short text) Int (integer), Float, Bool, Choice, MultiSelectChoice, EditableChoice (list of values, but a new one may be entered)
#function_mappings = {'checkbox_comor':checkbox_comor}

defaultInfoTable= ('currentDate',('Date2', QtCore.QDate.currentDate(),u'Date:'),
           'patientBirthday',('Date2',None , u'Date de Naissance'),
           'patientAge',('Line',u'0','Age patient'),
           'sex',('Choice',[u'Inconnu', u'Féminin', u'Masculin'], u'Sexe'),
           'handedness', ('Choice', ['Inconnu', 'Droitier', 'Gaucher', 'Ambidextre'], u'Latéralité manuelle'),
           'Consentement', ('Date2', None, u'Consentement Date Signature'),
           'ProtocolsRecherche',('Line',None,u'Protocol'),
           'ProtocolDate',('Date2',None, u"Date d'entrée dans le protocol"),
           'previousHistory', ('Choice2', [u'Inconnu', u'Aucun', u'Crises fébriles', u'Infection SNC', u'Souffrance néonatale', u'TC grave', u'Vasculaire'], u'Antécédents Personnels'),
           'previousFamilyHistory', ('Choice2', [u'Inconnu', u'Dépression', u'Epilepsie', u'Parkinson', u'TOC'], u'Antécédents Familiaux'),
           'causaleDisease', ('Choice2', [u'Inconnu', u'Génétique', u'Inflammatoire', u'Métabolique', u'Phacomatose'], u'Maladie neurologique causale'),
           'mriLesion', ('Choice2', [u'Inconnu', u'Aucune', u'Cicatrice', u'Malformation vasculaire', u'MCD', u'Sclérose hippocampique', u'Tumeur', u'Multiple'], u'Lésion IRM'),
           'mrilesionDroite',('ButtonBox',[u'T', u'F',u'P', u'O', u'Ins'], u'Droite'),
           'mrilesionGauche',('ButtonBox',[u'T', u'F',u'P', u'O', u'Ins'], u'Gauche'),
           'neurologiqueExam', ('Choice', [u'Inconnu', u'Normal', u'Anormal'],u'Examen Neurologique'),
           #'comorbiditylabel', ('Label',u'Comorbiditépouet'),
           'comoraucune', ('Checkbox',u'Aucune',u'Comorbidité','checkbox_comor'),
           'comorneuro', ('Choice2', [u'Inconnu', u'AVC', u'Dystonie', u'Epilepsie', u'Migraine', u'Parkinson'], u'comorbidité neurologique'),
           'comorpsy', ('Choice2', [u'Inconnu', u'Anxiété', u'Dépression', u'TOC'], u'comorbidité psychiatrique'),
           'comorautre', ('Choice2', [u'Inconnu', u'Parasomnies', u'SAS'], u'comorbidité autre')
           #'comorbidity', ('ChoiceSub', [u'Inconnu', u'Aucune', u'Neurologique', u'Psychiatrique', u'Autre'],u'Comorbidité'),
           #'medication', ('Compound',(
           #  'antiseizure',('Line', None, u'Antiépileptique'),
           #  'other',('Line', None, u'Autre'),), u'Traitement'),
           #'seizureFrequency', ('Choice', [' < 1/mois', 'mensuelle', 'hebdomadaire', 'quotidienne'], u'Fréquence des crises'),
           #'diseaseHistory': ('Text', u'Histoire de la maladie'),
           #'exams', ('Compound', (
           #    'checkExams', ('MultiSelectChoice', ['Evaluation neuropsy', 'Evaluation psy', 'Test de Wada', 'IRMf', 'SPECT', 'PET'], ''),
           #    'veegExam', ('Date', None, 'VEEG'),
           #    'megExam', ('Date', None, 'VEEG'),
           #    'fmriEegExam', ('Date', None, 'VEEG'),
           # ), u'Examens réalisés'),
           #'epilepsyType', ('Line',None, u'Type d\'épilepsie (avant SEEG)')
           )


def buildInt(parent, mini, maxi, label):
  lay = QtGui.QHBoxLayout()
  sp = QtGui.QSpinBox(parent)
  sp.setRange(mini, maxi)
  lay.addWidget(QtGui.QLabel(label, parent=parent))
  lay.addWidget(sp)
  return (lay, sp)

def buildChoice(parent, values, label):
  lay = QtGui.QHBoxLayout()
  combo = QtGui.QComboBox(parent)
  combo.clear()
  combo.addItems(values)
  if label is not None and len(label) > 0:
    lay.addWidget(QtGui.QLabel(label, parent=parent))
  lay.addWidget(combo)
  return (lay, combo)

def buildChoice2(parent, values, label):
  lay = QtGui.QHBoxLayout()
  combo = QtGui.QComboBox(parent)
  combo.clear()
  model = QtGui.QStandardItemModel(len(values), 1)# 5 rows, 1 col

  for i,area in enumerate(values):
    item = QtGui.QStandardItem(area)
    item.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled)
    item.setData(QtCore.Qt.Unchecked, QtCore.Qt.CheckStateRole)
    model.setItem(i+1, 0, item)

  combo.setModel(model)
  combo.setEditable(True)
  combo.activated.connect(lambda index, comb=combo: validatenewCombo(comb, index))
  #pdb.set_trace()

  #combo.addItems(values)
  if label is not None and len(label) > 0:
    lay.addWidget(QtGui.QLabel(label, parent=parent))
  lay.addWidget(combo)
  return (lay, combo)

def buildChoiceSub(parent, values, label):
  lay = QtGui.QHBoxLayout()

  combosub = QtGui.QMenu(parent)
  combosub.setTitle(u'Comorbitité menu')

  for i,area in enumerate(values):
    menu = QtGui.QMenu(area)

    combosub.addMenu(menu)

    #pdb.set_trace()
    if i<2:
      menu.addAction(QtGui.QAction(area, menu, checkable = True))
    else:
      #for ii, area2 in enumerate(subvalues)
      menu.addAction(QtGui.QAction("pouet", menu, checkable = True))



  if label is not None and len(label) > 0:
    lay.addWidget(QtGui.QLabel(label, parent=parent))
  lay.addWidget(combosub)
  return (lay, combosub)

def buildDate(parent, date, label):
  lay = QtGui.QHBoxLayout()
  qdate = QtGui.QDateEdit(parent)
  qdate.setDisplayFormat('dd/MM/yyyy')
  if date is None:
     qdate.setDate(QtCore.QDate.fromString('1900-01-01', 'yyyy-MM-dd'))
  else:
     qdate.setDate(date)
  lay.addWidget(QtGui.QLabel(label,parent=parent))
  lay.addWidget(qdate)
  return (lay, qdate)


def buildText(parent, label):
  lay = QtGui.QHBoxLayout()
  text = QtGui.QTextEdit(parent)
  lay.addWidget(QtGui.QLabel(label, parent=parent))
  lay.addWidget(text)
  return (lay, text)

def buildLine(parent, defaultValue, label):
  lay = QtGui.QHBoxLayout()
  line = QtGui.QLineEdit(parent)
  if label is not None and len(label) > 0:
    lay.addWidget(QtGui.QLabel(label, parent=parent))
  if defaultValue:
    line.setText(str(defaultValue))
  lay.addWidget(line)
  return (lay, line)

def buildCompound(parent, values, label):
  lay = QtGui.QHBoxLayout()
  frame = QtGui.QGroupBox(parent)
  frame.setLayout(lay)
  if label is not None:
    frame.setTitle(str(label))
  widgets = buildUI(frame, values)
  for wid in reversed(widgets.values()):
    lay.addWidget(wid)
  return (frame, [frame,]+widgets.values())

def buildButtonBox(parent,values,label):
   lay = QtGui.QHBoxLayout()
   buttbox = QtGui.QDialogButtonBox(parent)
   for i,area in enumerate(values):
     buttbox.addButton(area,3)

   for test in buttbox.findChildren(QtGui.QPushButton):
     test.setCheckable(True)

   buttbox.clicked.connect(lambda index, comb=buttbox: MRIlesionbox(comb, index))
   lay.addWidget(QtGui.QLabel(label, parent=parent))
   lay.addWidget(buttbox)
   return (lay,buttbox)

def buildLabel(parent,label):
  lay = QtGui.QHBoxLayout()
  labe = QtGui.QLabel(parent)
  labe.setText(label)

  lay.addWidget(labe)
  return (lay, labe)



def buildCheckbox(parent, label, label2, checkbox_connect):
  #function_mapping = {'checkbox_comor':checkbox_comor}
  lay = QtGui.QHBoxLayout()
  qcheck = QtGui.QCheckBox(parent)
  qcheck.setText(label)

  #qcheck.clicked.connect(lambda index, ChecK=qcheck: function_mapping[checkbox_connect](ChecK, index))

  lay.addWidget(QtGui.QLabel(label2, parent=parent))
  lay.addWidget(qcheck)
  return (lay, qcheck)



uiBuild = {'Int': buildInt,
           'Float':lambda:(None, None),# Not implemented
           'Line':buildLine,
           'Text':buildText,
           'Bool':lambda:(None, None),# Not implemented
           'Date':lambda p, value,label:(None, None),# Not implemented
           'Date2':buildDate,
           'Choice':buildChoice,
           'Choice2':buildChoice2,
           'ChoiceSub':buildChoiceSub,
           'ButtonBox':buildButtonBox,
           'MultiSelectChoice':lambda p, values, label:(None, None),# Not implemented
           'EditableChoice':lambda p,v,l : buildChoice(p,v,l,True),
           'Compound':buildCompound,
           'Label':buildLabel,
           'Checkbox':buildCheckbox,
          }

def buildUI(parent, infoTable=defaultInfoTable):
  """ Generates a UI from infoTable) and returns a dictionary with the generated widgets"""

  #si c est un choice2 c'est editable du coup faut aller chercher sur le json la derniere version de la liste
  wdi_global = ReadDiskItem('PatientInfoTemplate','Patient Template format')
  di_global = list(wdi_global.findValues({}, None, False))

  if len(di_global) > 0:
     if os.path.isfile(str(di_global[0])):
              print "read previous patienttemplate json"
              from codecs import open as opencodecs
              fin = opencodecs(str(di_global[0]),'r','latin1')
              info_dicti = json.loads(fin.read().decode('latin1'))
              fin.close()

              previous_lists_not_path = info_dicti['notPathoSpecific']
     else:
        previous_lists_not_path = {}
  else:
     previous_lists_not_path = {}

  lay = QtGui.QVBoxLayout(parent)
  lay.setSpacing(0)
  output = {}
  for i in range(int(len(infoTable)/2)):
    key = infoTable[2*i]
    value = infoTable[2*i+1]
    if key in previous_lists_not_path.keys():
       default_list = set(value[1])
       update_list = set(previous_lists_not_path[key])
       diff_list = update_list-default_list
       final_list = value[1] + list(diff_list)
       item1 = value[0]
       item3= value[2]
       value = (item1,final_list,item3)

    if value is not None:
      (l,wid) = uiBuild[value[0]](parent, *(value[1:]))
      if l is not None:
        if isinstance(l, QtGui.QLayout):
          lay.addLayout(l)
        elif isinstance(l, QtGui.QWidget):
          lay.addWidget(l)
        else:
          print "Neither a layout nor a widget : cannot add object to the layout !"
      if wid is not None:
        output[key]=wid

  return output

def dumpState(structuredWidgets):
  """Generates a dumpable object that represents the content of the patient info from the widgets"""
  pass

def validatenewCombo(combo_modified, selectedIndex):

  AllItems = [unicode(combo_modified.itemText(i)) for i in range(combo_modified.count()) if unicode(combo_modified.itemText(i)) != '']

  model = QtGui.QStandardItemModel(len(AllItems), 1)

  for i,area in enumerate(AllItems):
    item = QtGui.QStandardItem(area)
    item.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled)
    item.setData(QtCore.Qt.Unchecked, QtCore.Qt.CheckStateRole)
    model.setItem(i+1, 0, item)

  combo_modified.clear()
  combo_modified.setModel(model)
  combo_modified.setEditable(True)


def MRIlesionbox(buttbox_modiffied,selectedIndex):

   button_sender = buttbox_modiffied.sender().text()
   #print button_sender
   #buttbox_modiffied.sender().font().setBold(True)

   if buttbox_modiffied.sender().isChecked():
    buttbox_modiffied.sender().setStyleSheet("background-color: rgb(90, 255, 95);")
   else:
    buttbox_modiffied.sender().setStyleSheet("")


#def checkbox_comor(qcheck, qindex):
  #self.patientInfo['comoraucune'].isChecked()
  #pdb.set_trace()
