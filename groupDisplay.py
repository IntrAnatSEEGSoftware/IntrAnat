#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Application to view multi-subject data
#
# (c) Inserm U836 2012-2014 - Manik Bhattacharjee
#
# License GNU GPL v3
#
#
from soma.qt_gui.qt_backend import QtGui, QtCore, uic

import sys, pdb


from brainvisa import axon
from electrodedisplaywidget import ElectrodeDisplayWidget
from templatewidget import TemplateWidget
from patientdatafilterwidget import PatientDataFilterWidget


class GroupDisplayWidget(QtGui.QTabWidget):
  def __init__(self, app=None):
    QtGui.QTabWidget.__init__(self)
    # Load patient selection panel
    self.patientsPanel = PatientDataFilterWidget()
    # Load template panel (choose the common referential - MNI/Dartel custom template...)
    self.templatePanel = TemplateWidget() #Create a widget from uic.loadUi("groupPlots.ui", self)
    # Load plots panel
    self.plotsPanel = ElectrodeDisplayWidget(dataSubjects=self.patientsPanel.subjects) #Create a widget from uic.loadUi("groupPlots.ui", self)
    # Add them to the tabwidget itself
    self.addTab(self.patientsPanel, u"Subjects")
    self.addTab(self.templatePanel, u"Template")
    self.addTab(self.plotsPanel, u"Plots")
    self.setTabEnabled(1,False)
    # Connect tab change event to update the list of patients
    self.currentChanged.connect(self.tabSelected)

  def tabSelected(self, id):
    if id == 2:
      self.plotsPanel.setSubjects(self.patientsPanel.getSelectedPatientsNames(), self.patientsPanel.getSelectedPatients())




if __name__ == "__main__":
  app = QtGui.QApplication(sys.argv)
  axon.initializeProcesses()
  from brainvisa.data.readdiskitem import ReadDiskItem
  from brainvisa.data.writediskitem import WriteDiskItem
  QtCore.pyqtRemoveInputHook()
  window = GroupDisplayWidget()
  window.show()
  sys.exit(app.exec_())
