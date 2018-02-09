#! /usr/bin/env python
# -*- coding: utf-8 -*-

from anatomist import cpp as ana
from soma.qt_gui.qt_backend import QtGui, QtCore, uic
import pdb

class ControlFtract(ana.Control):
    
    def __init__(self, prio = 25):
        ana.Control.__init__(self,prio,'ControlFtract')
        
    def eventAutoSubscription(self,pool): #declare actions dispo et comment
        
        key = QtCore.Qt
        NoModifier = key.NoModifier
        ShiftModifier = key.ShiftModifier
        ControlModifier = key.ControlModifier
        AltModifier = key.AltModifier
        
        self.mousePressButtonEventSubscribe(key.LeftButton,NoModifier,pool.action('fTract_Action').show)
        
        
class StimulateResults(ana.Action):
    
    def show(self,x,y,globx,globy):
        
        print "Toto"        
        ww = self.view().aWindow()
        obj = ww.objectAtCursorPosition( x, y )
        poly = ww.polygonAtCursorPosition( x, y, obj )
        pdb.set_trace()
        if isinstance(obj,ana.MObject): #obj.size() == 2:
            vertexselected = obj[0].surface().polygon()[poly]
            #obj[1].texture1d()[0][2289]
            #obj[1].texture1d()[2289]
            #obj[1].texture1d()[0]
        #http://brainvisa.info/pyanatomist-4.5/sphinx/pyanatomist_examples.html
        
        otherwindow = [x for x in ana.Anatomist().getWindows() if x is not ww and x.parent() == ww.parent()]