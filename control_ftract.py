#!/usr/bin/env python
# -*- coding: utf-8 -*-
#  This software and supporting documentation are distributed by
#      Institut Federatif de Recherche 49
#      CEA/NeuroSpin, Batiment 145,
#      91191 Gif-sur-Yvette cedex
#      France
#
# This software is governed by the CeCILL license version 2 under
# French law and abiding by the rules of distribution of free software.
# You can  use, modify and/or redistribute the software under the 
# terms of the CeCILL license version 2 as circulated by CEA, CNRS
# and INRIA at the following URL "http://www.cecill.info". 
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license version 2 and that you accept its terms.
# -*- coding: utf-8 -*-
import anatomist.direct.api as anatomist
from soma import aims
from soma.aims import colormaphints
import sys, os, math
import pdb
from PyQt4 import QtCore, QtGui

sys.path.insert( 0, '.' )


# determine wheter we are using Qt4 or Qt3, and hack a little bit accordingly
# the boolean qt4 gloabl variable will tell it for later usage
#qt4 = False
#if sys.modules.has_key( 'PyQt4'):
  #qt4 = True
  #
  #qt = QtGui
  #from PyQt4.uic import loadUi
#else:
  #import qt, qtui
  #loadUi = qtui.QWidgetFactory.create

## do we have to run QApplication ?
#if qt.qApp.startingUp():
  #qapp = qt.QApplication( sys.argv )
  #runqt = True
#else:
  #runqt = False


class fTract_Action( anatomist.cpp.Action ):
  def name( self ):
    return 'fTract_Action'

  def resetRadius( self ):
    print 'reset radius to 1'
    s.setRadius( 1. )

  def startMoveRadius( self, x, y, globx, globy ):
    print 'start move radius', x, y
    self._initial = ( x, y )
    self._radius = s.radius()

  def endMoveRadius( self, x, y, globx, globy ):
    print 'end move radius', x, y

  def moveRadius( self, x, y, globx, globy ):
    print 'move radius', x, y
    s.setRadius( math.exp( 0.01 * ( self._initial[1] - y ) ) * self._radius )

  def takePolygon( self, x, y, globx, globy ):
    #print 'takePolygon', x, y

    print 'coucou', x,y

    #w = self.view().aWindow()
    #obj = w.objectAtCursorPosition( x, y )
    ##print 'object:', obj
    #if obj is not None:
      #print 'object:', obj, obj.name()
      #poly = w.polygonAtCursorPosition( x, y, obj )
      ##print 'polygon:', poly
      #mesh = anatomist.cpp.AObjectConverter.aims( obj )
      ##print 'mesh:', mesh
      #ppoly = mesh.polygon()[poly]
      #vert = mesh.vertex()
      ##print ppoly[0], ppoly[1], ppoly[2]
      ##print vert[ppoly[0]], vert[ppoly[1]], vert[ppoly[2]]
      #global selmesh, selanamesh
      #if selmesh is None:
        #selmesh = aims.AimsSurfaceTriangle()
      #selmesh.vertex().assign( [ vert[ppoly[0]], vert[ppoly[1]], vert[ppoly[2]] ] )
      #selmesh.polygon().assign( [ aims.AimsVector_U32_3( 0, 1, 2 ) ] )
      #if selanamesh is None:
        #selanamesh = anatomist.cpp.AObjectConverter.anatomist( selmesh )
        #a = anatomist.Anatomist()
        #a.execute( 'SetMaterial', objects=[selanamesh], diffuse=[0,0,1.,1.] )
        #a.execute( 'AddObject', objects=[selanamesh], windows=[w] )
      #selanamesh.setChanged()
      #selanamesh.notifyObservers()

class ftract_control( anatomist.cpp.Control ):
  def __init__( self, prio = 25 ):
    anatomist.cpp.Control.__init__( self, prio, 'ftract_control' )

  def eventAutoSubscription( self, pool ):

    key = QtCore.Qt
    NoModifier = key.NoModifier
    ShiftModifier = key.ShiftModifier
    ControlModifier = key.ControlModifier
    AltModifier = key.AltModifier

    pdb.set_trace()
    self.keyPressEventSubscribe( key.Key_C, NoModifier, 
                                 pool.action( 'fTract_Action' ).resetRadius )
    self.mouseLongEventSubscribe( key.LeftButton, NoModifier,
      pool.action( 'fTract_Action' ).startMoveRadius,
      pool.action( 'fTract_Action' ).moveRadius,
      pool.action( 'fTract_Action' ).endMoveRadius,
      False )
    self.mousePressButtonEventSubscribe(key.RightButton, NoModifier,
      pool.action( 'fTract_Action' ).takePolygon )

#a = anatomist.Anatomist()
#pix = qt.QPixmap('/home/b67-belledone/Desktop/epilepsie-manik/Logo-F-TRACT.xpm' )
#anatomist.cpp.IconDictionary.instance().addIcon( 'ftract_control',  pix )
#ad = anatomist.cpp.ActionDictionary.instance()
#ad.addAction( 'fTract_Action', lambda: fTract_Action() )
#cd = anatomist.cpp.ControlDictionary.instance()
#cd.addControl( 'ftract_control', lambda: ftract_control(), 25 )
#cm = anatomist.cpp.ControlManager.instance()
#cm.addControl( 'QAGLWidget3D', '', 'ftract_control' )

#s = aims.SurfaceGenerator.sphere(aims.Point3df(10, 10, 10), 4, 32)
##a.registerObject( s )
#aw = a.createWindow( '3D' )
#a.addObjects( [ s ], [ aw ] )

## run Qt
#if runqt:
  #if qt4:
    #qapp.exec_()
  #else:
    #qapp.exec_loop()