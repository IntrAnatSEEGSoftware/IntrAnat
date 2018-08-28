# -*- coding: utf-8 -*-
#  This software and supporting documentation are distributed by
#      Institut des Neurosciences de Grenoble
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

from brainvisa.processes import *
#import shfjGlobals
from brainvisa import anatomist
import glob
import brainvisa.registration as registration
#from editor import ElectrodeEditorDialog

name = 'Anatomist Show Electrode Model'
userLevel = 0
roles = ('viewer',)

def validation():
    anatomist.validation()

signature = Signature(
    'model', ReadDiskItem( 'Electrode Model', 'Electrode Model format' ),
)


def initialization( self ):
  pass

def execution( self, context ):
  a = anatomist.Anatomist()
  #elec = ElectrodeEditorDialog(a)
  #elec.open(self.model.fullPath())
  #meshes = elecDialog.getAnatomistObjects()
  w = a.createWindow('Axial')
  #a.addObjects(meshes, [w,])
  #return (w, elec, meshes)
  return (w,)

