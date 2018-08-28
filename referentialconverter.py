# -*- coding: utf-8 -*-
from soma import aims
from soma.aims import apctools
from numpy import *
from brainvisa.data.writediskitem import ReadDiskItem
import pdb

# Conversion de coordonnées d'un référentiel à un autre.
# Les transformations linéaires (matrices), le référentiel AC-PC,
#  le référentiel pseudo-Talairach de BrainVisa et le référentiel
#  Goetz peuvent être utilisés.
#
# (c) Inserm U836 2012-2014 - Manik Bhattacharjee
#
# License GNU GPL v3
class ReferentialConverter:
  """ This class allows to load linear transformations (matrices) from a referential to the "real" anatomist coordinates (image-based)
   All declared referentials/transformations can then be used to convert coordinates from one referential to another.
   Functions are provided to define the Goetz mesencephalon/PPN referential, the AC-PC referential (from the APC file in BrainVisa, brainvisa database access must be initialized by the caller),
  the BEN referential (Benabid normalization, AC-PC with height normalization using Thalamus height and AC-PC length, no Y normalization). Any matrix transformation can be loaded."""
  def __init__(self):
    self.availableRefs = {}
    self.withMatrixFromReal = {}
    self.withMatrixToReal = {}

  #def saveToFile
    #self.availableRefs, self.Hthal, self.Ac, self.Pc, self.Ih, self.withMatrixFromReal, self.withMatrixToReal
    #self.Oppn,self.XDppn, self.XGppn, self.Yppn, self.Zppn

  def availableReferentials(self):
    return self.availableRefs

  def isRefAvailable(self, ref):
    return self.availableRefs.has_key(ref)


  ############################ Generic transformation matrices ####################
  def setAnyRef2AnyRef(self, refFrom, refTo, transf, inverse = False):
      if not self.isRefAvailable(refFrom):
          if self.isRefAvailable(refTo): # If only the refTo is already known, just switch inputs
              tmpRef = refTo
              refTo = refFrom
              refFrom = tmpRef
              inverse = not inverse
          else:
              raise Exception("%s and %s are not known by referentialConverter ! Cannot set transform between them."%(str(refFrom), str(refTo)))
      transf = self.anyMatrix2AimsMotion(transf)
      if inverse:
          transf = transf.inverse()
      # Ok, now refFrom is a known referential, we just have to setup a T2 matrix such that fromRef --T1--> MainRef --T2-->toRef  = transf
      # T2 = transf*inv(T1)
      if refFrom in self.withMatrixToReal:
          t1inv = self.withMatrixToReal[refFrom].inverse()
      elif refFrom in self.withMatrixFromReal:
            t1inv = self.withMatrixFromReal[refFrom]
      else:
            raise Exception("Referential %s is not defined with a matrix !"%refFrom)
      t2 = transf*t1inv
      self.setTransformMatrix(refTo, t2, t2.inv())

  def anyMatrix2AimsMotion(self, mat):
      if type(mat) == aims.Motion:
          return mat
      elif type(mat) == list:
          return aims.Motion(mat)
      elif mat.__dict__.has_key('getInfos'):
          infos = anaTransf.getInfos()
          rot = infos['rotation_matrix']
          trans = infos['translation']
          return aims.Motion(rot[:3]+[trans[0]]+rot[3:6]+[trans[1]]+rot[6:]+[trans[2]]+[0,0,0,1])
      else:
          raise Exception("Unknown transformation type : %s"+type(mat))


  def setAnatomistTransform(self, refName, anaTransf, toRef=True):
    """ Adds a referential to the converter from an Anatomist Transformation object
      :arg refName name of the new referential
      :arg anaTransf anatomist transformation object (anatomist.Transformation)
      :arg toRef if true, the provided transformation goes from the native referential to refName referential. If false, the opposite.
    """
    infos = anaTransf.getInfos()
    rot = infos['rotation_matrix']
    trans = infos['translation']
    m = aims.Motion(rot[:3]+[trans[0]]+rot[3:6]+[trans[1]]+rot[6:]+[trans[2]]+[0,0,0,1])
    #pdb.set_trace()
    if toRef:
      self.setTransformMatrix(refName, m.inverse(), m)
    else:
      self.setTransformMatrix(refName, m, m.inverse())


  def setRawTransformMatrix(self, refName, matrix, toRef=True):
    """ Adds a referential to the converter from a raw transformation matrix ([1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1] for identity)"""
    m = aims.Motion(matrix)
    if toRef:
      self.setTransformMatrix(refName, m.inverse(), m)
    else:
      self.setTransformMatrix(refName, m, m.inverse())



  def setTransformMatrix(self, refName, matrixFromReal, matrixToReal):
    """ Adds a referential to the converter from an aims.Motion matrix """
    #pdb.set_trace()
    if matrixFromReal is not None:
      self.availableRefs[refName] = True
      self.withMatrixFromReal[refName] = matrixFromReal
    if matrixToReal is not None:
      self.availableRefs[refName] = True
      self.withMatrixToReal[refName] = matrixFromReal

  def applyMatrix(self, x, y, z, matrix):
    return matrix.transform( aims.Point3df(x,y,z) ).items()

  ############################ AC-PC intra-subject referential (no normalization)
  def loadACPC(self, image):
	""" Loads the .APC file linked to the provided Image Disk Item and sets AC-PC referential"""
	print "****refConv**** Loading ACPC referential"
	rdi = ReadDiskItem( 'Commissure coordinates', 'Commissure coordinates' )
	apcfile = rdi.findValue( image )
	points = apctools.apcRead(apcfile.fullPath())
	self.setACPC(points['acmm'], points['pcmm'], points['ihmm'])

  def setACPC(self, ac, pc, ih):
    """ Set AC, PC and InterHemispheric points to define the AC-PC referentials """
    self.Ac = ac
    self.Pc = pc
    self.Ih = ih
    self.availableRefs['AC-PC'] = True

  ####### AC-PC patient referential : no normalisation, origin at AC, X left to right, Y towards front, Z up
  def real2AcPc(self, x, y, z):
    """Convert 'real' native coordinates to AC-PC with origin at AC coordinates"""
    # Compute axes
    ac = array([self.Ac,])
    pc = array([self.Pc,])
    ih = array([self.Ih,])
    # Y axis is AC-PC axis towards the front
    nacpc = linalg.norm(ac - pc)
    taly = (ac-pc)/nacpc

    # X axis is perpendicular to the plane containing PC, AC and IH, so cross(Y axis, pseudoZ axis)
    talx = cross(taly, ih - pc)# x = y^z
    talx /= linalg.norm(talx)
    talz = cross(talx, taly)
    posT = array([[x, y, z]]) - ac
    return [dot(posT[0], talx[0]), dot(posT[0], taly[0]), dot(posT[0], talz[0])]


  def AcPc2Real(self, x, y, z):
    """Convert AC-PC with origin at AC coordinates to 'real' native coordinates """
    # Compute axes
    ac = array([self.Ac,])
    pc = array([self.Pc,])
    ih = array([self.Ih,])
    # Y axis is AC-PC axis towards the front
    nacpc = linalg.norm(ac - pc)
    taly = (ac-pc)/nacpc

    # X axis is perpendicular to the plane containing PC, AC and IH, so cross(Y axis, pseudoZ axis)
    talx = cross(taly, ih - pc)# x = y^z
    talx /= linalg.norm(talx)

    talz = cross(talx, taly)

    pos = -x*talx + y*taly + z*talz + ac
    print "     "+repr(pos.tolist()[0])
    return pos.tolist()[0]

  ############################### Anatomist Talairach Referential : bounding box normalizatiokn
  def loadTalairach(self, image):
    """From a DiskItem of an image, loads the pseudo-talairach transformation from BrainVisa database"""
    rdi = ReadDiskItem( 'Transform Raw T1 MRI to Talairach-AC/PC-Anatomist', 'Transformation matrix' )
    di = rdi.findValue(image)
    trans = aims.read( di.fullPath() ) # read a Motion object
    self.setTransformMatrix('Talairach', trans, trans.inverse())


  ##############################" Goetz Referential for the PPN region (pons/mesencephalon)
      # Conversion des coordonnées réelles vers les coordonnées Goetz : PROJECTION PERPENDICULAIRE AUX AXES (la projection utilisée pour les publis)
  def setGoetz(self, Oppn, XDppn, XGppn, Yppn, Zppn):
    """ Best referential for PPN studies according to Laurent Goetz ; PPNorthoGoetz -> orthoprojection on axes
        PPN parallel Goetz is also available (points are projected onto each axis parallel to other axes)
        Define the referential from the reference points : ponto-mesencephalon junction,
    """
    self.Oppn = Oppn
    self.XDppn = XDppn
    self.XGppn = XGppn
    self.Yppn = Yppn
    self.Zppn = Zppn
    self.availableRefs['PPNorthoGoetz'] = True
    self.availableRefs['PPNparaGoetz'] = True

  def real2Gorth(self, x, y, z):
    """ Convert native coordinates to normalized orthogonal Goetz referential"""
    u = array([x - self.Oppn[0], y - self.Oppn[1], z - self.Oppn[2]])
    # Vecteurs de la base
    oxg = array([self.XGppn[0] - self.Oppn[0], self.XGppn[1] - self.Oppn[1], self.XGppn[2] - self.Oppn[2]])
    oxd = array([self.XDppn[0] - self.Oppn[0], self.XDppn[1] - self.Oppn[1], self.XDppn[2] - self.Oppn[2]])
    oy = array([self.Yppn[0] - self.Oppn[0], self.Yppn[1] - self.Oppn[1], self.Yppn[2] - self.Oppn[2]])
    oz = array([self.Zppn[0] - self.Oppn[0], self.Zppn[1] - self.Oppn[1], self.Zppn[2] - self.Oppn[2]])
    # Normes des vecteurs de la base
    nxg = linalg.norm(oxg)
    nxd = linalg.norm(oxd)
    ny = linalg.norm(oy)
    nz = linalg.norm(oz)
    # vecteurs unitaires
    oxg /= nxg
    oxd /= nxd
    oy /= ny
    oz /= nz
    # Produit scalaire entre vecteurs unitaires des axes et la position u pour avoir les coordonnées orthogonale,
    #  division par la norme des vecteurs pour rapporter aux longueurs de reference
    if not self.isRightSideGoetz(u[0], u[1], u[2]): # Selon le coté
      return [dot(u, oxg)/nxg, dot(u, oy)/ny, dot(u, oz)/nz, 1]
    else: # Right side -> negative X
      return [dot(u, oxd)/nxd, dot(u, oy)/ny, dot(u, oz)/nz, -1]

  # Conversion des coordonnées réelles vers les coordonnées Goetz : PROJECTION PARALLELE AUX AXES : résolution d'un système d'équations (m*X+n*Y+o*Z = U)
  def real2G(self, x, y, z):
    """ Convert native coordinates to normalized parallel Goetz referential"""
    u = array([[x - self.Oppn[0], y - self.Oppn[1], z - self.Oppn[2]]])
    # Vecteurs de la base
    oxg = array([[self.XGppn[0] - self.Oppn[0], self.XGppn[1] - self.Oppn[1], self.XGppn[2] - self.Oppn[2]]])
    oxd = array([[self.XDppn[0] - self.Oppn[0], self.XDppn[1] - self.Oppn[1], self.XDppn[2] - self.Oppn[2]]])
    oy = array([[self.Yppn[0] - self.Oppn[0], self.Yppn[1] - self.Oppn[1], self.Yppn[2] - self.Oppn[2]]])
    oz = array([[self.Zppn[0] - self.Oppn[0], self.Zppn[1] - self.Oppn[1], self.Zppn[2] - self.Oppn[2]]])
    # Normes des vecteurs de la base
    nxg = linalg.norm(oxg)
    nxd = linalg.norm(oxd)
    ny = linalg.norm(oy)
    nz = linalg.norm(oz)
    # vecteurs unitaires
    oxg /= nxg
    oxd /= nxd
    oy /= ny
    oz /= nz

    if not self.isRightSideGoetz(u[0][0], u[0][1], u[0][2]): # Left side
      matrice = vstack((oxg, oy, oz))
      #print "MATRICE : " + repr(matrice)
      #print "SOLUTION : " + repr(linalg.solve(matrice, u[0]))
      #print "CHECK U : " + repr(dot(matrice, linalg.solve(matrice, u[0])))
      return (linalg.solve(matrice, u[0])/array([nxg, ny, nz])).tolist() + [1,]
    else:
      matrice = vstack((oxd, oy, oz))
      return (linalg.solve(matrice, u[0])/array([nxd, ny, nz])).tolist()+[-1,]



  # Conversion des coordonnées Goetz vers les coordonnées réelles : PROJECTION PARALLELE AUX AXES u=x_u*X+y_u*Y+z_u*Z
  def g2Real(self, xg, yg, zg, side):
    """ Convert normalized parallel Goetz referential to native coordinates"""
    if side == -1: # Right side
      return [self.Oppn[0] + xg*(self.XDppn[0] - self.Oppn[0]) + yg*(self.Yppn[0] - self.Oppn[0]) + zg*(self.Zppn[0] - self.Oppn[0]),\
	      self.Oppn[1] + xg*(self.XDppn[1] - self.Oppn[1]) + yg*(self.Yppn[1] - self.Oppn[1]) + zg*(self.Zppn[1] - self.Oppn[1]),\
	      self.Oppn[2] + xg*(self.XDppn[2] - self.Oppn[2]) + yg*(self.Yppn[2] - self.Oppn[2]) + zg*(self.Zppn[2] - self.Oppn[2])]
    elif side == 1: # left side
      return [self.Oppn[0] + xg*(self.XGppn[0] - self.Oppn[0]) + yg*(self.Yppn[0] - self.Oppn[0]) + zg*(self.Zppn[0] - self.Oppn[0]),\
	      self.Oppn[1] + xg*(self.XGppn[1] - self.Oppn[1]) + yg*(self.Yppn[1] - self.Oppn[1]) + zg*(self.Zppn[1] - self.Oppn[1]),\
	      self.Oppn[2] + xg*(self.XGppn[2] - self.Oppn[2]) + yg*(self.Yppn[2] - self.Oppn[2]) + zg*(self.Zppn[2] - self.Oppn[2])]
    else:
      print "ERROR : side value in invalid in g2Real : "+repr(side)

  # Est-on à droite de la ligne médiane ?
  def isRightSideGoetz(self, x, y, z):
    """ Check which side of the brain we are on, from PPN Goetz referential landmarks """
    oy = array([self.Yppn[0] - self.Oppn[0], self.Yppn[1] - self.Oppn[1], self.Yppn[2] - self.Oppn[2]])
    oz = array([self.Zppn[0] - self.Oppn[0], self.Zppn[1] - self.Oppn[1], self.Zppn[2] - self.Oppn[2]])
    pseudoX = cross(oy, oz)
    return dot(pseudoX, array([x,y,z])) > 0

  # Conversion des coordonnées Goetz vers les coordonnées réelles : PROJECTION PERPENDICULAIRE AUX AXES
  #  Résolution du système d'inconnu u tel que u.X = x_u, u.Y = y_u et u.Z = z_u ('.' est le produit scalaire)
  # OLD VERSION (slightly slower than the newer one using solve)
  def g2RealOrth2(self, xg, yg, zg, side):
    """ Convert normalized orthogonal Goetz referential to native coordinates (old, slow version)"""
    u = array([[xg, yg, zg]])
    # Vecteurs de la base
    oxg = array([[self.XGppn[0] - self.Oppn[0], self.XGppn[1] - self.Oppn[1], self.XGppn[2] - self.Oppn[2]]])
    oxd = array([[self.XDppn[0] - self.Oppn[0], self.XDppn[1] - self.Oppn[1], self.XDppn[2] - self.Oppn[2]]])
    oy = array([[self.Yppn[0] - self.Oppn[0], self.Yppn[1] - self.Oppn[1], self.Yppn[2] - self.Oppn[2]]])
    oz = array([[self.Zppn[0] - self.Oppn[0], self.Zppn[1] - self.Oppn[1], self.Zppn[2] - self.Oppn[2]]])
    # Normes des vecteurs de la base
    nxg = linalg.norm(oxg)
    nxd = linalg.norm(oxd)
    ny = linalg.norm(oy)
    nz = linalg.norm(oz)
    # vecteurs unitaires
    oxg /= nxg
    oxd /= nxd
    oy /= ny
    oz /= nz
    # MATLAB : ([norm(OXg), norm(OY), norm(OZ)].*[xg, yg, zg])*inv([OXg'./norm(OXg), OY'./norm(OY), OZ'./norm(OZ)]) + Oppn
    if side == -1: # Right side
      matrice = hstack((oxd.T, oy.T, oz.T))
      result = dot(array([[nxd, ny, nz]])*u, linalg.inv(matrice)) + array(self.Oppn)

    elif side == 1: # left side
      matrice = hstack((oxg.T, oy.T, oz.T))
      result = dot(array([[nxg, ny, nz]])*u, linalg.inv(matrice)) + array(self.Oppn)
    else:
      print "ERROR : side value is invalid in g2RealOrth : "+repr(side)

    return result.tolist()[0]


  # Conversion des coordonnées Goetz vers les coordonnées réelles : PROJECTION PERPENDICULAIRE AUX AXES
  #  Résolution du système d'inconnu u tel que u.X = x_u, u.Y = y_u et u.Z = z_u ('.' est le produit scalaire)
  def g2RealOrth(self, xg, yg, zg, side):
    """ Convert normalized orthogonal Goetz referential to native coordinates"""
    u = array([[xg, yg, zg]])
    # Vecteurs de la base
    oxg = array([[self.XGppn[0] - self.Oppn[0], self.XGppn[1] - self.Oppn[1], self.XGppn[2] - self.Oppn[2]]])
    oxd = array([[self.XDppn[0] - self.Oppn[0], self.XDppn[1] - self.Oppn[1], self.XDppn[2] - self.Oppn[2]]])
    oy = array([[self.Yppn[0] - self.Oppn[0], self.Yppn[1] - self.Oppn[1], self.Yppn[2] - self.Oppn[2]]])
    oz = array([[self.Zppn[0] - self.Oppn[0], self.Zppn[1] - self.Oppn[1], self.Zppn[2] - self.Oppn[2]]])
    # Normes des vecteurs de la base
    nxg = linalg.norm(oxg)
    nxd = linalg.norm(oxd)
    ny = linalg.norm(oy)
    nz = linalg.norm(oz)
    # vecteurs unitaires
    oxg /= nxg
    oxd /= nxd
    oy /= ny
    oz /= nz
    # MATLAB : ([norm(OXg), norm(OY), norm(OZ)].*[xg, yg, zg])*inv([OXg'./norm(OXg), OY'./norm(OY), OZ'./norm(OZ)]) + Oppn
    if side == -1: # Right side
      matrice = vstack((oxd, oy, oz))
      result = linalg.solve(matrice, array([nxd, ny, nz])*u[0]) + array(self.Oppn)
    elif side == 1: # left side
      matrice = vstack((oxg, oy, oz))
      result = linalg.solve(matrice, array([nxg, ny, nz])*u[0]) + array(self.Oppn)
    else:
      print "ERROR : side value is invalid in g2RealOrth : "+repr(side)
    return result.tolist()


  ############################### Ben's coordinates Benabid AC-PC with thalamus normalization
  ###### NORMALIZED BEN'S COORDINATES (normalized along Y and Z axis)
  def setBens(self, ac, pc, ih, hthal):
    """Set reference points for Ben's normalization : Anterior Commissure,
        Posterior Commissure, InterHemispheric point (above AC-PC),
        Thalamus height (top of the thalamus, bottom of lateral ventricle, up from AC-PC axial slice).
        Origin is at PC, no normalization on X axis,
        AC-PC normalization on Y (y/norm(ac-pc)*12.0), Thalamus Height normalization on Z (z/norm(thal-pc)*8.0)
    """
    self.setACPC(ac,pc,ih)
    self.Hthal = hthal
    self.availableRefs['Bens'] = True

  def real2Bens(self, x, y, z):
    """Convert native coordinates to normalized Ben's coordinates (IMPLEMENTATION NEVER USED/TESTED)"""
    # Compute axes :
    ac = array([self.Ac,])
    pc = array([self.Pc,])
    ih = array([self.Ih,])
    hthal = array([self.Hthal,])
    u = array([[x, y, z]])
    # Y axis is AC-PC axis towards the front
    nacpc = linalg.norm(ac - pc)
    thaly = (ac-pc)/nacpc
    # X axis is perpendicular to the plane containing PC, AC and IH, so cross(Y axis, pseudoZ axis)
    thalx = cross(thaly, ih - pc)# x = y^z
    thalx /= linalg.norm(thalx)
    thalz = cross(thalx, thaly)
    thalHeight = dot((hthal - pc)[0], thalz[0])
    return [dot(u, thalx), dot(u, thaly)*12.0/nacpc, dot(u, thalz)*8.0/thalHeight]

  def bens2Real(self, xb, yb, zb):
    """Convert normalized Ben's coordinates to native coordinates"""
    # x is not normalized, so the value is just the distance from the AC-PC line, with X<0 for the right hemisphere
    # Compute axes :
    ac = array([self.Ac,])
    pc = array([self.Pc,])
    ih = array([self.Ih,])
    hthal = array([self.Hthal,])

    # Y axis is AC-PC axis towards the front
    nacpc = linalg.norm(ac - pc)
    thaly = (ac-pc)/nacpc

    # X axis is perpendicular to the plane containing PC, AC and IH, so cross(Y axis, pseudoZ axis)
    thalx = cross(thaly, ih - pc)# x = y^z
    thalx /= linalg.norm(thalx)
    thalz = cross(thalx, thaly)
    thalHeight = dot((hthal - pc)[0], thalz[0])

    pos = (-xb * thalx)   +   (yb * nacpc / 12.0) * thaly   +   (zb * thalHeight / 8.0) * thalz + pc # -xb because the sign should ne changed for a direct referential

    return pos.tolist()[0]

  #############################" Generic Conversion ######################################
  # UNIVERSAL REFERENTIAL CONVERTER TO REAL MRI COORDS
  def anyRef2Real(self, coords, referential):
    """Convert coords [x,y,z] or [x,y,z,side] to native coordinates from any defined referential"""
    if not self.isRefAvailable(referential):
      print "Referential %s not available"%referential
      return None
    withSide = {'PPNparaGoetz':self.g2Real, 'PPNorthoGoetz':self.g2RealOrth}
    noSide = {'Bens':self.bens2Real, 'AC-PC':self.AcPc2Real,'real':lambda x,y,z:[x,y,z]}
    if referential in self.withMatrixToReal:
      return self.applyMatrix(coords[0], coords[1], coords[2], self.withMatrixToReal[referential])
    if referential in withSide and size(coords)>3:
	    return withSide[referential](coords[0], coords[1], coords[2], coords[3])
    elif referential in noSide:
	    return noSide[referential](coords[0], coords[1], coords[2])
    else:
	    print "anyRef2Real : no such referential or invalid coords : "+repr(referential)+" -> "+repr(coords)
	    return None

  # UNIVERSAL REFERENTIAL CONVERTER from REAL MRI COORDS
  def real2AnyRef(self, coords, referential):
    """Convert coords [x,y,z] or [x,y,z,side] to any defined referential from native coordinates"""
    if not self.isRefAvailable(referential):
      print "Referential %s not available"%referential
      return None
    withSide = {}
    noSide = {'Bens':self.real2Bens, 'AC-PC':self.real2AcPc,'real':lambda x,y,z:[x,y,z], 'PPNparaGoetz':self.real2G, 'PPNorthoGoetz':self.real2Gorth}
    #pdb.set_trace()
    if referential in self.withMatrixFromReal:
      return self.applyMatrix(coords[0], coords[1], coords[2], self.withMatrixFromReal[referential])
    elif referential in withSide and size(coords)>3:
	    return withSide[referential](coords[0], coords[1], coords[2], coords[3])
    elif referential in noSide:
	    return noSide[referential](coords[0], coords[1], coords[2])
    else:
	    print "anyRef2Real : no such referential or invalid coords : "+repr(referential)+" -> "+repr(coords)
	    return None

  # Converts coords from any ref to any other referential
  def anyRef2AnyRef(self, coords, referentialFrom, referentialTo):
	"""Convert coords [x,y,z] or [x,y,z,side] from any referential to any other referential"""
	reals = self.anyRef2Real(coords, referentialFrom)
	if reals is None:
	      print "Cannot convert to real !"
	      return None
	if size(coords) > 3: # There is side information in coords[3]
	  reals = reals[:3]+coords[3:] # Keep the side information in case it is needed by the destination referential
	return self.real2AnyRef(reals,referentialTo)
