Technical Documentation
***********************

.. contents::
   :depth: 3
..

IntrAnat

Technical documentation

Written : January 2014
Last updated : March 2021

List of software
===================

**IntrAnat** is a collection of software for Stereo-ElectroEncephaloGraphy (SEEG),
also called Intracerebral EEG and for Deep Brain Stimulation (DBS).

Main parts:

-  **ImageImport**: importing data into the IntrAnat/Brainvisa database

-  **locateElectrodes**: Display 3D imaging data and locate electrodes on
   images. Multimodal display (including atlases) and export electrode coordinates and parcels.

-  **groupDisplay**: Display electrode locations of multiple subjects in MNI space (**deprecated**)

-  **electrodeEditor**: creating a new model of electrode for *locateElectrodes*
   CombiView: visualisation de données de groupes de sujets

Software dependencies:
==========================

IntrAnat was developed and used mostly on Linux.

It should run on Windows and Mac OS X but bugs may be platform-specific.

Il uses multiple pieces of software :

-  **BrainVisa 4.6.1** (Started on 4.2) that should already be configured with a database.
   BrainVisa manages the database used by IntrAnat, provides many tools used by IntrAnat, notably
   Anatomist which displays 3D volumes and meshes.
-  **Qt4 and PyQt4** (included in BrainVisa)
-  **Python 2.6** (included in BrainVisa), and **Python 3**
-  **SPM12** (and **matlab**) fo intrasubject coregistration, and normalization to MNI referential.
-  **FreeSurfer** atlases can be imported and used to tag electrode contacts with the name of brain regions
-  **dcmtk** was used to import from DICOM servers but function was removed

| For Brainvisa 5 (released march 2021) code will be ported to Qt5 (see relevant git branch) e.g. rewriting
| *self.connect(self.electrodeList,
  QtCore.SIGNAL("currentRowChanged(int)"), self.electrodeSelect)*
| as
| *self.electrodeList.currentRowChanged[int].connect(self.electrodeSelect)*

**Matlab** is called through command line, so version changes should not be a problem
as long as Matlab binary is in the PATH.

For **SPM12**, functions called by IntrAnat were broken between versions in the past.
To update affected code, search for spm\_\* strings that contain spm calls in the code.

**Example** : spm_coregister ::

    spm_coregister = "try,VF=spm_vol(%s);VG=spm_vol(%s);x = spm_coreg(VF,VG);
    trm = spm_matrix(x(:)');trm = [trm(1:3,4)';trm(1:3,1:3)];
    dlmwrite(%s,trm, 'delimiter',' ','precision',16);
    catch, disp 'AN ERROR OCCURED'; end;quit;"

This string contains the matlab code using spm functions. It is then called
in *spmCoregister(self, image,target)*

How SPM coregister works : each Nifti image has an internal transform that goes to a
« scanner-based » referential. Sometimes this transform is not present in the file header
and SPM behavior is not always predictable (depends on versions and file header).

In principle, spm_coreg output transformation does not include scanner-based transforms.
So to get from an image to another, the transform should be::

    trm = VF.mat\spm_matrix(x(:)')*VG.mat

In IntrAnat, the matrix from one scanner-based ref to the other is directly stored and used by Anatomist
(Brainvisa's viewer)

Integration into BrainVisa
==========================

L'intégration comprend deux parties :

-  A toolbox « epilepsy » in BrainVisa **deprecated**
-  software using BrainVisa but that are called stand-alone (not from Brainvisa GUI)

Epilepsy toolbox
----------------

brainvisa-4.6.1/brainvisa/toolboxes/epilepsy/

This directory contains definitions of file types (electrode models, implantations
and so on), a description of how this data is organized in the database (inside
each subject directory), as well as a few « BrainVisa processes » that can be called from Brainvisa GUI

Definitions for the database
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Data stored for IntrAnat in BrainVisa database is defined in two steps:

-  in the **types** directory: epilepsy.py defines file formats (name
   and extension) and data types. Example for SEEG recordings using
   Micromed's TRC file format:
   Declare file format and extension::

      Format( 'EEG TRC format', 'f|*.trc' )

   Declare a data type called 'SEEG recording', which is not a subtype (so
   base type is 'Any Type')
   and which can be in multiple file formats::

      FileType( 'SEEG recording', 'Any Type', ['EEG TRC format', 'Elan EEG format'])#'ImaGIN matlab format'

   Now let's declare a type 'Raw SEEG recording', subtype of
   'SEEG recording', and which is stored in 'EEG TRC format'.::

      FileType( 'Raw SEEG recording', 'SEEG recording', 'EEG TRC format' )

   Various examples are visible in the main BrainVisa hierarchy
   (brainvisa-4.6.0/brainvisa/types/) but also in other toolboxes.
-  In *hierarchies* directory, multiple subdirectories, one for each version
   of the database structure.

   Currently, user database works with hierarchy brainvisa-3.2.0,
   in the same way as it is organized in brainvisa-4.6.0/brainvisa/hierarchies/.
   Directory *shared* corresponds to the internal database of Brainvisa (e.g. to
   store its templates, available for all studies).

   In brainvisa-4.6.0/brainvisa/toolboxes/epilepsy/hierarchies/brainvisa-3.2.0/,
   a few files are used to declare where IntrAnat-specific files will be inserted
   into the standard BrainVisa hierarchy.

   Example : in *images.py*, storing CT images in the database:
   Create a tuple *ct_content*, containing a string description representing
   the file name in the database. This name is an expression that will match
   a real directory in the database.

   Here, {acquisition} means that the directory name will be the same as the
   'acquisition' property of the same object. E.g. if a directory is called
   'postOp-2012-11-11', Brainvisa will set for each file in this directory
   a property 'acquisition' defined by the directory name.

   This property can be reused in directories and file names inside this directory.
   In this example, we will set a default value, and set this property as optional
   (e.g. we may want to store a subject's CT scan without setting an acquisition
   name if there will be only one::

      ct_content = (
      "{acquisition}", SetDefaultAttributeValue( 'acquisition',
      default_acquisition ), SetNonMandatoryKeyAttribute( 'acquisition' ),

   Then, declare the content of the directory::

      SetContent(
   A CT file (its type is declared just like 'SEEG recording') which name
   is subject's name, dash, and acquisition name. These values are properties
   that were declared previously and which values are known (e.g. declared
   through {acquisition})::

      "<subject>-<acquisition>", SetType( 'CT' ),

   A registration directory, to store referentials and geometrical transforms
   from CT imqge to other referentials::

      'registration', SetContent(
      'CT-<subject>_<acquisition>', SetType( 'Referential of CT' ),
      'CT-<subject>_<acquisition>_TO_Talairach-ACPC', SetType( 'Transform
      CT to Talairach-AC/PC-Anatomist' ),
      'CT-<subject>_<acquisition>_TO_Talairach-MNI', SetType( 'Transform CT
      to Talairach-MNI template-SPM'),
      'CT-<subject>_<acquisition>_TO_Scanner_Based', SetType(
      'Transformation to Scanner Based Referential' ),
   Here we are adding a new transform to another image of the same subject, with
   a specific modality and acquisition: those are new properties, declared with {}::

      'CT-<subject>_<acquisition>_TO_{modalityTarget}_{acquisitionTarget}',
      SetType( 'Transform CT to another image' ),
      'CT-<subject>_<acquisition>_Scanner_Based', SetType( 'Scanner Based
      Referential' ),
      ),
      )
      )
   Finally we insert this into the existing hierarchy:
   in directory '{protocol}/{subject}' we add a 'ct' directory, which has a 'modality'
   attribute with value 'ct'. We then add its content previously declared as ct_content.::

      insert( '{protocol}/{subject}',
      'ct', SetWeakAttr( 'modality', 'ct' ),
      apply( SetContent, ct_content)
      )
   Numerous other examples are available in main Brainvisa hierarchy and its toolboxes, e.g. ::

      brainvisa-4.3.0/brainvisa/hierarchies/brainvisa-3.1.0/base.py
      brainvisa-4.3.0/brainvisa/toolboxes/morphologist/hierarchies/brainvisa-3.1.0/anatomy.py

Patches
--------

Referential and transformation patches developed for IntrAnat were
integrated into BrainVisa from 4.5.0


Tools used for development
---------------------------

**Qt Designer / Qt Creator** to create Graphical User Interfaces.
It generates a .ui file that can be directly loaded from Python (e.g.
 \__init_\_ function in ImageImportWindow)::

   from PyQt4 import uic
   self.ui = uic.loadUi("epilepsie-electrodes.ui", self)

Programmed using Python/PyQt with bindings included in BrainVisa.

Editors: Kate (on KDE), PyCharm.

IntrAnat Librairies
--------------------

Python files with useful functions:

-  *electrode.py* is used to manage electrode models and their display in Anatomist

-  *dicomutilities.py* (DEPRECATED) functions to read and access DICOM files and servers

-  *externalprocesses.py* allows to call external software (synchronous or asynchronous,
   with callback functions), especially to run matlab code.

-  *referentialconverter.py* defines an object storing multiple referentials and can convert
   coordinates from one referential to another.


ImageImport
--------------

*ImageImport* is used to add patients in the BrainVisa database,
import images (MRI, CT scans, PET scans...), de register these images
to normalize them with SPM12 (to MNI referential), and to run the main
segmentation process of BrainVisa to get brain meshes and sulci.
It can also import Freesurfer segmentations.

UI is defined in ImageImportWindow.ui

Main code is contained in the class *ImageImportWindow* in ImageImportWindow.py
and ImageImport.py is used to run it.

Software structure :

-  Buttons and other UI elements are connected to functions in \__init_\_ function
   of *ImageImportWindow* class::

     self.connect(self.ui.regSubjectCombo,
     QtCore.SIGNAL('currentIndexChanged(QString)'),
     self.setCurrentSubject)
   For PyQt5, newer syntax must be used::

     self.ui.regSubjectCombo.currentIndexChanged[str].connect(self.setCurrentSubject)

   In this example *regSubjectCombo* (a combo box of subjects in registration tab) emitting
   *currentIndexChanged* calls *self.setCurrentSubject* function.

-  code comments should help understand the code.

LocateElectrodes
----------------

Ce logiciel permet de placer les modèles d'électrodes sur les images du
patient et d'obtenir les coordonnées des plots (fichiers PTS et .txt).
Il permet également d'afficher de nombreuses données du patient
(hémisphères cérébraux extraits de l'IRM T1, scanner CT, IRM T1, T2,
pre/post implantation/post-résection...) et des modèles d'électrodes
réalistes (ou bien les plots agrandis pour faciliter la visualisation.

L'interface est définie dans le fichier epilepsie-electrodes.ui

| Structure du logiciel :
| - quelques fonctions au début du fichier permettent de gérer plus
  facilement les électrodes (ajout, déplacement...)

- une classe principale définit les fonctions attachées aux éléments de l'interface graphique.


Formats de données
==================

**Images :** IRM, CT, PET : Nifti (.nii) ou nifti compressé (.nii.gz)

**SEEG**: TRC (micromed) .eeg (ELAN), .

| **Electrodes**: .elecmodel (variantes pickle et json)
| Pour l'instant, les fichiers elecmodel sont des dictionnaires python
  sauvegardés avec la librairie pickle de Python. A partir de brainvisa
  4.4, on pourra utiliser la librairie json qui utilise un format mieux
  défini et quasi universel. Ce n'est pas encore le cas.
| L'électrode est un ensemble de cylindres, représenté par un
  dictionnaire sous la forme {'Plot1', {...}, 'Plot2':{...}, 'Element
  1':{...}}

| Les éléments sont les morceaux non actifs du modèle d'électrode, les
  plots sont les contacts de l'électrode. L'électrode est définie dans
  un repère où l'extrémité (pointe) de l'électrode se situe en 0,0,0 et
  l'électrode est alignée avec l'axe Z. Plus on s'éloigne de l'extrémité
  plus la coordonnée z augmente.
| Chaque morceau de l'électrode est à son tour défini par un
  dictionnaire :

'Plot1': {'axis': 'Axe Z',

'diameter': 0.80000000000000004,

'length': 2.0,

'position': [0.0, 0.0, 0.0],

'type': 'Plot',

| 'vector': [0.0, 0.0, 1.0]}
| Comme on le voit, il y a la direction principale du cylindre qui
  compose l'élement, son diamètre en mm, sa longueur en mm, la position
  de son extrémité, son type (ici un Plot), et le vecteur qu'on ajout à
  la position pour trouver l'autre extrémité du cylindre.
| Pour lire ces fichiers depuis python :

| import pickle
| f=open('Dixi-D08-15BM.elecdef')
| d=pickle.load(f)

**Implantations d'électrodes** : .elecimplant (variantes pickle et
json), .pts, .txt

| Comme les fichiers elecmodel, les fichiers elecimplant sont des
  dictionnaires python sauvegardés avec la librairie pickle, et sont
  destinés à passer au format json.
| Il contiennent '2mni' qui devait servir à stocker la transformation
  linéaire vers le référentiel MNI (inutilisé), 'ReferentialUuid' qui
  est l'identifiant unique du référentiel utilisé pour les coordonnées
  des électrodes, electrodes qui contient la liste des électrodes
  implantées. Il s'agit du référentiel 'natif' Anatomist de l'IRM T1
  pré-opératoire.

Chaque électrode, élément dans la liste electrodes[] est un dictionnaire
qui contient les éléments suivants : 'entry' les coordonnées du point
d'entrée, 'model' le modèle d'électrode utilisé, 'name' le nom de
l'électrode, 'target' les coordonnées de la pointe de l'électrode.

{'2mni': None,

'ReferentialUuid': '2506a605-3d18-fa50-3557-a47922440c41',

'electrodes': [{'entry': [126.82000732421875,

121.16304779052734,

135.00004577636719],

'model': 'Dixi-D08-08AM',

'name': 'A',

'target': [101.98080444335938,

120.58821868896484,

132.00001525878906]},

{'entry': …....},]

}

**Guide** **développement**
===========================

Debugger
--------

Utiliser ipython -q4thread leFichier.py

import pdb;pdb.set_trace() permet de tomber dans le debugger à la ligne
ou ceci a été inséré dans le code.

Utiliser le Database Browser dans BrainVisa pour voir si les fichiers de
la base sont reconnus ou non. On peut aussi « Mettre à jour » la base de
données pour refaire les index s'ils ont été corrompus par un mauvais
fonctionnement du logiciel.

Ajouter un nouveau type de donnée dans la base
----------------------------------------------

(ATTENTION au changement de fonctionnement de la base dans la prochaine
version majeure de BrainVisa)

Comme explicité dans la partie III.1, pour ajouter un nouveau type de
données, il faut déclarer son type (s'il est nouveau), son emplacement
et son nom dans la base de données.

Par exemple, on souhaite ajouter un fichier par patient qui liste les
structures implantées.

Ce fichier se nommera structures_NomPatient.txt et sera ajouté dans le
répertoire « implantation » du répertoire du sujet dans la base. On a
besoin :

-  d'un **format de fichier**, ici le format txt qui est déjà déclaré
   dans la base. Si ce n'est pas le cas, il suffit de déclarer le format
   de fichier dans le fichier
   brainvisa/toolboxes/epilepsy/types/epilepsy.py sous la forme
   Format( 'PTS format', 'f|*.pts' ) # *Pour un simple fichier avec
   l'extension pts*
   Format ( 'Powerpoint file', ["f|*.ppt","f|*.pptx"] ) # *Pour un
   format avec plusieurs extensions
   *\ Les formats déjà connus peuvent avoir été déclarés dans
   brainvisa-4.3.0/brainvisa/types/*.py (c'est le cas pour les types
   basiques comme .txt) ou dans d'autres toolboxes :
   brainvisa-4.3.0/brainvisa/toolboxes/*/types/*.py
-  d'un **type de fichier** « Implanted Structures », qui n'est pas un
   sous-type d'un type existant (« Right Side Implanted Structures »
   pourrait être un sous-type d' « Implanted Structures »).
   On le déclare ainsi : FileType( '**Implanted Structures'**, 'Any
   Type', 'Text file' )
   dans le fichier brainvisa/toolboxes/epilepsy/types/epilepsy.py
-  d'une déclaration dans la hiérarchie de la base de données : on peut
   l'ajouter au fichier
   brainvisa-4.3.0/brainvisa/toolboxes/epilepsy/hierarchies/brainvisa-3.1.0/electrodes.py
   Dans ce fichier, on insère le dossier « implantation » dans le
   dossier du sujet, et on déclare son contenu. Il suffit donc d'ajouter
   une ligne de contenu sur le modèle des déclarations existantes :
   "structures_<subject>", SetType('**Implanted Structures**'),
   On déclare ainsi que dans le répertoire implantation, on peut avoir
   un fichier nommé structures_nomDuPatient.txt qui est du type
   'Implanted Structures'.

Accéder aux données depuis le code Python
-----------------------------------------

Pour accéder aux données, on utilise l'API python de BrainVisa.

On va donc importer les objets nécessaires et initialiser l'accès à la
base de données.

from brainvisa import axon

axon.initializeProcesses()

from brainvisa.data.readdiskitem import ReadDiskItem

from brainvisa.data.writediskitem import WriteDiskItem

On veut trouver tous les fichiers de type 'Implanted Structures' des
patients epileptiques. On va utiliser le type de fichier et les
attributs présents dans la base de données. La base de données contient
des dossiers de protocoles qui contiennent des dossiers de patients. Le
nom de ces dossiers correspond à un attribut défini pour toutes les
données contenues dans ces répertoires. En effet le nom du dossier de
protocole est déclaré comme '{protocol}' dans la hiérarchie, ce qui crée
un attribut 'protocol' contenant le nom réel du dossier.

rdi = ReadDiskItem( 'Implanted Structures', 'Text file' ,
requiredAttributes={'protocol':'Epilepsy'} )

Si on connaît le sujet, on peut ajouter une contrainte :

rdi2 = ReadDiskItem( 'Implanted Structures', 'Text file' ,
requiredAttributes={'protocol':'Epilepsy',
'subject' :'LYONNEURO_2013_DUPj'} )

On va obtenir la liste des résultats (il peut y en avoir un seul ou
plusieurs) :

implStructures = list( rdi._findValues( {}, None, False ) )

Oui on utilise une fonction interne (_findValues) parce qu'il n'existait
que cela quand j'ai posé la question. findValue existe mais ne retourne
qu'une seule valeur. On prend ensuite le premier objet retourné, qui est
un ReadDiskItem. Cet objet permet d'obtenir le chemin réel du fichier et
aussi ses attributs, par exemple le nom du sujet.

implS = implStructures[0]

| print 'Sujet '+implS.attributes()['subject']+'. Le fichier est là :'+
  implS.fullPath()
| Si on a un type de fichier qui peut être lu par Anatomist (ce n'est
  pas le cas ici), il suffit de faire :
| from brainvisa import anatomist
| anatomist.loadObject(implS)
| Sinon on utilise le chemin implS.fullPath() pour le lire.

| De même pour trouver le chemin d'un fichier que l'on veut écrire dans
  la base. C'est un petit peu plus complexe car on doit donner toutes
  les informations nécessaires pour que BrainVisa génère un nom de
  fichier : le type et les attributs. Comment savoir sinon dans quel
  protocole, dans quel sujet doit être placé le nouveau fichier ? Pour
  faciliter les choses, on peut fournir à BrainVisa le type de fichier
  et un autre diskItem dont les attributs sont suffisants pour trouver
  le nom de fichier. Par exemple, nous avons un ReadDiskItem de type IRM
  T1 (que nous avons ouvert précédemment) et nous voulons sauver le
  fichier 'Implanted Structures' qui correspond au même sujet. Il suffit
  alors de préciser le type de fichier et le diskItemT1. Exemples :
| wdi = WriteDiskItem( 'Implanted Structures', 'Text file' )
| di = wdi.findValue({'subject':'monSujet', 'protocol':'Epilepsy'} )

di2 = wdi.findValue(diskItemT1)

print 'Fichier de sortie : ' + di.fullPath()

Référentiels et Transformations géométriques
--------------------------------------------

Les voxels des images volumiques (IRM, CT, PET...) sont localisées dans
l'espace par rapport à un **référentiel géométrique**. La plupart des
logiciels utilisent en interne un **référentiel 'natif'**, par exemple
dans le cas d'Anatomist, le référentiel natif de l'image est défini par
une position 0,0,0 au centre du voxel le plus « en haut à droite au
fond ». Ensuite les coordonnées x,y,z sont la distance en mm le long des
axes de la matrice de voxels.

Malheureusement cette convention n'est pas la même selon les logiciels.

Le format DICOM définit en général une matrice de transformation qui
permet de calculer la position des voxels de l'image par rapport à un
référentiel de la machine. La conversion de DICOM vers Nifti conserve en
général cette transformation dans le header Nifti sous le nom
« \ **scanner-based**\ ». SPM utilise cette matrice pour les coordonnées
affichées lorsque l'on fait un « display » d'une image. On peut
l'afficher dans spm en chargeant une image et en regardant la matrice
disponible dans l'objet chargé avec
a=spm_vol('LYONNEURO_2013_AAAa.nii');a.mat

On peut également utiliser la commande « AimsFileInfo
LYONNEURO_2013_AAAa.nii » qui va afficher (entre autres) la liste des
transformations et des matrices stockées dans l'en-tête de l'image
Nifti :

'referentials' : [ 'Scanner-based anatomical coordinates' ],

'transformations' : [ [ -0.999992, 0, 0, 90.9604, 0, -1, 0, 134.016, 0,
0, -1, 121.85, 0, 0, 0, 1 ] ],

J'ai appelé ce référentiel **Scanner-based referential**.

BrainVisa peut stocker dans sa base de données des référentiels et des
transformations permettant de passer d'un référentiel à un autre. Un
référentiel est (comme tout les objets d'une base de données BrainVisa)
déclaré avec un identifiant unique (UUID). Les transformations sont
déclarées dans la base comme permettant de passer d'un UUID à un autre,
ce qui permet ensuite de les trouver automatiquement (cf
TransformationManager).

**ImageImport se charge donc de stocker le référentiel 'natif', le
référentiel scanner-based et la transformation correspondante pour
chaque image au moment de l'importation**. On peut trouver ces fichiers
(.referential, .trm) dans le répertoire registration de toute image
présente dans la base.

Si on entre CA-CP, des transformations vers le référentiel de Talairach
sont calculées et déclarées à leur tour.

Ca se complique lorsque l'on recale une IRM post-opératoire (ou un
scanner CT, PET...) vers une IRM T1 pré-opératoire. SPM coregister
calcule une matrice de transformation qui permet de passer du
référentiel scanner-based de l'IRM post vers le référentiel
scanner-based de l'IRM pre.

Pour convertir des coordonnées saisies dans le référentiel natif de la
t1post vers le référentiel natif de la t1pre, il faut donc appliquer :
natif post → scanner-based post → scanner-based pre → natif pre.

**Problème fréquent** : les images qui ont subit d'autres traitements
avec d'autres logiciels peuvent avoir perdu la matrice de transformation
vers le référentiel scanner-based. Dans ce cas, il est possible
qu'IntrAnat interprète de façon incorrecte les transformations des
headers.

En effet, il peut y avoir deux matrices différentes dans un header
nifti. Si aucune des deux n'est appelée scanner-based, IntrAnat ne sait
pas laquelle choisir comme base, et le choix retenu peut ne pas être le
même que celui que ferait SPM. Si on recale cette image, IntrAnat ne
saura pas à partir de quel référentiel SPM a calculé la matrice
coregister, et le recalage ne plantera pas, mais les résultats seront
faux.

**Le référentiel MNI**: la normalisation de SPM calcule à la fois une
transformation linéaire (une matrice, comme les .trm cités plus haut) et
une transformation non linéaire (tout est stocké dans le fichier
\_sn.mat).

BrainVisa ne gère pas les transformations non linéaires, donc IntrAnat
convertit les coordonnées des électrodes vers le référentiel MNI en
faisant appel à matlab et SPM. Les transformations \_sn.mat vont du
référentiel scanner-based de l'IRM vers le référentiel MNI.

L'API de BrainVisa contient depuis peu un TransformationManager, objet
qui permet de rechercher des référentiels et des transformations liées
aux objets. Cependant les fonctions en question sont incomplètes et
insuffisantes. Les patchs écrits pour BrainVisa et le fichier
referentialconverter.py permettent de contourner une partie de ces
limitations. En particulier, les patches permettent, lorsqu'on charge un
objet dans Anatomist, de rechercher les transformations qui le lient aux
référentiels déjà présents en cherchant un « chemin » de
transformations. La version standard du code trouvera par exemple une
transformation qui lie T1pre natif à T1 post natif, mais pas le chemin
pre natif-> pre scanner-based → post scanner-based → post natif, ce que
fera le patch. Ainsi, les images pre, post, CT, PET, T2 etc sont
synchronisées dans l'affichage d'Anatomist.

ReferentialConverter permet quant à lui de déclarer un ensemble de
transformations (CA-CP, Talairach, Goetz ou toute autre transformation
linéaire...) et de convertir des coordonnées de points d'un référentiel
à un autre.

Calcul des coordonnées de plots des électrodes :

-  les coordonnées des électrodes saisies dans IntrAnat locateElectrodes
   sont enregistrées dans le référentiel natif de l'IRM T1.
-  Pour l'exportation en .pts, ces coordonnées sont converties à l'aide
   de la transformation natif → scanner-based pour les exprimer dans le
   référentiel du patient (grâce aux fonctions de ReferentialConverter).
-  Les coordonnées scanner-based sont sauvées dans un fichier
   temporaire, et un code matlab est exécuté pour convertir la
   transformation sn.mat en champ de vecteurs (y_field.nii)
-  Ce champ de vecteurs est utilisé pour convertir les coordonnées, qui
   sont sauvées dans le fichier temporaire.
-  A la fin de l'exécution de matlab, le fichier de sortie est relu
   depuis python, puis le fichier PTS est sauvé dans le référentiel MNI.

Ajouter un processus BrainVisa
------------------------------

Il suffit d'ajouter un fichier python dans le répertoire
brainvisa-4.3.0/brainvisa/toolboxes/epilepsy/processes.

Ce fichier doit suivre le modèle standard des processus BrainVisa :
importer quelques fichiers, déclarer une signature (de quoi a-t-il
besoin comme paramètres), quelques variables (nom...), une fonction
d'initialisation et une fonction d'exécution qui correspond au bouton
'run' dans l'interface BrainVisa.

Voici un exemple très simple. Pour les processes plus complexes, il faut
s'inspirer des processes existants dans les toolboxes.

| from neuroProcesses import \*
| import shfjGlobals
| from brainvisa import anatomist
| import glob, registration
| name = 'Anatomist Show Electrode Model' **# Nom du processus dans
  l'interface BrainVisa**
| userLevel = 0 **# niveau 0 accessible à tous, niveau 1 : utilisateurs
  avancés, niveau 2 experts**
| roles = ('viewer',) **# certains processus ont des rôles particuliers.
  Celui-ci est un viewer pour le type de données déclaré en premier dans
  la signature, ce qui signifie qu'il sera utilisé dès que l'on clique
  sur l'icône œil à côté d'un fichier de ce type dans BrainVisa. Il
  existe d'autres rôles, comme converter qui permet de transformer un
  format de fichier en un autre.**
| def validation(): **# Permet de vérifier que les paramètres sont
  corrects**
| anatomist.validation()
| **# Ici on veut un seul paramètre, un fichier Electrode Model en
  lecture.**
| signature = Signature(
| 'model', ReadDiskItem( 'Electrode Model', 'Electrode Model format' ),
| )
| **# Si on souhaite préremplir certains paramètres au lancement du
  processus**
| def initialization( self ):
| pass
| **# La fonction qui sera exécuté quand on appuie sur le bouton
  « run »**
| def execution( self, context ):
| a = anatomist.Anatomist()
| elec = ElectrodeEditorDialog(a)
| elec.open(self.model.fullPath()) **# on accède au paramètre déclaré
  dans la signature**
| meshes = elecDialog.getAnatomistObjects()
| w = a.createWindow('Axial')
| a.addObjects(meshes, [w,])
| return (w, elec, meshes) **# On renvoie tous les objets qui ne doivent
  pas être détruit à la fin de l'exécution de la fonction, ici, les
  objets 3D à afficher qui seront détruits plus tard par BrainVisa.**

API Brainvisa
-------------

| BrainVisa est développé principalement au CEA Neurospin par Denis
  Rivière, Yann Cointepas et Isabelle Denghien.
| Il y a beaucoup de documentation en ligne, à partir de :

http://brainvisa.info/doc/cartointernet/cartointernet_pg/en/html/index.html

En particulier dans l'API, le plus important est pyaims (les bindings
Python d'aims qui gère les coordonnées 3D, les maillages, les
images...), pyanatomist (le contrôle d'Anatomist depuis Python, charger
et afficher des images et des maillages 3D).

| L'API BrainVisa est documentée ici :
  http://brainvisa.info/doc/axon-4.4/sphinx/index.html
| On peut aussi regarder directement dans le code source
  brainvisa-4.3.0/brainvisa, par exemple le fichier registration.py qui
  définit les fonctions du transformationManager.

Sinon, on peut poser des questions sur les forums (il y a des forums
pour les développeurs qui sont accessibles sur demande)
http://brainvisa.info/forum/

Interactions avec Matlab
------------------------

Les fonctions présentes dans externalprocesses.py facilitent les appels
à matlab.

| Le plus simple est d'écrire le code matlab dans une chaîne python avec
  des %s à la place des paramètres, et qui se termine par « quit ; »
  comme dans ImageImport.py :
| from externalprocesses import \*

spm_coregister = "VF=spm_vol(%s);VG=spm_vol(%s);x = spm_coreg(VF, VG);\\

trm = spm_matrix(x(:)');trm = [trm(1:3,4)';trm(1:3,1:3)];\\

dlmwrite(%s,trm, 'delimiter',' ','precision',16); \\

quit;"

On remplit les paramètres :

| call = spm_coregister%("'monFichier.img,1'", "'AutreFichier.img,1'",
  "'fichierOutput'")
| On lance l'exécution
| matlabRun(call)

**Cette fonction est bloquante**, donc le logiciel va être bloqué
pendant toute l'exécution du code matlab.

On peut aussi utiliser un **appel non-bloquant** qui crée un objet
Qthread, et connecter son signal de fin d'exécution à une fonction. On
stocke l'objet pour que la thread ne soit pas détruite à la fin de la
fonction qui la crée, et on lance l'execution avec start :

| thr = matlabRunNB(call)
| thr.finished.connect(lambda:self.taskfinished(u"SPM Coregister
  terminé", thr))
| self.threads.append(thr)
| thr.start()

| Si on a besoin de créer un **fichier temporaire** dans lequel la
  fonction matlab va écrire, on peut en créer avec
| tempfile = getTmpFilePath('txt')

**Il faudra l'effacer** à la main ensuite.

Execution distribuée
--------------------

Avec BrainVisa, une librairie est fournie pour l'exécution distribuée
(multi-core sur une même machine, ou sur un cluster de machines). C'est
bien plus complet et puissant qu'externalprocesses donc à utiliser dans
le futur.

http://brainvisa.info/soma/soma-workflow/index.html
