# -*- coding: utf-8 -*-
#
# This software and supporting documentation are distributed by
#      Institut des Neurosciences de Grenoble - INSERM U836
#      France
#
# License : GNU General Public License version 3
#
# Upload/anonymizing easy-to-use software for F-TRACT DB
# Anonymizer functions
#

import os, tempfile, shutil, struct, datetime
import traceback
from struct import pack, unpack
import codecs

############################### DICOM ###############################################

import dicom
from anonymizeDicom import anonymize as anonymizePyDicomFile
from dicom.filereader import InvalidDicomError

def anonymizeDicomFile(path, outputPath, newName="anonymous"):
  oldsize = os.path.getsize(path)
  try: 
    anonymizePyDicomFile(path, outputPath, new_person_name=newName, new_patient_id=newName)
  except InvalidDicomError: # Not dicom file -> just copy the file
    print "Not a dicom file, no anonymization done."
    shutil.copy(path, outputPath)
    if os.path.exists(outputPath):
        newsize = os.path.getsize(outputPath)
        if newsize != oldsize:
            print('ERROR non-DICOM file copy failed  ! '+repr(outputPath))
            raise Exception('ERROR non-DICOM file copy failed  ! '+repr(outputPath))
    else:
        print('ERROR non-DICOM file copy failed (not created) ! '+repr(outputPath))
        raise Exception('ERROR non-DICOM file copy failed (not created) ! '+repr(outputPath))
        
                
  except Exception, err:
    print ("Other error while anonymizing DICOM file %s : %s"%(path, repr(err)))
    traceback.print_exc()
  newsize = os.path.getsize(outputPath)
  if oldsize > 0 and newsize == 0:
      print('ERROR : Dicom file size is 0 after anonymization' + repr(path))
      raise Exception('DICOM anonymization failed : file size is 0 : ' + repr(path))
  return outputPath


def anonymizeDicomDir(path,newName="anonymous", outputPath=None):
  """Anonymizes a directory containing dicom files.
     Non-Dicom files are copied along with the dicom files.
     Warning ! This creates a directory that must be erased after use (if outputPath is empty, of course) !
     @return Path of the temporary directory that contains an anonymous copy of the directory
  """
  if not os.path.isdir(path):
    print "AnonymizeDicom : %s is not a valid directory"%path
    return None
  if outputPath is None:
    tmp = tempfile.mkdtemp(prefix='ftract')
  else:
    tmp = outputPath
  print ('AnonymizeDicomDir will use output directory '+repr(tmp))
  for r,d,s in os.walk(path):
    dOut = os.path.relpath(r, path)
    if dOut != '.':
        print('AnonymizeDicomDir newSubdir')
        try:
            os.mkdir(os.path.join(tmp, dOut))
            print('Making directory '+repr(os.path.join(tmp, dOut)))
        except:
            print('Failed to create directory '+repr(os.path.join(tmp, dOut)))
    for fileIn in s:
        print('Anonymizing file '+repr(os.path.join(r, fileIn)) + ' into '+ os.path.join(tmp, dOut, fileIn) + ' with name = '+newName)
        anonymizeDicomFile(os.path.join(r, fileIn), os.path.join(tmp, dOut, fileIn), newName)
      
  return tmp
  
def anonymizeDicomPaths(paths, newName="anonymous"):
  out = {}
  for p in paths:
    if os.path.isdir(p):
        print('Anonymizing dicom dir '+repr(p))
        out[p]=tempfile.mkdtemp(prefix='ftract')
        anonymizeDicomDir(p, newName, out[p])
        print('Got '+repr(out[p]) + ' from anonymizeDicomDir')
    elif os.path.isfile(p):
        print('Anonymizing dicom file '+repr(p))
        out[p] = tempfile.mkstemp(prefix='ftract')[1]
        anonymizeDicomFile(p, out[p], newName)
        print('Got '+repr(out[p]) + ' from anonymizeDicomFile')
      
  return out


############################### SEEG ###############################################

def anonymizeSeeg(path, newName = "anonymous", filetype = "Micromed"):
  if filetype == "Micromed":
    return anonymizeTrcPaths(path, newName)
  elif filetype == "Nihon Kohden":
    return anonymizeNihonKohdenPaths(path, newName)
  elif filetype == "EDF":
    return anonymizeEdfPaths(path, newName)
  elif filetype == "Nicolet":
    return anonymizeNicoletEPaths(path, newName)
  elif filetype == "xltek" or filetype == "Deltamed":
    print "ANONYMIZATION DOES NOT WORK FOR XLTEK/DELTAMED SEEG DATA !"
    return dict([(p,p) for p in path])
  else:
    print "UNKNOWN FORMAT (should be Micromed, xltek, Nicolet (EDF), or Nihon Kohden)"
    return dict([(p, p) for p in path])
  
def anonPaths(paths, newName="anonymous", dirAnonymizer=None, fileAnonymizer = None):
  """ A function to anonymize files and dirs in a list of paths using fileAnonymizer and dirAnonymizer functions"""
  out = {}
  if type(paths) == str:
      paths = [paths,]
  if dirAnonymizer is None:
    dirAnonymizer = lambda p,n:anonDir(path=p, newName = n, fileAnonymizer = fileAnonymizer)
  if fileAnonymizer is None:
    print "No anonymizer !!!!"
    return None
  for p in paths:
    if os.path.isdir(p):
      out[p] = dirAnonymizer(p, newName)
    elif os.path.isfile(p):
      fileout = tempfile.mkstemp(prefix='ftract')[1]
      shutil.copy (p, fileout)
      fileAnonymizer(fileout, newName)
      out[p] = fileout
      
  return out
  
def anonDir(path, newName = "anonymous", fileAnonymizer = None):
  """ A function to anonymize a full directory by running fileAnonymizer on each file on a copy of the directory"""
  if not os.path.isdir(path):
    print "AnonymizeDir : %s is not a valid directory"%path
    return None
  tmp = tempfile.mkdtemp(prefix='ftract')
  # 
  for r,d,s in os.walk(path):
    dOut = os.path.relpath(r, path)
    if dOut != '.':
      os.mkdir(os.path.join(tmp, dOut))
    for fileIn in s:
      fileout = os.path.join(tmp, dOut, fileIn)
      shutil.copy(os.path.join(r, fileIn), fileout)
      fileAnonymizer(fileout, newName)
  return tmp

##########################################" Nihon Kohden format
  

anonymizeNihonKohdenPaths = lambda paths, newName="anonymous":anonPaths(paths, newName, None, anonymizeNihonKohdenFile)    


def anonymizeNihonKohdenFile(filepath, newName = "anonymous"):
  nihonKohdenAnon = {'CND':anonNihonKohdenDeviceBlock,
                   'CN2':anonNihonKohdenDeviceBlock,
                   'PNT':anonNihonKohdenPnt,
                   'LOG':anonNihonKohdenDeviceBlock,
                   'EEG':anonNihonKohdenDeviceBlock,
                   'EGF': anonNihonKohdenEGF,
                   'CMT':anonNihonKohdenDeviceBlock,}
           # Ignore other files without patient info (.21E, .11D, .REG...)  
  # Check this is a Nihon Kohden file
  extension = os.path.splitext(filepath)[1][1:].upper()
  if extension not in nihonKohdenAnon.keys():
      print "Not a Nihon Kohden file with patient name (extension unknown - not anonymized) : "+extension
      return filepath  
  # Better expensive checks : read some headers and anonymize using the appropriate function
  nihonKohdenAnon[extension](filepath, newName)
  return filepath


def anonNihonKohdenDeviceBlock(filepath, newName):
  print "Anon Nihon Device block for "+filepath
  f = codecs.open(filepath, "r+b", encoding='cp1252') # Open read-write binary mode
  # should check for a header !
  deviceType = f.read(16)
  filePrev = f.read(16)
  fileNext = f.read(16)
  name = os.path.splitext(os.path.basename(filepath))[0]
  prev = os.path.splitext("".join([a for a in filePrev if a!='\x00']))[0]
  next = os.path.splitext("".join([a for a in fileNext if a!='\x00']))[0]
  if (prev != '' and prev != name) or (next != '' and next != name) or not deviceType.startswith('EEG-'):
    print "Invalid previous/next file name in header or device type is not EEG :\n   %s\n  %s\n  %s"%(repr(filePrev), repr(fileNext), repr(deviceType))
    return
   
  idNumber = f.read(16)
  dateTime = f.read(14)
  charType = chr(struct.unpack('B',f.read(1))[0]) # Should be '0' (48) for ascii, '1' (49) for japanese 
  if charType not in ['0','1'] or len([a for a in dateTime if not  a.isdigit()]) != 0:
    print "Nihon Kohden -> invalid char type"
    return
  patientName = f.read(32)
  reserved = f.read(12) # Should be all 0x00
  additionalVersion = f.read(5)
  if len([r for r in reserved if r != '\x00']) != 0 or [a for a in idNumber if not  a.isdigit() and a != '\x00']:
    print "reserved should be all 0x00"
    return
  # Overwrite patient name
  f.seek(79)
  f.write(newName.ljust(32,'\x00')[:32]) 
  f.close()

def anonNihonKohdenPnt(filepath, newName):
  print "Anon PNT "+filepath
  # First anonymize the device block
  anonNihonKohdenDeviceBlock(filepath, newName)
  # Go for specific patient block
  f=codecs.open(filepath, "r+b", encoding='cp1252')
  f.seek(146)
  val = struct.unpack('<L', f.read(4))[0]

  if val != 1024:
    print "Strange Patient block in data -> offset in not 1024 !"
  f.seek(val)
  if f.read(1) != '\x03':
    #print "Invalid value at offset 1024"
    return
  device = f.read(16) # -> device model "EEG-1100A V01.00"
  if f.read(1) != '\x7b':
    return
  # Information about the rest of the header : [[name length, value length, value type, unit],...[]] - 4 bytes, 123=0x7b times
  datainfo = [struct.unpack('BBBB',f.read(4)) for j in xrange(123)]
  # Read the data according to the datainfo
  data = {}
  for i in xrange(len(datainfo)):
    k = f.read(datainfo[i][0]).replace('\x00','').lower()
    data[k] = (f.tell(), f.read(datainfo[i][1]).replace('\x00',''), datainfo[i][1]) # length of the name of the data type
  # Did we find a name ? Overwrite it
  if 'name' in data:
    f.seek(data['name'][0])
    f.write(newName.ljust(data['name'][2],'\x00')[:data['name'][2]])
  # Did we find a date of birth ? Overwrite it
  if 'date of birth' in data:
    f.seek(data['date of birth'][0] + 8) # 8 because we leave as is the first 8 chars "YYYY/MM/" and just set the day to 01
    f.write('01')
  f.close()
  
  
def anonNihonKohdenEGF(filepath, newName):
  print "Anon EGF "+filepath
  #PATNAME=Dupont Jean
  #BIRTHD=2003/01/01
  f = codecs.open(filepath,'rb', encoding='cp1252')
  lines = f.read().split(u"\r\n")
  f.close()
  updateNeeded = False
  for i in xrange(len(lines)):
    if lines[i].startswith('PATNAME='):
      lines[i] = 'PATNAME='+newName
      updateNeeded = True
    elif lines[i].startswith('BIRTHD='):
      lines[i] = lines[i][:15]+'01'#Just sets the day to 01 : BIRTHD=2010/05/18 -> BIRTHD=2010/05/01
      updateNeeded = True
  if updateNeeded:
    import pdb; pdb.set_trace()
    f = codecs.open(filepath, 'wb', encoding='cp1252')
    f.write(u"\r\n".join(lines))
    f.close()
  

############### Nicolet One .e files ############ EDF files ############
class NicoletFile:
    LABELSIZE = 32
    TSLABELSIZE = 64
    UNITSIZE = 16
    ITEMNAMESIZE  = 64
    tagnames = {'ExtraDataTags': 'ExtraDataTags',
          'SegmentStream': 'SegmentStream',
          'DataStream': 'DataStream',
          'InfoChangeStream': 'InfoChangeStream',
          'InfoGuids': 'InfoGuids',
          '{A271CCCB-515D-4590-B6A1-DC170C8D6EE2}': 'TSGUID',
          '{8A19AA48-BEA0-40D5-B89F-667FC578D635}': 'DERIVATIONGUID',
          '{F824D60C-995E-4D94-9578-893C755ECB99}': 'FILTERGUID',
          '{02950361-35BB-4A22-9F0B-C78AAA5DB094}': 'DISPLAYGUID',
          '{8E94EF21-70F5-11D3-8F72-00105A9AFD56}': 'FILEINFOGUID',
          '{E4138BC0-7733-11D3-8685-0050044DAAB1}': 'SRINFOGUID',
          '{C728E565-E5A0-4419-93D2-F6CFC69F3B8F}': 'EVENTTYPEINFOGUID',
          '{D01B34A0-9DBD-11D3-93D3-00500400C148}': 'AUDIOINFOGUID',
          '{BF7C95EF-6C3B-4E70-9E11-779BFFF58EA7}': 'CHANNELGUID',
          '{2DEB82A1-D15F-4770-A4A4-CF03815F52DE}': 'INPUTGUID',
          '{5B036022-2EDC-465F-86EC-C0A4AB1A7A91}': 'INPUTSETTINGSGUID',
          '{99A636F2-51F7-4B9D-9569-C7D45058431A}': 'PHOTICGUID',
          '{55C5E044-5541-4594-9E35-5B3004EF7647}': 'ERRORGUID',
          '{223A3CA0-B5AC-43FB-B0A8-74CF8752BDBE}': 'VIDEOGUID',
          '{0623B545-38BE-4939-B9D0-55F5E241278D}': 'DETECTIONPARAMSGUID',
          '{CE06297D-D9D6-4E4B-8EAC-305EA1243EAB}': 'PAGEGUID',
          '{782B34E8-8E51-4BB9-9701-3227BB882A23}': 'ACCINFOGUID',
          '{3A6E8546-D144-4B55-A2C7-40DF579ED11E}': 'RECCTRLGUID',
          '{D046F2B0-5130-41B1-ABD7-38C12B32FAC3}': 'GUID TRENDINFOGUID',
          '{CBEBA8E6-1CDA-4509-B6C2-6AC2EA7DB8F8}': 'HWINFOGUID',
          '{E11C4CBA-0753-4655-A1E9-2B2309D1545B}': 'VIDEOSYNCGUID',
          '{B9344241-7AC1-42B5-BE9B-B7AFA16CBFA5}': 'SLEEPSCOREINFOGUID',
          '{15B41C32-0294-440E-ADFF-DD8B61C8B5AE}': 'FOURIERSETTINGSGUID',
          '{024FA81F-6A83-43C8-8C82-241A5501F0A1}': 'SPECTRUMGUID',
          '{8032E68A-EA3E-42E8-893E-6E93C59ED515}': 'SIGNALINFOGUID',
          '{30950D98-C39C-4352-AF3E-CB17D5B93DED}': 'SENSORINFOGUID',
          '{F5D39CD3-A340-4172-A1A3-78B2CDBCCB9F}': 'DERIVEDSIGNALINFOGUID',
          '{969FBB89-EE8E-4501-AD40-FB5A448BC4F9}': 'ARTIFACTINFOGUID',
          '{02948284-17EC-4538-A7FA-8E18BD65E167}': 'STUDYINFOGUID',
          '{D0B3FD0B-49D9-4BF0-8929-296DE5A55910}': 'PATIENTINFOGUID',
          '{7842FEF5-A686-459D-8196-769FC0AD99B3}': 'DOCUMENTINFOGUID',
          '{BCDAEE87-2496-4DF4-B07C-8B4E31E3C495}': 'USERSINFOGUID',
          '{B799F680-72A4-11D3-93D3-00500400C148}': 'EVENTGUID',
          '{AF2B3281-7FCE-11D2-B2DE-00104B6FC652}': 'SHORTSAMPLESGUID',
          '{89A091B3-972E-4DA2-9266-261B186302A9}': 'DELAYLINESAMPLESGUID',
          '{291E2381-B3B4-44D1-BB77-8CF5C24420D7}': 'GENERALSAMPLESGUID',
          '{5F11C628-FCCC-4FDD-B429-5EC94CB3AFEB}': 'FILTERSAMPLESGUID',
          '{728087F8-73E1-44D1-8882-C770976478A2}': 'DATEXDATAGUID',
          '{35F356D9-0F1C-4DFE-8286-D3DB3346FD75}': 'TESTINFOGUID'}
    def __init__(self, filename, anonymizeName=None):
        self.fileName = filename
        self.patientInfo = None
        self.segments = None
        self.sections = None
        self.index = None
        self.sigInfo = None
        self.tsInfo = None
        self.chInfo = None
        self.notchFreq = None
        self.Qi = None
        self.Qii = None
        self.allIndexIDs = None
        self.useTSinfoIdx = 1
        
        if os.path.splitext(filename)[1][1:].lower()  != "e":
            print "ERROR : File extension must be .e. Filename is :"+str(filename)
            return
        self.filename = os.path.realpath(filename)
        filename = self.filename
        h = open(filename,'r+b')
        #Get init 
        misc1 = unpack('<5I', h.read(20))#ok<NASGU>
        unknown =unpack('<I',  h.read(4))[0] #ok<NASGU>
        indexIdx = unpack('<I',  h.read(4))[0]
        
        #Get TAGS structure and Channel IDS
        h.seek(172)
        nrTags =  unpack('<I',  h.read(4))[0]
        print nrTags," TAGS"
        Tags = [None] * nrTags # preallocate
        for i in xrange(nrTags):
            Tags[i] = {}
            Tags[i]['tag'] = "".join(map(chr, unpack('<40H', h.read(80)))).strip().rstrip('\x00')
            Tags[i]['index'] = unpack("<I", h.read(4))[0]
            if Tags[i]['tag']  in self.tagnames:
                Tags[i]['IDStr'] =  self.tagnames[Tags[i]['tag'] ]
            elif Tags[i]['tag'].isdigit():
                Tags[i]['IDStr'] = int(Tags[i]['tag'])
            else:
                Tags[i]['IDStr'] ="UNKNOWN"
    
        self.sections = Tags;
        
        # QI index
        h.seek(172208)
        self.Qi={}
        self.Qi['nrEntries'] = unpack("<I", h.read(4))[0]
        self.Qi['misc1'] = unpack("<I", h.read(4))[0]
        self.Qi['indexIdx'] = unpack("<I", h.read(4))[0]
        self.Qi['misc3'] = unpack("<I", h.read(4))[0]
        self.Qi['LQi'] = unpack("<Q", h.read(8))[0]
        self.Qi['firstIdx'] = unpack("<%dQ"%nrTags, h.read(nrTags*8)); # Read nrTags 64 bits integers

        #% Don't know what this index is for... Not required to get data and
        #% can be huge...
        h.seek(188664)
        self.Qindex  = [None]*self.Qi['LQi']
        for i in xrange(self.Qi['LQi']):
          self.Qindex[i] = {
              'ftel':h.tell(),
              'index': unpack("<2H", h.read(4)),
              'misc1': unpack("<I", h.read(4))[0],
              'indexIdx': unpack("<I", h.read(4))[0],
              'misc2': unpack("<3I", h.read(12)),
              'sectionIdx': unpack("<I", h.read(4))[0],
              'misc3': unpack("<I", h.read(4))[0],
              'offset': unpack("<Q", h.read(8))[0],
              'blockL': unpack("<I", h.read(4))[0],
              'dataL': unpack("<I", h.read(4))[0],
              }
          
        self.Qi['index'] = self.Qindex

        ## Get Main Index: 
        # Index consists of multiple blocks, after each block is the pointer
        # to the next block. Total number of entries is in obj.Qi.nrEntries
      
        Index = []
        curIdx = 0
        nextIndexPointer = indexIdx;
        print 'Parsing index '
        curIdx2 = 1
        while curIdx < self.Qi['nrEntries']:
            if curIdx2 % 20 == 0:
                print '.'
            else:
                print "\n."
        
            h.seek(nextIndexPointer);
            nrIdx = unpack('Q', h.read(8))[0]
            for i in xrange (nrIdx):
                Index.append ( {} ) 
                Index[curIdx + i]['sectionIdx'] = unpack('<Q', h.read(8))[0]
                Index[curIdx + i]['offset'] = unpack('<Q', h.read(8))[0]
                Index[curIdx + i]['blockL'] =  unpack("<I", h.read(4))[0]
                Index[curIdx + i]['sectionL'] =  unpack("<I", h.read(4))[0]
            
            nextIndexPointer =unpack('<Q', h.read(8))[0]
            curIdx = curIdx + nrIdx;
            curIdx2=curIdx2+1;
        
        print 'done'
        self.index = Index; 
        self.allIndexIDs = [idx['sectionIdx'] for idx in self.index]
      
      # Get PatientGUID
        info ={}
        
        infoProps = [ 'patientID', 'firstName','middleName','lastName',
            'altID','mothersMaidenName','DOB','DOD','street','sexID','phone',
            'notes','dominance','siteID','suffix','prefix','degree','apartment',
            'city','state','country','language','height','weight','race','religion',
            'maritalStatus'];
        ifnoIdx = [(i,tag) for i, tag in enumerate(Tags) if tag['IDStr'] == 'PATIENTINFOGUID']
        indexInstance = [idx for idx in Index if idx['sectionIdx'] == ifnoIdx[0][0]]
        h.seek(indexInstance[0]['offset'])
        guid = unpack('<16B', h.read(16))
        lSection = unpack('<Q', h.read(8))[0]
        #reserved = unpack('<3H', h.read(6))
        nrValues = unpack('<Q', h.read(8))[0]
        nrBstr = unpack('<Q', h.read(8))[0]
        self.segments = []
        for i in xrange(nrValues):
            id = unpack('<Q', h.read(8))[0]
            print "ID "+str(id)
            self.segments.append({})
            if id in [7,8]:
                unix_time = unpack('<d', h.read(8))[0]
                self.segments[i]['dateStr'] =unix_time
                unix_time = (unix_time*(3600*24)) - 2209161600 # The date is stored as a double in days from 1/1/1900 -> convert to seconds and remove duration until 1/1/1970 (unix base time)
                #thetime = datetime.datetime.fromtimestamp(unix_time) # This can crash, and we don't care about dates anyway
                #self.segments[i]['dateStr'] = thetime.strftime('%Y-%m-%d %H:%M:%S')
                #datestr(unix_time/86400 + datenum(1970,1,1));
                #value = [thetime.day, thetime.month, thetime.year] # Day month year
                value = unix_time
            if id in [23,24]:
                value = unpack('<d', h.read(8))[0]
            else:
                value = 0

            info[infoProps[id]] = value;  
        
      
        strSetup = unpack('<'+str(2*nrBstr)+'Q', h.read(8*2*nrBstr))
        toWrite = {}
        for i in xrange (0, nrBstr*2, 2):
            id  = strSetup[i]
            print "INFO : ", infoProps[id], " AT ", h.tell()
            if anonymizeName is not None and  infoProps[id] in [ 'firstName','middleName','lastName','mothersMaidenName',  'patientID',  'altID']:
                beginCurrent = h.tell()
                nbCharU16 = strSetup[i+1]
                newVal = anonymizeName[:nbCharU16].ljust(nbCharU16, ' ')+'\0'
                #newVal = '\0'.join([nv for nv in newVal])+'\0'
                print("Will write %s at %s with length %s"%(str(newVal), str(beginCurrent), str(nbCharU16)))
                #h.write(pack('<'+str(nbCharU16 + 1)+'H', *[ord(c) for c in newVal]))
                toWrite[beginCurrent] = pack('<'+str(nbCharU16 + 1)+'H', *[ord(c) for c in newVal])
                h.seek(beginCurrent)
                print( " Going back to "+str(beginCurrent))
            print ("Reading now ->")
            value = "".join(map (chr, unpack('<'+str( strSetup[i+1] + 1)+'H', h.read(2*(strSetup[i+1] + 1))))).rstrip(' \t\r\n\0')
            info[infoProps[id]] = value;
            print "        "+str(value)
        self.patientInfo = info;
        # Close and reopen for write, because writing directly fails in windows XP with no error message...
        h.close()
        print("Reopening Nicolet file  to write inside")
        h = open(filename,'r+b')
        for offset in toWrite:
          h.seek(offset)
          h.write(toWrite[offset])
          print("Writing at " + str(offset))
        h.close()
        
        
def anonymizeNicoletE(filepath, newName = "anonymous"):
    """Anonymizes a .e file (from Nicolet)
    WARNING overwrites the file ! """
    try:
      nf = NicoletFile(filepath, anonymizeName = newName)
    except Exception as e:
        print "Could not anonymize "+ filepath+" : NicoletFile exception when reading file !"+str(e)
        raise

    
    
def anonymizeNicoletEPaths(paths, newName="anonymous"):
  out = {}
  for p in paths:
    if os.path.isdir(p):
      out[p] = anonymizeNicoletEDir(p, newName)
    elif os.path.isfile(p):
      fileout = tempfile.mkstemp(prefix='ftract')[1]
      shutil.copy (p, fileout)
      anonymizeNicoletE(fileout, newName)
      out[p] = fileout    
  return out
  
def anonymizeNicoletEDir(path, newName = "anonymous"):
  """Anonymizes a directory containing TRC files.
     Non-TRC files are copied along with the TRC files.
     Warning ! This creates a directory that must be erased after use !
     @return Path of the temporary directory that contains an anonymous copy of the directory
  """
  if not os.path.isdir(path):
    print "AnonymizeNicoletE : %s is not a valid directory"%path
    return None
  tmp = tempfile.mkdtemp(prefix='ftract')
  # 
  for r,d,s in os.walk(path):
    dOut = os.path.relpath(r, path)
    if dOut != '.':
      os.mkdir(os.path.join(tmp, dOut))
    for fileIn in s:
      fileout = os.path.join(tmp, dOut, fileIn)
      shutil.copy(os.path.join(r, fileIn), fileout)
      anonymizeNicoletE(fileout, newName)
  return tmp


###############       EDF files     ############
def anonymizeEdf(filepath, newName = "anonymous"):
  """Anonymizes a EDF file (from Nicolet One)
    WARNING overwrites the file !
  """
  f = open(filepath, 'r+b')
  f.seek(88)
  if f.read(9) != 'Startdate':
    print "No Startdate in file : probably not EDF -> no anonymization !"
    return
  # reserved header : if it starts with EDF+ this is an EDF+ file
  f.seek(192)
  isEdfPlus = False
  if f.read(4) == 'EDF+':
      isEdfPlus = True
  # Go to the patient data and overwrite it using EDF+ format (space separated subfields, X for no value)
  # patient code (no spaces, use _ instead), sex (F or M),  dd-MMM-yyyy birthdate, patient name
  f.seek(8)
  patientInfo = f.read(80)
  if isEdfPlus:
      patientInfo = patientInfo.split(u' ')
      patId = patientInfo[0]
      patId = newName.replace(u' ', u'_')[:30]
      patSex = patientInfo[1]
      patBirthdate = patientInfo[2]
      patName = patientInfo[3]
      patName = newName.replace(u' ', u'_')[:30]
      if len(patBirthdate) >10: # Should be 21-AUG-2015
          patBirthdate ='01'+patBirthdate[2:]
      patientInfo = ' '.join([patId, patSex, patBirthdate, patName])
  else: # EDF
      patientInfo = ' '.join([newName.replace(u' ', u'_')[:30], u'X', u'X',  newName.replace(u' ', u'_')[:30]])
  f.seek(8)
  f.write(patientInfo[:80].ljust(80,' '))
  f.close()

def anonymizeEdfPaths(paths, newName=u"anonymous"):
  out = {}
  for p in paths:
    if os.path.isdir(p):
      out[p] = anonymizeEdfDir(p, newName)
    elif os.path.isfile(p):
      fileout = tempfile.mkstemp(prefix='ftract')[1]
      shutil.copy (p, fileout)
      anonymizeEdf(fileout, newName)
      out[p] = fileout    
  return out
  
def anonymizeEdfDir(path, newName = u"anonymous"):
  """Anonymizes a directory containing TRC files.
     Non-TRC files are copied along with the TRC files.
     Warning ! This creates a directory that must be erased after use !
     @return Path of the temporary directory that contains an anonymous copy of the directory
  """
  if not os.path.isdir(path):
    print "AnonymizeEDF : %s is not a valid directory"%path
    return None
  tmp = tempfile.mkdtemp(prefix='ftract')
  # 
  for r,d,s in os.walk(path):
    dOut = os.path.relpath(r, path)
    if dOut != '.':
      os.mkdir(os.path.join(tmp, dOut))
    for fileIn in s:
      fileout = os.path.join(tmp, dOut, fileIn)
      shutil.copy(os.path.join(r, fileIn), fileout)
      anonymizeEdf(fileout, newName)
  return tmp

#################### MICROMED TRC FORMAT System 98 type 4
def anonymizeTrcDir(path, newName = "anonymous"):
  """Anonymizes a directory containing TRC files.
     Non-TRC files are copied along with the TRC files.
     Warning ! This creates a directory that must be erased after use !
     @return Path of the temporary directory that contains an anonymous copy of the directory
  """
  if not os.path.isdir(path):
    print "AnonymizeTRC : %s is not a valid directory"%path
    return None
  tmp = tempfile.mkdtemp(prefix='ftract')
  # 
  for r,d,s in os.walk(path):
    dOut = os.path.relpath(r, path)
    if dOut != '.':
      os.mkdir(os.path.join(tmp, dOut))
    for fileIn in s:
      fileout = os.path.join(tmp, dOut, fileIn)
      shutil.copy(os.path.join(r, fileIn), fileout)
      anonymizeTRC_Sys98_t4(fileout, newName, newName)
  return tmp
  
def anonymizeTrcPaths(paths, newName="anonymous"):
  out = {}
  for p in paths:
    if os.path.isdir(p):
      out[p] = anonymizeTrcDir(p, newName)
    elif os.path.isfile(p):
      fileout = tempfile.mkstemp(prefix='ftract')[1]
      shutil.copy (p, fileout)
      anonymizeTRC_Sys98_t4(fileout, firstname=newName, lastname=newName)
      out[p] = fileout
  return out

def anonymizeTRC_Sys98_t4(filepath, firstname="No", lastname="Name", nowrite=False, overwriteMontage=True, overwriteUndocumentedElectrode = True): # Anonymize Micromed's System 98 type 4 TRC files
  """ Anonymize Micromed's System 98 type 4 TRC files """
  f=open(filepath, "r+b") # Open read-write binary mode
  # should check for a TRC header !
  f.seek(2)
  headerTitle = f.read(26)
  f.seek(175)
  headerType = f.read(1)
  trcVersion =  struct.unpack('B',headerType)[0]
  if headerTitle != "MICROMED  Brain-Quick file" or (trcVersion != 4 and trcVersion != 3):
    print "Not a MICROMED System 98 Type 3/4 TRC file -> ignoring"
    f.close()
    return False
  f.seek(64) # go to patient data offset in header
  # get a 22-char string, padding with spaces, convert to integers and pack as 22 unsigned chars in little endian
  if nowrite:
    f.close()
    return True
  f.write(struct.pack('<22B',*[ord(a) for a in lastname[:22].ljust(22,' ')])) 
  # Same with 20 chars 
  f.write(struct.pack('<20B',*[ord(a) for a in firstname[:20].ljust(20,' ')]))
  # Same with date (3 unsigned chars, for example [10, 05, 72] for october 5th 1972
  f.seek(107)
  f.write(struct.pack('<B',1))
  
  if overwriteMontage:
    try:
      if trcVersion == 3:
          sizeMontage = 3072
      elif trcVersion == 4:
          sizeMontage = 4096
      # Read unsigned short int number of montages
      f.seek(152)
      nbMontages = struct.unpack('H',f.read(2))[0] 
      # Reading Montage header to get offset of montage data
      f.seek(288)
      if f.read(8) != 'MONTAGE ':
        raise Exception("Incorrect MONTAGE header")
      montageOffs = struct.unpack('I', f.read(4))[0] 
      # Montage name
      for i in range(nbMontages):
        newMontageName = "Montage "+str(i)
        f.seek(montageOffs + sizeMontage*i+264)
        # To print description string : desc = f.read(64); print desc[:desc.find('\x00')]
        f.write(struct.pack('<64B',*[ord(a) for a in newMontageName.ljust(64,'\x00')]))
        
        
      # Montage name in HISTORY : find history offset, 
      f.seek(336)
      if f.read(8) != 'HISTORY ':
          raise Exception ("Incorrect HISTORY header")
      historyOffs = struct.unpack('I', f.read(4))[0] 
      # There is first an area unsigned long int[MAX_SAMPLE] --> 128* 4-bytes values (the sample where themontage was changed when viewing)
      # We can skip that, then MAX_HISTORY = 128  "specific montage" structures which are identical to the montages above.
      # Description string starts at offset 264 of each montage structure and is 64 bytes long,, just like above.
      for i in range(30):
          newMontageName = "Montage History "+str(i)
          f.seek(historyOffs + 128*4 + sizeMontage*i+264)
          desc = f.read(64)
          # To print description string : desc = f.read(64); print desc[:desc.find('\x00')]
          if "".join(["\x00" for j in range(64)]) != desc:
              f.seek(historyOffs + 128*4 + sizeMontage*i+264)
              f.write(struct.pack('<64B',*[ord(a) for a in newMontageName.ljust(64,'\x00')]))
              
      # undocumented string in place of the last electrode in  name in HISTORY : find history offset, 
      f.seek(192)
      if f.read(8) != 'LABCOD  ':
          raise Exception ("Incorrect LABCOD header")
      elecOffs = struct.unpack('I', f.read(4))[0]
      # There should be 640 electrode structures of 128 bytes. But the last one is not an electrode structure. It seems to be a 32 bits integer and a 64 bytes string that contains a montage name...
      if overwriteUndocumentedElectrode:
          f.seek(elecOffs + 639*128 + 4)
          f.write(struct.pack('<64B',*[ord(a) for a in "undocumented montage".ljust(64,'\x00')]))
    
    except:
      print "Could not overwrite Montage name"
  
  f.close()
  return True
