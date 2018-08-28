# -*- coding: utf-8 -*-
""" Function to anonymize patient data in TRC files 
    (c) 2014 - INSERM U836 - Manik Bhattacharjee
    License GNU GPL v3  
"""
import struct

def anonymizeTRC(filepath, firstname="No", lastname="Name", nowrite=False, overwriteMontage=False):
  """ Anonymize patient data in TRC files """
  
  f=open(filepath, "r+b") # Open read-write binary mode
  # should check for a TRC header !
  try:
    f.seek(2)
    headerTitle = f.read(8)
    f.seek(175)
    version = struct.unpack('B', f.read(1))[0]
    if headerTitle != "MICROMED" or version != 4:
      print "Not a MICROMED TRC System 98 v4 file -> ignoring"
      f.close()
      return False
  except:
    # Could not read header info -> wrong file type
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
  f.write(struct.pack('<3B',1,1,1))
  
  if overwriteMontage:
    try:
      # Read unsigned short int number of montages
      f.seek(152)
      nbMontages = struct.unpack('H',f.read(2))[0] 
      # Reading Montage header to get offset of montage data
      f.seek(288)
      if f.read(8) != 'MONTAGE ':
        raise Exception("Incorrect MONTAGE header")
      montageOffs = struct.unpack('I', f.read(4))[0]
      
      for i in range(nbMontages):
        newMontageName = "Montage "+str(i)
        f.seek(montageOffs + 4096*i+264)
        # To print description string : desc = f.read(64); print desc[:desc.find('\x00')]
        f.write(struct.pack('<64B',*[ord(a) for a in newMontageName.ljust(64,'\x00')]))
    except:
      print "Could not overwrite Montage name"
  
  f.close()
    
  return True

