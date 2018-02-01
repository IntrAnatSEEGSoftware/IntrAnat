#! /usr/bin/env python
# -*- coding: utf-8 -*-

import csv, sys, os
import pdb

def readFunctionalTractography(functional_result_file):

  if os.path.isfile(functional_result_file):
      
      if functional_result_file.split(".")[-1] == 'csv':
          
          full_dictionnary = {}
          lastParam = None
          with open(functional_result_file, 'r') as csvfile:
            
            spamreader = csv.reader(csvfile, delimiter='\t')
            for row in spamreader:
              if len(row)<1:
                  continue
              if row[0] == 'Patient' or row[0] == 'Atlas':
                  full_dictionnary.update({row[0]:row[1]})
              elif row[0] == 'Param':
                  lastParam = row[1]
                  full_dictionnary.update({row[1]:{}})
              elif row[0] == 'Parcel_name':
                  import copy
                  row_parcel_names = copy.deepcopy(row[1:])
                  row_parcel_names = [x.strip() for x in row_parcel_names]
                  #for i_parcels in range(len(row)-1):
                      #full_dictionnary[lastParam].update({row[i_parcels+1].strip():{}})
              elif row[0] == 'EndParcelsNames':
                  lastParam = None
              else:
                 if lastParam is not None:
                     if row[0].strip() not in full_dictionnary[lastParam].keys():
                       full_dictionnary[lastParam].update({row[0].strip():{}})
                     else:
                         pdb.set_trace()
                     #if row[0].strip() in full_dictionnary[lastParam].keys():
                     for i_parcels_bis in range(len(row_parcel_names)):
                        full_dictionnary[lastParam][row[0].strip()].update({row_parcel_names[i_parcels_bis]:row[i_parcels_bis+1]})
                             
                             
            return full_dictionnary
          
      else:
          print("error, wrong file extension")
      
  else:
      print("error, the file doesn't exist")

