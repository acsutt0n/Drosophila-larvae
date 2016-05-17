"""
staAnalysis --

Clampfit can spit out these .sta files containing a few measurements
that may indicate the health of the cell, etc, such as Ra, Rm, ...

Currently, this should be used interactively. Eventually I may make
a version that can work with an argument (.sta file name) passed
from the command line.
"""


# Imports

import sys
if sys.version_info.major > 2:
  print('Need to run with python 2.7 (for PyMC capability)!!')
else:
  import pymc as pm
  
  
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import subprocess
import collections



#########################################################################
# Analyzing stats files

def loadStats(stafil, returnmeta=False):
  """
  Rudimentary parsing of a sta file.
  """
  meta = collections.namedtuple('meta', 'comment, starttime, startdate, stopwatchtime, exported')
  meta_ = {}
  data = None
  
  with open(stafil, 'r') as fIn:
    start = False
    for line in fIn:
      if line:
        lline = line.strip()
        try:
          splitLine = lline.split(None)
        
          # I don't feel great about hard-coding this, but time *should* be first...
          if len(splitLine) > 0:
            if splitLine[0] == '"Time' and start is False:
              # print(lline)
              data = {lline.split('"')[u]: [] for u in [1,3,5,7,9,11]}
              order = [lline.split('"')[u] for u in [1,3,5,7,9,11] ]
              start = True
            elif len(splitLine) > 3 and start is True:
              for j in range(len(splitLine)):
                data[order[j]].append(splitLine[j])
        except:
          pass
    # Now the whole dictionary should be populated
  
  if data is None:
    print('Nothing found for %s' %stafil)
    return None
  
  if returnmeta:
    return pd.DataFrame(data=data), meta_
  return pd.DataFrame(data=data,)



def statsByCell(df, pathlistfile=None):
  """
  Return two stats snapshots per sta file belonging to each cell.
  """
  # get the attributes first
  if pathlistfile is None:
    pathlistfile = '~/data/misc/maria/temp.txt'
    stas = findExtInDF(df, ext='sta', refCol=0, outfile=pathlistfile)
  paths = []
  with open(pathlistfile, 'r') as fIn:
    for line in fIn:
      paths.append(line.strip().split(','))
  
  stas = findExtInDF(df, ext='sta', refCol=0, outfile=None)
  # For each cell, get the sta files associated with it
  cells = {}
  # Get the conditions:
  conds = allConditions(df, ) # Can add options here
  for s in range(len(paths)):
    try:
      #print(paths[s])#, getFullPath(paths[s][1])[0])
      temp_ = loadStats(getFullPath(paths[s][1])[0]) # temp dataframe
      # print(temp_.columns)
      if stas[s][0] not in cells.keys():
        cells[stas[s][0]] = {k: [] for k in temp_.columns[0:]}
      for c in range(0,temp_.shape[1]):
        try:
          cells[stas[s][0]][temp_.columns[c]].append(temp_.ix[0,c])
          cells[stas[s][0]][temp_.columns[c]].append(temp_.ix[temp_.shape[0]-1,c])
        except:
          print('no home for %s' %temp_.columns[c])
    except:
      pass
  
  # Clean the conditions a little
  for k in conds.keys():
    cnt = 0
    while cnt < len(conds[k]):
      if 'Notes' in conds[k][cnt]:
        conds[k].pop(cnt)
      else:
        cnt += 1
  
  # Now process each cell by cell type
  for c in cells.keys(): # For each cell
    rowid = list(df.ix[:,0].values).index(c)
    cells[c]['type'] = df.ix[rowid,5]
    cells[c]['genotype'] = df.ix[rowid, 4] # And get genotype
    try:
      cells[c]['tx'] = conds[c]
    except:
      pass
      
  return cells



def simpleProp(cell_, where=0):
  # Pass a sub-dict of cells
  proto = {}
  for p in ['Rm', 'Ra', 'Cm']:
    for k in cell_.keys():
      if p in k:
        proto[p] = cell_[k][where] # This is a list
  if len(proto) < 3:
    print('Found fewer than 3 properties!')
    print(proto)
  return proto

  

def baselineStats(cells, show=True):
  """
  baseline stats -- pass a cell dict from statsByCell
  """
  # All by genotype
  props, glist, clist = {'Rm': {}, 'Cm': {}, 'Ra' : {}}, [], []
  for c in cells.keys():
    gen_, typ_ = cells[c]['genotype'], cells[c]['type']
    glist.append(gen_)
    clist.append(typ_)
    for pr in props.keys():
      if gen_ not in props[pr].keys():
        props[pr][gen_] = {}
      if typ_ not in props[pr][gen_].keys():
        props[pr][gen_][typ_] = []
      
    # Now that the dict is prepared, add the info
    temp_ = simpleProp(cells[c], where=-1)
    for k in temp_.keys():
      try:
        props[k][gen_][typ_].append(float(temp_[k]))
      except:
        print(k, gen_, typ_)
  
  # Round out the dictionaries
  for p in props.keys():
    for g in glist:
      for c in clist:
        try:
          _ = props[p][g][c]
        except:
          props[p][g][c] = []
  
  # Some basic plotting
  if show:
    genoByCell(props)
  return props



def deltaStats(cells, show=True):
  """
  Show the change in stats from the first timepoint to the last,
  assume chonological order of sta files.
  """
  # Make the dictionary
  props, glist, clist = {'Rm': {}, 'Cm': {}, 'Ra' : {}}, [], []
  for c in cells.keys():
    gen_, typ_ = cells[c]['genotype'], cells[c]['type']
    glist.append(gen_)
    clist.append(typ_)
    for pr in props.keys():
      if gen_ not in props[pr].keys():
        props[pr][gen_] = {}
      if typ_ not in props[pr][gen_].keys():
        props[pr][gen_][typ_] = []
      
    # Then calculate the % change
    temp_0 = simpleProp(cells[c], where=0)
    temp_1 = simpleProp(cells[c], where=-1)
    for k in temp_0.keys():
      try:
        props[k][gen_][typ_].append(float(temp_1[k])/float(temp_0[k]))
      except:
        print(k, gen_, typ_)
    
  # Round out the dictionaries
  for p in props.keys():
    for g in glist:
      for c in clist:
        try:
          _ = props[p][g][c]
        except:
          props[p][g][c] = []
  
  # show the plot
  if show:
    genoByCell(props)
  return props

  

#

def minTreatments(cells, show=True):
  """
  cells is a dict of cells, each cell has a tx:[] member.
  """
  # Minimum treatments
  txlist = []
  for c in cells.keys():
    for t in cells[c]['tx']:
      txlist.append(t)
  txcnt = [txlist.count(t) for t in list(set(txlist))]
  mintx = [t for t in list(set(txlist)) if txlist.count(t)==max(txcnt)]
  
  # Make dict
  props, glist, clist = {'Rm': {}, 'Cm': {}, 'Ra' : {}}, [], []
  for c in cells.keys():
    gen_, typ_ = cells[c]['genotype'], cells[c]['type']
    glist.append(gen_)
    clist.append(typ_)
    for pr in props.keys():
      if gen_ not in props[pr].keys():
        props[pr][gen_] = {}
      if typ_ not in props[pr][gen_].keys():
        props[pr][gen_][typ_] = []
      
    # Only add %-change for minimum tx (no tz)
    if set(cells[c]['tx']) == set(mintx):
      temp_0, temp_1 = simpleProp(cells[c], where=0), simpleProp(cells[c], where=-1)
      for k in temp_0.keys():
        try:
          props[k][gen_][typ_].append(float(temp_1[k])/float(temp_0[k]))
        except:
          print(k, gen_, typ_)
  
  # Round out the dictionaries
  for p in props.keys():
    for g in glist:
      for c in clist:
        try:
          _ = props[p][g][c]
        except:
          props[p][g][c] = []
  
  # Some basic plotting
  if show:
    genoByCell(props)
  return props







########################################################################

if __name__ == "__main__":
  args = sys.argv
  print('Module should be used interactively')














