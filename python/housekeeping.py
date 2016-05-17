"""
Housekeeping --

These functions help find cells in data frames associated with .abf 
files, or certain treatments, can help identify sample sizes,
can get full paths of cells where only a file name is known, etc.

These are generally very useful but take a few minutes to figure out
how they work (and what kind of input they prefer).

These functions will not do anything when called from the command line;
they need to be used interactively.
"""


# Make sure we're in python 2.7
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




def writeNewDF(df, original_name, new_name):
  """
  Write the new df (with the clust_inds column) to a new csv file,
  but use all the columns from the original csv.
  """
  try:
    orig = pd.read_csv(original_name)
  except:
    print('Could not find or load %s' %original_name)
  
  newcol = []
  orig_tms, new_tms = list(orig.times.values), list(df.times.values)
  for i in range(orig.shape[0]):
    if orig_tms[i] in new_tms:
      newcol.append(df.clust_inds.values[new_tms.index(orig_tms[i])])
    else:
      newcol.append(np.nan)
  
  orig['clust_inds'] = newcol
  print('Assigned %i values, the rest (%i) were nan' 
        %(len([i for i in newcol if i != np.nan]), newcol.count(np.nan)))
  
  orig.to_csv(new_name)
  return
  


def findExtInDF(df, ext='abf', labrow=1, refCol=0, outfile=None):
  """
  Find all filenames in a dataframe (from csv, excel) that match
  a certain extension -- used to find all .abf files in Maria's stuff.
  refCol can be any column number, i.e.: where cell type is or associated stats files.
  """
  fils = []
  for i in range(df.shape[0]): # For every row
    for c in df.columns:
      try:
        temp_ = df.iloc[i][c]
      except:
        print(i, c)
      try:
        if temp_.split('.')[-1] == ext:
          if labrow is not None:
            col = df[c].values[labrow]
          else:
            col = c
          try: # Could be multiple files
            temp_ = ''.join([u for u in temp_ if u != ' '])
            csplit = temp_.split(',')
            for ent in csplit:
              fils.append([df.iloc[i][df.columns[refCol]], ent, col])
          except:
            fils.append([df.iloc[i][df.columns[refCol]], temp_, col]) # there was only one file
      except:
        pass
  print('Found %i hits' %len(fils))
  
  # This part will save only the filenames in a txt doc
  # which can be passed to getPaths.sh to return the full paths to each file.
  if outfile is not None:
    with open(outfile, 'w') as fOut:
      for f in fils:
        fOut.write(','.join([str(i) for i in f]))
        fOut.write('\n')
  return fils



def burstersFromCSV(csvfile, whichCol=1, includes=['B'], outfile=None):
  """
  This opens a csvfile and looks in _whichCol_ for any string that
  includes _includes_ and keeps the path for analysis. Can write path.
  """
  df = pd.read_csv(csvfile)
  includes = [o.lower() for o in includes]
  paths = []
  for u in range(df.shape[0]):
    for i in includes:
      try:
        this = df.ix[u,whichCol].lower()
      except: # Like for NaN
        this = ''
      if i in this and df.ix[u,0] not in paths:
        paths.append(df.ix[u,0])
  
  if outfile is not None:
    with open(outfile, 'w') as fOut:
      for p in paths:
        fOut.write('%s\n' %p)
    return
  return paths

#


def rmPath(infile, rmpath='/Users/mariagenco/Documents/',
           newpath='/media/alex/BACKUP/mcgenco/', outfile=None):
  """
  Removes rmpath from each path in infile and replaces with newpath.
  """
  newpaths = []
  if '/' in infile:
    paths = []
    with open(infile, 'r') as fIn:
      for line in fIn:
        paths.append(line.strip())
  else:
    paths = infile
  
  for p in paths:
    if rmpath in p:
      t_p = p.split(rmpath)[1]
      newpaths.append(newpath+t_p)
  if outfile is not None:
    with open(outfile, 'w') as fOut:
      for p in newpaths:
        write('%s\n' %p)
    return
  return newpaths
  
#

  
def getPathfromCSV(flist, full_csv_list, outfile=None):
  """
  Given a txt file of filenames, this returns a path for each filename
  (if it exists).
  """
  # Get the file list
  fils, csvfils = [], []
  with open(flist, 'r') as fIn:
    for line in fIn:
      if line:
        fils.append(line.split('.')[0].strip())
  with open(full_csv_list, 'r') as fIn:
    for line in fIn:
      if line:
        csvfils.append([line.split('/')[-1].split('_')[0].strip(), # Filename only
                        line.strip()]) # File path only
  
  # replace it with the path list
  paths = []
  for f in fils:
    if f in [o[0] for o in csvfils]:
      idx = [o[0] for o in csvfils].index(f)
      paths.append(csvfils[idx][1])
    else:
      print('Could not find: %s' %f)
  
  print('Retrieved %i paths (of %i)' %(len(paths), len(fils)))
  if outfile is not None:
    with open(outfile, 'w') as fOut:
      for p in paths:
        fOut.write(p)
        fOut.write('\n')
  
  return paths
  
#


def getFullPath(pathfil, direct='/media/alex/BACKUP/mcgenco/', outfile=None):
  """
  Give a hint of the directory to look in, this finds the full path.
  """
  if '/' in pathfil:
    pshorts = []
    with open(pathfil, 'r') as fIn:
      for line in fIn:
        pshorts.append(line.strip())
  else:
    pshorts = [pathfil]
  
  # Get the full paths
  paths = []
  for nam in pshorts:
    p = subprocess.Popen(['find', direct, '-name', nam],
          stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    out = str(p.stdout.read())
    out = '/'.join(out.split('/')[1:])
    out = '/'+out.split('\\')[0]
    paths.append(str(out.strip()))
  
  if outfile is not None:
    with open(outfile, 'w') as fOut:
      for p in paths:
        fOut.write(p)
        fOut.write('\n')
    print('%s written' %outfile)
    return
  return paths

#


def allConditions(df, startCol=9, labrow=1, refCol=0, outfile=None):
  """
  Note all the treatments/conditions that a cell has been through,
  indicated by any non-empty df element.
  """
  # Start at row labrow+1 and column 
  cells = {}
  for c in range(labrow+1, df.shape[0]):
    txs = []
    for tx in range(startCol, df.shape[1]):
      if not pd.isnull(df.ix[c][tx]):
        txs.append(df.ix[labrow][tx])
    cells[df.ix[c][refCol]] = txs
  
  if outfile is not None:
    with open(outfile, 'w') as fOut:
      for ckey in cells.keys():
        fOut.write(','.join(cells[ckey]))
        fOut.write('\n')
    print('%s written' %outfile)
    return
  return cells

#

  









#########################################################################

if __name__ == "__main__":
  print('Module should be used interactively.')






