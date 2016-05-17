"""
activityAnalysis --

This is best used interactively with either a list of spike-sorted
DataFrames (pandas) or a single DF. Much of this is used to explore
the data and analyses are currently ad hoc.

The first half of the software is ANALYSIS code used to classify
bursting, tonic and silent time and compare these factors across groups,
genotypes, cell types, etc.

The second half is PLOTTING oriented, but since almost every analysis
function has an associated plotting method it doesn't make sense to
separate them.
"""


# Imports

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




def outlier(arr, as_nan=True, thresh=0.05, show=False, report=False):
  """
  Return nan instead (more robust) of nothing (loss of index parity).
  Median is more robust than mean.
  """
  if len(arr) < 3:
    return arr
  if show:
    plt.subplot(1,2,1) # Plot part 1 first
    plt.plot(np.random.random(len(arr)), thing1, 'o', color='blue',
             markeredgecolor='none', alpha=0.4)
    plt.title('With outliers')
  
  med_res = [(np.median(arr)-i)**2 for i in arr] 
  med_res_ix = [u for u in med_res] # Create index
  arr_copy = [u for u in arr] # The copy will be edited first
  stds = []
  med_res.sort(reverse=True) # Largest to smallest
  # print(med_res[:10])
  numPts = max([int(len(arr)*thresh), 2])
  # print('Testing largest %i residuals' %numPts)
  
  # Pretend to remove 10% of points
  for i in range(numPts): #for i in range(int(len(arr)*.1)): #
    stds.append(np.std(arr_copy))
    rm_ix = med_res_ix.index(med_res[i])
    try:
      rm = arr[rm_ix]
    except:
      print('tried to remove ix %i but arr is len %i'
              %(rm_ix, len(arr)))
    try:      
      arr_copy.pop(arr_copy.index(rm))
    except:
      print('tried to remove %f but not in arr_copy' %rm)
  
  # Find the greatest d(std)
  dstd = np.diff(stds)
  dstd = [abs(i) for i in dstd]
  rm_to = list(dstd).index(max(dstd))+1 # len(diff) = len(arr)-1

  #print('Mean d(std): %.3f, removing all above %.3f (%i pts)'
        # %(np.mean(dstd), dstd[rm_to-1], rm_to))
  
  for i in range(rm_to):
    arr[med_res_ix.index(med_res[i])] = np.nan
    
  if show: # Show
    plt.subplot(1,2,2)
    plt.plot(np.random.random(len(arr)), arr, 'o',
             color='red', markeredgecolor='none', alpha=0.4)
    plt.title('Without outliers')
    plt.show()
  if as_nan:
    return arr
  return [i for i in arr if not pd.isnull(i)] # Else just eliminate it.



def batchAnalysis(groupfil):
  """
  Compile lists of the features listed below grouped by the groupfil.
  A line of groupfil is /path/to/file.csv,group name
  """
  groups = []
  with open(groupfil, 'r') as fIn:
    for line in fIn:
      groups.append(line.strip().split(','))
  
  checks = ['maxV', 'maxDerivV', 'maxDerivdV', 'minDerivV',
        'minDerivdV', 'preMinV', 'postMinV', 'preMaxCurveV',
        'preMaxCurveK', 'postMaxCurveV', 'postMaxCurveK',
        'height', 'repolarizationV', 'intervals', 'frequencies']
  props = {ch: {gr: {} for gr in list(set([g[1] for g in groups]))}
                for ch in checks} # A dict of dicts
  #  props [properties] [group name] [cell name]
  cells = [f[0].split('/')[-1].split('_')[0] for f in groups]
  
  # Add a few more keys
  props['activity'] = {gr: {} for gr in list(set([g[1] for g in groups]))}
  
  # Assign all the properties to the props dict
  for g in groups:
    df = pd.read_csv(g[0])
    df = df.drop('Unnamed: 33', 1) # Garbage
    df = df.drop('freq', 1) # These are downsampled
    df = df.dropna() # Dropna
    
    # If there are multiple clusters, add them in order
    if max(df.clust_inds) == 1: # Two clusters
      numClusts = int(max(df.clust_inds)+1)
      for ch in checks:
        for clust in range(numClusts):
          try:
            props[ch][g[1]][cells[groups.index(g)]].append(df[df['clust_inds']==clust][ch].dropna().values)
          except:
            props[ch][g[1]][cells[groups.index(g)]] = [df[df['clust_inds']==clust][ch].dropna().values]
    else: # Just one cluster
      for ch in checks:
        props[ch][g[1]][cells[groups.index(g)]] = [df[ch].dropna().values]
    # Get activity profile
    tIn, cBouts = timeInClusters(df)
    props['activity'][g[1]][cells[groups.index(g)]] = [tIn, cBouts]
  
  return props



def cellAnalysis(celltypelist, fullcsvpaths):
  """
  Organize by cell type first, then treatment 2nd.
  celltypelist: 'unknown,15827010.abf,F/I steps (sometimes ramp)'
  Generated from 
  """
  typelist, paths = [], []
  with open(celltypelist, 'r') as fIn:
    for line in fIn:
      typelist.append(line.strip().split(','))
  with open(fullcsvpaths, 'r') as fIn:
    for line in fIn:
      paths.append(line.strip())
  
  # Create the default dicts
  types = list(set([p[0] for p in typelist]))
  groups = list(set([p[2] for p in typelist]))
  checks = ['maxV', 'maxDerivV', 'maxDerivdV', 'minDerivV',
        'minDerivdV', 'preMinV', 'postMinV', 'preMaxCurveV',
        'preMaxCurveK', 'postMaxCurveV', 'postMaxCurveK',
        'height', 'repolarizationV', 'intervals', 'frequencies']
  props = {typ: {ch: {gr: {} for gr in groups} for ch in checks} for typ in types}
  # Add a few more keys
  for typ in types:
    props[typ]['activity'] = {gr: {} for gr in groups}
    props[typ]['duration'] = {gr: {} for gr in groups}
  
  # Find the matching csv files
  paths = [p for p in paths if p.split('_')[-1]=='clusters.csv'] # If it's a clusters file
  reffils = [f.split('/')[-1].split('_')[0].split('.')[0] for f in paths] # ref to cluster file
  typepaths = []
  #print(
  
  for fil in typelist:
    t_ = fil[1].split('.')[0]
    if t_ in reffils:
      typepaths.append(paths[reffils.index(t_)])
    else:
      typepaths.append('none')
  
  # Populate the dictionary
  fail, success = [], []
  print('%i (of %i) files seem to be present' %(len(typepaths)-typepaths.count('none'),
                                                len(typepaths)))
  for g in range(len(typepaths)): # This retains the order of typelist
    try:
      df = pd.read_csv(typepaths[g])
      df = df.drop('Unnamed: 33', 1) # Garbage
      df = df.drop('freq', 1) # These are downsampled
      df = df.dropna() # Dropna
      
      # If there are multiple clusters, add them in order
      if max(df.clust_inds) == 1: # Two clusters
        numClusts = int(max(df.clust_inds)+1)
        for ch in checks:
          type_ = typelist[g][0]
          group_ = typelist[g][2]
          cell_ = typelist[g][1].split('.')[0]
          for clust in range(numClusts):
            props[type_][ch][group_][cell_].append(df[df['clust_inds']==clust][ch].dropna().values)
      else: # Just one cluster
        for ch in checks:
          props[type_][ch][group_][cell_] = [df[ch].dropna().values]
      
      # Get activity profile
      tIn, cBouts = timeInClusters(df)
      props[type_]['activity'][group_][cell_] = [tIn, cBouts]
      props[type_]['duration'][group_][cell_] = df.times.iloc[-1]
      success.append(typelist[g])
    
    except:
      fail.append(typelist[g])
  
  #print(failed)
  return props, success, fail



def activityProps(pdict):
  """
  Analysis of activity (burst, tonic, etc). Props comes from 
  """
  # For baseline, establish the precentage of time each cell of each type
  act = {gr: {'burst': [], 'tonic': [], 'silent': [],
              'burstLoc': [], 'tonicLoc': []} for gr in pdict.keys()}
  
  # Populate the dict
  for group in pdict.keys():
    for cell in pdict[group]['intervals']['GapFree I=0 / Baseline recording'].keys():
      inters, timeSpent = [], []
      for clust in range(len(pdict[group]['intervals']['GapFree I=0 / Baseline recording'][cell])):
          inters.append(np.mean(pdict[group]['intervals']['GapFree I=0 / Baseline recording'][cell][clust]))
          timeSpent.append(np.mean(pdict[group]['activity']['GapFree I=0 / Baseline recording'][cell][0][clust]))
      
      # Add these percentages
      maxT = pdict[group]['duration']['GapFree I=0 / Baseline recording'][cell]
      if len(inters) > 1:
        time_sort =[x for (y,x) in sorted(zip(inters, timeSpent))]
        inter_sort = [i for i in sorted(inters)]
        act[group]['burst'].append(time_sort[0]/maxT)
        act[group]['tonic'].append(time_sort[1]/maxT)
        act[group]['silent'].append(1-(time_sort[0]+time_sort[1])/maxT)
        act[group]['burstLoc'].append(inter_sort[0])
        act[group]['tonicLoc'].append(inter_sort[1])
      else:
        act[group]['tonic'].append(timeSpent[0]/maxT)
        act[group]['tonicLoc'].append(inters[0])
        act[group]['silent'].append(1-(timeSpent[0]/maxT))
      
      # Each cell done
    # Group done
  # All groups done
  return act



  
# Make a list of props dicts for each treatment combination ??
def actSnapshot(df, where=-1, T=20000, var='intervals'):
  """
  Get a snapshot of activity for the first (where=0) or last (-1) t-milliseconds.
  """
  if type(df) is str:
    df = pd.read_csv(df)
  if var in ['counts', 'freq', 'count']:
    var_ = 'intervals'
  else:
    var_ = var
  keep, su_, g = [], 0., 0
  inters = df[var_].dropna().values
  go = list(range(len(inters)))
  if where == -1: # Reverse order, otherwise leave it
    go = go[::-1]
  while su_ < T and g < len(go):
    t_ = df[var_].values[go[g]]
    if not pd.isnull(t_):
      keep.append(t_)
      su_ = su_ + t_
    g += 1
  if len(keep) < 1:
    return None
  if var == 'intervals':
    return 1./(np.mean(keep)*1000)
  elif var in ['count', 'counts']:
    return len(keep)
  elif var == 'freq':
    return len(keep)/(float(T)/T)
  return np.mean(keep)





def byTreatment(df, keep=['GapFree', 'Pilo', 'CCh', 'ModA', 'Washout', 'MCA', 'Nico'],
                goi='OK371-GFP-Gal4', var='intervals'):
  """
  Show treatments by cell type if they contain any of 'keep' list.
  In 'keep' list, *baseline must be first!*.
  Show %-change from baseline.
  """
  linked = findExtInDF(df, ext='abf', labrow=1, refCol=0, outfile=None)
  # Make the cell dict first -- easier this way
  props, glist, clist, newl = {}, [], [], []
  for l in linked: # linked: ['2015_05_11_c2', '15511005.abf', '3.3mM calcium ModA'],
    for k in keep:
      if k in l[2] and l not in newl:
        if l[2].count('.') > 1: # Multiple files, probably 
          l[2] = '.'.join
        newl.append(l)
  print('Keeping %i (of %i) files' %(len(newl), len(linked)))
  
  # Organize by cell & treatment
  cell_tx = {}
  for l in newl:
    #print(l[2])
    if l[0] not in cell_tx.keys(): # Add new cell
      cell_tx[l[0]] = {'type': df.ix[ list(df.ix[:,0].values).index(l[0]), 5 ],
                       'genotype': df.ix[ list(df.ix[:,0].values).index(l[0]), 4 ],
                       }
    t_path = getFullPath(l[1].split('.')[0]+'_props.csv')[0]
    df_t = None
    try:
      df_t = pd.read_csv(t_path)
    except:
      print('could not load %s' %t_path)
    if df_t is not None:
      if keep[0] in l[2]: # This is baseline
        cell_tx[l[0]]['baseline'] = actSnapshot(df_t, where=0, var=var) # From the beginning
      # Check the other treatments
      for k in keep:
        if k in l[2]:
          snap = actSnapshot(df_t, where=-1, var=var) # From the end
          if snap is not None:
            cell_tx[l[0]][l[2]] = snap
          else:
            print('Could not add %s (%s)' %(l[2], t_path))
  
  for l in newl:
    if l[2] not in props.keys(): # Check each treatment
      props[l[2]] = {}
    gen_ = df.ix[ list(df.ix[:,0].values).index(l[0]), 4 ] # Get the genotype
    if gen_ not in props[l[2]].keys(): # Check each genotype
      props[l[2]][gen_] = {}
      glist.append(gen_)
    cell_ = df.ix[ list(df.ix[:,0].values).index(l[0]), 5 ] # Get the celltype
    if cell_ not in props[l[2]][gen_].keys():
      props[l[2]][gen_][cell_] = []
      clist.append(cell_)
  # return props
  
  clist, glist = list(set(clist)), list(set(glist))
  print(props.keys(), glist, clist)
  # Round out the dictionaries
  for p in props.keys():
    for g in glist:
      if g not in props[p].keys():
        props[p][g] = {}
      for c in clist:
        if c not in props[p][g].keys():
          props[p][g][c] = []
  
  # Populate the props dictionary with the property, here activity
  for ck in cell_tx.keys():
    cel_ = cell_tx[ck]
    for tx in cel_.keys():
      if tx not in ['type', 'genotype', 'baseline']:
        if var == 'count':
          try:
            props[tx] [cel_['genotype']] [cel_['type']].append(float(cel_[tx])-float(cel_['baseline']))
          except:
            props[tx] [cel_['genotype']] [cel_['type']].append(float(cel_[tx]))
        else:
          try:
            props[tx] [cel_['genotype']] [cel_['type']].append(float(cel_[tx])/float(cel_['baseline']))
          except:
            props[tx] [cel_['genotype']] [cel_['type']].append(1.)
          
  for i in list(props.keys()):
    if goi is not None:
      genoByCell({i: {goi: props[i][goi]}})
  return props



  
def burstActivity(csvfile, tryburst=True, show=True):
  """
  Show the bursting activity, or whatever.
  """
  # Clean the filename (if as abf)
  if '.abf' in csvfile:
    try: # Try for clusters first
      df = loadCSVfeatures(csvfile.split('.')[0]+'_props_clusters.csv')
      tryburst = False
    except:
      df = loadCSVfeatures(csvfile.split('.')[0]+'_props.csv')
  else:
    try: # Assume it's already a csv
      df = loadCSVfeatures(csvfile)
      if 'clusters' in csvfile:
        tryburst = False
    except: # It must be a df
      df = csvfile
      if 'clust_inds' in df.columns:
        tryburst = False
  
  # Get the bursts (if needed)
  if tryburst:
    print('Trying to find the bursts!')
    cents = getCenters(df, show=False)
    clusts = runMCMC(df, cents, show=False)
    df = assignSpikes(clusts, df, show=False, force=True)
    timeIn, cluster_bouts = timeInClusters(df, thresh=2000., show=False)
    
  # Show the bursting activity
  if 'clust_inds' not in df.columns:
    print('Could not segregate clusters!')
    return df
    
  if show:
    # First is the simple color-by-burst plot
    collist = ['blue', 'red', 'forestgreen', 'goldenrod', 'purple', 'yellowgreen',
             'skyblue', 'tomato', 'darkgray']
    for i in range(df.shape[0]):
      plt.plot([df.ix[i].times, df.ix[i].times], 
               [df.ix[i].clust_inds-1, df.ix[i].clust_inds], 
               color=collist[df.ix[i].clust_inds], linewidth=1.)
    patches = []
    for u in range(int(max(df.clust_inds))):
      patches.append(mpatches.Patch(color=collist[u],
                                    label='Cluster %i' %u))
    plt.legend(handles=patches)
    
    # Next is the burst activity patterns
    plt.figure()
    checks = ['maxDerivV', 'maxDerivdV',
              'minDerivdV', 'preMaxCurveK', 'postMaxCurveK',
              'height', 'repolarizationV', 'intervals', 'frequencies']
   # for ch in range(len(checks)):
   #   plt.subplot(2,int(len(checks)/2 +1), checks.index(ch)+1)
   #   for clust in range(max(df.clust_inds)+1): # For each cluster
   #     plotthis = df[df.clust_inds==clust][ch]
   #     plt.plot([i+clust for i in np.random.random(len(
  return




def assignToBurst(abfroot, burst, show=True, rmoutliers=True):
  """
  Show bursting activity by cell. abfroot should be without .abf ext.
  Burst should be either a df (bursttms) or a path to that df.
  """
  # Find out type of input
  if type(abfroot) is str:
    if '.' in abfroot:
      abfroot = abfroot.split('.')[0]
    if '/' not in abfroot:
      try:
        dfroot = getFullPath(abfroot+'_props_clusters.csv')[0]
        df = pd.read_csv(dfroot)
      except:
        dfroot = getFullPath(abfroot+'_props.csv')[0]
        df = pd.read_csv(dfroot)
    else:
      dfroot = abfroot
    # %print('Trying to load %s ....' %dfroot)
      df = pd.read_csv(dfroot)
  if type(burst) is str:
    if '/' not in burst:
      try:
        burst = getFullPath(burst, '/home/alex/data/misc')[0]
      except:
        burst = getFullPath(burst)[0]
      burst = pd.read_csv(burst)
  
  # Now have both as data frames
  bs_cells = [i.split('s')[0].split('_')[1] for i in burst.columns] # Make sure it's in burst df
  if abfroot not in bs_cells:
    df['in_burst'] = [False for f in range(df.shape[0])]
    return df # No bursts, just return the df
  cell_id = 'id_'+ abfroot
  start = burst[cell_id+'start'].dropna().values
  stop = burst[cell_id+'stop'].dropna().values
  
  in_burst = [] # Check if each spike belongs to a burst
  for i in range(df.shape[0]):
    t_ = df.ix[i]['times']/1000. # For each spike time
    ibs = False
    for bur in range(len(start)): # Check if it fits inside a burst!
      if start[bur] < t_ < stop[bur]:
        ibs = True
    in_burst.append(int(ibs))
  df['in_burst'] = in_burst
  
  # Now do all the plotting!
  if show:
    # First is the simple color-by-burst plot
    collist = ['blue', 'red', 'forestgreen', 'goldenrod', 'purple', 'yellowgreen',
             'skyblue', 'tomato', 'darkgray']
    for i in range(df.shape[0]):
      plt.plot([df.ix[i].times, df.ix[i].times], 
               [df.ix[i].in_burst-1, df.ix[i].in_burst], 
               color=collist[int(df.ix[i].in_burst)], linewidth=1.)
    patches = []
    labs = ['Tonic', 'Burst']
    for u in range(int(max(df.in_burst)+1)):
      patches.append(mpatches.Patch(color=collist[u],
                                    label=labs[u]))
    plt.legend(handles=patches)
    plt.ylim([-1.5, max(df.in_burst)+.5])
    
    # Next is the burst activity patterns
    plt.figure()
    checks = ['maxDerivV', 'maxDerivdV',
              'minDerivdV', 'preMaxCurveK', 'postMaxCurveK',
              'height', 'repolarizationV', 'intervals', 'frequencies']
    for ch in range(len(checks)):
      plt.subplot(2,int(len(checks)/2 +1), ch+1)
      labels = ['Tonic', 'Burst']
      
      for clust in range(max(df.in_burst)+1): # For each cluster
        plotthis = df[df.in_burst==clust][checks[ch]].values
        if rmoutliers:
          plotthis = outlier(plotthis, as_nan=False)
          plotthis = outlier(plotthis, as_nan=False)
        plt.plot([i*0.2+clust for i in np.random.random(len(plotthis))],
                 plotthis, 'o', color=collist[clust], markeredgecolor='none',
                 alpha=0.3)
        plt.plot([clust, clust+.2], [np.mean(plotthis), np.mean(plotthis)],
                 color='black', lw=2)
        plt.plot([clust+.1, clust+.1], 
                 [np.percentile(plotthis, 25), np.percentile(plotthis, 75)],
                 color='black', lw=2)
      labels = [labels[i] for i in range(max(df.in_burst)+1)]
      poses = [i+.1 for i in range(len(labels))]
      plt.xticks(poses, labels, rotation=45)
      plt.xlim([-.1, max(df.in_burst)+.3])
      plt.title(checks[ch])
    plt.show()
  
  return df

#



def burstDFhelper(tdf, temp, bs, cell_id):
  """
  Populate the temp dictionary.
  """
  def ibi_cv(bstart, bstop):
    """
    Calculate inter-burst interval coefficient of variation.
    """
    ibis = []
    for b in range(len(bstart)-1):
      if bstart[b+1] > bstop[b]: # ortho, correct
        ibis.append(bstart[b+1] - bstop[b])
      else:
        print('    In %s, %.2f starts before burst ends at %.2f' 
              %(cell_id, bstart[b+1], bstop[b]))
    return np.mean(ibis), np.std(ibis)/np.mean(ibis)
  
  def spikesperburst(tdf, bstart, bstop):
    """
    Count spikes per burst and spikes/burst CV.
    """
    tms = list(tdf.times.dropna().values)
    bursts = [[tms[u] for u in range(len(tms)) if bstart[k]<(tms[u]/1000.)<bstop[k] ]
              for k in range(len(bstart))]
    bursts = [len(i) for i in bursts]
    return np.mean(bursts), np.std(bursts)/np.mean(bursts)
  
  def burst_time(temp, bstart, bstop):
    """
    Make sure bstop[i] is always after bstart[i]; also burst length
    """
    to_sum = []
    for b in range(len(bstart)):
      if bstop[b]-bstart[b] >= 0:
        to_sum.append(bstop[b]-bstart[b])
      elif bstop[b]-bstart[b] < 0 and b == len(bstop)+1: # Make it go to end
        to_sum.append(temp['length']/1000.-bstart[b])
      else:
        pass
    return np.mean(to_sum), np.std(to_sum)/np.mean(to_sum), sum(to_sum)/(temp['length']/1000.)
  
  bs_cells = [i.split('s')[0].split('_')[1] for i in bs.columns]
  #print(cell_id, bs_cells)
  if cell_id in bs_cells:
    
    bstart = bs['id_'+cell_id+'start'].dropna().values
    bstop = bs['id_'+cell_id+'stop'].dropna().values
    temp['numbursts'] = len(bstart) # Number of bursts
    print('  --> Found %i bursts ' %temp['numbursts'])
    temp['burst_length'], temp['burst_length_cv'], \
         temp['burst'] = burst_time(temp, bstart, bstop)
    temp['spikespburst'], temp['spikespburst_cv'] = \
                                      spikesperburst(tdf, bstart, bstop)
    if temp['burst'] < 0:
      print('      Warning! Found %.4f burst time for %s!' 
            %(temp['burst'], temp['file']))
      temp['burst'] = 0.
    else:
      temp['burst'] = temp['burst']/(temp['length']/1000.) # Burst time in s!!!
    temp['ibi_length'], temp['ibi_cv'] = ibi_cv(bstart, bstop)
  else:  # Else, it doesn't burst
    temp['burst'], temp['burst_length_cv'], temp['ibi_cv'] = 0., np.nan, np.nan
  temp['tonic'] = sum(tdf[tdf.in_burst==0]['intervals'].dropna().values)/temp['length']
  temp['silent'] = 1. - (temp['burst']+temp['tonic'])
  
  return temp




def tonicDFhelper(tdf, temp, bs, cell_id, burstpresent=True):
  """
  Tonic isi and tonic isi cv
  """
  if type(tdf) is str: # Load the data frame
    if '.' in tdf:
      tdf = tdf.split('.')[0]
    if '/' not in tdf:
      try:
        dfroot = getFullPath(tdf+'_props_clusters.csv')[0]
        tdf = pd.read_csv(dfroot)
      except:
        dfroot = getFullPath(tdf+'_props.csv')[0]
        tdf = pd.read_csv(dfroot)
    else:
      dfroot = tdf
    # %print('Trying to load %s ....' %dfroot)
      tdf = pd.read_csv(dfroot)
  
  if burstpresent:
    bstart = bs['id_'+cell_id+'start'].dropna().values
    bstop = bs['id_'+cell_id+'stop'].dropna().values
  else:
    bstart, bstop = [0.], [0.000001]
  isis = []
  for t in range(tdf.shape[0]):
    tms = t/1000.
    add = True
    for b in range(len(bstart)):
      if bstart[b] < bstop[b]:
        if bstart[b] < tms < bstop[b]: # time is in burst
          add = False
      elif bstop[b] == 0:
        if bstart[b] < tms < temp['length']/1000.:
          add = False
    if add:
      isis.append(tms)
  return np.mean(isis), np.std(isis)/np.mean(isis)
  



def burstVtonic(filelengths, autodf, burstdf='/home/alex/data/misc/maria/bursts/bursts_times.csv',
                refRow=1, refCol=5):
  """
  autodf is a cell list where we'll find all abf files for GapFree (refRow)
  and matching whatever is in refCol (cell type, genotype, etc).
  burstdf is start/stop 
  """
  if type(burstdf) is str:
    bs = pd.read_csv(burstdf) # bs is the bursting data
  else:
    bs = burstdf
  if type(autodf) is str:
    auto = pd.read_csv(autodf)
  else:
    auto = autodf
  
  # Get the cell types for baseline props
  cellfiles = findExtInDF(auto, refCol=5)
  calc = pd.DataFrame(index=[f.split('/')[-1].split('.')[0] 
                             for f in filelengths.keys()],
                      columns=['file', 'cell', 'length', 'burst', 'tonic',
                               'silent', 'numbursts', 'burst_length_cv', 'ibi_cv',
                               'spikespburst', 'spikespburst_cv',
                               'burst_length', 'ibi_length',
                               'tonic_isi', 'tonic_isi_cv'])
  
  for c in cellfiles:
    if 'GapFree' in c[2]:
      for k in filelengths.keys():
        if c[1] in k:# A hit
          print('Getting df for %s ...' %c[1])
          temp = {}
          
          try:
            tdf = assignToBurst(c[1], bs, rmoutliers=True, show=False)
            # print(tdf.shape)
            for f in filelengths.keys(): # First get file length in ms
              if c[1] in f:
                temp['length'] = float(filelengths[f])
            print('Creating the dictionary for %s' %c[1])
            temp = burstDFhelper(tdf, temp, bs, c[1].split('.')[0])
            temp['tonic_isi'], temp['tonic_isi_cv'] = \
                   tonicDFhelper(tdf, temp, bs, c[1].split('.')[0])
          
          except:
            try: # Just get tonic stuff
              temp['tonic_isi'], temp['tonic_isi_cv'] = \
                       tonicDFhelper(c[1].split('.')[0], temp, 
                                     bs, c[1].split('.')[0], False)
            except:
              print('Could not create dict for %s' %c[1])
                     
          temp['cell'] = c[0]
          temp['file'] = c[1].split('.')[0]
          for k in temp.keys():
            calc[k][temp['file']] = temp[k]
  # With df made, can now analyze
  return calc



# Plot calc stuff
def plot_props(df, p1, p2, size=None, colors='cell', cinv=True, 
               axes=None, sizefactor=100, title=None):
  contcols = ['lightskyblue', 'brown', 'orange', 'springgreen',
            'fuchsia', 'tomato', 'gold', 'indigo',
            'darkslateblue', 'black', 'darkgreen', 'aqua',
            'darkorchid', 'grey', 'salmon', 'plum',
            'coral', 'sienna', 'darkkhaki', 'yellowgreen',
            'deeppink', 'ivory', 'orchid', 'lightsteelblue']
  cellcolors = [contcols[list(set(df[colors].values)).index(u)]
                for u in df[colors].values]
  if size is None:
    size=20.
  else:
    if cinv: # Invert from 0
      size=[(1-i)*sizefactor for i in df[size]]
    else:
      size=[(i)*sizefactor for i in df[size]]
  plt.scatter(df[p1].values, df[p2].values,
              color=cellcolors, s=size, alpha=0.5)
  patches = []
  for u in range(len(list(set(df[colors].values)))):
    patches.append(mpatches.Patch(color=contcols[u],
                                  label=list(set(df[colors].values))[u]))
  plt.legend(handles=patches, fontsize=15)
  if axes is not None:
    plt.xlabel(axes[0])
    plt.ylabel(axes[1])
  else:
    plt.xlabel(p1)
    plt.ylabel(p2)
  if title is not None:
    plt.title(title)
  plt.show()
  return




def df_pca(df_in, keep=None, expvar=False, rmoutliers=True, show=True,
           colorcol=None):
  """
  Run a simple PCA on the df features of keep.
  If expvar is True, a plot of explained variance is also shown.
  Heavily inspired by http://sebastianraschka.com/Articles/2015_pca_in_3_steps.html
  """
  from sklearn.preprocessing import StandardScaler
  if keep is None:
    keep = ['maxV', 'maxDerivV', 'maxDerivdV', 'minDerivV', 
            'minDerivdV', 'preMinV', 'postMinV', 'preMaxCurveK', 
            'postMaxCurveK', 'postMaxCurveV', 'preMaxCurveV', 'height', 
            'repolarizationV', 'intervals']
  # Clean the data frame
  df = df_in.copy()
  for col in df.columns:
    if col not in keep:
      df = df.drop(col, 1)
    else:
      if col != colorcol:
        df[col] = outlier(df[col].values)
  df = df.dropna()
  if colorcol is not None:
    colors = df[colorcol].values
    df = df.drop(colorcol, 1)
  # Make into np.array
  data = []
  for col in df.columns:
    temp_ = df[col]
    data.append(temp_)
  data = np.array(data).T # Make as array and transpose
  data = StandardScaler().fit_transform(data) # Standardize data
  
  # run pca (svd)
  u, eigvals, eigvecs = np.linalg.svd(data, full_matrices=False)
  eigpairs = [(np.abs(eigvals[i]), eigvecs[:,i])
              for i in range(len(eigvals))]
  eigpairs.sort()
  eigpairs.reverse()
  mat_w = np.hstack((eigpairs[0][1].reshape(eigvals.shape[0],1),
                      eigpairs[1][1].reshape(eigvals.shape[0],1)))
  Y = data.dot(mat_w) # Re-transform by matrix
  
  # Plot these data
  if show:
    contcols = ['lightskyblue', 'brown', 'orange', 'springgreen',
            'fuchsia', 'tomato', 'gold', 'indigo',
            'darkslateblue', 'black', 'darkgreen', 'aqua',
            'darkorchid', 'grey', 'salmon', 'plum',
            'coral', 'sienna', 'darkkhaki', 'yellowgreen',
            'deeppink', 'ivory', 'orchid', 'lightsteelblue']
    plt.figure()
    if colorcol is not None:
      try:
        colors = [contcols[list(set(colors)).index(u)] for u in colors]
      except:
        colors = 'blue'
    else:
      colors='blue'
    plt.scatter(Y[:,0], Y[:,1], color=colors, edgecolor='none',
                alpha=0.7)
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.tight_layout()
  
    # Explained variance
    if expvar: # eigvals come pre-sorted
      var_exp = [i/sum(eigvals)*100. for i in eigvals]
      cum_var_exp = np.cumsum(var_exp)
      #with plt.style.context('seaborn_whitegrid'):
      plt.figure()
      plt.bar(range(len(var_exp)), var_exp, alpha=0.5, align='center',
              label='individual explained variance')
      plt.step(range(len(cum_var_exp)), cum_var_exp, where='mid',
               label='cumulative explained variance')
      plt.xlabel('Principal components')
      plt.ylabel('Explained variance (\%100)') # \\%
      plt.legend(loc='best')
      plt.tight_layout()
    
    plt.show() # Show the plots
  return Y


########################################################################
# Plotting




def scatterdf(df, xfactor, yfactor, showcross=True):
  """
  xfactor is the grouping, yfactor is the measurement being plotted.
  """
  labels = list(set(df[xfactor]))
  collist = ['blue', 'red', 'forestgreen', 'goldenrod', 'purple', 'yellowgreen',
             'skyblue', 'tomato', 'darkgray']
  for i in range(df.shape[0]):
    xval = labels.index(df.ix[i][xfactor])
    plt.plot(xval+np.random.random()*.2, df.ix[i][yfactor], 'o', markersize=10,
             color=collist[xval], alpha=0.4, markeredgecolor='white')
  if showcross:
    for lab in labels:
      mn = np.mean(calc[calc[xfactor]==lab][yfactor].values)
      iqr = np.percentile(calc[calc[xfactor]==lab][yfactor].values, [25,75])
      print(iqr)
      plt.plot([labels.index(lab)-0.1, labels.index(lab)+0.3],
               [mn, mn], lw=2, color='black')#collist[labels.index(lab)]) # Mean
      plt.plot([labels.index(lab)+0.1, labels.index(lab)+0.1],
               [iqr[0], iqr[1]], lw=2, color='black')#collist[labels.index(lab)]) # IQR
  
  plt.xlim([-1, len(labels)])
  plt.xticks(range(len(labels)), labels, rotation=45)
  plt.xlabel(xfactor)
  plt.ylabel(yfactor)
  plt.show()
  return





def stackeddf(df, xfactor, yfactors, axes=None):
  """
  xfactor is grouping (i.e.: celltype), yfactors is a list of measurements.
  """
  labels = list(set(df[xfactor]))
  collist = ['blue', 'red', 'forestgreen', 'goldenrod', 'purple', 'yellowgreen',
             'skyblue', 'tomato', 'darkgray']

  for lab in labels:
    bots = 0.
    for yf in yfactors:
      mn = np.mean(df[df[xfactor]==lab][yf])
      plt.bar(labels.index(lab), mn, bottom=bots, color=collist[yfactors.index(yf)],
              edgecolor='white', alpha=0.7)
      bots = bots+mn
  
  # With all the xfactors and yfactors plotted and stacked, do cosmetics
  patches = []
  for yf in yfactors:
    patches.append(mpatches.Patch(color=collist[yfactors.index(yf)],
                                    label=yf))
  plt.legend(handles=patches)
  plt.xticks([i+0.2 for i in range(len(labels))], labels, rotation=45)
  if axes is not None:
    plt.xlabel(axes[0])
    plt.ylabel(axes[1])
  plt.show()
  return





def baselineProp(pdict, prop):
  return  {k: pdict[k][prop]['GapFree I=0 / Baseline recording'] 
                   for k in list(pdict.keys())}



def plotProp(pdict, title=None, sameax=True, showmean=True, 
             bounds=[None,None]):
  """
  Some basic plotting stuff. Collapses multiple spike types into single.
  """
  try:
    pdict.pop('all stats')
  except:
    pass
  spk, groups = [], list(pdict.keys())
  fig = plt.figure()
  c_colors = {}
  
  if sameax:
    ax = fig.add_subplot(111)
    for g in range(len(groups)):
      sofar = []
      for cell in pdict[groups[g]].keys():
        if cell not in c_colors.keys():
          c_colors[cell] = np.random.random(3)
        this = [u for u in pdict[groups[g]][cell][0]]
        if len(pdict[groups[g]][cell]) > 1:
          for sp in pdict[groups[g]][cell][1]:
            this.append(sp)
        ax.plot([i for i in np.random.normal(loc=g, scale=0.1, size=len(this))], this, 'o',
                 color=c_colors[cell], label=groups[g], alpha=0.3,
                 markeredgecolor='none', markersize=1)
        for t in this:
          sofar.append(t)
      if showmean:
        ax.plot([g-.5,g+.5], [np.mean(sofar), np.mean(sofar)],
                '--', color='black', lw=2)
    # Cosmetics
    plt.xticks(range(len(groups)), groups, rotation=30)
    plt.ylim([bounds[0], bounds[1]])
    
  else:
    plots = [fig.add_subplot(1, len(groups)+1, p) for p in range(len(groups))]
    for g in range(len(groups)):
      for cell in pdict[groups[g]].keys():
        if cell not in c_colors.keys():
          c_colors[cell] = np.random.random(3)
        this = [u for u in pdict[groups[g]][cell][0]]
        if len(pdict[groups[g]][cell]) > 1:
          for sp in pdict[groups[g]][cell][1]:
            this.append(sp)
        plots[g].plot([i+g for i in np.random.random(len(this))], this, 'o',
                 color=c_colors[cell], label=groups[g], alpha=0.3,
                 markeredgecolor='none')
  
  if title:
    plt.title(title)
  plt.show()
  return

#


def activityPlot(act):
  """
  Plots activity by cell type.
  """
  # Plot 1 is simple stacked bar
  plt.figure(figsize=(9,4), dpi=100)
  ax1 = plt.subplot(1,2,1)
  labels = [gr for gr in act.keys()]
  poses = [i+.5 for i in range(len(labels))]
  # b_means, b_stds, t_means, t_stds, s_means, s_stds = [], [], [], [], [], []
  stat = {'b_means': [], 'b_stds': [], 't_means': [], 't_stds': [],'s_means': [], 's_stds': []}
  grkey = {'b_means': 'burst', 'b_stds': 'burst', 't_means': 'tonic', 't_stds': 'tonic','s_means': 'silent', 's_stds': 'silent'}
  fnkey = {'b_means': np.mean, 'b_stds': np.std, 't_means': np.mean, 't_stds': np.std,'s_means': np.mean, 's_stds': np.std}
  
  
  for gr in labels:
    for k in stat.keys():
      try:
        temp_ = fnkey[k](act[gr][grkey[k]])
        if str(temp_) == 'nan':
          stat[k].append(0.)
        else:
          stat[k].append(temp_)
      except:
        stat[k].append(0.)
    
  p_b = ax1.bar(poses, stat['b_means'], color='blue', alpha=0.6, 
                yerr=stat['b_stds'], edgecolor='white')
  p_t = ax1.bar(poses, stat['t_means'], bottom=stat['b_means'], color='red', alpha=0.6, 
                yerr=stat['t_stds'], edgecolor='white')
  p_s = ax1.bar(poses, stat['s_means'], bottom=[stat['b_means'][i]+\
                stat['t_means'][i] for i in range(len(stat['b_means']))],
                color='purple', alpha=0.6, yerr=stat['s_stds'],
                edgecolor='white')
  # Cosmetics
  plt.xticks(poses, labels, rotation=30)
  plt.legend((p_b[0], p_t[0], p_s[0]), ('Burst', 'Tonic', 'Silent'))
  
  # Plot 2 is complex
  # ax2 = plt.subplot2grid((1,3), (0,1), colspan=2)
  ax2 = plt.subplot(1,2,2)
  for gr in range(len(labels)):
    ax2.plot(np.random.normal(loc=poses[gr], scale=.1, size=len(act[labels[gr]]['burstLoc'])), 
             act[labels[gr]]['burstLoc'], 'o', color='blue', alpha=0.6,
             markeredgecolor='none')
    ax2.plot(np.random.normal(loc=poses[gr], scale=.1, size=len(act[labels[gr]]['tonicLoc'])), 
             act[labels[gr]]['tonicLoc'], 'o', color='red', alpha=0.6,
             markeredgecolor='none')
  
  # Cosmetics
  plt.xticks(poses, labels, rotation=30)
  print(stat)
  plt.show()
  return



def dfActivity(df):
  inters = df.intervals.dropna().values
  tonic = sum([i for i in inters if i>80.])/max(df.times)
  burst = sum([i for i in inters if i<=80.])/max(df.times)
  silent = 1.-(tonic+burst)
  thing = {'burst': burst, 'tonic': tonic, 'silent': silent} 
  # print(thing)
  return thing



def ActivityPlot2(df, goi='OK371-GFP-Gal4', hidefliers=True):
  """
  sdf.
  """
  linked = findExtInDF(df, ext='abf', labrow=1, refCol=0, outfile=None)
  # Make the cell dict first -- easier this way
  glist, clist, newl = [], [], []
  for l in linked: # linked: ['2015_05_11_c2', '15511005.abf', '3.3mM calcium ModA'],
    for k in keep:
      if k in l[2] and l not in newl:
        if l[2].count('.') > 1: # Multiple files, probably 
          l[2] = '.'.join
        newl.append(l)
  print('Keeping %i (of %i) files' %(len(newl), len(linked)))
  
  # Organize by cell & treatment
  cell_tx = {}
  for l in newl:
    #print(l[2])
    if l[0] not in cell_tx.keys(): # Add new cell
      cell_tx[l[0]] = {'type': df.ix[ list(df.ix[:,0].values).index(l[0]), 5 ],
                       'genotype': df.ix[ list(df.ix[:,0].values).index(l[0]), 4 ],
                       }
    t_path = getFullPath(l[1].split('.')[0]+'_props.csv')[0]
    df_t = None
    try:
      df_t = pd.read_csv(t_path)
    except:
      print('could not load %s' %t_path)
      df_t = None
    if df_t is not None:
      # Check the other treatments
      if 'GapFree' in l[2]:
        snap = dfActivity(df_t) # From the end
        if snap is not None:
          cell_tx[l[0]][l[2]] = snap
        else:
          print('Could not add %s (%s)' %(l[2], t_path))
  
  props = {}
  for l in newl:
    if l[2] not in props.keys(): # Check each treatment
      props[l[2]] = {}
    gen_ = df.ix[ list(df.ix[:,0].values).index(l[0]), 4 ] # Get the genotype
    if gen_ not in props[l[2]].keys(): # Check each genotype
      props[l[2]][gen_] = {}
      glist.append(gen_)
    cell_ = df.ix[ list(df.ix[:,0].values).index(l[0]), 5 ] # Get the celltype
    if cell_ not in props[l[2]][gen_].keys():
      props[l[2]][gen_][cell_] = []
      clist.append(cell_)
  # return props
  
  clist, glist = list(set(clist)), list(set(glist))
  print(props.keys(), glist, clist)
  # Round out the dictionaries
  for p in props.keys():
    for g in glist:
      if g not in props[p].keys():
        props[p][g] = {}
      for c in clist:
        if c not in props[p][g].keys():
          props[p][g][c] = []
  
  # Populate the props dictionary with the property, here activity
  for ck in cell_tx.keys():
    cel_ = cell_tx[ck]
    props['GapFree I=0 / Baseline recording'][cel_['genotype']] [cel_['type']].append(cel_['GapFree I=0 / Baseline recording'])
  print(props)
  # print(props['GapFree I=0 / Baseline recording'])
  genoBar({'Baseline': props['GapFree I=0 / Baseline recording']})
  return props

#


def simpleActivityPlot(df, goi='OK371-GFP-Gal4'):
  """
  Given only a single condition (gap free) (genotypes) (cell types)
  """
  linked = findExtInDF(df, ext='abf', labrow=1, refCol=0, outfile=None)
  # Make the cell dict first -- easier this way
  props, glist, clist, newl = {}, [], [], []
  for l in linked: # linked: ['2015_05_11_c2', '15511005.abf', '3.3mM calcium ModA'],
    if 'GapFree' in l[2] and l not in newl:
      if l[2].count('.') > 1: # Multiple files, probably 
        l[2] = '.'.join
      newl.append(l)
  print('Keeping %i (of %i) files' %(len(newl), len(linked)))
  
  # Organize by cell & treatment
  cell_tx = {}
  for l in newl:
    if l[0] not in cell_tx.keys(): # Add new cell

      cell_tx[l[0]] = {'type': df.ix[ list(df.ix[:,0].values).index(l[0]), 5 ],
                       'genotype': df.ix[ list(df.ix[:,0].values).index(l[0]), 4 ],
                       }
    df_t, path_t= None, getFullPath(l[1].split('.')[0]+'_props.csv')[0]
    try:
      df_t = pd.read_csv(path_t)
    except:
      print('could not load %s' %path_t)
    if df_t is not None:
      inters = df_t.intervals.dropna().values
      tonic = sum([i for i in inters if i>80.])/max(df_t.times)
      burst = sum([i for i in inters if i<=80.])/max(df_t.times)
      silent = 1.-(tonic+burst)
      cell_tx[l[0]]['tonic'] = tonic
      cell_tx[l[0]]['burst'] = burst
      cell_tx[l[0]]['silent'] = silent
    else:
      print('no data found for %s' %path_t)
  
  # Condition the plots
  glist = list(set([cell_tx[k]['genotype'] for k in cell_tx.keys()]))
  clist = list(set([cell_tx[k]['type'] for k in cell_tx.keys()]))
  bygene = {g: {c: {'burst': [], 'tonic': [], 'silent': []} for c in clist} for g in glist}
  for c in cell_tx.keys():
    cel_ = cell_tx[c]
    bygene[cel_['genotype']][cel_['type']]['burst'].append(cel_['burst'])
    bygene[cel_['genotype']][cel_['type']]['tonic'].append(cel_['tonic'])
    bygene[cel_['genotype']][cel_['type']]['silent'].append(cel_['silent'])
  
  # Now plot these as a box
  wid = 0.5
  ind = [[i+wid*k for i in np.arange(len(clist))] for k in range(len(glist))]
  collist = ['blue', 'red', 'forestgreen', 'goldenrod', 'purple', 'yellowgreen',
             'skyblue', 'tomato', 'darkgray']
  print(glist)
  for g in [2]: #len(glist)
    # Burst
    plt.bar(ind[g], [np.mean([bygene[glist[g]][c]['burst']]) for c in clist],
            color=collist[0], edgecolor='white', width=wid, alpha=0.6,
            yerr=[np.std([bygene[glist[g]][c]['burst']]) for c in clist],)
    plt.bar(ind[g], [np.mean([bygene[glist[g]][c]['tonic']]) for c in clist],
            color=collist[1], edgecolor='white', width=wid, alpha=0.6,
            yerr=[np.std([bygene[glist[g]][c]['tonic']]) for c in clist],
            bottom=[np.mean([bygene[glist[g]][c]['burst']]) for c in clist])
    plt.bar(ind[g], [np.mean([bygene[glist[g]][c]['silent']]) for c in clist],
            color=collist[2], edgecolor='white', width=wid, alpha=0.6,
            yerr=[np.std([bygene[glist[g]][c]['silent']]) for c in clist],
            bottom=[x+y for x,y in zip([np.mean([bygene[glist[g]][c]['burst']]) for c in clist],
                                        [np.mean([bygene[glist[g]][c]['tonic']]) for c in clist])],
            )
  
  # Cosmetics
  plt.xticks(ind[0], clist, rotation=30)
  plt.title('Baseline activity')
  plt.ylabel('Percent activity')
  patches = []
  for u in [['burst', 'blue'], ['tonic', 'red'], ['silent', 'green']]:
    patches.append(mpatches.Patch(color=u[1],
                                    label=u[0]))
  plt.legend(handles=patches)
  plt.show()
  return

#


def simpleBar(props, goi=['OK371-GFP-Gal4']):
  """
  """
  cells = {g: {} for g in goi}
  actions = ['burst', 'tonic', 'silent']
  for g in goi:
    for celltype in props[g].keys():
      cells[g][celltype] = {'burst': [], 'tonic': [], 'silent': []}
      for c_ in props[g][celltype]: # List of dicts
        for act in actions:
          try:
            cells[g][celltype][act].append(c_[act])
          except:
            pass
  
  # Plotting
  collist = ['blue', 'red', 'forestgreen', 'goldenrod', 'purple', 'yellowgreen',
             'skyblue', 'tomato', 'darkgray']
  wid=0.3
  ind = [[i+wid*k for i in list(range(len(cells[goi[0]].keys())))] for k in range(len(goi))]
  for g in goi:
    for c in cells[g].keys():
      sofar = 0.
      for act in range(len(actions)):
        this = [i for i in cells[g][c][actions[act]] if not pd.isnull(i)]
        plt.bar(ind[goi.index(g)][list(cells[g].keys()).index(c)], np.mean(this),
                #yerr=np.std(this), 
                bottom=sofar, width=wid-.05,
                color=collist[act], edgecolor='white')
        sofar = sofar + np.mean(this)
      
      # If multiple genes, indicate with dashed lines
      if len(goi) > 1:
        plt.axvline(ind[goi.index(g)][list(cells[g].keys()).index(c)]+wid/2.,
                    0, 1, linestyle='--', lw=2., color=collist[len(collist)-1-goi.index(g)],
                    alpha=0.5)
  
  patches = []
  for u in [['burst', 'blue'], ['tonic', 'red'], ['silent', 'green']]:
    patches.append(mpatches.Patch(color=u[1],
                                    label=u[0]))
  if len(goi) > 1:
    for g in goi:
      patches.append(mpatches.Patch(color=collist[len(collist)-1-goi.index(g)],
                                    label=g))
  plt.legend(handles=patches, fontsize=15)
  plt.xticks([i+.5 for i in ind[0]], list(cells[g].keys()), rotation=30, fontsize=15)
  plt.ylabel('% activity', fontsize=15)
  plt.ylim([0,1.2])
  plt.show()
  return

#


def genoBar(props, hidefliers=True):
  """
  Bar.
  """
  pkeys = list(props.keys())
  gkeys = list(props[pkeys[0]].keys())
  ckeys = list(props[pkeys[0]][gkeys[0]].keys())
  actions = ['burst', 'tonic', 'silent']
  fig = plt.figure()
  axs = [fig.add_subplot(len(props.keys()), 1, p+1) for p in range(len(pkeys))]
  wid = 0.5
  ind = [[i+wid*k for i in np.arange(len(ckeys))] for k in range(len(gkeys))] # Plot specifics
  print(gkeys)
  collist = ['blue', 'red', 'forestgreen', 'goldenrod', 'purple', 'yellowgreen',
             'skyblue', 'tomato', 'darkgray']
  
  for p in range(len(pkeys)): # For each property
    lns = []
    for g in range(len(gkeys)): # For each genotype
      #bars.append(axs[p].boxplot(ind[g], [np.mean(props[pkeys[p]][gkeys[g]][c]) for c in ckeys],
      #                       # color=['blue', 'red', 'forestgreen'][g], # wid,
      #                       #yerr=[np.std(props[pkeys[p]][gkeys[g]][c]) for c in ckeys],
      #                       ))#alpha=0.7))
      for c in range(len(ckeys)):
        try:
          for act in range(len(actions)):
            this = [props[pkeys[p]][gkeys[g]][ckeys[c]][i][actions[act]]
                    for i in range(len(props[pkeys[p]][gkeys[g]][ckeys[c]]))]
            print(this)
            if hidefliers:
              this = outlier(this, as_nan=False)
            else:
              pass
            sofar = 0.
            axs[p].bar(ind[g], np.mean(this), yerr=np.std(this),
                        color=collist[act], bottom=sofar,
                        edgecolor='none', alpha=0.5)
            sofar = sofar + np.mean(this)
        except:
          pass
    
    # Plot cosmetics
    plt.ylabel(pkeys[p], fontsize=15)
    plt.xlim([-.25, len(ckeys)+.25])
    plt.sca(axs[p])
    plt.xticks(ind[g], ckeys, rotation=30)
    print(ind[g])
  
  if p == len(pkeys)-1:
    patches = []
    for g in range(len(actions)):
      patches.append(mpatches.Patch(color=collist[g],
                                    label=actions[g]))
    plt.legend(handles=patches) #lns, gkeys
  
  plt.show()
  return 



def genoByCell(props, hidefliers=True, shaders=None):
  """
  Grouped scatter plot.
  """
  pkeys = list(props.keys())
  gkeys = list(props[pkeys[0]].keys())
  ckeys = list(props[pkeys[0]][gkeys[0]].keys())
  fig = plt.figure()
  axs = [fig.add_subplot(len(props.keys()), 1, p+1) for p in range(len(pkeys))]
  wid, ind = 0.25, [[i+.25*k for i in np.arange(len(ckeys))] for k in range(len(gkeys))] # Plot specifics
  print(gkeys)
  collist = ['blue', 'red', 'forestgreen', 'goldenrod', 'purple', 'yellowgreen',
             'skyblue', 'tomato', 'darkgray']
  
  for p in range(len(pkeys)): # For each property
    lns = []
    for g in range(len(gkeys)): # For each genotype
      #bars.append(axs[p].boxplot(ind[g], [np.mean(props[pkeys[p]][gkeys[g]][c]) for c in ckeys],
      #                       # color=['blue', 'red', 'forestgreen'][g], # wid,
      #                       #yerr=[np.std(props[pkeys[p]][gkeys[g]][c]) for c in ckeys],
      #                       ))#alpha=0.7))
      for c in range(len(ckeys)):
        if hidefliers:
          this = outlier(props[pkeys[p]][gkeys[g]][ckeys[c]], as_nan=False)
        else:
          this = props[pkeys[p]][gkeys[g]][ckeys[c]]
        axs[p].plot(np.random.normal(ind[g][c],0.05, len(this)),
                    this, 'o', 
                    color=collist[g],
                    markeredgecolor='none', alpha=0.5)
        axs[p].plot([ind[g][c]-0.1, ind[g][c]+.1],
                    [np.mean(props[pkeys[p]][gkeys[g]][ckeys[c]]),
                    np.mean(props[pkeys[p]][gkeys[g]][ckeys[c]])],
                    lw=1, alpha=0.6, 
                    color=collist[g],
                    label=gkeys[g])
        if len(props[pkeys[p]][gkeys[g]][ckeys[c]]) == 0:
          maxH = 20
        else:
          maxH = max(props[pkeys[p]][gkeys[g]][ckeys[c]])
        axs[p].axvline(ind[g][c], 0, maxH, alpha=0.2,
                       linestyle='--', color=collist[g])
    
    # Plot cosmetics
    axs[p].set_ylabel(pkeys[p], fontsize=15)
    axs[p].set_xlim([-.25, len(ckeys)+.25])
    plt.sca(axs[p])
    plt.xticks(ind[g], ckeys, rotation=30)
    print(ind[g])
  
  if p == len(pkeys)-1:
    patches = []
    for g in range(len(gkeys)):
      patches.append(mpatches.Patch(color=collist[g],
                                    label=gkeys[g]))
    plt.legend(handles=patches) #lns, gkeys
  
  plt.show()
  return 







########################################################################
# Program control sequence


# If getting as a csv
def csvSpikes(fname, show=False):
  """
  Run most of the above spike sorting stuff for a csv file.
  """
  df = loadCSVfeatures(fname, rmoutliers=True)
  cents = getCenters(df, show=show)
  clusts = runMCMC(df, cents, show=show)
  df = assignSpikes(clusts, df, show=show, force=True)
  timeincluster, _ = timeInClusters(df, show=show)
  # Give it a new name to differentiate the sorted spikes
  newn = fname.split('.')[0]+'_clusters.csv'
  writeNewDF(df, fname, newn)
  return df
  






########################################################################

if __name__ == "__main__":
  args = sys.argv
  print('Module should be used interactively')










