# Stuff for Maria


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




def loadCSVfeatures(fname, rmoutliers=False):
  """
  Load a csv, keep only the features (and times), dropna, return df.
  """
  df = pd.read_csv(fname)
  checks = ['maxV', 'maxDerivV', 'maxDerivdV', 'minDerivV',
          'minDerivdV', 'preMinV', 'postMinV', 'preMaxCurveV',
          'preMaxCurveK', 'postMaxCurveV', 'postMaxCurveK', 'times',
          'height', 'repolarizationV', 'intervals', 'frequencies',
          'clust_inds', 'mslength']
  for col in df.columns:
    if col not in checks:
      df = df.drop(col, 1)
  df = df.dropna()
  
  if rmoutliers: # Treat outliers
    for col in checks:
      if col != 'times':
        df[col] = outlier(df[col])
  
  df = df.dropna()
  return df
  



def containsSpikes(filelist):
  """
  Given a txt file which is a list of .csv files, this returns the
  files that have valid spikes (>= 1 spike) and the files that do not.
  """
  fnames, numlines = [], []
  # Load the filenames
  with open(filelist, 'r') as fIn:
    for line in fIn:
      if line:
        fnames.append(line.split(None)[0])
  
  # Run through the filenames and get the number of lines
  for nam in fnames:
    # runstr = 'wc -l %s' %nam
    p = subprocess.Popen(['wc', '-l', nam],
        stdout=PIPE, stderr=PIPE, stdin=PIPE)
    out = p.stdout.read()
    try:
      numlines.append(int(out.split(None)[0]))
    except:
      numlines.append(0)
  # print(numlines[:10])
  
  # Find the min (as long as it's not zero)
  minlines = np.inf
  for n in range(len(numlines)):
    if numlines[n] < minlines and numlines[n] != 0:
      minlines = numlines[n]
  # And replace that and min with nan, export the resulting list
  newlines = [1 if i > minlines else 0 for i in numlines]
  # print(newlines[:10])
  
  return [fnames[u] for u in range(len(newlines)) if newlines[u] > 0]



def showprofile(csvfile, color='rand'):
  """
  Show all properties of the data frame. The fields listed below
  are ignored.     
  """
  ignore = ['n1List', 'n2List', 'maxVtms', 'maxVinds', 'maxDerivtms', 
     'maxDerivinds', 'minDerivtms', 'minDerivinds', 'preMintms', 'preMininds',
     'postMintms', 'postMininds', 'preMaxCurvetms', 'preMaxCurveinds',
     'postMaxCurvetms', 'postMaxCurveinds', 'times']
  f = pd.read_csv(csvfile)
  columns = [col for col in f.columns if col not in ignore and len(col.split(None)) == 1]
  ncol = int(len(columns)/2.) + 1
  if color is 'rand':
    color = np.random.random(3)
  
  # Plot these mofos
  fig = plt.figure()
  plots = [fig.add_subplot(2, ncol, i+1) for i in range(len(columns))]
  for col in range(len(columns)):
    try:
      plots[col].hist(f[columns[col]].dropna(), bins=50, facecolor=color,
                      edgecolor='none', alpha=0.5)
    except:
      print(columns[col])
    plots[col].set_title(columns[col])
    plots[col].set_ylim()
  plt.show()
  return



def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    from math import factorial
    
    try:
      window_size = np.abs(np.int(window_size))
      order = np.abs(np.int(order))
    except ValueError: #, msg:
      raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
      raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
      raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')



def getCenters(df, show=False):
  """
  Get the bin centers 
  """
  def whichPeaks(trace):
    """Find the peaks for the dist."""
    peaks = []
    df = np.diff(trace)
    for t in range(len(df)-4):
        if df[t] > 0 and df[t+1] > 0:
            if df[t+2] < 0 and df[t+3] < 0: # Potential peak
                if trace[t+2] > np.mean(trace):
                    peaks.append([t+2, trace[t+2]])
    return peaks
  
  # Get the interval data and bin it
  int_data = df.intervals
  hist, bin_e = np.histogram(int_data, bins=50)
  bin_cents = (bin_e[:-1]+bin_e[1:])*.5
  
  # Smooth the data and identify the probable peaks
  histhat = savitzky_golay(hist, 11, 3)
  pks = whichPeaks(histhat)
  print('Found peaks at :')
  cents = [bin_cents[p[0]] for p in pks]
  # print(cents)
  
  if show:
    plt.bar(bin_cents, hist, color='blue', edgecolor='white',
            alpha=0.4)
    plt.plot(bin_cents, histhat, color='blue', lw=2)
    for c in cents:
      plt.axvline(c, 0, max(hist), '--', color='red', lw=1)
    plt.show()
  
  return cents
  


def runMCMC(df, cents, show=False):
  """
  Run the MCMC algo for as many centers as needed
  """
  if type(cents) is not list:
    cents = [cents]
  numCents = len(cents)
  p = None
  
  # Tau = the precision of the normal distribution (of the above peaks)
  taus = 1. / pm.Uniform('stds', 0, 100, size=numCents)**2 # tau = 1/sigma**2
  centers = pm.Normal('centers', cents, [0.0025 for i in cents],
                      size=numCents)
  
  if numCents == 2: # Assignment probability
    p = pm.Uniform('p', 0, 1)
    assignment = pm.Categorical('asisgnment', [p, 1-p],
                                size=len(df.intervals))
    @pm.deterministic
    def center_i(assignment=assignment, centers=centers):
      return centers[assignment]
    @pm.deterministic
    def tau_i(assignment=assignment, taus=taus):
      return taus[assignment]
    observations = pm.Normal('obs', center_i, tau_i, value=df.intervals,
                             observed=True)
    # Create the model 2 peaks
    mcmc = pm.MCMC([p, assignment, observations, taus, centers])
    
  else:
    observations = pm.Normal('obs', value=df.intervals, observed=True)
    mcmc = pm.MCMC([observations, taus, centers]) # Create model, 1 peak
  
  # Run the model
  mcmc.sample(50000)
  center_trace = mcmc.trace("centers")[:]
  try:
    clusts = [center_trace[:,i] for i in range(numCents)]
  except:
    clusts = [center_trace]
  
  if show:
    for i in range(numCents):
      plt.hist(center_trace[:,i], bins=50, histtype='stepfilled',
               color=['blue', 'red'][i], alpha=0.7)
    plt.show()
  
  print('Evolved clusters at:')
  print([np.mean(c) for c in clusts])
  return clusts



def assignSpikes(clusts, df, show=False, force=True):
  """
  Assign each spiking event to either cluster 1 or cluster 2.
  """
  if 'clust_inds' in df.columns and force is False:
    print('Data frame already contains clust_inds')
    return
    
  def assignTms(clusts, tms):
    # Assign a delta_tms to cluster1 or cluster2
    assns = [abs(np.mean(clusts[c])-tms) for c in range(len(clusts))]
    return assns.index(min(assns))
  
  # Assign each spike time to a cluster
  clust_tms = [ [] for c in clusts]
  for t in range(len(df.times)-1):
    t_clust = assignTms(clusts, df.times.values[t+1]-df.times.values[t])
    clust_tms[t_clust].append(df.times.values[t])
  
  # Group spikes from same spike type together
  type_tms = []
  for c in range(len(clust_tms)):
    for t in clust_tms[c]:
      type_tms.append([t, c]) # [spk tms, clust index]
  
  # Group these together 
  clust_id = []
  for i in range(df.shape[0]):
    if df.iloc[i].times in [k[0] for k in type_tms]:
      clust_id.append(type_tms[[k[0] for k in type_tms].index(df.iloc[i].times)][1])
    else: # Not matching spike found -- happens w/ isolated spikes
      clust_id.append(np.nan)
  
  df['clust_inds'] = clust_id
  print([clust_id.count(j) for j in list(set(clust_id))], list(set(clust_id)))
  if show: # Show the cluter spikes
    for c in range(max(clust_id)+1): # Plot cluster spikes individually
      temp_spikes = df[df['clust_inds']==c]['times']
      plt.plot(temp_spikes, [c+1 for i in temp_spikes], '|', 
               color=['blue', 'red'][c])
    plt.ylim([0,3])
    plt.show()
  
  return df



def timeInClusters(df, thresh=2000., show=False):
  """
  thresh (ms): how much time should elapse before a bout is 
  automatically ended.
  """
  clust_tms = [ list(df[df.clust_inds==i]['times']) for i in 
                range(int(max(df.clust_inds))+1) ]
  cluster_bouts = [ [] for c in clust_tms ]
  for c in range(len(clust_tms)): # Each new cluster
    on = False
    for t in range(len(clust_tms[c])-1):
      if clust_tms[c][t+1] - clust_tms[c][t] > thresh: # 2-s threshold
        if on: # If this is a bout, end it
          cluster_bouts[c].append(clust_tms[c][t])
          on = False
      # Else, ignore it (too sparse for tonic)
      else: # Close enough to be a continuation or new bout
        if on is False: # Not an active bout, so start it
          cluster_bouts[c].append(clust_tms[c][t])
          on = True
          # Else, ignore it (continue the bout)
    if on: # If a bout is active at the end, end it with
        cluster_bouts[c].append(clust_tms[c][-1])

  # Calculate time spent in each cluster (assumes unit is ms)
  timeIn = [sum([cluster_bouts[c][2*i+1]-cluster_bouts[c][2*i]
                 for i in range(len(cluster_bouts[c])/2)])
            for c in range(len(cluster_bouts))]
  percentIn = [i/max([max(cl) for cl in clust_tms]) for i in timeIn]
  for t in timeIn: # Time is dependent on dT
    print('Time (percent) spent in cluster %i: %.3f s (%.2f)'
          %(timeIn.index(t), t/10000., percentIn[timeIn.index(t)]))
  
  
  if show: # Show plots
    for clust in range(len(cluster_bouts)):
      for b in range(len(cluster_bouts[clust])/2):
        plt.plot([cluster_bouts[clust][b*2], cluster_bouts[clust][b*2+1]],
                 [clust+1, clust+1], lw=5, color=['blue', 'red'][clust])
    plt.ylim([0,len(cluster_bouts)+1])
    plt.show()
  
  print(timeIn, cluster_bouts)
  return timeIn, cluster_bouts
  


#########################################################################
# Feature-specific stuff, outliers

def showVs(df, feat1, feat2):
  """
  Show some sample features.
  """
  colors = ['blue', 'red', 'green', 'coral']
  for u in range(len(cBouts)):
    plt.plot(f[f['clust_ind'] == u][feat1],
             f[f['clust_ind'] == u][feat2], 'o', color=colors[u],
             alpha=0.6, markeredgecolor='none')
    plt.xlabel(feat1)
    plt.ylabel(feat2)
  plt.show()
  return



def discrim(a, b):
  return (np.mean(a)-np.mean(b))/ \
         np.sqrt(0.5*(np.var(a)**2 + np.var(b)**2))



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


#










        
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




  




#########################################################################
# Rudimentary plotting




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





#

#########################################################################
# Scripts not related to spikes/clustering, just housekeeping

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
          try: # Could be multiple files
            temp_ = ''.join([u for u in temp_ if u != ' '])
            csplit = temp_.split(',')
            for ent in csplit:
              fils.append([df.iloc[i][df.columns[refCol]], ent, df[c].values[labrow]])
          except:
            fils.append([df.iloc[i][df.columns[refCol]], temp_, df[c].values[labrow]]) # there was only one file
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



### Sample usage:
"""
In [103]: list(calc.burst).index(max(calc.burst))
Out[103]: 58

In [104]: calc.ix[58]
Out[104]: 
file             15722003
cell                   1b
length             120000
burst        0.0009435492
tonic          0.02137228
silent          0.9776842
numbursts              17
burstfreq       0.1416667
Name: 15722003, dtype: object

In [105]: tdf = assignToBurst(calc.file['15722003'], bs, True, True)

"""

