"""
Burst/tonic spike sorting algorithm.

The algorithm is a Bayesian burst/tonic distinguisher that makes use of
interspike/interburst interval clustering. It works best when both
tonic and bursting events are present. If only tonic OR bursting
events are present, it is better to simply assign every event to
cluster 0 (bursts) or cluster 1 (tonic).

ISIs are binned, smoothed using a Savitzky-Golay filter, and
then Markov chain Monte Carlo sampling is done with the original
intervals to assign them to one cluster or another.

The input is a csv file that will be a pandas data frame.

At the moment (5/2016) this hasn't been adapted for Python 3.x, but
that shouldn't be too hard.
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



def sort_spikes(csvfile):
  """
  Sort the spikes using MCMC sampling.
  """
  df = loadCSVfeatures(csvfile) # Load the csv file's features
  cents = getCenters(df)
  clusters = runMCMC(df, cents) # Don't show the results
  newdf = assignSpikes(clusters, df)
  return newdf



def write_csv(df, infile=None, outfile=None):
  """
  Write the new file to the same location but a slightly different name
  so we know the spike sorting already occured.
  """
  if outfile is None:
    if infile is None:
      outCsv = "temp.csv"
    else:
      outCsv = infile.split('.')[0]+'_clusters.csv'
  else:
    outCsv = outfile
  df.to_csv(outCsv, index=False)
  print('DataFrame written to %s' %outCsv)
  return





#########################################################################

if __name__ == "__main__":
  args = sys.argv
  if len(args) < 2:
    print("Need a csv file, else run interactively")
  elif len(args) >= 2:
    csvIn = args[1]
  if len(args) > 2:
    csvOut = args[2]
  else:
    csvOut = None
  try:
    dfOut = sort_spikes(csvIn)
    write_csv(dfOut, csvIn, csvOut)
  except:
    print("Could not sort spikes for some reason")







