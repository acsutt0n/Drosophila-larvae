#################
# dfTrace

# Return a string of properties from a given data frame file


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from maria import *


"""
Many features of a trace are returned from the data frame.
The length of the trace, in sec, must be passed as well as the df file
if the length of the trace is not also a column in the df (lab='length')
"""
class dfTrace():
  
  def __init__(self, dfile, tracelength=None):
    # Initialize variables
    self.df = pd.read_csv(dfile)
    if tracelength is not None:
      self.length = float(tracelength)
    else:
      try:
        self.length = np.mean(df.length.dropna().values())
      except:
        print('No column named length in data frame and no tracelength given!')
        return
    self.items = {'burst': None, 'tonic': None, 'silent': None,
                  'numbursts': None, 'ibi_cv': None, 'freq': None,
                  
    # Make sure bursts are identified
    
    # Run operations
    
    
  def percentActivity(self):
    





