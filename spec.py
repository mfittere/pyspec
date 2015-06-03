import numpy as _np
from doros import doros

class spec():
  """class to do spectral analysis with wrapper to read in:
  - getdoros: orbit data from DOROS BPMs
  - getadt:   orbit data from ADT
  """
  def __init__(self,data):
    self.data=_np.array(data)
    
