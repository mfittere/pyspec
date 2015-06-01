import numpy as _np
import doros

class spec():
  """class to do spectral analysis with wrapper to read in:
  - getdoros: orbit data from DOROS BPMs
  - getadt:   orbit data from ADT
  """
  def __init__(self,b1h=[],b1v=[],b2h=[],b2v=[]):
    self.data={'b1h':_np.array(b1h),'b1v':_np.array(b1v),'b2h':_np.array(b2h),'b2v':_np.array(b2v)}
  @classmethod
  def getdoros(cls,fn=''):
   return doros.getdata(fn)
