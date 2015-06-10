"""module with various fft functions
"""
import numpy as _np
import matplotlib.pyplot as _pl

def rfft_avg(x,NFFT=256,NAVG=None,window=None,noverlap=0):
  """average the FFT over several samples
  noverlap - not yet implemented -> noverlap=0

  Call signature::

    rfft_avg(x,NFFT=256,NAVG=None,window=None,noverlap=0)

  Compute and plot the average fft of data in *x*. Data are split into
  *NFFT* length segments and the spectrum of each section is
  computed.  The windowing function *window* is applied to each
  segment, and the amount of overlap of each segment is
  specified with *noverlap*. The average over *NAVG* spectra is then
  computed. If *NAVG* is None, the maximum number of 
  averages = len(*x*)/NFFT is computed.

  *x*: 1-D array or sequence
       Array or sequence containing the data
  *window*: callable or ndarray
    A function or a vector of length *NFFT*.
    default is None, then no window function is applied
    If a function is passed as the
    argument, it must take a data segment as an argument and
    return the windowed version of the segment.
  """
  navgmax=(len(x)-NFFT)/(NFFT-noverlap)
  if(navgmax<NAVG):
    print 'WARNING: maximum number of averages:   %4.0f'%(navgmax)
    print '         number of averages requestes: %4.0f'%(NAVG)
    print '         -> compute %4.0f avgerages'%(navgmax)
    NAVG=navgmax
  if NAVG==None:
    NAVG=navgmax
    print 'NAVG= %4.0f'%(NAVG)
  lfft=[]
  if(window==None):#no window
    ww=_np.ones(NFFT)
  elif(type(window)==_np.ndarray):#self defined window given as array
    ww=window
  else:#window function
    ww=window(NFFT)
  for idx in range(NAVG):
    lfft.append(_np.fft.rfft(x[idx*(NFFT-noverlap):idx*(NFFT-noverlap)+NFFT]*ww))
  return _np.mean(_np.array(lfft),axis=0)

