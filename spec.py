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

def psd(data,nfft=None,n0=0,window=_np.hanning,fs=11245.0):
    """calculate the PSD [m**2/Hz]. For correct normalisation
    see CERN-THESIS-2013-028 and scipy.signal.welch, which
    normalizes by the window function:
      psd=1/(fs*sum(win)**2)*abs(fft)**2
    window = window function
    nfft   = number of data points used for FFT
    n0     = use data[n0:n0+nfft]
    """
    if(nfft==None):
      nfft=len(data)-n0
    ff=_np.arange(nfft/2+1)*fs/(nfft-1)
    if window==None:
      data=data[n0:n0+nfft]
      scale=nfft**2
    else:
      win=window(nfft)
      data=data[n0:n0+nfft]*win
      scale=win.sum()**2
    psd=1/(fs*scale)*_np.abs(_np.fft.rfft(data,n=nfft,axis=0))**2
    return ff,psd
