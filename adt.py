"""class to import data from ADT"""
import numpy as _np
import matplotlib.pyplot as _pl
import scipy.io as sio
import spec as spec
#import os
#import cPickle as pickle
#import glob as glob
#from scipy.signal import welch

def get_fn_data(fn):
  """get time stamp, plane and bunch number
  from filename fn"""
  fn=(fn.split('/')[-1]).split('.mat')[0]#remove .mat file ending
  idxbunch=int(fn.split('bunch')[1])
  #convert time stamp to timber input format yyyy-mm-dd hh:mm:ss.sss
  date=fn.split('-')[0]
  time=fn.split('-')[1]
  date='-'.join([date[0:4],date[4:6],date[6:8]])
  time=':'.join([time[0:2],time[3:5],time[6:8]+'.000'])
  timestamp=date+' '+time
  #get position and plane
  planepos=fn.split('-')[-2]#e.g. 'ADTBposHorQ9LB1'
  pos=planepos[-5:-2]#position of ADT: Q9[LR],Q7[LR]
  beam=planepos[-2:].lower()#beam
  dic={'Hor':'h','Ver':'v'}
  plane=beam+dic[planepos.split('pos')[1].split('Q')[0]] 
  return idxbunch,timestamp,pos,plane

class adt():
  """class to import data from ADT"""
  def __init__(self,idxbunch,timestamp,pos,plane,fs=11245.0,betastar=0,dd=[]):
    self.idxbunch  = idxbunch
    self.timestamp = timestamp
    self.pos       = pos
    self.plane     = plane
    self.fs        = fs
    self.betastar  = betastar
    self.data={plane:_np.array(dd)}
  @classmethod
  def getdata(cls,fn):
    idxbunch,timestamp,pos,plane = get_fn_data(fn)
    dd   = sio.loadmat(fn)
    data = ((sio.loadmat(fn))['Mem1']).flatten()
    return cls(idxbunch,timestamp,pos,plane,data)
  def psd(self,nfft=None,n0=0,window=_np.hanning):
    """calculate the PSD [m**2/Hz]. For correct normalisation
    see CERN-THESIS-2013-028
    window = window function
    nfft   = number of data points used for FFT
    n0     = use data[n0:n0+nfft]
    for the psd the average value is subtracted from the data:
      data=data-mean(data)
    """
    xx=self.data[self.plane]
    xx=xx-_np.mean(xx)#substract average value from data
    return spec.psd(data=xx,nfft=nfft,n0=n0,window=window,fs=self.fs)
  def psd_welch(self,n0=0,n1=None,window=_np.hanning,nperseg=4096,noverlap=None):
    """calculate the PSD [m**2/Hz] using the Welche method.
    For more information of input parameters see
    scipy.signal.welch.
    For correct normalisation see CERN-THESIS-2013-028
    n0,n1     = use data[n0:n1]
    """
    xx=self.data[self.plane]
    xx=xx-mean(xx)
    if(n1==None):
      n1=len(xx)-n0
    xx=xx[n0:n1]
    ff,psd=welch(xx,fs=self.fs,window=window(nperseg),nperseg=nperseg,noverlap=noverlap,nfft=None,detrend=False,return_onesided=True,scaling='density')
    return ff,psd
  def opt_plot_psd(self,xlog,ylog):
    _pl.xlim(1,100)
    _pl.xlabel(r'f [Hz]')
    _pl.ylabel(r'PSD [$\mathrm{\mu m}^2$/Hz]')
    _pl.grid(which='both')
    if(xlog):
      _pl.xscale('log')
    if(ylog):
      _pl.yscale('log')
  def plot_psd(self,nfft=None,n0=0,window=None,scale=1,lbl=None,color='b',linestyle='-',xlog=True,ylog=True):
    """plot the PSD spectrum in m**2/Hz
    window = window function
    nfft   = number of data points used for FFT
    n0     = use data[n0:n0+nfft]
    offset = curve is shifted by offset in plot to distinguish lines
             which would coincide
    """
    if lbl==None:
      lbl='%s %s'%(self.pos,self.plane)
    ff,psd=self.psd(nfft=nfft,window=window,n0=n0,fs=self.fs)
    _pl.plot(ff[1:],scale*psd[1:],color=color,linestyle=linestyle,label=lbl)#do not plot DC offset
    self.opt_plot_psd(xlog,ylog)
  def plot_psd_welch(self,n0=0,n1=None,window=_np.hanning,nperseg=4096,noverlap=None,scale=1,lbl=None,color='b',linestyle='-',xlog=True,ylog=True):
    """plot the PSD spectrum in m**2/Hz using the welch method
    window = window function (default: hanning)
    n0     = use data[n0:n1]
    offset = curve is shifted by offset in plot to distinguish lines
             which would coincide
    """
    if lbl==None:
      lbl='%s %s'%(self.pos,self.plane)
    ff,psd=self.psd_welch(bb=bb,n0=n0,n1=n1,window=window,nperseg=nperseg,noverlap=noverlap)
    _pl.plot(ff[1:],scale*psd[1:],color=color,linestyle=linestyle,label=lbl)#do not plot DC offset
    self.opt_plot_psd(xlog,ylog)