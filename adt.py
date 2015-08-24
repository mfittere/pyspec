"""class to import data from ADT"""
import numpy as _np
import matplotlib.pyplot as _pl
import scipy.io as sio
import spec as spec
from scipy.signal import welch
from rdmstores import *
import glob as glob
import os as os
import cPickle as pickle

def get_timestamp(fn):
  if '/' in fn:
    fn = fn.split('/')[-1]
  date=fn.split('-')[0]
  time=fn.split('-')[1]
  date='-'.join([date[0:4],date[4:6],date[6:8]])
  time=':'.join([time[0:2],time[3:5],time[6:8]+'.000'])
  return date+' '+time
def get_fn_data(fn):
  """get time stamp, plane and bunch number
  from filename fn"""
  fn=(fn.split('/')[-1]).split('.mat')[0]#remove .mat file ending
  idxbunch=int(fn.split('bunch')[1])
  #convert time stamp to timber input format yyyy-mm-dd hh:mm:ss.sss
  timestamp=get_timestamp(fn)
  #get position and plane
  planepos=fn.split('-')[-2]#e.g. 'ADTBposHorQ9LB1'
  pos=planepos[-5:-2]#position of ADT: Q9[LR],Q7[LR]
  beam=planepos[-2:].lower()#beam
  dic={'Hor':'h','Ver':'v'}
  plane=beam+dic[planepos.split('pos')[1].split('Q')[0]] 
  return idxbunch,timestamp,pos,plane

def getbeta(dn,force=False):
  """get beta function for all files in folder *dn*,
  data is saved in a pickled file with the normal
  format as from rdmstores/timber"""
  files=glob.glob(dn+'/*.mat')
  if( force==False and os.path.isfile(dn+'/betastar.p')):
    beta = pickle.load(open(dn+'/betastar.p',"rb"))
    print '%s found!'%(dn+'/betastar.p')
  else:
    #get earliest and latest timestamp
    fmt='%Y-%m-%d %H:%M:%S.SSS'
    fn=files[0]
    ts = _np.array([ strpunix(get_timestamp(fn),fmt) for fn in files ])
    #add 60 min at beginning and end`
    start = dumpdate(addtimedelta(ts.min(),-60*60))
    end   = dumpdate(addtimedelta(ts.max(),60*60) )
    try:
      beta  = logdb.get(['HX:BETASTAR_IP1','HX:BETASTAR_IP2','HX:BETASTAR_IP5','HX:BETASTAR_IP8'],start,end)
      pickle.dump(beta,open(dn+'/betastar.p',"wb"))
    except IOError:
      print 'ERROR: logdb can not be accessed!'
      beta=None
  return beta

class adt():
  """class to import data from ADT"""
  bitstomum = 1/3.0 #3 bits per mum for conversion of the signal
  def __init__(self,idxbunch=0,timestamp='2008-09-10 08:30:00.000',pos='Q9L',plane='h',data=[],fs=11245.0,betastar=0):
    """idxbunch:  bunch number
       timestamp: timestamp of measurement (starting time)
       pos:       ADT pickup (Q9[LR],Q7[LR] beam 1 or beam2)
       plane:     plane (horizontal or vetical)
       fs:        revolution frequency
       betastar:  betastar during *timestamp*
       data:      orbit in mum (without mean value subtracted"""
    self.idxbunch  = idxbunch
    self.timestamp = timestamp #timestamp in local time
    self.t0        = strpunix(timestamp)#timestamp in unix time 
    self.pos       = pos
    self.plane     = plane
    self.fs        = fs
    self.betastar  = betastar
    self.data={plane:_np.array(data)}
  @classmethod
  def getdata(cls,fn,force=False):
    """load data and return orbit in mum"""
    idxbunch,timestamp,pos,plane = get_fn_data(fn)
    dd   = sio.loadmat(fn)
    data = ((sio.loadmat(fn))['Mem1']).flatten()
    data = cls.bitstomum*data #convert bits to mum
    betastar = getbeta(dn,force)  
    fs = 11245.0 #revolution frequency
    return cls(idxbunch,timestamp,pos,plane,data,fs,betastar)
  def orb(self):
    """subtract mean value from orbit"""
    xx=self.data[self.plane]
    xx=(xx-_np.mean(xx))#subtract average value from data
    return xx
  def psd(self,nfft=None,n0=0,window=_np.hanning,scale=1.0):
    """calculate the PSD [m**2/Hz]. For correct normalisation
    see CERN-THESIS-2013-028
    window = window function
    nfft   = number of data points used for FFT
    n0,n1     = use data[n0:n0+nfft]
    scale  = scale input data by *scale*
    for the psd the average value is subtracted from the data:
      data=data-mean(data)
    """
    xx=scale*self.orb()
    return spec.psd(data=xx,nfft=nfft,n0=n0,window=window,fs=self.fs)
  def psd_welch(self,n0=0,n1=None,window=_np.hanning,nperseg=4096,noverlap=None):
    """calculate the PSD [m**2/Hz] using the Welche method.
    For more information of input parameters see
    scipy.signal.welch.
    For correct normalisation see CERN-THESIS-2013-028
    n0,n1     = use data[n0:n1]
    """
    xx=self.orb()
    if(n1==None):
      n1=len(xx)-n0
    xx=xx[n0:n1]
    ff,psd=welch(xx,window=window(nperseg),nperseg=nperseg,noverlap=noverlap,nfft=None,detrend=False,return_onesided=True,scaling='density')
    return ff,psd
  def plot_psd(self,nfft=None,n0=0,window=None,scale=1.0,lbl=None,color='b',linestyle='-',xlog=True,ylog=True):
    """plot the PSD spectrum in m**2/Hz
    window = window function
    nfft   = number of data points used for FFT
    n0     = use data[n0:n0+nfft]
    scale  = scale input data by *scale*
    """
    if lbl==None:
      lbl='%s %s'%(self.pos,self.plane)
    ff,psd=self.psd(nfft=nfft,window=window,n0=n0,fs=self.fs,scale=scale)
    _pl.plot(ff[1:],psd[1:],color=color,linestyle=linestyle,label=lbl)#do not plot DC offset
    self.opt_plot_psd(xlog,ylog)
  def plot_orb_psd(self,nfft=None,window=None,n0=0,scale=1.0,lbl=None,color='b',linestyle='-'):
    """plot the orbit [mum] and PSD spectrum in mum**2/Hz
    window = window function
    nfft   = number of data points used for FFT
    n0     = use data[n0:n0+nfft]
    scale  = scale input data by *scale*
    """
    if lbl==None:
      lbl='%s %s'%(self.pos,self.plane)
    xx     = self.orb()
    ff,psd=self.psd(nfft=nfft,window=window,n0=n0,scale=scale)
    _pl.clf()
    _pl.subplot(211)
    _pl.plot(xx,label=lbl)
    self.opt_plot_orb(ylim=(-30,30))
    _pl.legend(loc='lower left')
    _pl.subplot(212)
    _pl.plot(ff[1:],psd[1:],color=color,linestyle=linestyle,label='%s, scale=%4.2f'%(lbl,scale))#do not plot DC offset
    self.opt_plot_psd(True,True)
    _pl.legend(loc='lower left')
    _pl.ylim(1.e-12,1.e-2)
  def plot_psd_welch(self,n0=0,n1=None,window=_np.hanning,nperseg=4096,noverlap=None,scale=1,lbl=None,color='b',linestyle='-',xlog=True,ylog=True):
    """plot the PSD spectrum in m**2/Hz using the welch method
    window: window function (default: hanning)
    n0,n1    : use data[n0:n1]
    nperseg : int, optional
    Length of each segment.  Defaults to 4096.
    noverlap: int, optional
    Number of points to overlap between segments. If None,
    """
    if lbl==None:
      lbl='%s %s'%(self.pos,self.plane)
    ff,psd=self.psd_welch(n0=n0,n1=n1,window=window,nperseg=nperseg,noverlap=noverlap)
    _pl.plot(ff[1:],scale*psd[1:],color=color,linestyle=linestyle,label=lbl)#do not plot DC offset
    self.opt_plot_psd(xlog,ylog)
  def opt_plot_orb(self,ylim=(-30,30)):
    _pl.ylim(ylim)
    _pl.grid(which='both')
    _pl.xlabel('number of turns')
    _pl.ylabel(r'z [$\mu$m]')
  def opt_plot_psd(self,xlog,ylog):
    _pl.xlim(1,100)
    _pl.xlabel(r'f [Hz]')
    _pl.ylabel(r'PSD [$\mathrm{\mu m}^2$/Hz]')
    _pl.grid(which='both')
    _pl.title(r'$\beta*_{\rm IP1}=$%4.2f'%self.betastar)
    if(xlog):
      _pl.xscale('log')
    if(ylog):
      _pl.yscale('log')
