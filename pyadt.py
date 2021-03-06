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
import timbertools as timb
#from timbertools import *


def get_timestamp_adt(fn):
  if '/' in fn:
    fn = fn.split('/')[-1]
  date=fn.split('-')[0]
  time=fn.split('-')[1]
  date='-'.join([date[0:4],date[4:6],date[6:8]])
  time=':'.join([time[0:2],time[3:5],time[6:8]+'.000'])
  return date+' '+time
def sort_files(files):
  return (timb.mk_sort_files(get_timestamp_adt))(files)
def get_fn_data(fn):
  """get time stamp, plane and bunch number
  from filename fn"""
  fn=(fn.split('/')[-1]).split('.mat')[0]#remove .mat file ending
  idxbunch=int(fn.split('bunch')[1])
  #convert time stamp to timber input format yyyy-mm-dd hh:mm:ss.sss
  timestamp=get_timestamp_adt(fn)
  #get position and plane
  planepos=fn.split('-')[-2]#e.g. 'ADTBposHorQ9LB1'
  pos=planepos[-5:-2]#position of ADT: Q9[LR],Q7[LR]
  beam=planepos[-2:].lower()#beam
  dic={'Hor':'h','Ver':'v'}
  plane=beam+dic[planepos.split('pos')[1].split('Q')[0]] 
  return idxbunch,timestamp,pos,plane,beam
def getbeta_adt(dn,force=False,verbose=False):
  """get beta function for all files in folder *dn*,
  data is saved in a pickled file with the normal
  format as from rdmstores/timber"""
  files=glob.glob(dn+'/*.mat')
  if( force==False and os.path.isfile(dn+'/betastar.p')):
    beta = pickle.load(open(dn+'/betastar.p',"rb"))
    if verbose == True: print '%s found!'%(dn+'/betastar.p')
  else:
    #get earliest and latest timestamp
    fmt='%Y-%m-%d %H:%M:%S.SSS'
    fn=files[0]
    ts = _np.array([ strpunix(get_timestamp_adt(fn),fmt) for fn in files ])
    #add 60 min at beginning and end`
    start = dumpdate(addtimedelta(ts.min(),-60*60))
    end   = dumpdate(addtimedelta(ts.max(),60*60) )
    try:
      beta  = logdb.get(['HX:BETASTAR_IP1','HX:BETASTAR_IP2','HX:BETASTAR_IP5','HX:BETASTAR_IP8'],start,end)
      pickle.dump(beta,open(dn+'/betastar.p',"wb"))
    except IOError:
      print 'ERROR: logdb can not be accessed! No beta* value obtained.'
      beta=None
  return beta

class adt():
  """class to import data from ADT"""
  bitstomum = 1/3.0 #3 bits per mum for conversion of the signal
  beta_adt={'Q9Lb1h':130.580,'Q7Lb1h':131.252,'Q7Rb1v':114.997,'Q9Rb1v':125.984,'Q7Rb2h':201.770,'Q9Rb2h':113.417,'Q7Lb2v':167.025,'Q9Lb2v':138.484}#beta functions [m] at ADT pickups
  def __init__(self,idxbunch=0,timestamp='2008-09-10 08:30:00.000',pos='Q9L',plane='b1h',beam='B1',data=[],fs=11245.0,beta=0,fn='',betadt=0):
    """idxbunch:  bunch number
       filename:  source data file
       timestamp: timestamp of measurement (starting time)
       pos:       ADT pickup (Q9[LR],Q7[LR] beam 1 or beam2)
       plane:     plane (horizontal or vetical)
       fs:        revolution frequency
       beta:      beta at IP1/2/5/8 during *timestamp*
       data:      orbit in mum (without mean value subtracted"""
    self.idxbunch  = idxbunch
    self.timestamp = timestamp #timestamp in local time
    self.t0        = strpunix(timestamp)#timestamp in unix time 
    self.pos       = pos
    self.plane     = plane
    self.beam      = beam
    self.fs        = fs
    self.beta      = beta
    self.betadt    = betadt
    self.data={plane:_np.array(data)}
    self.filename  = fn  
  @classmethod
  def getdata(cls,fn,force=False):
    """load data and return orbit in mum"""
    idxbunch,timestamp,pos,plane,beam = get_fn_data(fn)
    dd   = sio.loadmat(fn)
    data = ((sio.loadmat(fn))['Mem1']).flatten()
    data = cls.bitstomum*data #convert bits to mum
    dn = os.path.split(os.path.abspath(fn))[0]#directory of file
    beta = getbeta_adt(dn,force) 
    mtime = strpunix(get_timestamp_adt(fn),'%Y-%m-%d %H:%M:%S.SSS') #get timestamp
    betasample=timb.getbetasample(beta,mtime)
    fs = 11245.0 #revolution frequency
    betadt=cls.beta_adt[pos+plane]
    return cls(idxbunch,timestamp,pos,plane,beam,data,fs,betasample,fn,betadt)
  def orb(self):
    """subtract mean value from orbit"""
    xx=self.data[self.plane]
    xx=(xx-_np.mean(xx))#subtract average value from data
    return xx
  def fft(self,nfft=None,n0=0,window=_np.hanning,scale=1.0):
    """calculate the real fft (calls np.fft.rfft)
    window = window function
    nfft   = number of data points used for FFT
    n0,n1     = use data[n0:n0+nfft]
    scale  = scale input data by *scale*
    for the fft the average value is subtracted from the data:
      data=data-mean(data)
    """
    xx=scale*self.orb()
    if(nfft==None):
      nfft=len(xx)-n0
    if window==None:
      xx=xx[n0:n0+nfft]
    else:
      win=window(nfft)
      xx=xx[n0:n0+nfft]*win
    ff =_np.arange(nfft/2+1)*self.fs/(nfft-1)
    fft=_np.fft.rfft(xx)
    return ff,fft
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
    """calculate the PSD [mum**2/Hz] using the Welche method.
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
  def fft_bartlett(self,nfft=None,aavg=20,nslice=None,window=_np.hanning,scale=1.0):
    directory=os.path.dirname(self.filename)
    filenames='%s/*%s%s*%s.mat'%(directory,self.pos,self.plane[:-1].upper(),self.idxbunch)
    print filenames
    files=sort_files(glob.glob(filenames))
    start= files.index(self.filename)
    if nslice == None or nslice==0 or nslice==1:
      n0s = [0]
      nfft = None #then takes the full length
    else:
      nfft = int(len(self.data[self.plane]-7)/nslice)#if length not given take full length minus 7 points margin in case points are missing
      n0s   = [int(n) for n in _np.arange(nslice)*nfft ]
      aavg = int(round(aavg/nslice))
    end  = start+aavg
    if end > len(files):#in case not enough files are there for averaging
      end = len(files)
      aavg = end - start
      print 'average only over %s files'%aavg
    lfft={'f':[],'fft':[]}
    for ii in start+_np.arange(aavg):#take average over aavg files with larger timestamp
      print ii
      fn=files[ii]
      dd=adt.getdata(fn)
      for n0 in n0s:
        f_fft,fft=dd.fft(nfft=nfft,n0=n0,window=window,scale=scale)
        lfft['f'].append(f_fft)
        lfft['fft'].append(_np.abs(fft))
    f_fft_avg=_np.mean(lfft['f'],axis=0)
    fft_avg=_np.mean(lfft['fft'],axis=0)
    return f_fft_avg,fft_avg
  def psd_bartlett(self,nfft=None,aavg=20,nslice=None,window=_np.hanning,scale=1.0):
    directory=os.path.dirname(self.filename)
    filenames='%s/*%s%s*%s.mat'%(directory,self.pos,self.plane[:-1].upper(),self.idxbunch)
    files=sort_files(glob.glob(filenames))
    start= files.index(self.filename)
    if nslice == None or nslice==0 or nslice==1:
      n0s = [0]
      nfft = None #then takes the full length
    else:
      nfft = int(len(self.data[self.plane]-7)/nslice)#if length not given take full length minus 7 points margin in case points are missing
      n0s   = [int(n) for n in _np.arange(nslice)*nfft ]
      aavg = int(round(aavg/nslice))
    end  = start+aavg
    if end > len(files):#in case not enough files are there for averaging
      end = len(files)
      aavg = end - start
      print 'average only over %s files'%aavg
    lpsd={'f':[],'psd':[]}
    for ii in start+_np.arange(aavg):#take average over aavg files with larger timestamp
      fn=files[ii]
      dd=adt.getdata(fn)
      for n0 in n0s:
        f_psd,psd=dd.psd(nfft=nfft,n0=n0,window=window,scale=scale)
        lpsd['f'].append(f_psd)
        lpsd['psd'].append(_np.abs(psd))
    f_psd_avg=_np.mean(lpsd['f'],axis=0)
    psd_avg=_np.mean(lpsd['psd'],axis=0)
    return f_psd_avg,psd_avg
  def plot_psd_bartlett(self,nfft=None,aavg=20,nslice=None,window=None,scale=1.0,lbl=None,color='b',linestyle='-',xlog=True,ylog=True):
    """plot the PSD spectrum in mum**2/Hz
    averaged over *aavg* samples (bartlett
    method)
    window = window function
    nfft   = number of data points used for FFT
    n0     = use data[n0:n0+nfft]
    scale  = scale input data by *scale*
    """
    if lbl==None:
      lbl='%s %s'%(self.pos,self.plane)
    ff,psd=self.psd_bartlett(nfft=nfft,aavg=aavg,nslice=nslice,window=window,scale=scale)
    _pl.plot(ff[1:],psd[1:],color=color,linestyle=linestyle,label=lbl)#do not plot DC offset
    self.opt_plot_psd(xlog,ylog)
    #if scale !=1 values are scales with 1/sqrt(beta) -> units are mu 1.e-12*m
    if abs(scale-1.0)>1.e-6:  _pl.ylabel(r'PSD [$p\mathrm{m}$/Hz]',fontsize=14)
  def plot_fft_bartlett(self,nfft=None,aavg=20,nslice=None,window=None,scale=1.0,lbl=None,color='b',linestyle='-',xlog=True,ylog=True):
    """plot the fft spectrum in mum
    where the fft is averaged with the bartlett method
    window = window function
    nfft   = number of data points used for FFT
    n0     = use data[n0:n0+nfft]
    scale  = scale input data by *scale*
    """
    if lbl==None:
      lbl='%s %s'%(self.pos,self.plane)
    ff,fft=self.fft_bartlett(nfft=nfft,aavg=aavg,nslice=nslice,window=window,scale=scale)
    _pl.plot(ff[1:],_np.abs(fft[1:]),color=color,linestyle=linestyle,label=lbl)#do not plot DC offset
    self.opt_plot_fft(xlog,ylog)
    #if scale !=1 values are scales with 1/sqrt(beta) -> units are mu sqrt(m)
    if abs(scale-1.0)>1.e-6: _pl.ylabel(r'amplitude [$\mu\sqrt{\rm m}$]')
  def plot_fft(self,nfft=None,n0=0,window=None,scale=1.0,lbl=None,color='b',linestyle='-',xlog=True,ylog=True):
    """plot the fft spectrum in mum
    window = window function
    nfft   = number of data points used for FFT
    n0     = use data[n0:n0+nfft]
    scale  = scale input data by *scale*
    """
    if lbl==None:
      lbl='%s %s'%(self.pos,self.plane)
    ff,fft=self.fft(nfft=nfft,window=window,n0=n0,scale=scale)
    _pl.plot(ff[1:],_np.abs(fft[1:]),color=color,linestyle=linestyle,label=lbl)#do not plot DC offset
    self.opt_plot_fft(xlog,ylog)
    #if scale !=1 values are scales with 1/sqrt(beta) -> units are mu sqrt(m)
    if abs(scale-1.0)>1.e-6: _pl.ylabel(r'amplitude [$\mu\sqrt{\rm m}$]')
  def plot_orb_fft(self,nfft=None,window=None,n0=0,scale=1.0,lbl=None,color='b',linestyle='-'):
    """plot the orbit [mum] and fft spectrum in mum
    window = window function
    nfft   = number of data points used for FFT
    n0     = use data[n0:n0+nfft]
    scale  = scale input data by *scale*
    """
    if lbl==None:
      lbl='%s %s'%(self.pos,self.plane)
    xx     = self.orb()
    if(nfft==None):
      nfft=len(xx)-n0
    xx     = xx[n0:n0+nfft] 
    _pl.clf()
    _pl.gcf().set_tight_layout(True)
    _pl.subplot(211)
    _pl.plot(xx,label='%s %s'%(self.pos,self.plane),color=color,linestyle=linestyle)#remove scale from lbl
    self.opt_plot_orb()
    _pl.legend(loc='lower left')
    _pl.subplot(212)
    self.plot_fft(nfft=nfft,window=window,n0=n0,scale=scale,lbl=lbl,color=color,linestyle=linestyle)
    if abs(scale-1.0)>1.e-6: _pl.ylabel(r'amplitude [$\mu\sqrt{\rm m}$]')
    _pl.legend(loc='lower left')
    _pl.ylim(1.e1*scale,1.e6*scale)
    _pl.title('')
  def plot_psd(self,nfft=None,n0=0,window=None,scale=1.0,lbl=None,color='b',linestyle='-',xlog=True,ylog=True):
    """plot the PSD spectrum in mum**2/Hz
    window = window function
    nfft   = number of data points used for FFT
    n0     = use data[n0:n0+nfft]
    scale  = scale input data by *scale*
    """
    if lbl==None:
      lbl='%s %s'%(self.pos,self.plane)
    ff,psd=self.psd(nfft=nfft,window=window,n0=n0,scale=scale)
    _pl.plot(ff[1:],psd[1:],color=color,linestyle=linestyle,label=lbl)#do not plot DC offset
    self.opt_plot_psd(xlog,ylog)
    if abs(scale-1.0)>1.e-6:  _pl.ylabel(r'PSD [$p\mathrm{m}$/Hz]',fontsize=14)
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
    if(nfft==None):
      nfft=len(xx)-n0
    xx     = xx[n0:n0+nfft] 
    _pl.clf()
    _pl.gcf().set_tight_layout(True)
    _pl.subplot(211)
    _pl.plot(xx,label='%s %s'%(self.pos,self.plane),color=color,linestyle=linestyle)#remove scale from lbl
    self.opt_plot_orb()
    _pl.legend(loc='lower left')
    _pl.subplot(212)
    self.plot_psd(nfft=nfft,n0=n0,window=window,scale=scale,color=color,linestyle=linestyle,lbl='%s, scale=%4.2f'%(lbl,scale))
    #if scale !=1 values are scales with 1/sqrt(beta) -> units are mu 1.e-12*m
    if abs(scale-1.0)>1.e-6:  _pl.ylabel(r'PSD [$p\mathrm{m}$/Hz]',fontsize=14)
    _pl.legend(loc='lower left')
    _pl.ylim(1.e-12*scale**2,1.e-2*scale**2)
    _pl.title('')
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
  def opt_plot_orb(self,ylim=(-20,20)):
    _pl.ylim(ylim)
    _pl.grid(which='both')
    _pl.xlabel('number of turns')
    _pl.ylabel(r'z [$\mu$m]')
    _pl.title(r'$\beta_{\rm IP1}=$%4.2f cm, %s'%(self.beta['betaIP1'],dumpdate(self.t0,fmt='%Y-%m-%d %H:%M:%S')))
  def opt_plot_fft(self,xlog,ylog):
    _pl.xlim(1,100)
    _pl.xlabel(r'f [Hz]')
    _pl.ylabel(r'amplitude [$\mathrm{\mu m}$]')#abs(fft)
    _pl.grid(which='both')
    _pl.title(r'$\beta_{\rm IP1}=$%4.2f cm, %s'%(self.beta['betaIP1'],dumpdate(self.t0,fmt='%Y-%m-%d %H:%M:%S')))
    if(xlog):
      _pl.xscale('log')
    if(ylog):
      _pl.yscale('log')
  def opt_plot_psd(self,xlog,ylog):
    _pl.xlim(1,100)
    _pl.xlabel(r'f [Hz]')
    _pl.ylabel(r'PSD [$\mathrm{\mu m}^2$/Hz]')
    _pl.grid(which='both')
    _pl.title(r'$\beta_{\rm IP1}=$%4.2f cm, %s'%(self.beta['betaIP1'],dumpdate(self.t0,fmt='%Y-%m-%d %H:%M:%S')))
    if(xlog):
      _pl.xscale('log')
    if(ylog):
      _pl.yscale('log')
