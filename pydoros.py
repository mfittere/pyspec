"""class to import data from DOROS BPMs"""
import numpy as _np
import matplotlib.pyplot as _pl
import os,time
import cPickle as pickle
import glob as glob
import spec as spec
from scipy.signal import welch
from rdmstores import *
import timbertools as timb
import datetime
import time
def ipaddress_to_irlr(ipaddress):
  if(ipaddress=='172_18_66_233'):
    ip='1'
    side='R'
    front='frontend 1'
  if(ipaddress=='172_18_66_234'):
    ip='1'
    side='L'
    front='frontend 2'
  if(ipaddress=='172_18_41_214'):
    ip='5'
    side='R'
    front='frontend 1'
  if(ipaddress=='172_18_53_135'):
    ip='5'
    side='L'
    front='frontend 2'
  return ip,side,front
def getbeta_doros(dn,force=False):
  files=glob.glob(dn+'/*.bin')
  if( force==False and os.path.isfile(dn+'/betastar.p')):
    beta = pickle.load(open(dn+'/betastar.p',"rb"))
    print '%s found!'%(dn+'/betastar.p')
  else:
    ts = [ os.path.getmtime(fn) for fn in files ]#get timestamps
    tmin =_np.min(ts) 
    tmax =_np.max(ts) 
    #add 60 min at beginning and end`
    start = dumpdate(addtimedelta(tmin,-60*60))
    end   = dumpdate(addtimedelta(tmax,60*60) )
    try:
      beta  = logdb.get(['HX:BETASTAR_IP1','HX:BETASTAR_IP2','HX:BETASTAR_IP5','HX:BETASTAR_IP8'],start,end)
      pickle.dump(beta,open(dn+'/betastar.p',"wb"))
    except IOError:
      print 'ERROR: measurement database can not be accessed!'
      beta=None
  return beta
def _get_data_2(fn,beta=None):
  HeaderOffset =2
  CaptureDecimationFactor = 20 
  ff = open(fn+'.bin','rb')
  data=_np.fromfile(ff,dtype='uint32')#binary file written in uint32
  ip,side,front=ipaddress_to_irlr(fn.split('IP')[1])
  if beta==None: beta=getbeta_doros(os.path.dirname(fn))
  mtime=os.path.getmtime(fn+'.bin')
  betasample=timb.getbetasample(beta,mtime)
  B1NumofChan=data[0]& (2**16-1)#ersten 16 bit von rechts gelesen
  B2NumofChan=(data[0]&((2**16-1)<<16))>>16
  Numofsamples = data[1]
  ADCB1dorCHdata,ADCB2dorCHdata=decode_channel_2(data,B1NumofChan,B2NumofChan,Numofsamples,HeaderOffset)
  b1h1,b1h2,b1v1,b1v2=ADCB1dorCHdata[:4]
  b2h1,b2h2,b2v1,b2v2=ADCB2dorCHdata[:4]
  return [b1h1,b1h2,b1v1,b1v2,b2h1,b2h2,b2v1,b2v2,mtime,betasample,ip]   
def decode_channel_2(data,n1,n2,nsample,header):
  """decode binary files for files older than
  30/08/2015"""
  imax32bit=2**32-1
  if n1==0 and n2==0: d1,d2=[],[]
  if n1 >0 and n2==0:
    d = _np.array(data[header:header+n1*nsample],dtype=_np.float64)/imax32bit
    d1,d2 = d.reshape(n1,nsample),[]
  if n1==0 and n2 > 0:
    d= _np.array(data[header:header+n2*nsample],dtype=_np.float64)/imax32bit
    d1,d2 = [],d.reshape(n2,nsample)
  if n1 >0 and n2 > 0:
    d1aux= _np.array(data[header:header+n1*nsample],dtype=_np.float64)/imax32bit
    d1=d1aux.reshape(n1,nsample)
    d2aux= _np.array(data[header+n1*nsample:header+n1*nsample+n2*nsample],dtype=_np.float64)/imax32bit
    d2=d2aux.reshape(n2,nsample)
  return d1,d2

def _get_data_1(fn,beta=None):
  """function to get data before 30/08/2015"""
  #---- define parameters used for import the datafile ----
  #ADC=Analog to Digital Converter parameters
  #Decoding parameters
  nChan = 8  # number of channel 
  nByte = 4  # number of bytes per sample
  nFramePerADCb = 15#number of frames per ADC buffer
  nADCBuff = 2#number of ADC buffers (or ADCs)
  nADCchan = 8#Number of ADC channels
  LheaderData = 8# Data header length
  LheaderUdp = 20#UDP header length
  LADCbBuff = nFramePerADCb*nChan*nByte
  SampleDecimation = 1#data decimation factor
  #UDP=User Datagram Protocol structure addressing
  ADCb1nFrameAddr = 17
  ADCb2nFrameAddr = 18
  FIFOovfladdr = 19
  SeqNumberOffset = 1
  #UDP per frame
  LdataUdp = nByte*nChan*nFramePerADCb*nADCBuff #number of bytes of UDP data
  Ludp = LheaderUdp + LdataUdp #UDP header + data length
  data,udpcheck,ip=readbinary(fn+'.bin',LheaderData,LheaderUdp,LdataUdp)
  if(udpcheck):#only process data if udp check passed
    if beta==None: beta=getbeta_doros(os.path.dirname(fn))
    ADC1chanTable=decode_channel_1(data,ADCb1nFrameAddr,nADCchan,nByte,LADCbBuff,nChan,SampleDecimation,beam='b1')
    ADC2chanTable=decode_channel_1(data,ADCb2nFrameAddr,nADCchan,nByte,LADCbBuff,nChan,SampleDecimation,beam='b2')
    #number of valid frames per udp ADC1 buffer
    [b1h1,b1h2,b1v1,b1v2]=ADC1chanTable[0:4]/(2**24-1)
    [b2h1,b2h2,b2v1,b2v2]=ADC2chanTable[0:4]/(2**24-1)
    mtime=os.path.getmtime(fn+'.bin')#ctime gives back the timestamp when the file was created for MAC, instead use mtime, which seems to give back the correct timestamp
    betasample=timb.getbetasample(beta,mtime)
    #store already processed orbit data in *.p
    pickle.dump([b1h1,b1h2,b1v1,b1v2,b2h1,b2h2,b2v1,b2v2,mtime,betasample,ip],open(fn+'.p',"wb"))  
    print '... store b1h1,b1h2,b2h1,b2h2 etc. in file %s.p for faster reload'%(fn.split('/')[-1])
  else:
    ip = '1'
    betasample = None
    mtime=0
    [b1h1,b1h2,b1v1,b1v2,b2h1,b2h2,b2v1,b2v2]=[[] for x in range(8)]
  return [b1h1,b1h2,b1v1,b1v2,b2h1,b2h2,b2v1,b2v2,mtime,betasample,ip]
def decode_channel_1(data,ADCbnFrameAddr,nADCchan,nByte,LADCbBuff,nChan,SampleDecimation,beam='b1'):
  #extract ADCb[12] buffer
  ADCbnFrames = data['header'][:,ADCbnFrameAddr-1]
  nTotalADCbBuffDataFrames = sum(ADCbnFrames)
  if(beam=='b1'):
    ADCbData = _np.hstack(_np.array([ data['data'][idx,0:ADCbnFrames[idx]*nADCchan*nByte] for idx in range(len(ADCbnFrames)) ]))
  if(beam=='b2'):
    ADCbData = _np.hstack(_np.array([ data['data'][idx,LADCbBuff:LADCbBuff+ADCbnFrames[idx]*nADCchan*nByte] for idx in range(len(ADCbnFrames)) ]))
  #decode channel data
  ADCbdataDec=_np.array([ (2**24*ADCbData[4*idx]+2**16*ADCbData[4*idx+1]+2**8*ADCbData[4*idx+2]+ADCbData[4*idx+3]+2**23) & (2**24-1) for idx in range(nTotalADCbBuffDataFrames*nChan) ],dtype=float)
  ADCchanTable = _np.array([ ADCbdataDec[idx:idx+nTotalADCbBuffDataFrames*nChan:SampleDecimation*nChan] for idx in range(nChan) ])
  return ADCchanTable
def readbinary(fn,LheaderData=8,LheaderUdp=20,LdataUdp=980):
  """read in the binary file with filename fn
  and return the header in data['header']
  and the data in data['data']
  file format:
    Dataheader (LheaderData)
    UDPheader (LheaderUdp)
    Data
    (...)   
    UDPheader (LheaderUdp)
    Data"""
  # dataAll=_np.fromfile(file=fn,dtype=_np.uint8,count=-1,sep="") #get the full data
#    adc=ADCb1Data.view('uint32')
  ptostr={'L':'left','R':'right'}
  ip,side,front=ipaddress_to_irlr(fn.split('IP')[1].split('_Data')[0])
  print '... read data from %s (IR%s %s)'%(front,ip,ptostr[side])
  ff = open(fn,'rb')
  header = map(ord,ff.read(LheaderData))
  NumOfFileDumps   = header[0]+header[1]*2**8+header[2]*2**16+header[3]*2**24 #number of files dumped
  NumOfUDPsPerDump = header[4]+header[5]*2**8+header[6]*2**16+header[7]*2**24 #number of UDPs per file
  ff.seek(LheaderData)# remove the Dataheader
  ftype=_np.dtype([('header','%su1'%(LheaderUdp)),('data','%su1'%(LdataUdp))])
  data=_np.fromfile(file=ff,dtype=ftype,count=-1,sep="") #get the full data
  udpcheck=True
  if(NumOfUDPsPerDump!=len(data['header'])):
    udpcheck=False
    print 'WARNING: file corrupted as number of UDP dumps != number of UDPs in file'
    print 'number of UDP dumps:    %s'%(NumOfUDPsPerDump)
    print 'number of UDPs in file: %s'%(len(data['header']))
    print 'lenght of file in bytes:%s'%(os.path.getsize(fn))
#    print '... remove zeros from end of file ...'
#    Ludp = LheaderUdp + LdataUdp
#    TotalLenOfFile = NumOfUDPsPerDump*NumOfFileDumps*Ludp+LheaderData
#    nUdp = NumOfUDPsPerDump*NumOfFileDumps
#    data['header']=data[LheaderData+Ludp:LheaderData + Ludp*(nUdp-1) +LheaderUdp:Ludp]
#  return data[0:TotalLenOfFile],udpcheck
  return data,udpcheck,ip

class doros():
  """class to import data from DOROS BPMS"""
  sIP=21.475 #distance from IP [m]
  def __init__(self,b1h1=[],b1h2=[],b1v1=[],b1v2=[],b2h1=[],b2h2=[],b2v1=[],b2v2=[],mtime=None,beta=None,ip='1',filename=[]):
    self.filename = filename
    self.ip    = ip
    self.beta  = beta 
    self.mtime = mtime
    self.data  = {'b1h1':_np.array(b1h1),'b1h2':_np.array(b1h2),'b1v1':_np.array(b1v1),'b1v2':_np.array(b1v2),'b2h1':_np.array(b2h1),'b2h2':_np.array(b2h2),'b2v1':_np.array(b2v1),'b2v2':_np.array(b2v2)}
    #times for which binary file format changed
    self.t1=time.mktime(datetime.date(2015,8,30).timetuple())
    #pick - up slope in um, theory is 1/4 of the electrode distance [d],
    self.elecdist = 61000#electrode distance
    self.PUgain=self.elecdist/4
    if mtime < self.t1:
      self.FsamADC = (40*1.e6)/3564.0#sampling frequency
    else:
      self.FsamADC = 11245.55
  @classmethod
  def getdata(cls,fn,beta=None,force=False):
    """get the data from the binary file fn.bin and save it
    it in a pickle file fn.p. Add timestamp from creation
    time of binary file. If pickle file already exists, just
    load the pickle file.
    beta is the beta* extracted from timber. read in file
    with beta=pickle.load(open(fnbeta,"rb")).
    
    Note
    ---------
    uses:
        _get_data_1 for timestamps ... -> 30/08/2015
        _get_data_2 for timestamps 30/08/2015 -> ..."""
    #times for which binary file format changed
    cls.t1=time.mktime(datetime.date(2015,8,30).timetuple())
    if(fn.split('.')[-1]=='bin'):
      fn=fn[0:-4]
    if(fn.split('.')[-1]=='p'):
      fn=fn[0:-2]
    if(os.path.isfile(fn+'.p') and force==False):
      print 'from *.p file'
      b1h1,b1h2,b1v1,b1v2,b2h1,b2h2,b2v1,b2v2,mtime,betasample,ip=pickle.load(open(fn+'.p',"rb"))
    else:
      mtime=os.path.getmtime(fn+'.bin')#ctime gives back the timestamp when the file was created for MAC, instead use mtime, which seems to give back the correct timestamp
      if mtime<cls.t1: [b1h1,b1h2,b1v1,b1v2,b2h1,b2h2,b2v1,b2v2,mtime,betasample,ip]=_get_data_1(fn,beta)
      if mtime>cls.t1: [b1h1,b1h2,b1v1,b1v2,b2h1,b2h2,b2v1,b2v2,mtime,betasample,ip]=_get_data_2(fn,beta)
    return cls(b1h1,b1h2,b1v1,b1v2,b2h1,b2h2,b2v1,b2v2,mtime,betasample,ip,fn)
  def checkdata(self,fn):
    """check for errors in data acquisition for
    files with timestamps ... -> 30/08/2015"""
    if os.path.getmtime(fn+'.bin')<self.t1:
      LeAll = os.path.getsize(fn)#total file byte length
      nUdp = (LeAll - self.LheaderData)/self.Ludp #number of UDPs
      dataAll=_np.fromfile(file=fn,dtype=_np.uint8,count=-1,sep="")
      #Extract how many UDPs were lost from sequential number
      idx=self.LheaderData + self.SeqNumberOffset-1
      iUdpFirst = 2**24*dataAll[idx] + 2**16*dataAll[idx + 1] + 2**8*dataAll[idx + 2] + dataAll[idx + 3]
      iUdpLast =2**24*dataAll[idx + self.Ludp*(nUdp - 1)] + 2**16*dataAll[idx + 1 + self.Ludp*(nUdp - 1)] + 2**8*dataAll[idx + 2 + self.Ludp*(nUdp - 1)] + dataAll[idx + 3 + self.Ludp*(nUdp - 1)]
      nUdpLost=iUdpLast - iUdpFirst - nUdp + 1 # number of lost UDPs
      print "iUdpFirst=%d"%(iUdpFirst)
      print "iUdpLast =%d"%(iUdpLast)
      print "nUdpLost =%d"%(nUdpLost)
      #check FIFO overflow
      FIFOovfl = sum(dataAll[self.LheaderData + self.FIFOovfladdr-1:self.LheaderData + self.FIFOovfladdr+self.Ludp*(nUdp-1)-1:self.Ludp])
      print "FIFOovfl =%d"%(FIFOovfl)
    else:
      print "WARNING: check only necessary for files with
      timestampes <  30/08/2015"""
  @classmethod
  def process_dir(cls,dn,force=False):
    """process orbit data in directory dn and
    store it in *.p files"""
    files = glob.glob(dn+'/*.bin')
    beta  = getbeta_doros(dn,force=True) #rdmstore object with beta* values
    for fn in files:
      if(os.path.isfile(fn[0:-4]+'.p') and force==False):
        print fn+' is already processed'
      else:
        print 'processing '+fn
        cls.getdata(fn,beta=beta,force=force)
  def acqutime(self):
    """returns the acquisition time in s assuming a revolution frequency of doros.fsamADC"""
    ll=len(self.data['b1h1'])
    lflag=True
    for bb in 'b1','b2':
      for pp in 'h1','v1','h2','v2':
        if((len(self.data[bb+pp])-ll)>1.e-3):
          lflag=False
          print '%s: data acquisition of %4.0f samples over %4.2f seconds'%(bb+pp,len(self.data[bb+pp]),len(self.data[bb+pp])/self.FsamADC)
    if(lflag):
      print 'data acquisition over of %4.0f samples over %4.2f seconds'%(ll,ll/self.FsamADC)
    else:
      print 'WARNING: not all b[12][hv][12] arrays have the same length!' 
    return ll,ll/self.FsamADC
  def dumpdate(self,fmt='%Y-%m-%d %H:%M:%S.SSS'):
    return dumpdate(self.mtime,fmt=fmt) 
  def betabpm(self):
    """return beta at DOROS BPMs for beta*=beta
    beta_doros=beta+s**2/beta
    as only a drift space is between the IP and
    the DOROS BPMs"""
    #self.beta['betaIP'+ip]=0 if no value is saved
    if abs(self.beta['betaIP'+self.ip]) > 1.e-8:
      return self.beta['betaIP'+self.ip]*1.e-2+self.sIP**2/(self.beta['betaIP'+self.ip]*1.e-2) 
    else:
      return 1.0 #return 1.0, so that psd normalization is not influenced
  def orb(self):
    """returns a dictionary with the orbit b1h,b1v,b2h,b2v in mum,
    where the orbit is defined over the signals of the two electrodes:
    b1haux=PUgain*(b1h1-b1h2)/(b1h1+b1h2)
    b1h=b1haux/mean(b1haux) etc."""
    dic={}
    for bb in 'b1','b2':
      for pp in 'h','v':
        dic[bb+pp]=self.PUgain*(self.data[bb+pp+'1']-self.data[bb+pp+'2'])/(self.data[bb+pp+'1']+self.data[bb+pp+'2'])
        #subtract average orbit
        dic[bb+pp]=dic[bb+pp]-_np.mean(dic[bb+pp])
    return dic
  def psd(self,bb='b1h',nfft=None,n0=0,window=_np.hanning,scale=1.0):
    """calculate the PSD [mum**2/Hz]. For correct normalisation
    see CERN-THESIS-2013-028. This function calls spec.psd().
    window = window function
    nfft   = number of data points used for FFT
    n0     = use data[n0:n0+nfft]
    scale  = scale orbit by *scale*, e.g. 1/self.betabpm()
    """
    xx=scale*self.orb()[bb]
    return spec.psd(data=xx,nfft=nfft,n0=n0,window=window,fs=self.FsamADC)
  def psd_welch(self,bb='b1h',n0=0,n1=None,window=_np.hanning,nperseg=4096,noverlap=None,scale=1.0):
    """calculate the PSD [mum**2/Hz] using the Welche method.
    For more information of input parameters see 
    scipy.signal.welch.
    For correct normalisation see CERN-THESIS-2013-028
    n0,n1: use data[n0:n1]
    """
    xx=scale*self.orb()[bb]
    if(n1==None):
      n1=len(xx)-n0
    xx=xx[n0:n1]
    ff,psd=welch(xx,fs=self.FsamADC,window=window(nperseg),nperseg=nperseg,noverlap=noverlap,nfft=None,detrend=False,return_onesided=True,scaling='density')
    return ff,psd
  def plot_psd(self,bb='b1h',nfft=None,window=None,n0=0,scale=1.0,lbl='',color='b',linestyle='-',xlog=True,ylog=True):
    """plot the PSD spectrum in mum**2/Hz
    window = window function
    nfft   = number of data points used for FFT
    n0     = use data[n0:n0+nfft]
    scale  = scale input data by *scale*
    """
    ff,psd=self.psd(bb=bb,nfft=nfft,window=window,n0=n0,scale=scale)
    _pl.plot(ff[1:],psd[1:],color=color,linestyle=linestyle,label=lbl)#do not plot DC offset
    self.opt_plot_psd(xlog,ylog)
    #if scale !=1 values are scales with 1/sqrt(beta) -> units are mu 1.e-12*m
    if abs(scale-1.0)>1.e-6:  _pl.ylabel(r'PSD [$p\mathrm{m}$/Hz]',fontsize=14)
  def plot_orb_fft(self,bb='b1h',nfft=None,window=None,n0=0,scale=1.0,lbl=None,color='b',linestyle='-'):
    """plot the orbit [mum] and FFT spectrum in mum
    window = window function
    nfft   = number of data points used for FFT
    n0     = use data[n0:n0+nfft] (to be implemented)
    scale  = scale input data by *scale*
    """
    if lbl == None:
      lbl=bb.upper()
    xx     = self.orb()[bb]
    _pl.clf()
    _pl.subplot(211)
    _pl.plot(xx,label=lbl)
    self.opt_plot_orb(ylim=(-30,30))
    _pl.legend(loc='lower left')
    _pl.subplot(212)
    ff,orbfft = self.fft(bb=bb,nfft=nfft,window=window,scale=scale)
    _pl.plot(ff[1:],_np.abs(orbfft[1:]),color=color,linestyle=linestyle,label='%s, scale=%4.2f'%(lbl,scale))#do not plot DC offset
    self.opt_plot_fft(True,True)
    #if scale !=1 values are scales with 1/sqrt(beta) -> units are mu sqrt(m)
    if abs(scale-1.0)>1.e-6: _pl.ylabel(r'amplitude [$\mu\sqrt{\rm m}$]')
    _pl.title('')
    _pl.legend(loc='lower left')
    _pl.ylim(1.e-4,1.e5)
    _pl.tight_layout()
  def plot_orb_psd(self,bb='b1h',nfft=None,window=None,n0=0,scale=1.0,lbl=None,color='b',linestyle='-'):
    """plot the orbit [mum] and PSD spectrum in mum**2/Hz
    window = window function
    nfft   = number of data points used for FFT
    n0     = use data[n0:n0+nfft]
    scale  = scale input data by *scale*
    """
    if lbl == None:
      lbl=bb.upper()
    xx     = self.orb()[bb]
    ff,psd = self.psd(bb=bb,nfft=nfft,window=window,n0=n0,scale=scale)
    _pl.clf()
    _pl.subplot(211)
    _pl.plot(xx,label=lbl)
    self.opt_plot_orb(ylim=(-30,30)) 
    _pl.legend(loc='lower left')
    _pl.subplot(212)
    _pl.plot(ff[1:],psd[1:],color=color,linestyle=linestyle,label='%s, scale=%4.2f'%(lbl,scale))#do not plot DC offset
    self.opt_plot_psd(True,True)
    #if scale !=1 values are scales with 1/sqrt(beta) -> units are mu 1.e-12*m
    if abs(scale-1.0)>1.e-6:  _pl.ylabel(r'PSD [$p\mathrm{m}$/Hz]',fontsize=14)
    _pl.title('')
    _pl.legend(loc='lower left')
    _pl.ylim(1.e-16,1.e-6)
    _pl.tight_layout()
  def plot_orb_all(self):
    """plot orbit position change in mum"""
    _pl.clf()
    for bb in [1,2]:
      for pp in 'h','v':
        _pl.subplot(int(210+bb))
        _pl.plot(self.orb()['b'+str(bb)+pp],label=('b'+str(bb)+pp).upper()) 
        self.opt_plot_orb()
    _pl.legend()
  def plot_psd_welch(self,bb='b1h',n0=0,n1=None,window=_np.hanning,nperseg=4096,noverlap=None,scale=1.0,lbl='',color='b',linestyle='-',xlog=True,ylog=True):
    """plot the PSD spectrum in mum**2/Hz using the welch method
    window = window function (default: hanning)
    n0     = use data[n0:n1]
    offset = curve is shifted by offset in plot to distinguish lines
             which would coincide
    scale  = scale input data by *scale*
    """
    ff,psd=self.psd_welch(bb=bb,n0=n0,n1=n1,window=window,nperseg=nperseg,noverlap=noverlap,scale=scale)
    _pl.plot(ff[1:],psd[1:],color=color,linestyle=linestyle,label=lbl)#do not plot DC offset
    self.opt_plot_psd(xlog,ylog)
    #if scale !=1 values are scales with 1/sqrt(beta) -> units are mu 1.e-12*m
    if abs(scale-1.0)>1.e-6:  _pl.ylabel(r'PSD [$p\mathrm{m}$/Hz]',fontsize=14)
  def opt_plot_orb(self,ylim=(-80,80)):
    _pl.ylim(ylim)
    _pl.grid(which='both')
    _pl.xlabel('number of turns') 
    _pl.ylabel(r'z [$\mu$m]') 
    _pl.title(r'$\beta_{IP%s}=%s $ cm, %s'%(self.ip,self.beta['betaIP'+self.ip],dumpdate(self.mtime,fmt='%Y-%m-%d %H:%M:%S')))
  def opt_plot_psd(self,xlog,ylog):
    _pl.xlim(1,100)
    _pl.xlabel(r'f [Hz]',fontsize=16)
    _pl.ylabel(r'PSD [$\mathrm{\mu m}^2$/Hz]',fontsize=14)
    _pl.grid(which='both')
    _pl.title(r'$\beta_{IP%s}=%s $ cm, %s'%(self.ip,self.beta['betaIP'+self.ip],dumpdate(self.mtime,fmt='%Y-%m-%d %H:%M:%S')))
    if(xlog):
      _pl.xscale('log')
    if(ylog):
      _pl.yscale('log')
  def abs_fft_dB(self,bb='b1h',nfft=None,window=None,n0=0,unit='dB'):
    """calculate the amplitude of the FFT spectrum in dB, 
    where the amplitude is normalized
    in respect of the maximum value:
    fft_1=2*abs(fft(b1h1-b1h2))
    fft=20*log10(fft_1/max(fft_1))
    window = window function
    nfft   = number of data points used for FFT
    n0     = use data[n0:n0+nfft]
    """
    bz1=bb+'1'
    bz2=bb+'2'
    if(nfft==None):
      if(len(self.data[bz1])==len(self.data[bz2])):
        nfft=len(self.data[bz1]) #take all datapoints for the FFT 
        print 'nfft='+str(nfft)
      else:
        print 'WARNING: length of '+bz1+' and '+bz2+'differ! Use len('+bz1+') for FFT'
    if(window==None):
      p1=self.data[bz1][n0:n0+nfft]
      p2=self.data[bz2][n0:n0+nfft]
      pdiff=p1-p2
    else:
      p1=self.data[bz1][n0:n0+nfft]
      p2=self.data[bz2][n0:n0+nfft]
      pdiff=window(nfft)*(p1-p2)
    orbfft=2*_np.abs(_np.fft.rfft(pdiff))
    #normalize to maximum signal
    orbfft=20*_np.log10(orbfft/(orbfft.max()))
    #scale from data points to frequency [Hz]
    FFtscale=self.FsamADC/(nfft-1)
    ff=_np.arange(nfft/2+1)*FFtscale
    return ff,orbfft
  def fft(self,bb='b1h',nfft=None,window=_np.hanning,scale=1.0):
    """calculate the fft, orbit in mum with the mean
    value subtracted (see self.orb) 
    window = window function
    nfft   = number of data points used for FFT
    scale  = scale orbit by *scale*, e.g. 1/self.betabpm()
    """
    xx=scale*self.orb()[bb]
    if nfft==None: nfft=len(xx)
    xx=xx[0:nfft]
    if(window is not None):
      xx=window(len(xx))*xx
    ff=_np.arange(nfft/2+1)*self.FsamADC/(nfft-1)
    orbfft=_np.fft.rfft(xx)
    return ff,orbfft
  def opt_plot_fft(self,xlog,ylog):
    _pl.xlim(1,100)
    _pl.xlabel(r'f [Hz]',fontsize=16)
    _pl.ylabel(r'amplitude [$\mu$m]')
    _pl.grid(which='both')
    _pl.title(r'$\beta_{IP%s}=%s $ cm, %s'%(self.ip,self.beta['betaIP'+self.ip],dumpdate(self.mtime,fmt='%Y-%m-%d %H:%M:%S')))
    if(xlog):
      _pl.xscale('log')
    if(ylog):
      _pl.yscale('log')
  def opt_plot_dB(self,log):
    _pl.xlim(0,100)
    _pl.xlabel('f [Hz]')
    _pl.ylabel('dB')
    _pl.grid(which='both')
    if(log):
      _pl.xscale('log')
  def plot_fft(self,bb='b1h',nfft=None,window=None,n0=0,offset=0,lbl='',color='b',linestyle='-',xlog=True,ylog=True,scale=1.0):
    """plot the FFT spectrum , where the orbit data is normalized
    with the beta at the bpms:
    window = window function
    nfft   = number of data points used for FFT
    n0     = use data[n0:n0+nfft]
    offset = curve is shifted by offset in plot to distinguish lines
             which would coincide
    """
    ff,orbfft=self.fft(bb=bb,nfft=nfft,window=window,scale=scale)
    _pl.plot(ff[1:],offset+_np.abs(orbfft)[1:],color=color,linestyle=linestyle,label=lbl)#don't plot the DC part
    self.opt_plot_fft(xlog,ylog)
    #if scale !=1 values are scales with 1/sqrt(beta) -> units are mu sqrt(m)
    if abs(scale-1.0)>1.e-6: _pl.ylabel(r'z [$\mu\sqrt{\rm m}$]')
  def plot_fft_dB(self,bb='b1h',nfft=None,window=None,n0=0,offset=0,lbl='',color='b',linestyle='-',log=True):
    """plot the FFT spectrum in dB, where the amplitude is normalized
    in respect of the maximum value:
    fft_1=2*abs(fft(b1h1-b1h2))
    fft=20*log10(fft_1/max(fft_1))
    window = window function
    nfft   = number of data points used for FFT
    n0     = use data[n0:n0+nfft]
    offset = curve is shifted by offset in plot to distinguish lines
             which would coincide
    """
    ff,orbfft=self.abs_fft_dB(bb=bb,nfft=nfft,window=window,n0=n0)
    _pl.plot(ff,offset+orbfft,color=color,linestyle=linestyle,label=lbl)
    self.opt_plot_dB(log)
