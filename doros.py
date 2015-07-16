"""class to import data from DOROS BPMs"""
import numpy as _np
import matplotlib.pyplot as _pl
import os,time
import cPickle as pickle
import glob as glob
import spec as spec
from scipy.signal import welch
import localdate

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
  if(fn.split('IP')[1].split('_Data')[0]=='172_18_66_233'):
    print '... read data from frontend 1 (IR1 right)'
  if(fn.split('IP')[1].split('_Data')[0]=='172_18_66_234'):
    print '... read data from frontend 2 (IR1 left)'
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
  return data,udpcheck

def decode_chan(data,ADCbnFrameAddr,nADCchan,nByte,LADCbBuff,nChan,SampleDecimation,beam='b1'):
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

class doros():
  """class to import data from DOROS BPMS"""
  #---- define parameters used for import the datafile ----
  ImSize = 1000
  #pick - up slope in um, theory is 1/4 of the electrode distance [d],
  #for BPMC (d = 49) the database says 12.68 mm, 
  #for BPMSW (d = 61) 15.25 mm
  PUgain=15250
  #ADC=Analog to Digital Converter parameters
  FsamADC = (40*1.e6)/3564.0#sampling frequency
  Vref = 2.5
  FullScale = 2^24
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
  def __init__(self,b1h1=[],b1h2=[],b1v1=[],b1v2=[],b2h1=[],b2h2=[],b2v1=[],b2v2=[],mtime=None):
    self.mtime=mtime
    self.data={'b1h1':_np.array(b1h1),'b1h2':_np.array(b1h2),'b1v1':_np.array(b1v1),'b1v2':_np.array(b1v2),'b2h1':_np.array(b2h1),'b2h2':_np.array(b2h2),'b2v1':_np.array(b2v1),'b2v2':_np.array(b2v2)}
  @classmethod
  def getdata(cls,fn,force=False):
    """get the data from the binary file fn.bin and save it
    it in a pickle file fn.p. Add timestamp from creation
    time of binary file. If pickle file already exists, just
    load the pickle file"""
    if(fn.split('.')[-1]=='bin'):
      fn=fn[0:-4]
    if(fn.split('.')[-1]=='p'):
      fn=fn[0:-2]
    if(os.path.isfile(fn+'.p') and force==False):
      b1h1,b1h2,b1v1,b1v2,b2h1,b2h2,b2v1,b2v2,mtime=pickle.load(open(fn+'.p',"rb"))
    else:
      mtime=os.path.getmtime(fn+'.bin')#ctime gives back the timestamp when the file was created for MAC, instead use mtime, which seems to give back the correct timestamp
      data,udpcheck=readbinary(fn+'.bin',cls.LheaderData,cls.LheaderUdp,cls.LdataUdp)
      if(udpcheck):#only process data if udp check passed
        ADC1chanTable=decode_chan(data,cls.ADCb1nFrameAddr,cls.nADCchan,cls.nByte,cls.LADCbBuff,cls.nChan,cls.SampleDecimation,beam='b1')
        ADC2chanTable=decode_chan(data,cls.ADCb2nFrameAddr,cls.nADCchan,cls.nByte,cls.LADCbBuff,cls.nChan,cls.SampleDecimation,beam='b2')
        #number of valid frames per udp ADC1 buffer
        [b1h1,b1h2,b1v1,b1v2]=ADC1chanTable[0:4]/(2**24-1)
        [b2h1,b2h2,b2v1,b2v2]=ADC2chanTable[0:4]/(2**24-1)
        #store already processed orbit data in *.p
        pickle.dump([b1h1,b1h2,b1v1,b1v2,b2h1,b2h2,b2v1,b2v2,mtime],open(fn+'.p',"wb"))  
        print '... store b1h1,b1h2,b2h1,b2h2 etc. in file %s.p for faster reload'%(fn.split('/')[-1])
      else:
        mtime=0
        [b1h1,b1h2,b1v1,b1v2,b2h1,b2h2,b2v1,b2v2]=[[] for x in range(8)]
    return cls(b1h1,b1h2,b1v1,b1v2,b2h1,b2h2,b2v1,b2v2,mtime)
  def acqutime(self,frev=11245.0):
    """returns the acquisition time in s assuming a revolution frequency of frev"""
    ll=len(self.data['b1h1'])
    lflag=True
    for bb in 'b1','b2':
      for pp in 'h1','v1','h2','v2':
        if((len(self.data[bb+pp])-ll)>1.e-3):
          lflag=False
          print '%s: data acquisition of %4.0f samples over %4.2f seconds'%(bb+pp,len(self.data[bb+pp]),len(self.data[bb+pp])/frev)
    if(lflag):
      print 'data acquisition over of %4.0f samples over %4.2f seconds'%(ll,ll/frev)
    else:
      print 'WARNING: not all b[12][hv][12] arrays have the same length!' 
  def checkdata(self,fn):
    """check for errors in data acquisition"""
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
  @classmethod
  def process_dir(cls,dn,force=False):
    """process orbit data in directory dn and
    store it in *.p files"""
    for fn in glob.glob(dn+'/*.bin'):
      if(os.path.isfile(fn[0:-4]+'.p') and force==False):
        print fn+' is already processed'
      else:
        print 'processing '+fn
        cls.getdata(fn,force=force)
  def dumpdate(self,fmt='%Y-%m-%d %H:%M:%S.SSS'):
    return localdate.dumpdate(self.mtime,fmt=fmt) 
  def orb(self):
    """returns a dictionary with the orbit b1h,b1v,b2h,b2v in mum,
    where the orbit is defined by the signals by:
    b1h=(b1h1-b1h2)/(b1h1+b1h2) etc."""
    #number of valid frames per udp ADC1 buffer
    dic={}
    for bb in 'b1','b2':
      for pp in 'h','v':
        dic[bb+pp]=(self.data[bb+pp+'1']-self.data[bb+pp+'2'])/(self.data[bb+pp+'1']+self.data[bb+pp+'2'])
    return dic
  def plot_orb(self,refidx=None):
    """plot orbit position change in mum
    if refidx=None: plot e.g. data['b1h']
    if refidx=7: plot e.g. data['b1h']-data['b1h'][7]
    """
    f,ax=_pl.subplots(2)
    for bb in [1,2]:
      for pp in 'h','v':
        if(refidx==None):
          ax[bb-1].plot(self.PUgain*self.orb()['b'+str(bb)+pp],label=('b'+str(bb)+pp).upper()) 
        else:
          ax[bb-1].plot(self.PUgain*(self.orb()['b'+str(bb)+pp]-self.orb()['b'+str(bb)+pp][refidx]),label=('b'+str(bb)+pp).upper())
        ax[bb-1].set_xlabel('number of turns') 
        ax[bb-1].set_ylabel(r'z [$\mu$m]') 
        ax[bb-1].legend()
  def psd(self,bb='b1h',nfft=None,n0=0,window=_np.hanning):
    """calculate the PSD [m**2/Hz]. For correct normalisation
    see CERN-THESIS-2013-028. This function calls spec.psd().
    window = window function
    nfft   = number of data points used for FFT
    n0     = use data[n0:n0+nfft]
    for the psd the average value is subtracted from the data:
      data=data-mean(data)
    """
    xx=self.PUgain*(self.data[bb+'1']-self.data[bb+'2'])/(self.data[bb+'1']+self.data[bb+'2']) #orbit in mum
    xx=xx-_np.mean(xx)#substract average value from data
    return spec.psd(data=xx,nfft=nfft,n0=n0,window=window,fs=self.FsamADC)
  def psd_welch(self,bb='b1h',n0=0,n1=None,window=_np.hanning,nperseg=4096,noverlap=None):
    """calculate the PSD [m**2/Hz] using the Welche method.
    For more information of input parameters see 
    scipy.signal.welch.
    For correct normalisation see CERN-THESIS-2013-028
    n0,n1: use data[n0:n1]
    """
    xx=self.PUgain*(self.data[bb+'1']-self.data[bb+'2'])/(self.data[bb+'1']+self.data[bb+'2']) #orbit in mum
    if(n1==None):
      n1=len(xx)-n0
    xx=xx[n0:n1]
    ff,psd=welch(xx,fs=self.FsamADC,window=window(nperseg),nperseg=nperseg,noverlap=noverlap,nfft=None,detrend=False,return_onesided=True,scaling='density')
    return ff,psd
  def plot_psd_welch(self,bb='b1h',n0=0,n1=None,window=_np.hanning,nperseg=4096,noverlap=None,scale=1,lbl='',color='b',linestyle='-',xlog=True,ylog=True):
    """plot the PSD spectrum in m**2/Hz using the welch method
    window = window function (default: hanning)
    n0     = use data[n0:n1]
    offset = curve is shifted by offset in plot to distinguish lines
             which would coincide
    """
    ff,psd=self.psd_welch(bb=bb,n0=n0,n1=n1,window=window,nperseg=nperseg,noverlap=noverlap)
    _pl.plot(ff[1:],scale*psd[1:],color=color,linestyle=linestyle,label=lbl)#do not plot DC offset
    self.opt_plot_psd(xlog,ylog)
  def opt_plot_psd(self,xlog,ylog):
    _pl.xlim(1,100)
    _pl.xlabel(r'f [Hz]',fontsize=16)
    _pl.ylabel(r'PSD [$\mathrm{\mu m}^2$/Hz]',fontsize=14)
    _pl.grid(which='both')
    if(xlog):
      _pl.xscale('log')
    if(ylog):
      _pl.yscale('log')
  def plot_psd(self,bb='b1h',nfft=None,window=None,n0=0,scale=1,lbl='',color='b',linestyle='-',xlog=True,ylog=True):
    """plot the PSD spectrum in m**2/Hz
    window = window function
    nfft   = number of data points used for FFT
    n0     = use data[n0:n0+nfft]
    offset = curve is shifted by offset in plot to distinguish lines
             which would coincide
    """
    ff,psd=self.psd(bb=bb,nfft=nfft,window=window,n0=n0)
    _pl.plot(ff[1:],scale*psd[1:],color=color,linestyle=linestyle,label=lbl)#do not plot DC offset
    self.opt_plot_psd(xlog,ylog)
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
  def opt_plot(self,log):
    _pl.xlim(0,100)
    _pl.xlabel('f [Hz]')
    _pl.ylabel('dB')
    _pl.grid(which='both')
    if(log):
      _pl.xscale('log')
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
    self.opt_plot(log)
