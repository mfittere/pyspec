"""class to import data from DOROS BPMs"""
import numpy as _np
import matplotlib.pyplot as _pl
import os
import cPickle as pickle
import glob as glob

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
  FsamADC = (40*1.e6)/3564.0
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
  def __init__(self,b1h1=[],b1h2=[],b1v1=[],b1v2=[],b2h1=[],b2h2=[],b2v1=[],b2v2=[]):
    self.data={'b1h1':_np.array(b1h1),'b1h2':_np.array(b1h2),'b1v1':_np.array(b1v1),'b1v2':_np.array(b1v2),'b2h1':_np.array(b2h1),'b2h2':_np.array(b2h2),'b2v1':_np.array(b2v1),'b2v2':_np.array(b2v2)}
  @classmethod
  def getdata(cls,fn,force=False):
    if(fn.split('.')[-1]=='bin'):
      fn=fn[0:-4]
    if(fn.split('.')[-1]=='p'):
      fn=fn[0:-2]
    if(os.path.isfile(fn+'.p') and force==False):
      b1h1,b1h2,b1v1,b1v2,b2h1,b2h2,b2v1,b2v2=pickle.load(open(fn+'.p',"rb"))
    else:
      data,udpcheck=readbinary(fn+'.bin',cls.LheaderData,cls.LheaderUdp,cls.LdataUdp)
      if(udpcheck):#only process data if udp check passed
        ADC1chanTable=decode_chan(data,cls.ADCb1nFrameAddr,cls.nADCchan,cls.nByte,cls.LADCbBuff,cls.nChan,cls.SampleDecimation,beam='b1')
        ADC2chanTable=decode_chan(data,cls.ADCb2nFrameAddr,cls.nADCchan,cls.nByte,cls.LADCbBuff,cls.nChan,cls.SampleDecimation,beam='b2')
        #number of valid frames per udp ADC1 buffer
        [b1h1,b1h2,b1v1,b1v2]=ADC1chanTable[0:4]/(2**24-1)
        [b2h1,b2h2,b2v1,b2v2]=ADC2chanTable[0:4]/(2**24-1)
        #store already processed orbit data in *.p
        pickle.dump([b1h1,b1h2,b1v1,b1v2,b2h1,b2h2,b2v1,b2v2],open(fn+'.p',"wb"))  
        print '... store b1h1,b1h2,b2h1,b2h2 etc. in file %s.p for faster reload'%(fn.split('/')[-1])
      else:
        [b1h1,b1h2,b1v1,b1v2,b2h1,b2h2,b2v1,b2v2]=[[] for x in range(8)]
    return cls(b1h1,b1h2,b1v1,b1v2,b2h1,b2h2,b2v1,b2v2)
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
  def orb(self):
    """returns a dictionary with the orbit b1h,b1v,b2h,b2v in mum,
    where the orbit is defined by the signlas by e.g.:
    b1h=(b1h1-b1h2)/(b1h1+b1h2)"""
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
  def plot_fft_max(self,bb='b1h',NFFT=None):
    """plot the fft spectrum in dB, where the amplitude is normalized
    in respect of the maximum value:
    fft_1=2*abs(fft(b1h1-b1h2))
    fft=20*log10(fft_1/max(fft_1))
    """
    bz1=bb+'1'
    bz2=bb+'2'
    if(NFFT==None):
      if(len(self.data[bz1])==len(self.data[bz2])):
        NFFT=len(self.data[bz1]) #take all datapoints for the FFT 
        print 'NFFT='+str(NFFT)
      else:
        print 'WARNING: length of '+bz1+' and '+bz2+'differ! Use len('+bz1+') for FFT'
    p1=self.data[bz1][0:NFFT]
    p2=self.data[bz2][0:NFFT]
    orbfft=2*_np.abs(_np.fft.rfft(p1-p2))
    #normalize to maximum signal
    orbfft=20*_np.log10(orbfft/(orbfft.max()))
    #scale from data points to frequency [Hz]
    FFtscale=self.FsamADC/(NFFT-1)
    ff=_np.arange(NFFT/2+1)*FFtscale
    _pl.plot(ff,orbfft)
    _pl.xlabel('f [Hz]')
    _pl.ylabel('dB')
    _pl.grid()
