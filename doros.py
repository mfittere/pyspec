"""class to import data from DOROS BPMs"""
import numpy as _np
import os

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
  ff = open(fn,'rb')
  header = map(ord,ff.read(LheaderData))
  NumOfFileDumps   = header[0]+header[1]*2**8+header[2]*2**16+header[3]*2**24 #number of files dumped
  NumOfUDPsPerDump = header[4]+header[5]*2**8+header[6]*2**16+header[7]*2**24 #number of UDPs per file
  ff.seek(LheaderData)# remove the Dataheader
  ftype=_np.dtype([('header','%su1'%(LheaderUdp)),('data','%su1'%(LdataUdp))])
  data=_np.fromfile(file=ff,dtype=ftype,count=-1,sep="") #get the full data
  if(NumOfUDPsPerDump!=len(data['header'])):
    print 'WARNING: file corrupted as number of UDP dumps != number of UDPs in file'
    print 'number of UDP dumps:    %s'%(NumOfUDPsPerDump)
    print 'number of UDPs in file: %s'%(len(data['header']))
  #cut away 0 at end of file
  TotalLenOfFile = NumOfUDPsPerDump*NumOfFileDumps*LdataUdp
  return data[0:TotalLenOfFile]

def chan_to_orb(ADCchanTable):
  [bh1,bh2,bv1,bv2]=ADCchanTable[0:4]/(2**24-1)
  bh=(bh1-bh2)/(bh1+bh2)
  bv=(bv1-bv2)/(bv1+bv2)
  return bh,bv

class doros():
  """class to import data from DOROS BPMS"""
  #---- define parameters used for import the datafile ----
  ImSize = 1000
  #pick - up slope in um, theory is 1/4 of the electrode distance [d],
  #for BPMC (d = 49) the database says 12.68 mm, 
  #for BPMSW (d = 61) 15.25 mm
  PUgain=15250
  #ADC=Analog to Digital Converter parameters
  FsamADC = (40*10^6)/3564.0
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
  def __init__(self,b1h=[],b1v=[],b2h=[],b2v=[]):
    self.data={'b1h':_np.array(b1h),'b1v':_np.array(b1v),'b2h':_np.array(b2h),'b2v':_np.array(b2v)}
  @classmethod
  def getdata(cls,fn):
    data=readbinary(fn,cls.LheaderData,cls.LheaderUdp,cls.LdataUdp)
    #extract ADCb[12] buffer
    ADCb1nFrames = data['header'][:,cls.ADCb1nFrameAddr-1]
    ADCb2nFrames = data['header'][:,cls.ADCb2nFrameAddr-1]
    nTotalADCb1BuffDataFrames = sum(ADCb1nFrames)
    nTotalADCb2BuffDataFrames = sum(ADCb2nFrames)
    #number of valid frames per udp ADC1 buffer
    ADCb1Data = _np.hstack(_np.array([ data['data'][idx,0:ADCb1nFrames[idx]*cls.nADCchan*cls.nByte] for idx in range(len(ADCb1nFrames)) ]))
    ADCb2Data = _np.hstack(_np.array([ data['data'][idx,cls.LADCbBuff:cls.LADCbBuff+ADCb2nFrames[idx]*cls.nADCchan*cls.nByte] for idx in range(len(ADCb2nFrames)) ]))
    #decode channel data
    ADCb1dataDec=_np.array([ (2**24*ADCb1Data[4*idx]+2**16*ADCb1Data[4*idx+1]+2**8*ADCb1Data[4*idx+2]+ADCb1Data[4*idx+3]+2**23) & (2**24-1) for idx in range(nTotalADCb1BuffDataFrames*cls.nChan) ],dtype=float)
    ADCb2dataDec=_np.array([ (2**24*ADCb2Data[4*idx]+2**16*ADCb2Data[4*idx+1]+2**8*ADCb2Data[4*idx+2]+ADCb2Data[4*idx+3]+2**23) & (2**24-1) for idx in range(nTotalADCb2BuffDataFrames*cls.nChan) ], dtype=float)
    ADC1chanTable = _np.array([ ADCb1dataDec[idx:idx+nTotalADCb1BuffDataFrames*cls.nChan:cls.SampleDecimation*cls.nChan] for idx in range(cls.nChan) ])
    ADC2chanTable = _np.array([ ADCb2dataDec[idx:idx+nTotalADCb2BuffDataFrames*cls.nChan:cls.SampleDecimation*cls.nChan] for idx in range(cls.nChan) ])
    b1h,b1v=chan_to_orb(ADC1chanTable)
    b2h,b2v=chan_to_orb(ADC2chanTable)
    return cls(b1h,b1v,b2h,b2v)
#    adc=ADCb1Data.view('uint32')
#    LeAll = os.path.getsize(fn)#total file byte length
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
   
#  def getdatadir(dirname=''):
#    """get data from a directory"""
#  
#    return dirname
