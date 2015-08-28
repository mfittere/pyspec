from rdmstores import *
import os as os
import cPickle as pickle
import glob as glob
from adt import getbeta_adt
from doros import getbeta_doros
import matplotlib.pyplot as _pl

def getbeta(dndata,force=False):
  """get a file with the beta* values and
  if it doesn't exist get the data from timber
  and save the values in dndata/betastar.p.
  This function just calls the correct function
  for DOROS/ADT.
  """
  if len(glob.glob(dndata+'/*.mat'))>0:#for ADT data
    beta=getbeta_adt(dndata,force=force)
  if len(glob.glob(dndata+'/*.bin'))>0:#for DOROS data
    beta=getbeta_doros(dndata,force=force)
  return beta
def plot_beta_intensity(dndata,dnpic,force=False):
  """plot beta and intensity, save the beta*
  in dndata/betastar.p and the intensity in
  dndata/beam_intensity.p.
  dndata: folder with ADT/DOROS data files
  dnpic:  folder to store th plots
  """
  #- plot beta*
  beta=getbeta(dndata,force=force)
  start= dumpdate(beta.t1)
  end  = dumpdate(beta.t2)
  print '%s <-> %s'%(start,end)
  _pl.clf()
  beta.plot_2d(vscale=1)
  _pl.savefig(dnpic+'/betastar.png')
  _pl.figure()
  if( force==False and os.path.isfile(dndata+'/beam_intensity.p')):
    bint = pickle.load(open(dndata+'/beam_intensity.p',"rb"))
    print '%s found!'%(dndata+'/beam_intensity.p')
  else:
    try:
      bint=logdb.get(['LHC.BCTFR.A6R4.B1:BEAM_INTENSITY','LHC.BCTFR.A6R4.B2:BEAM_INTENSITY'],start,end)
      pickle.dump(bint,open(dndata+'/beam_intensity.p',"wb"))
      print 'save beam intensity in %s!'%(dndata+'/beam_intensity.p')
    except IOError:
      pass
  if(os.path.isfile(dndata+'/beam_intensity.p')):
    bint.plot_2d()
    _pl.savefig(dnpic+'/beam_intensity.png')
  else:
    print 'WARNING: beam_intensity.p not found and timber not available'
    _pl.close(gcf())
def get_bunch_pattern(dndata,force=False):
  """print the bunch pattern and save the 
  bunch intensities in bunch_intensity.p
  dndata: folder with ADT/DOROS data files
  """
  #- plot beta*
  beta=getbeta(dndata,force=force)
  start= dumpdate(beta.t1)
  end  = dumpdate(beta.t2)
  print '%s <-> %s'%(start,end)
  if( force==False and os.path.isfile(dndata+'/bunch_intensity.p')):
    bdat = pickle.load(open(dndata+'/bunch_intensity.p',"rb"))
    print '%s found!'%(dndata+'/bunch_intensity.p')
  else:
    try:
      bdat=logdb.get(['LHC.BCTFR.A6R4.B1:BUNCH_INTENSITY','LHC.BCTFR.A6R4.B2:BUNCH_INTENSITY'],start,end)
      pickle.dump(bdat,open(dndata+'/bunch_intensity.p',"wb"))
    except IOError:
      pass
  if(os.path.isfile(dndata+'/bunch_intensity.p')):
    for bb in '12':
      bintall=(bdat.data['LHC.BCTFR.A6R4.B'+bb+':BUNCH_INTENSITY'][1]).max(axis=0)
      bint=bintall[np.nonzero(bintall)]
      print bint
      bintmin=bint.min()
      bintmax=bint.max()
      print 'number of bunches b%s      : %4.0f'%(bb,len(bint))
      print 'min/max bunch intensity b%s: %4.2e/%4.2e'%(bb,bintmin,bintmax)
  else:
    print 'WARNING: bunch_intensity.p not found and timber not available'

#---- beta*
def getbetasample(beta,mtime):
  """get beta* at time *mtime*"""
  betasample={}
  if beta == None:
    for ii in ['1','2','5','8']:
      betasample['betaIP'+ii] = 0.0
  else:
    for ii in ['1','2','5','8']:
      idxaux = argmtime(beta.data['HX:BETASTAR_IP'+ii][0],t0=mtime)
      betasample['betaIP'+ii] = beta.data['HX:BETASTAR_IP'+ii][1][idxaux]
  return betasample

#---- sorting files
def mk_cmp_timestamp(get_timestamp=os.path.getmtime):
  def cmp_timestamp(fn1,fn2):
    """compares the timestamp of the files with
    filename *fn1* and *fn2*"""
    if type(get_timestamp(fn1))==str:
      d1=int(1.e3*strpunix(get_timestamp(fn1)))
      d2=int(1.e3*strpunix(get_timestamp(fn2)))
    elif type(get_timestamp(fn1))==float:
      d1=int(1.e3*get_timestamp(fn1))
      d2=int(1.e3*get_timestamp(fn2))
    else:
      print 'get_timestamp must return int or float'
    return d1-d2
  return cmp_timestamp
def mk_sort_files(get_timestamp=os.path.getmtime):
  """returns function *sort_files* which sorts the files
  after the timestamp using the function
  *get_timestamp* to obtain the timestamp"""
  def sort_files(files):
    return sorted(files,cmp=mk_cmp_timestamp(get_timestamp))
  return sort_files
