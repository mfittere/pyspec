from rdmstores import *
import os as os

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
