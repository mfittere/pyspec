from rdmstores import *

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
