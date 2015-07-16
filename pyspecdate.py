import time,datetime

def dumpdate(t,fmt='%Y-%m-%d %H:%M:%S.SSS'):
  """converts unix time to locale time"""
  ti=int(t)
  tf=t-ti
  s=time.strftime(fmt,time.localtime(t))
  if 'SSS' in s:
    s=s.replace('SSS','%03d'%(tf*1000))
  return s

def addtimedelta(t,s,fmt='%Y-%m-%d %H:%M:%S.SSS'):
  """add *m* minutes and *s* seconds to 
  time t: tnew = t + m + s"""
  fmtdatetime='%Y-%m-%d %H:%M:%S.%f'
  dt=datetime.timedelta(seconds=s)
  s = datetime.datetime.strptime(t,fmt)
  return (datetime.datetime.strptime(start,fmtdatetime)).strftime(fmt)
