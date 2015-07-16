import time,datetime

def dumpdate(t,fmt='%Y-%m-%d %H:%M:%S.SSS'):
  """converts unix time to locale time"""
  ti=int(t)
  tf=t-ti
  s=time.strftime(fmt,time.localtime(t))
  if 'SSS' in s:
    s=s.replace('SSS','%03d'%(tf*1000))
  return s

def addtimedelta(t,dt):
  """add *dt* seconds to t in unix time"""
  dt=datetime.timedelta(seconds=dt)
  t = datetime.datetime.fromtimestamp(t)+dt
  return time.mktime(t.timetuple())
