from imports import *
from util import *

n1 = 1001
ri = (n1-1)/2
fpeak = 0.1
shiftMin = -15
shiftMax = 15
sampling = Sampling(n1,0.002,0.0)

def main(args):
  trc = Synthetic.getTrace(sampling,fpeak,200)
  s = 2 # resample synthetic trace by this factor
  o = 11 # offset for synthetic trace to warp
  r = 100 # resample factor for DynamicWarpingSlopes
  pt = True   # plot traces
  pe = False  # plot errors
  pa = False  # plot accumulated errors
  pw = False  # plot warped traces
  #doTest(trc,s,o,r,pt,pe,pa,pw)
  myTest()

  #pv = PointsView(sa);
  #pv.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
  #pp = PlotPanel();
  #pp.addTiledView(pv);
  #pf = PlotFrame(pp);
  #pf.setVisible(true);

  pd = True  # plot differences
  #doTest2d(trc,s,o,r,pt,pe,pa,pd)

def myTest():
  n_samplings = 101
  shift_val = 2
  s = Sampling(1001,0.002,0.0)
  fpeak = 30
  seed = 200
  st = zerofloat(n_samplings,1001)
  f = Synthetic.getTrace(s,fpeak,seed)
  st[0] = f;
  for i2 in range(1,n_samplings-1):
    Synthetic.shift(st[i2-1],st[i2],shift_val)
  myPlot(st)

def doTest(trc,s,o,r,pt,pe,pa,pw):
  dw  = DynamicWarping(shiftMin,shiftMax)
  #dws = DynamicWarpingSlopes(shiftMin,shiftMax,r)
  #dws.setStrainMax(0.1)
  trcS = shift(trc,o)           # shift original
  #trcR = resample(trc,s)        # resample original
  #trcSR = resample(trcS,s)      # resample shifted trace
  #trcRI = interpolate(trcR,r)   # interpolate resampled trace
  #trcSRI = interpolate(trcSR,r) # interpolate shifted/resampled trace
  
  # Compute Errors
  e1 = dw.computeErrors(trc,trcS)
  #e2 = dw.computeErrors(trcR,trcSR)
  #es = dws.computeErrors(trcR,trcSR)
  
  # Accumulate Errors
  d1 = dw.accumulateForward(e1)
  #d2 = dw.accumulateForward(e2)
  #ds = dws.accumulateForward(es)

  # Backtrack to get shifts
  u1 = dw.backtrackReverse(d1,e1)
  #u2 = dw.backtrackReverse(d2,e2)
  #us = dws.backtrackReverse(ds,es)
  
  # Apply shifts
  h1 = dw.applyShifts(u1,trcS)
  #h2 = dw.applyShifts(u2,trcSR)
  #hs = dws.applyShifts(us,trcSR)

  si = ri/s
  #print "Original trace shifted by %i, computed resampled shift = %f" % (o,(o/float(s)))
  #print """shift at center sample index:
  #         index= %i: Original Sampling              u1=%f
  #                    Resampled DynamicWarping       u2=%f
  #         index=%i:  Resampled DynamicWarpingSlopes us=%f"""
  #         %(si,u1[si],u2[si],si*r)#(us[si*r])/r)
  if pt:
    plot(trc,   "f"  )
    plot(trcS,  "fs" )
    #plot(trcR,  "fr" )
    #plot(trcRI, "fri")
    #plot(trcSR, "fsr")
    #plot(trcSRI,"fsri")
  if pe:
    plotError(e1,u1,title="Original Sampling")
    #plotError(e2,u2,title="Resampled (%i) DW"%(s))
    #plotError(es,us,title="Resampled(%i) DWS"%(s))
  if pa:
    plotError(d1,u1,title="Accumulated - Original Sampling")
    #plotError(d2,u2,title="Accumulated - Resampled DW")
    #plotError(ds,us,title="Accumulated - Resampled DWS")
  if pw:
    plot(trc,"DW orig",h1)
    #plot(trcR,"DW rsmpl",h2)
    #plot(trcR,"DWS rsmpl",hs)

def doTest2d(trc,s,o,r,pt,pe,pa,pd):
  dw  = DynamicWarping(shiftMin,shiftMax)
  #dws = DynamicWarpingSlopes(shiftMin,shiftMax,r)
  f = make2d(trc,o)
  g = shiftImage(f)
  d = subtract(f,g)
  fr = resample2(f,s)
  gr = resample2(g,s)
  fri = interpolate2(fr,r)
  gri = interpolate2(gr,r)
  
  # Compute Errors
  e1 = dw.computeErrors(f,g)
  e2 = dw.computeErrors(fr,gr)
  #es = dws.computeErrors(fr,gr)

  # Accumulate errors
  a = dws.accumulateForward(es)
  
  # Get shifts
  u1 = dw.findShifts(f,g)
  u2 = dw.findShifts(fr,gr)
  us = dws.backtrackReverse(a,es)

  # Apply shifts and diff them
  h1 = dw.applyShifts(u1,g)
  h2 = dw.applyShifts(u2,gr)
  hs = dws.applyShifts(us,gr)
  d1 = subtract(f,h1)
  d2 = subtract(fr,h2)
  ds = subtract(fr,hs)

  si = ri/s+o
  print "Original trace shifted by %i, actual shift = %f" % (o,(o/float(s)))
  print "shift at trace %i, si %i: u1=%f, u2=%f, us=%f"%(5,si,u1[5][si],u2[5][si],us[5][si*r]/r)
  print "Original min=%f, max=%f" % (min(d1),max(d1))
  print "Resampled DW min=%f, max=%f" % (min(d2),max(d2))
  print "Resampled DWS min=%f, max=%f" % (min(ds),max(ds))
  if pt:
    plotPixels(f,"f")
    plotPixels(g,"g")
    plotPixels(d,"d")
    plotPixels(fr,"fr")
    plotPixels(gr,"gr")
    plotPixels(fri,"fri")
    plotPixels(gri,"gri")
  if pe:
    plotError(e1[5],u=u1[5],title="Original Sampling")
    plotError(e2[5],u=u2[5],title="Resampled - DynamicWarping")
    plotError(es[5],u=us[5],title="Resampled - DynamicWarpingSlopes")
  if pa:
    plotError(a[5],u=us[5],title="Accumulated - DWS")
  if pd:
    #plotPixels(h1,"h1")
    #plotPixels(h2,"h2")
    #plotPixels(hs,"hs")
    #plotPixels(subtract(f,h1),"f-h1",cmin=-1.0,cmax=1.0)
    plotPixels(subtract(fr,h2),"f-h2",cmin=-1.0,cmax=1.0)
    plotPixels(subtract(fr,hs),"f-hs",cmin=-1.0,cmax=1.0)
  
def make2d(trc,offset):
  f = zerofloat(n1,10)
  f[0] = trc
  for i1 in range(1,10):
    f[i1] = shift(f[i1-1],offset)
  return f

def shift(trc,shift):
  n1 = len(trc)
  st = zerofloat(n1)
  for i1 in range(n1):
    si = i1-shift
    if si < n1-1:
      st[i1] = trc[si]
  return st

def shiftImage(f):
  n1,n2 = len(f[0]),len(f)
  g = zerofloat(n1,n2)
  for i2 in range(n2):
    ish = i2+1
    if ish < n2-1:
      g[i2] = f[ish]
    else:
      g[i2] = f[i2]
  return g

def subtract(f,g):
  return sub(f,g)

def interpolate2(f,r):
  n2 = len(f)
  fi = []
  for i2 in range(n2):
    fi.append(interpolate(f[i2],r))
  return fi
    
def interpolate(trc,r):
  n1 = len(trc)
  si = SincInterpolator()
  #si.setUniformSampling(n1,1.0,0.0)
  #si.setUniformSamples(trc)
  nx = n1*r
  dx = 1.0/r
  x = rampfloat(0.0,dx,nx)
  y = zerofloat(nx)
  si.interpolate(n1,1.0,0.0,input array,nx,x,y)
  return y

def resample(f,dx):
  n1 = len(f)
  n1r = int(ceil(float(n1)/dx))
  #print "resample length is",n1r
  g = zerofloat(n1r)
  for i1 in range(n1r):
    g[i1] = f[i1*dx]
  return g

def resample2(f,dx):
  n1,n2 = len(f[0]),len(f)
  g = zerofloat(n1,n2)
  for i2 in range(n2):
    g[i2] = resample(f[i2],dx)
  return g
  
######################################################################
# Plotting

def plot(trc,title,trc2=None):
  pp = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  pp.setTitle(title)
  pp.addPoints(trc)
  pp.addColorBar()
  pp.setColorBarWidthMinimum(100)
  if trc2:
    pv = pp.addPoints(trc2)
    pv.setLineColor(Color.RED)
  pf = PlotFrame(pp)
  pf.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE)
  pf.setSize(200,1200)
  pf.setVisible(True)

def myPlot(test):
  pv = PixelsView(test)
  pv.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
  pp = PlotPanel()
  pp.addTiledView(pv)
  pp.addColorBar()
  pp.setColorBarWidthMinimum(100)
  pf = PlotFrame(pp)
  pf.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE)
  pf.setVisible(True)

def plotPixels(f,title,cmin=None,cmax=None):
  pp = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  pp.setTitle(title)
  pv = pp.addPixels(f)
  if cmin and cmax:
    pv.setClips(cmin,cmax)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pp.addColorBar()
  pf = PlotFrame(pp)
  pf.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE)
  pf.setSize(300,1200)
  pf.setVisible(True)

def plotError(e,u=None,slopeError=None,title=None):
  et = DynamicWarping.transposeLag(e)
  nl = len(et)
  fl = -(nl-1)/2
  dl = 1.0
  if slopeError:
    fl = fl/10
    dl = 0.1
  s2 = Sampling(nl,dl,fl)
  s1 = Sampling(len(et[0]),1.0,0)
  pp = PlotPanel()
  if title:
    pp.setTitle(title)
  pv = pp.addPixels(s1,s2,et)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(ColorMap.JET)
  pp.addColorBar()
  if u:
    ptsv = pp.addPoints(u)
    ptsv.setLineColor(Color.WHITE)
  pf = PlotFrame(pp)
  pf.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE)
  pf.setVisible(True)

######################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())


