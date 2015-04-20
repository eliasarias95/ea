import sys

from java.awt import *
from java.io import *
from java.lang import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.ogl.Gl import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from dnp import *
from slopes import *
from util import *

#pngDir = None
pngDir = "/Users/earias/Documents/research/figures/slides/"
    
seismicDir = "/Users/earias/Home/git/ea/bench/src/util/data/tp/"
#s1 = Sampling(401,0.002,0)
#s2 = Sampling(357,0.025,0)
#s3 = Sampling(161,0.025,0)
s1 = Sampling(251,.002,0)
s2 = Sampling(357,.025,0)
s3 = Sampling(161,.025,0)
ss1 = Sampling(251,1,0)
ss2 = Sampling(357,1,0)
ss3 = Sampling(161,1,0)

pmax = 2.0
n1,n2,n3 = s1.count,s2.count,s3.count
d1,d2,d3 = s1.delta,s2.delta,s3.delta

k1,k2,k3 = 250,0,50; azimuth=240; elevation=20 # for 3D views
fmin,fmax = -5.5,5.5

#fxfile = "tpsz_subz_401_4_400"
fxfile = "tpst_subt_251_4_500"
gtfile = "_flattened"
p2file = "_p2"
p3file = "_p3"
epfile = "_ep"

def main(args):
  #goSlopeLSF()
  #goSlopePWD()
  #goSlopeSDW()
  goFlattenLSF()
  #goFlattenPWD()
  #goFlattenSDW()

def goSlopeLSF():
  f = readImage(fxfile)
  m = readImage("tpmt_subt_251_4_500")
  p2 = copy(f)
  p3 = copy(f)
  ep = copy(f)
  sigma1 = 8.0
  sigma2 = 2.0
  sigma3 = 2.0
  lsf = LocalSlopeFinder(sigma1,sigma2,sigma3,pmax)
  lsf.findSlopesE(f,p2,p3,ep);
  #ep = zerofloat(n1,n2,n3) #comment these lines out for real
  #fill(1.0,ep)             #planarities from structure tensors
  zm = ZeroMask(m)
  zero = 0.00;
  tiny = 0.01;
  zm.apply(zero,p2);
  zm.apply(zero,p3);
  zm.apply(tiny,ep);
  writeImage("lsf"+p2file,p2)
  writeImage("lsf"+p3file,p3)
  writeImage("lsf"+epfile,ep)

def goSlopePWD():
  f = readImage(fxfile)
  p2 = copy(f)
  p3 = copy(f)
  sd = Sfdip(-pmax,pmax)
  rect1 = 55
  rect2 = 6 
  rect3 = 6
  sd.setRect(rect1,rect2,rect3)
  sd.setOrder(2)
  sd.setN4(2)
  sd.findSlopes(s1,s2,s3,f,p2,p3);
  ep = zerofloat(n1,n2,n3)
  fill(1.0,ep)
  writeImage("pwd"+p2file,p2)
  writeImage("pwd"+p3file,p3)
  writeImage("pwd"+epfile,ep)

def goSlopeSDW():
  f = readImage(fxfile)
  m = readImage("tpmt_subt_251_4_500")
  p2 = copy(f)
  p3 = copy(f)
  k  = 10
  r1 = 0.1
  r2 = 0.3
  r3 = 0.3
  h1 = 50.0
  h2 =  9.0
  h3 =  9.0
  ss1 = Sampling(n1);
  ss2 = Sampling(n2);
  ss3 = Sampling(n3);  
  sdw = DynamicWarpingSlopes(k,pmax,h1,h2,h3,r1,r2,r3,ss1,ss2,ss3)
  sdw.findSmoothSlopes(ss1,f,p2,p3)
  ep = zerofloat(n1,n2,n3)
  fill(1.0,ep)
  zm = ZeroMask(m)
  zero = 0.00;
  tiny = 0.01;
  zm.apply(zero,p2);
  zm.apply(zero,p3);
  zm.apply(tiny,ep);
  writeImage("sdw"+p2file,p2)
  writeImage("sdw"+p3file,p3)
  writeImage("sdw"+epfile,ep)

def goFlattenLSF():
  f = readImage(fxfile)
  p2 = readImage("lsf"+p2file)
  p3 = readImage("lsf"+p3file)
  #ep = readImage("lsf"+epfile)
  ep = readImage("lsf_ep_orig")
  p2 = mul(d1/d2,p2)
  p3 = mul(d1/d3,p3)
  ep = pow(ep,6.0)
  s = sub(1.0,ep)
  #fl = Flattener3()
  #fl.setIterations(0.01,100)
  #fm = fl.getMappingsFromSlopes(s1,s2,s3,p2,p3,ep)
  #gt = fm.flatten(f)
  #writeImage("lsf"+gtfile,gt)
  gt = readImage("lsf"+gtfile)
  #plot3(f,g=s,clab="Smoothing",cint=0.2)
  plot3(f,g=ep,clab="Planarity",cint=0.2)#,png="planarity")
  #plot3(gt,clab="Amplitude",cint=2)
  #plot3s(f,clabel="Amplitude",png="before_flattening")
  #plot3s(gt,clabel="Amplitude",png="after_flattenning_lsf_ep=1")
  #plot3s(gt,clabel="Amplitude")

def goFlattenPWD():
  f = readImage(fxfile)
  m = readImage("tpmt_subt_251_4_500")
  p2 = readImage("pwd"+p2file)
  p3 = readImage("pwd"+p3file)
  ep = readImage("pwd"+epfile)
  p2 = mul(d1/d2,p2)
  p3 = mul(d1/d3,p3)
  #fl = Flattener3()
  #fl.setIterations(0.01,100)
  #fm = fl.getMappingsFromSlopes(s1,s2,s3,p2,p3,ep)
  #gt = fm.flatten(f)
  #writeImage("pwd"+gtfile,gt)
  gt = readImage("pwd"+gtfile)
  plot3(gt,clab="Amplitude",cint=2)
  #plot3s(gt,clabel="Amplitude",png="after_flattenning_pwd")
  plot3s(gt,clabel="Amplitude")

def goFlattenSDW():
  f = readImage(fxfile)
  m = readImage("tpmt_subt_251_4_500")
  p2 = readImage("sdw"+p2file)
  p3 = readImage("sdw"+p3file)
  ep = readImage("sdw"+epfile)
  p2 = mul(d1/d2,p2)
  p3 = mul(d1/d3,p3)
  #fl = Flattener3()
  #fl.setIterations(0.01,100)
  #fm = fl.getMappingsFromSlopes(s1,s2,s3,p2,p3,ep)
  #gt = fm.flatten(f)
  #writeImage("sdw"+gtfile,gt)
  gt = readImage("sdw"+gtfile)
  plot3(gt,clab="Amplitude",cint=2)
  #plot3s(gt,clabel="Amplitude",png="after_flattenning_sdw")
  plot3s(gt,clabel="Amplitude")

#############################################################################
# read/write files

def readImage(name):
  fileName = seismicDir+name+".dat"
  n1,n2,n3 = s1.count,s2.count,s3.count
  image = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def writeImage(name,image):
  fileName = seismicDir+name+".dat"
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(image)
  aos.close()
  return image

def readSlice3(name):
  fileName = seismicDir+name+".dat"
  n1,n2 = s1.count,s2.count
  image = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def slice12(k3,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n1,n2)
  SimpleFloat3(f).get12(n1,n2,0,0,k3,s)
  return s

def slice13(k2,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n1,n3)
  SimpleFloat3(f).get13(n1,n3,0,k2,0,s)
  return s

def slice23(k1,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n2,n3)
  SimpleFloat3(f).get23(n2,n3,k1,0,0,s)
  return s

#############################################################################
# graphics

gray = ColorMap.GRAY
jet = ColorMap.JET
def jetFill(alpha):
  return ColorMap.setAlpha(ColorMap.JET,alpha)

def addColorBar(frame, clab=None, cint=None):
  cbar = ColorBar(clab)
  if cint:
    cbar.setInterval(cint)
  cbar.setFont(Font("Arial", Font.PLAIN, 32))  # size by experimenting
  cbar.setWidthMinimum
  cbar.setBackground(Color.WHITE)
  frame.add(cbar, BorderLayout.EAST)
  return cbar

def plot3(f,g=None,cmin=None,cmax=None,
          cmap=None,clab=None,cint=None,png=None):
  n1,n2,n3 = s1.count,s2.count,s3.count
  d1,d2,d3 = s1.delta,s2.delta,s3.delta
  f1,f2,f3 = s1.first,s2.first,s3.first
  l1,l2,l3 = s1.last,s2.last,s3.last
  sf = SimpleFrame(AxesOrientation.XRIGHT_YOUT_ZDOWN)
  cbar = None
  if g==None:
    ipg = sf.addImagePanels(s1,s2,s3,f)
    if cmap!=None:
      ipg.setColorModel(cmap)
    if cmin!=None and cmax!=None:
      ipg.setClips(cmin,cmax)
    else:
      ipg.setClips(fmin,fmax)
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMapListener(cbar)
  else:
    ipg = ImagePanelGroup2(s1,s2,s3,f,g)
    ipg.setClips1(fmin,fmax)
    if cmin!=None and cmax!=None:
      ipg.setClips2(cmin,cmax)
    if cmap==None:
      cmap = jetFill(0.4)
    ipg.setColorModel2(cmap)
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMap2Listener(cbar)
    sf.world.addChild(ipg)
  if cbar:
    cbar.setWidthMinimum(120)
  ipg.setSlices(k1,k2,k3)
  if cbar:
    sf.setSize(1137,900)
  else:
    sf.setSize(1000,900)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)

  ov = sf.getOrbitView()
  zscale = 0.75*max(n2*d2,n3*d3)/(n1*d1)
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setScale(1.0)
  ov.setAzimuth(azimuth)
  ov.setElevation(elevation)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(137,1,pngDir+png+"cbar.png")

def plot3s(s,c=None,clabel="",cmin=0,cmax=0,png=None):
  st = zerofloat(n1,n3,n2)
  Util.transpose23(s,st);
  pp = PlotPanelPixels3(
    PlotPanelPixels3.Orientation.X1DOWN_X2RIGHT,
    PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM,
    ss1,ss3,ss2,st)
  pp.setSlices(k1,k3,k2)
  pp.setLabel1("Samples")
  pp.setLabel2("Crossline (traces)")
  pp.setLabel3("Inline (traces)")
  pp.setClips(fmin,fmax)
  if c:
    cb = pp.addColorBar(clabel)
    #cb.setInterval(1.0)
    pp.setColorBarWidthMinimum(140)
    pp.setLineColor(Color.BLACK)
  else:
    pp.setLineColor(Color.YELLOW)
    cb = pp.addColorBar("Amplitude")
    cb.setInterval(2.0)
  pp.setInterval1(100)
  pp.setInterval2(100)
  pp.setInterval3(100)
  #pp.mosaic.setHeightElastic(1,80)
  #pp.mosaic.setWidthElastic(1,150)
  if c:
    pv12 = PixelsView(ss1,ss2,slice12(k3,c))
    pv12.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    pv13 = PixelsView(ss1,ss3,slice13(k2,c))
    pv13.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    pv23 = PixelsView(ss2,ss3,slice23(k1,c))
    pv23.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP)
    for pv in [pv12,pv13,pv23]:
      pv.setColorModel(ColorMap.getJet(0.5))
      if cmin!=cmax:
        pv.setClips(cmin,cmax)
    pp.pixelsView12.tile.addTiledView(pv12)
    pp.pixelsView13.tile.addTiledView(pv13)
    pp.pixelsView23.tile.addTiledView(pv23)
  pf = PlotFrame(pp)
  pf.setFontSizeForSlide(1.0,0.95,16.0/9.0)
  #pf.setFontSizeForSlide(0.7,0.8,16.0/9.0)
  pf.setSize(960,800)
  pf.setVisible(True)
  if png and pngDir:
    pf.paintToPng(300,6,pngDir+png+".png")

#############################################################################
# Run the function main on the Swing thread
import sys
class _RunMain(Runnable):
  def __init__(self,main):
    self.main = main
  def run(self):
    self.main(sys.argv)
def run(main):
  SwingUtilities.invokeLater(_RunMain(main)) 
run(main)
