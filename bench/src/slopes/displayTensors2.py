import sys

"""
Computes structure eigen-tensors.
"""
from java.awt import *
from java.io import *
from java.lang import *
from javax.swing import *

from edu.mines.jtk.io import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.awt import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util.ArrayMath import *
from slopes import *

k1,k2,k3 = 250,0,32
n1 = 251
n2 = 357
n3 = 161
s1 = Sampling(n1)
s2 = Sampling(n2)
s3 = Sampling(n3)
PATH = "/Users/earias/Home/git/ea/bench/src/util/"
efile = "data/tp/tpet.dat" # eigen-tensors

def main(args):
  makeTensors2()
  display2()
  #makeTensors3()
  #display3()

def makeTensors2():
  f = read2(n1,n2,"data/tp/tp73.dat")
  sigma1 = 8.0
  sigma2 = 2.0
  lof = LocalOrientFilter(sigma1,sigma2)
  e = lof.applyForTensors(f)
  e.invertStructure(0.0,1.0)
  writeTensors(efile,e)

def display2():
  f = read2(n1,n2,"data/tp/tp73.dat")
  et = readTensors(efile)
  plot(s1,s2,f,et)

def makeTensors3():
  f = read3(n1,n2,n3,"data/tp/tpst_subt_251_4_500.dat")
  m = read3(n1,n2,n3,"data/tp/tpmt_subt_251_4_500.dat")
  sigma1 = 8.0
  sigma2 = 2.0
  sigma3 = 2.0
  mask = ZeroMask(m)
  lof = LocalOrientFilter(sigma1,sigma2,sigma3)
  e = lof.applyForTensors(f)
  e.invertStructure(0.0,2.0)
  mask.apply((1.0,0.0,0.0,1.0,0.0,1.0),e)
  writeTensors(efile,e)

def display3():
  f = read3(n1,n2,n3,"data/tp/tpst_subt_251_4_500.dat")
  et = readTensors(efile)
  world = World()
  ipg = addImageToWorld(world,f)
  ipg.setClips(-5.5,5.5)
  ipg.setSlices(k1,k2,k3)
  addTensorsInImage(ipg.getImagePanel(Axis.X),et,40)
  addTensorsInImage(ipg.getImagePanel(Axis.Y),et,40)
  addTensorsInImage(ipg.getImagePanel(Axis.Z),et,40)
  frame = makeFrame(world)
  frame.setSize(1460,980)
  frame.orbitView.setAzimuth(240.0)
  frame.orbitView.setElevation(20.0)
  background = Color(254,254,254)
  frame.viewCanvas.setBackground(background)

def read2(n1,n2,name):
  """ 
  Reads an image from a file with specified name.
  name: base name of image file; e.g., "tpsz"
  """
  fileName = PATH+name
  image = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def read3(n1,n2,n3,name):
  """ 
  Reads an image from a file with specified name.
  name: base name of image file; e.g., "tpsz"
  """
  fileName = PATH+name
  image = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

from org.python.util import PythonObjectInputStream
def readTensors(name):
  """
  Reads tensors from file with specified basename; e.g., "tpet".
  """
  fis = FileInputStream(PATH+name)
  ois = PythonObjectInputStream(fis)
  tensors = ois.readObject()
  fis.close()
  return tensors

def writeTensors(name,tensors):
  """
  Writes tensors to file with specified basename; e.g., "tpet".
  """
  fos = FileOutputStream(PATH+name)
  oos = ObjectOutputStream(fos)
  oos.writeObject(tensors)
  fos.close()

#############################################################################
# graphics

def plot(s1,s2,f,et): 
  pp = PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT)
  pv1 = pp.addPixels(s1,s2,f)
  pv2 = TensorsView(s1,s2,et)
  pv2.setOrientation(TensorsView.Orientation.X1DOWN_X2RIGHT)
  pv2.setLineColor(Color.CYAN)
  pv2.setLineWidth(3)
  pv2.setScale(0.7)
  pv2.setEllipsesDisplayed(12)
  pv1.setInterpolation(PixelsView.Interpolation.NEAREST)
  pp.addTiledView(pv1)
  pp.addTiledView(pv2)
  pp.setHLabel("Traces")
  pp.setVLabel("Samples")
  pp.addColorBar("Amplitude")
  pp.setColorBarWidthMinimum(120)
  pf = PlotFrame(pp)
  pf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE)
  pf.setSize(850,600)
  pf.setVisible(True)
  dpi = 720
  pf.setFontSizeForPrint(8.0,222.0)
  pf.paintToPng(720,3.08,"/Users/earias/Documents/tthesis/figures/tensors2")

def addImageToWorld(world,image):
  ipg = ImagePanelGroup(s1,s2,s3,image)
  world.addChild(ipg)
  return ipg

def addImage2ToWorld(world,image1,image2):
  ipg = ImagePanelGroup2(s1,s2,s3,image1,image2)
  ipg.setColorModel1(ColorMap.getGray())
  ipg.setColorModel2(ColorMap.getJet(0.3))
  world.addChild(ipg)
  return ipg

def addTensorsInImage(ip,et,esize):
  tp = TensorsPanel(s1,s2,s3,et)
  tp.setEllipsoidSize(esize)
  ip.getFrame().addChild(tp)
  return tp

def makeFrame(world):
  n1,n2,n3 = s1.count,s2.count,s3.count
  d1,d2,d3 = s1.delta,s2.delta,s3.delta
  f1,f2,f3 = s1.first,s2.first,s3.first
  l1,l2,l3 = s1.last,s2.last,s3.last
  frame = SimpleFrame(world)
  view = frame.getOrbitView()
  zscale = 0.75*max(n2*d2,n3*d3)/(n1*d1)
  view.setAxesScale(1.0,1.0,zscale)
  view.setScale(1.3)
  #view.setAzimuth(75.0)
  #view.setAzimuth(-75.0)
  view.setAzimuth(-65.0)
  view.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  frame.viewCanvas.setBackground(frame.getBackground())
  frame.setSize(1250,900)
  frame.setVisible(True)
  return frame

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

#############################################################################
#run(main)
