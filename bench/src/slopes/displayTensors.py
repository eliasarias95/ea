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
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util.ArrayMath import *
from slopes import *

#setupForSubset("subz_401_4_600")
k3,k2,k1 = 67,320,210
n1 = 251
n2 = 357
n3 = 161
s1 = Sampling(n1)
s2 = Sampling(n2)
s3 = Sampling(n3)
PATH = "/Users/earias/Home/git/ea/bench/src/util/"

efile = "data/tp/tpet.dat" # eigen-tensors

def main(args):
  #makeTensors()
  display()

def makeTensors():
  f = readImage(n1,n2,n3,"data/tp/tpst_subt_251_4_500.dat")
  m = readImage(n1,n2,n3,"data/tp/tpmt_subt_251_4_500.dat")
  sigma1 = 8.0
  sigma2 = 2.0
  sigma3 = 2.0
  mask = ZeroMask(m)
  lof = LocalOrientFilter(sigma1,sigma2,sigma3)
  e = lof.applyForTensors(f)
  e.invertStructure(0.0,1.0,1.0)
  mask.apply((1.0,0.0,0.0,0.0,1.0,0.0),e)
  writeTensors(efile,e)

def display():
  f = readImage(n1,n2,n3,"data/tp/tpst_subt_251_4_500.dat")
  et = readTensors(efile)
  world = World()
  ipg = addImageToWorld(world,f)
  ipg.setClips(-5.5,5.5)
  ipg.setSlices(k1,k2,k3)
  addTensorsInImage(ipg.getImagePanel(Axis.X),et,30)
  addTensorsInImage(ipg.getImagePanel(Axis.Y),et,30)
  addTensorsInImage(ipg.getImagePanel(Axis.Z),et,30)
  frame = makeFrame(world)
  frame.setSize(850,700)
  frame.orbitView.setAzimuth(220.0)
  frame.orbitView.setElevation(20.0)
  frame.orbitView.setScale(1.2)
  background = Color(254,254,254)
  frame.viewCanvas.setBackground(background)
  frame.paintToFile("/Users/earias/Documents/tthesis/figures/tensors3.png")

def readImage(n1,n2,n3,name):
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

def writeImage(name,image):
  """ 
  Writes an image to a file with specified name.
  name: base name of image file; e.g., "tpgp"
  image: the image
  """
  fileName = PATH+name
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(image)
  aos.close()
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
