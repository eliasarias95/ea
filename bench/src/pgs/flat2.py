import sys

from java.nio import *
from java.lang import *
from javax.swing import *

from edu.mines.jtk.dsp import Sampling
from edu.mines.jtk.io import *
from edu.mines.jtk.util.ArrayMath import *

from dnp import Flattener2
from slopes import *
from util import *

pngDir = System.getProperty("user.home")+"/Documents/research/figures/testing/"
seismicDir = System.getProperty("user.home")+"/Home/git/ea/bench/src/pgs/data/"

#n1,n2 = 2001,38 #offset_gather_Triton
#n1,n2 = 1001,81 #999
n1,n2 = 1601,41 #000
#n1,n2 = 1601,46 #EK
s1,s2 = Sampling(n1),Sampling(n2)

pmax = 4.0
T = True
F = False
paint = F
fw = 0.8
fh = 0.8

#fxfile = "offset_gather_Triton"
#fxfile = "bad_gather_999"
#fxfile = "good_gather_999"
fxfile = "good_gather_000"
#fxfile = "EK_Triton_azang_nears_test_1ang_az4"
#fxfile = "EK_Triton_ang_nears_test_1ang"
gtfile = "_flattened"
pfile = "_p"
elfile = "_el"

def main(args):
  goSlopeSDW()
  goFlattenSDW()

def goSlopeSDW():
  f = readImage(fxfile)
  p  = copy(f)
  k  = 10
  r1 = 0.1
  r2 = 0.5
  h1 = 5.0
  h2 = 2.0
  sdw = DynamicWarpingSlopes(k,pmax,h1,h2,r1,r2,s1,s2)
  sdw.findSmoothSlopes(s1,f,p)
  el = zerofloat(n1,n2)
  fill(1.0,el)
  writeImage("sdw"+pfile,p)
  writeImage("sdw"+elfile,el)
  cm = 2.5
  plot2(f,p,cm,"sdw_slope")

def goFlattenSDW():
  f = readImage(fxfile)
  p  = readImage("sdw"+pfile)
  el = readImage("sdw"+elfile)
  fl = Flattener2()
  fl.setIterations(0.01,100)
  fm = fl.getMappingsFromSlopes(s1,s2,p,el)
  gt = fm.flatten(f)
  writeImage("sdw"+gtfile,gt)
  gt = readImage("sdw"+gtfile)
  cm = 40000000
  cm = 3000
  plot1(f,cm,"sdw_before_flattening")
  plot1(gt,cm,"sdw_after_flattening")

#############################################################################
# read/write files

def readImage(name):
  bo = ByteOrder.LITTLE_ENDIAN
  fileName = seismicDir+name+".dat"
  image = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName,bo)
  ais.readFloats(image)
  ais.close()
  return image

def writeImage(name,image):
  bo = ByteOrder.LITTLE_ENDIAN
  fileName = seismicDir+name+".dat"
  aos = ArrayOutputStream(fileName,bo)
  aos.writeFloats(image)
  aos.close()
  return image

#############################################################################
# graphics

def plot1(f,cm,ttl):
  cbl = "Amplitude" #colorbar label
  #clip,interp,ttl,paint,color,slide,one
  Plot.plot(s1,s2,f,ttl,cbl,fw,fh,-cm,cm,T,F,F,paint,F,T,T);

def plot2(f,p,cm,ttl):
  cbl = "slope (samples/trace)" #colorbar label
  #clip, title, paint, slide, no. columns
  Plot.plot(s1,s2,f,p,ttl,cbl,fw,fh,-cm,cm,T,F,paint,T,T)

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
