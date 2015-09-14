import sys

"""
Compares slope methods for given real and synthetic seismic images.
"""
from java.lang import *
from javax.swing import *

from edu.mines.jtk.dsp import *
from edu.mines.jtk.util.ArrayMath import *
from util import *

T = True
F = False
paint = F
fw = 0.8
fh = 0.8
PATH = "/Users/earias/Home/git/ea/bench/src/test/"

def main(args):
  #triton3D()
  santos3D()

def santos3D():
  n1 = 501
  n2 = 201
  n3 = 201
  cm = 2.0 #colorbar limits
  f = Util.readImageL(n1,n2,n3,PATH+"Santos_small_3D.dat")
  p2_sdw = Util.readImageL(n1,n2,n3,PATH+"Santos_small_3D_sdw2_dzdy.dat")
  p2_sxy = Util.readImageL(n1,n2,n3,PATH+"Santos_small_3D_sxy_dzdy.dat")

  plot3Ds(n1,n2,n3,f,p2_sdw,p2_sxy,cm,"Santos_sub_3D_sdw") # sdw

def triton3D():
  n1 = 1201
  n2 = 1561
  n3 = 251
  cm = 2.0 #colorbar limits
  f = Util.readImageL(n1,n2,n3,PATH+"Triton_small_3D.dat")
  p2_sdw = Util.readImageL(n1,n2,n3,PATH+"Triton_small_3D_sdw1_dzdx.dat")
  p2_sxy = Util.readImageL(n1,n2,n3,PATH+"Triton_small_3D_sxy_dzdx.dat")

  plot3Ds(n1,n2,n3,f,p2_sdw,p2_sxy,cm,"Triton_small_3D")

def plot3Ds(n1,n2,n3,f,p2,p3,cm,ttl):
  s1 = Sampling(n1)
  s2 = Sampling(n2)
  s3 = Sampling(n3)
  Plot.plot(s1,s2,s3,f,p2,ttl+"_p2_slices",-cm,cm,paint);
  Plot.plot(s1,s2,s3,f,p3,ttl+"_p3_slices",-cm,cm,paint);
#############Run the function main on the Swing thread#############
import sys
class _RunMain(Runnable):
  def __init__(self,main):
    self.main = main
  def run(self):
    self.main(sys.argv)
def run(main):
  SwingUtilities.invokeLater(_RunMain(main)) 
run(main)
