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
one = T
paint = F
fw = 0.8
fh = 0.8
PATH = "/users/elias.arias/Home/git/ea/bench/src"

n1 = 501
n2 = 501
n3 = 1
s1 = Sampling(n1)
s2 = Sampling(n2)
s3 = Sampling(n3)

def main(args):
  #constant2D()
  complex2D()

def constant2D():
  cm = 1.0 #colorbar limits
  f = Util.readImageL(n1,n2,PATH+"/pgs/data/constant2D.dat")
  p_sdw = Util.readImageL(n1,n2,PATH+"/pgs/data/constant2D_sdw_out.dat")
  p_sxy = Util.readImageL(n1,n2,PATH+"/pgs/data/constant2D_sxy_out_dzdx.dat")
  plot2D(f,p_sdw,cm,"constant2D_sdwc") # sdw_c
  plot2D(f,p_sxy,cm,"constant2D_slopexy") # slopexy
  print "mean slopexy=", sum(p_sxy)/(n1*n2)

def complex2D():
  cm = 4.0 #colorbar limits
  f = Util.readImageL(n1,n2,PATH+"/pgs/data/complex2D.dat")
  p_sdwc = Util.readImageL(n1,n2,PATH+"/pgs/data/complex2D_sdw_out.dat")
  p_sdwj = Util.readImage(n1,n2,PATH+"/pgs/data/sdw_half_complex2D.dat")
  p_sxy  = Util.readImageL(n1,n2,PATH+"/pgs/data/complex2D_sxy_out_dzdx.dat")
  plot2D(f,p_sdwc,cm,"complex2D_sdwc") # sdwc
  plot2D(f,p_sdwj,cm,"complex2D_sdwj") # sdwj
  plot2D(f,p_sxy,cm,"complex2D_slopexy") # slopexy


def plot2D(f,p,cm,ttl):
  cbl = "slope (samples/trace)" #colorbar label
  #clip, title, paint, slide, no. columns
  Plot.plot(s1,s2,f,p,ttl,cbl,fw,fh,-cm,cm,T,F,paint,T,one)

def plot3D(f,p,cm,ttl):
  Plot.plot(s1,s2,s3,f,p,"slope (samples/trace)",ttl,fw,fh,-cm,cm,paint,T,
      one);
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
