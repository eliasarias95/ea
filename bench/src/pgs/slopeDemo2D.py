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
cm = 4.0 #colorbar limits
fw = 0.8
fh = 0.8
PATH = "/users/elias.arias/Home/git/ea/bench/src"

n1 = 501
n2 = 501
n3 = 1
#f = Util.readImageL(n1,n2,PATH+"/pgs/data/constant2D.dat")
#p_sdw = Util.readImageL(n1,n2,PATH+"/pgs/data/constant2D_sdw_out.dat")
#p_sxy = Util.readImageL(n1,n2,PATH+"/pgs/data/constant2D_sxy_out_dzdx.dat")
f = Util.readImageL(n1,n2,PATH+"/pgs/data/complex2D.dat")
p_sdw = Util.readImageL(n1,n2,PATH+"/pgs/data/complex2D_sdw_out.dat")
p_sxy = Util.readImageL(n1,n2,PATH+"/pgs/data/complex2D_sxy_out_dzdx.dat")
s1 = Sampling(n1)
s2 = Sampling(n2)
s3 = Sampling(n3)

def main(args):
  plot2D(f,p_sdw) # sdw_c
  plot2D(f,p_sxy) # slopexy
  #print "mean slopexy=", sum(p_sxy)/(n1*n2)

def plot2D(f,p):
  cbl = "slope (samples/trace)" #colorbar label
  #clip, title, paint, slide, no. columns
  Plot.plot(s1,s2,f,p,"2D plot",cbl,fw,fh,-cm,cm,T,F,F,T,one)

def plot3D(f,p):
  Plot.plot(s1,s2,s3,f,p,"slope (samples/trace)","p_panels",fw,fh,-cm,cm,F,T,
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
