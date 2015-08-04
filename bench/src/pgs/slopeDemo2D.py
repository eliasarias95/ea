import sys

"""
Compares slope methods for given real and synthetic seismic images.
"""
from java.lang import *
from javax.swing import *

from edu.mines.jtk.dsp import *
from edu.mines.jtk.util.ArrayMath import *
from util import *

cmax = 1.0
PATH = "/users/elias.arias/Home/git/ea/bench/src"

n1 = 21
n2 = 15
n3 = 1
f = Util.readImageL(n1,n2,n3,PATH+"/pgs/data/constant2D.dat")
p = Util.readImageL(n1,n2,n3,PATH+"/pgs/data/constant2D_out.dat")
s1 = Sampling(n1)
s2 = Sampling(n2)
s3 = Sampling(n3)

def main(args):
  plot3D(f,p)

def plot3D(f,p):
  one = True;
  Plot.plot(s1,s2,s3,f,p,"slope (samples/trace)","p_panels",
      0.6,0.8,-cmax,cmax,False,True,one);

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
