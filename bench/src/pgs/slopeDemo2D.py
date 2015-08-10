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
paint = T
fw = 0.8
fh = 0.8
PATH = "/users/elias.arias/Home/git/ea/bench/src"

def main(args):
  #constant2D()
  #complex2D()
  santos2D()
  #triton2D()
  #complex3D()

def constant2D():
  n1 = 501
  n2 = 501
  cm = 1.0 #colorbar limits
  f = Util.readImageL(n1,n2,PATH+"/pgs/data/constant2D.dat")
  p_sdw = Util.readImageL(n1,n2,PATH+"/pgs/data/constant2D_sdw.dat")
  p_sxy = Util.readImageL(n1,n2,PATH+"/pgs/data/constant2D_sxy_dzdx.dat")

  plot2D(n1,n2,f,p_sdw,cm,"constant2D_sdwc") # sdw_c
  plot2D(n1,n2,f,p_sxy,cm,"constant2D_slopexy") # slopexy
  print "mean slopexy=", sum(p_sxy)/(n1*n2)

def complex2D():
  n1 = 501
  n2 = 501
  cm = 4.0 #colorbar limits
  f = Util.readImageL(n1,n2,PATH+"/pgs/data/complex2D.dat")
  p = Util.readImageL(n1,n2,PATH+"/pgs/data/complex2D_p.dat")
  p_sdwc = Util.readImageL(n1,n2,PATH+"/pgs/data/complex2D_sdw.dat")
  p_sdwj = Util.readImage(n1,n2,PATH+"/pgs/data/sdw_zero_complex2D.dat")
  p_sxy  = Util.readImageL(n1,n2,PATH+"/pgs/data/complex2D_sxy_dzdx.dat")

  print "SDW C++ RMS error:"
  e_sdwc = Util.rmsError(p_sdwc,p,T)
  print "SDW Java RMS error:"
  e_sdwj = Util.rmsError(p_sdwj,p,T)
  print "Slopexy RMS error:"
  e_sxy = Util.rmsError(p_sxy,p,T)

  plot2D(n1,n2,f,p_sdwc,cm,"complex2D_sdwc") # sdwc
  #plot2D(n1,n2,f,p_sdwj,cm,"complex2D_sdwj") # sdwj
  #plot2D(n1,n2,f,p_sxy,cm,"complex2D_slopexy") # slopexy

def santos2D():
  n1 = 2001
  n2 = 809
  cm = 4.0 #colorbar limits
  f = Util.readImageL(n1,n2,PATH+"/pgs/data/Santos_2D.dat")
  #p_sdwc = Util.readImageL(n1,n2,PATH+"/pgs/data/Santos_2D_sdw.dat")
  #p_sdwj = Util.readImage(n1,n2,PATH+"/pgs/data/Santos_2D_java.dat")
  p_sxy  = Util.readImageL(n1,n2,PATH+"/pgs/data/Santos_2D_sxy_dzdx.dat")

  #plot2D(n1,n2,f,p_sdwc,cm,"santos2D_sdwc") # sdwc
  #plot2D(n1,n2,f,p_sdwj,cm,"santos2D_sdwj") # sdwj
  plot2D(n1,n2,f,p_sxy,cm,"santos2D_slopexy") # slopexy

def triton2D():
  n1 = 2001
  n2 = 2921
  cm = 4.0 #colorbar limits
  f = Util.readImageL(n1,n2,PATH+"/pgs/data/Triton_2D.dat")
  p_sdwc = Util.readImageL(n1,n2,PATH+"/pgs/data/Triton_2D_sdw.dat")
  #p_sdwj = Util.readImage(n1,n2,PATH+"/pgs/data/Triton_2D_java.dat")
  p_sxy  = Util.readImageL(n1,n2,PATH+"/pgs/data/Triton_2D_sxy_dzdx.dat")

  plot2D(n1,n2,f,p_sdwc,cm,"triton2D_sdwc") # sdwc
  #plot2D(n1,n2,f,p_sdwj,cm,"triton2D_sdwj") # sdwj
  plot2D(n1,n2,f,p_sxy,cm,"triton2D_slopexy") # slopexy

def complex3D():
  n1 = 101
  n2 = 102
  n3 = 103
  cm = 1.5 #colorbar limits
  f = Util.readImageL(n1,n2,n3,PATH+"/pgs/data/complex3D.dat")
  p2 = Util.readImageL(n1,n2,n3,PATH+"/pgs/data/complex3D_p2.dat")
  p3 = Util.readImageL(n1,n2,n3,PATH+"/pgs/data/complex3D_p3.dat")
  p2_sdwc = Util.readImageL(n1,n2,n3,PATH+"/pgs/data/complex3D_sdw_dzdx.dat")
  #p2_sdwj = Util.readImage(n1,n2,n3,PATH+"/pgs/data/sdw_half_complex2D.dat")
  p2_sxy  = Util.readImageL(n1,n2,n3,PATH+"/pgs/data/complex3D_sxy_dzdx.dat")

  p3_sdwc = Util.readImageL(n1,n2,n3,PATH+"/pgs/data/complex3D_sdw_dzdy.dat")
  #p3_sdwj = Util.readImage(n1,n2,n3,PATH+"/pgs/data/sdw_half_complex2D.dat")
  p3_sxy  = Util.readImageL(n1,n2,n3,PATH+"/pgs/data/complex3D_sxy_dzdy.dat")

  print "Xline slope:"
  print "SDW C++ RMS error:"
  e2_sdwc = Util.rmsError(p2_sdwc,p2,T)
  print "SDW Java RMS error:"
  #e2_sdwj = Util.rmsError(p2_sdwj,p2,T)
  print "Slopexy RMS error:"
  e2_sxy = Util.rmsError(p2_sxy,p2,T)

  print "Subline slope:"
  print "SDW C++ RMS error:"
  e3_sdwc = Util.rmsError(p3_sdwc,p3,T)
  print "SDW Java RMS error:"
  #e3_sdwj = Util.rmsError(p3_sdwj,p3,T)
  print "Slopexy RMS error:"
  e3_sxy = Util.rmsError(p3_sxy,p3,T)

  plot3D(n1,n2,n3,f,p2_sdwc,p3_sdwc,cm,"complex3D_sdwc") # sdwc
  #plot3D(n1,n2,n3,f,p2_sdwj,p3_sdwj,cm,"complex3D_sdwj") # sdwj
  plot3D(n1,n2,n3,f,p2_sxy,p3_sxy,cm,"complex3D_sxy") # slopexy

def plot2D(n1,n2,f,p,cm,ttl):
  s1 = Sampling(n1)
  s2 = Sampling(n2)
  cbl = "slope (samples/trace)" #colorbar label
  #clip, title, paint, slide, no. columns
  Plot.plot(s1,s2,f,p,ttl,cbl,fw,fh,-cm,cm,T,F,paint,T,T)

def plot3D(n1,n2,n3,f,p2,p3,cm,ttl):
  s1 = Sampling(n1)
  s2 = Sampling(n2)
  s3 = Sampling(n3)
  Plot.plot(s1,s2,s3,f,p2,"slope (samples/trace)",ttl+"_p2",fw,fh,-cm,cm,paint,
      T,T);
  Plot.plot(s1,s2,s3,f,p3,"slope (samples/trace)",ttl+"_p3",fw,fh,-cm,cm,paint,
      T,T);
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
