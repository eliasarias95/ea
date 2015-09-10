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
PATH = "/datadb/eariasp/data/"
nodes = 2
h3 = 10

def main(args):
  #subsample()
  #subcomplex()
  #constant2D()
  #complex2D()
  #santos2D()
  #triton2D()
  #complex3D()
  santos3D()
  #triton3D()

def subsample():
  n1 = 101
  n2 = 102
  n  = 103
  n3 = n-(n/nodes*(nodes-1)+n%nodes)+h3
  print n3
  #f = Util.readImageL(n1,n2,n,PATH+"complex3D.dat")
  #f = Util.readImageL(n1,n2,n,PATH+"first100_complex3D.dat")
  #g = zerofloat(n1,n2,n3)
  p2 = Util.readImageL(n1,n2,n,PATH+"complex3D_sdw_dzdx.dat")
  q2 = zerofloat(n1,n2,n3)
  p3 = Util.readImageL(n1,n2,n,PATH+"complex3D_sdw_dzdy.dat")
  q3 = zerofloat(n1,n2,n3)

  #r2 = Util.readImageL(n1,n2,n3+h3,PATH+"test61_complex3D_sdw_dzdx.dat")
  #s2 = zerofloat(n1,n2,n3)
  #r3 = Util.readImageL(n1,n2,n3+h3,PATH+"test61_complex3D_sdw_dzdy.dat")
  #s3 = zerofloat(n1,n2,n3)
  #p2_sdwc = Util.readImageL(n1,n2,n,PATH+"complex3D_sdw_dzdx.dat")
  #p3_sdwc = Util.readImageL(n1,n2,n,PATH+"complex3D_sdw_dzdy.dat")
  for i3 in range(n3):
    for i2 in range(n2):
      for i1 in range(n1):
        #g[i3][i2][i1]  = f [i3+(n-n3)][i2][i1]
        q2[i3][i2][i1] = p2[i3+(n-n3)][i2][i1]
        q3[i3][i2][i1] = p3[i3+(n-n3)][i2][i1]
        #s2[i3][i2][i1] = r2[i3+h3][i2][i1]
        #s3[i3][i2][i1] = r3[i3+h3][i2][i1]
        #s2[i3][i2][i1] = r2[i3+(n-n3)][i2][i1]
        #s3[i3][i2][i1] = r3[i3+(n-n3)][i2][i1]
  #Util.writeBinaryL(g,"test61_complex3D.dat");
  Util.writeBinaryL(q2,"test61_complex3D_sdw_dzdx_multi.dat");
  Util.writeBinaryL(q3,"test61_complex3D_sdw_dzdy_multi.dat");
  #Util.writeBinaryL(s2,"test51_complex3D_sdw_dzdx_single.dat");
  #Util.writeBinaryL(s3,"test51_complex3D_sdw_dzdy_single.dat");
  #Util.writeBinaryL(s2,"test_complex3D_sdw_dzdx.dat");
  #Util.writeBinaryL(s3,"test_complex3D_sdw_dzdy.dat");

def constant2D():
  n1 = 501
  n2 = 501
  cm = 1.0 #colorbar limits
  f = Util.readImageL(n1,n2,PATH+"constant2D.dat")
  p_sdw = Util.readImageL(n1,n2,PATH+"constant2D_sdw.dat")
  p_sxy = Util.readImageL(n1,n2,PATH+"constant2D_sxy_dzdx.dat")

  plot2D(n1,n2,f,p_sdw,cm,"constant2D_sdwc") # sdw_c
  plot2D(n1,n2,f,p_sxy,cm,"constant2D_slopexy") # slopexy
  print "mean slopexy=", sum(p_sxy)/(n1*n2)

def complex2D():
  n1 = 501
  n2 = 501
  cm = 4.0 #colorbar limits
  f = Util.readImageL(n1,n2,PATH+"complex2D.dat")
  p = Util.readImageL(n1,n2,PATH+"complex2D_p.dat")
  p_sdwc = Util.readImageL(n1,n2,PATH+"complex2D_sdw.dat")
  #p_sdwj = Util.readImage(n1,n2,PATH+"sdw_zero_complex2D.dat")
  #p_sdwj = Util.readImage(n1,n2,PATH+"sdw_complex2D_p.dat")
  #p_sxy  = Util.readImageL(n1,n2,PATH+"complex2D_sxy_dzdx.dat")

  print "SDW C++ & SDW Java RMS error:"
  #e_sdw = Util.rmsError(p_sdwc,p_sdwj,T)
  print "SDW C++ RMS error:"
  e_sdwc = Util.rmsError(p_sdwc,p,T)
  #print "SDW Java RMS error:"
  #e_sdwj = Util.rmsError(p_sdwj,p,T)
  #print "Slopexy RMS error:"
  #e_sxy = Util.rmsError(p_sxy,p,T)

  plot2D(n1,n2,f,p_sdwc,cm,"complex2D_sdwc") # sdwc
  #plot2D(n1,n2,f,p_sdwj,cm,"complex2D_sdwj") # sdwj
  #plot2D(n1,n2,f,p_sxy,cm,"complex2D_slopexy") # slopexy

def santos2D():
  n1 = 2001
  n2 = 3809
  cm = 4.0 #colorbar limits
  f = Util.readImageL(n1,n2,PATH+"Santos_2D.dat")
  p_sdw = Util.readImageL(n1,n2,PATH+"Santos_2D_sdw.dat")
  p_sxy = Util.readImageL(n1,n2,PATH+"Santos_2D_sxy.dat")

  plot2D(n1,n2,f,p_sdw,cm,"santos2D_sdw") # sdw
  #plot2D(n1,n2,f,p_sxy,cm,"santos2D_slopexy") # slopexy

def santos3D():
  #n1 = 1001
  #n2 = 1905
  #n3 = 1913
  n1 = 501
  n2 = 201
  n3 = 201
  cm = 2.0 #colorbar limits
  #f = Util.readImageL(n1,n2,n3,PATH+"Santos_sub_3D.dat")
  #p2_sdwc = Util.readImageL(n1,n2,n3,PATH+"Santos_sub_3D_sdw_dzdx.dat")
  #p3_sdwc = Util.readImageL(n1,n2,n3,PATH+"Santos_sub_3D_sdw_dzdy.dat")
  f = Util.readImageL(n1,n2,n3,PATH+"Santos_small_3D.dat")
  p2_sdwc = Util.readImageL(n1,n2,n3,PATH+"Santos_small_3D_sdw_dzdx.dat")
  p3_sdwc = Util.readImageL(n1,n2,n3,PATH+"Santos_small_3D_sdw_dzdy.dat")
  #p2_sxy  = Util.readImageL(n1,n2,n3,PATH+"Santos_small_3D_sxy_dzdx.dat")
  #p3_sxy  = Util.readImageL(n1,n2,n3,PATH+"Santos_small_3D_sxy_dzdy.dat")
  #p2_sdwc2 = Util.readImageL(n1,n2,n3,PATH+"Santos_small_3D_sdw_dzdx2.dat")
  #p3_sdwc2 = Util.readImageL(n1,n2,n3,PATH+"Santos_small_3D_sdw_dzdy2.dat")

  plot3Ds(n1,n2,n3,f,p2_sdwc,p3_sdwc,cm,"Santos_sub_3D_sdw") # sdw
  #plot3Ds(n1,n2,n3,f,p2_sxy,p3_sxy,cm,"Santos_sub_3D_sxy") # sxy
  #plot3Dp(n1,n2,n3,f,p2_sdwc,p3_sdwc,cm,"Santos_sub_3D_sdw") # sdw
  #plot3Dp(n1,n2,n3,f,p2_sxy,p3_sxy,cm,"Santos_sub_3D_sxy") # sxy

def triton3D():
  n1 = 1201
  n2 = 1561
  n3 = 251
  cm = 2.0 #colorbar limits
  f = Util.readImageL(n1,n2,n3,PATH+"Triton_small_3D.dat")
  p2_sdw = Util.readImageL(n1,n2,n3,PATH+"Triton_small_3D_sdw_dzdx.dat")
  p3_sdw = Util.readImageL(n1,n2,n3,PATH+"Triton_small_3D_sdw_dzdy.dat")
  #p2_sxy = Util.readImageL(n1,n2,n3,PATH+"Triton_small_3D_sxy_dzdx.dat")
  #p3_sxy = Util.readImageL(n1,n2,n3,PATH+"Triton_small_3D_sxy_dzdy.dat")

  #plot3Ds(n1,n2,n3,f,p2_sdw,p3_sdw,cm,"Santos_sub_3D_sdw") # sdw
  #plot3Ds(n1,n2,n3,f,p2_sxy,p3_sxy,cm,"Santos_sub_3D_sxy") # sxy
  plot3Dp(n1,n2,n3,f,p2_sdw,p3_sdw,cm,"Santos_sub_3D_sdw") # sdw
  #plot3Dp(n1,n2,n3,f,p2_sxy,p3_sxy,cm,"Santos_sub_3D_sxy") # sxy

def triton2D():
  n1 = 2001
  n2 = 2921
  cm = 4.0 #colorbar limits
  f = Util.readImageL(n1,n2,PATH+"Triton_2D.dat")
  p_sdw = Util.readImageL(n1,n2,PATH+"Triton_2D_sdw.dat")
  p_sxy = Util.readImageL(n1,n2,PATH+"Triton_2D_sxy.dat")

  plot2D(n1,n2,f,p_sdw,cm,"triton2D_sdw") # sdw
  #plot2D(n1,n2,f,p_sxy,cm,"triton2D_slopexy") # slopexy

def complex3D():
  n1 = 101
  n2 = 102
  n3 = 103
  cm = 1.5 #colorbar limits
  f = Util.readImageL(n1,n2,n3,PATH+"complex3D.dat")
  p2 = Util.readImageL(n1,n2,n3,PATH+"complex3D_p2.dat")
  p3 = Util.readImageL(n1,n2,n3,PATH+"complex3D_p3.dat")
  p2_sdwc = Util.readImageL(n1,n2,n3,PATH+"complex3D_sdw_dzdx.dat")
  p3_sdwc = Util.readImageL(n1,n2,n3,PATH+"complex3D_sdw_dzdy.dat")
  p2_sdwc2 = Util.readImageL(n1,n2,n3,PATH+"complex3D_sdw_dzdx2.dat")
  p3_sdwc2 = Util.readImageL(n1,n2,n3,PATH+"complex3D_sdw_dzdy2.dat")

  print "Xline slope:"
  print "SDW C++ RMS error:"
  e2_sdwc = Util.rmsError(p2_sdwc,p2,T)
  e2_sdwc2 = Util.rmsError(p2_sdwc2,p2,T)

  print "Subline slope:"
  print "SDW C++ RMS error:"
  e3_sdwc = Util.rmsError(p3_sdwc,p3,T)
  e3_sdwc2 = Util.rmsError(p3_sdwc2,p3,T)

  plot3Ds(n1,n2,n3,f,p2_sdwc,p3_sdwc,cm,"multi_node") # sdwc
  plot3Ds(n1,n2,n3,f,p2_sdwc2,p3_sdwc2,cm,"single_node") # sdwc

def subcomplex():
  n1 = 101
  n2 = 102
  n  = 103
  n3 = n-(n/nodes*(nodes-1)+n%nodes)+h3
  print n3
  f  = Util.readImageL(n1,n2,n3,PATH+"test61_complex3D.dat")
  p2 = Util.readImageL(n1,n2,n3,PATH+"test61_complex3D_sdw_dzdx_multi.dat")
  p3 = Util.readImageL(n1,n2,n3,PATH+"test61_complex3D_sdw_dzdy_multi.dat")

  q2 = Util.readImageL(n1,n2,n3,PATH+"test61_complex3D_sdw_dzdx.dat")
  q3 = Util.readImageL(n1,n2,n3,PATH+"test61_complex3D_sdw_dzdy.dat")

  print "Xline slope:"
  print "SDW C++ RMS error:"
  e2 = Util.rmsError(q2,p2,T)

  print "Subline slope:"
  print "SDW C++ RMS error:"
  e3 = Util.rmsError(q3,p3,T)

  cm = 0.5 #colorbar limits
  plot3Ds(n1,n2,n3,f,p2,p3,cm,"full") # multi node
  plot3Ds(n1,n2,n3,f,q2,q3,cm,"chunk") # single node

def plot2D(n1,n2,f,p,cm,ttl):
  s1 = Sampling(n1)
  s2 = Sampling(n2)
  cbl = "slope (samples/trace)" #colorbar label
  #clip, title, paint, slide, no. columns
  Plot.plot(s1,s2,f,p,ttl,cbl,fw,fh,-cm,cm,T,F,paint,T,T)

def plot3Dp(n1,n2,n3,f,p2,p3,cm,ttl):
  s1 = Sampling(n1)
  s2 = Sampling(n2)
  s3 = Sampling(n3)
  Plot.plot(s1,s2,s3,f,p2,"slope (samples/trace)",ttl+"_p2",fw,fh,-cm,cm,paint,
      T,T);
  Plot.plot(s1,s2,s3,f,p3,"slope (samples/trace)",ttl+"_p3",fw,fh,-cm,cm,paint,
      T,T);

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
