import sys

"""
Compares slope methods for given real and synthetic seismic images.
"""
from java.lang import *
from javax.swing import *

from edu.mines.jtk.dsp import *
from edu.mines.jtk.util.ArrayMath import *
from slopes import *
from util import *

noise = 0.0
pmax = 4.0
PATH = "/users/elias.arias/Home/git/ea/bench/src"

n1 = 101
n2 = 102
n3 = 103
f  = zerofloat(n1,n2,n3)
p2 = zerofloat(n1,n2,n3)
p3 = zerofloat(n1,n2,n3)
s1 = Sampling(n1)
s2 = Sampling(n2)
s3 = Sampling(n3)
s3D = Slopes(noise,pmax,s1,s2,s3)
lsf_title = "lsf_complex3D"
pwd_title = "pwd_complex3D"
sdw_title = "sdw_complex3D"
pgs_title = "pgs_complex3D"
known_title = "known_complex3D"
pgs_p2 = PATH+"/pgs/pgs_slope_out_dzdy.bin"
pgs_p3 = PATH+"/pgs/pgs_slope_out_dzdx.bin"
Slopes.makeSyntheticComplex(noise,f,p2,p3)

def main(args):
  #littleToBig(PATH+"/pgs/synth3D_testbin",PATH+"/pgs/synth3D_testbinB")
  #testLsf()
  #testPwd()
  #testSdw()
  testPgs()
  #testAll()

def testLsf():
  s3D.estimateLSF(f,p2,p3,lsf_title)
  s3D.plot3D(f,lsf_title)

def testPwd():
  s3D.estimatePWDM(f,p2,p3,pwd_title)
  s3D.plot3D(f,pwd_title)

def testSdw():
  s3D.estimateSDW(k,f,p2,p3,sdw_title)
  s3D.plot3D(f,sdw_title)

def testPgs():
  littleToBig(pgs_p2,PATH+"/util/data/pgs_complex3D_p2.dat")
  littleToBig(pgs_p3,PATH+"/util/data/pgs_complex3D_p3.dat")
  s3D.plot3D(f,pgs_title)

def testAll():
  s3D.estimateLSF(f,p2,p3,lsf_title)
  s3D.plot3D(f,lsf_title)
  s3D.estimatePWDM(f,p2,p3,pwd_title)
  s3D.plot3D(f,pwd_title)
  s3D.estimateSDW(k,f,p2,p3,sdw_title)
  s3D.plot3D(f,sdw_title)
  littleToBig(pgs_p2,PATH+"/util/data/pgs_complex3D_p2.dat")
  littleToBig(pgs_p3,PATH+"/util/data/pgs_complex3D_p3.dat")
  s3D.plot3D(f,pgs_title)

def littleToBig(str_in,str_out):
  Util.writeBinary(Util.readImageL(n1,n2,n3,str_in),str_out)

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
