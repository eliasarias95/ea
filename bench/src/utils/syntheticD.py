from imports import *
from utils import *
from math import pi
from pwd import *

T = True
F = False 

def main(args):
  test()
  PWD.goSfdipSynth()

def test():  
  nx        = 501
  nt        = 501
  dx        = 1
  dt        = 1
  fx        = 0
  ft        = 0
  noise     = 0.5

  sx = Sampling(nx,dx,fx)
  st = Sampling(nt,dt,ft)

  synthAndSlope  = zerofloat(nt,nx,2)
  synthAndSlope  = FakeData.seismicAndSlopes2d2014B(noise,False)
  synth          = zerofloat(nt,nx)
  slope          = zerofloat(nt,nx)
  synth          = synthAndSlope[0]
  slope          = synthAndSlope[1]

  Synthetic.writeBinary(synth,
          "/Users/earias/Home/git/ea/bench/src/slopes/data/synth.dat");
  Synthetic.writeBinary(slope,
          "/Users/earias/Home/git/ea/bench/src/slopes/data/exact_slope.dat");
  Synthetic.plot(st,sx,synth,"Synthetic Seismic",F,F,F,F)
  Synthetic.plot(st,sx,slope,"Exact Slope Values",F,F,T,T)

######################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())

