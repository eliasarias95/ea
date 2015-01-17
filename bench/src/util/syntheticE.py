from imports import *
from utils import *
from math import pi

T = True
F = False 

def main(args):
  test()
  sweep(1,5,0.5,1000)

def sweep(f_start, f_end, interval, n_steps):
  f = zerofloat(n_steps)
  for i in range(n_steps):
    delta = i / float(n_steps)
    t = interval * delta
    f[i] = t
  #Synthetic.plot(f,sin(Synthetic.sweep(f_start,f_end,interval,n_steps)),
  #        "Sweep")

def test():  
  #shift_val = 0.01
  A         = 0.25
  nx        = 1000
  nt        = 501
  dx        = 0.0016
  dt        = 0.002
  fx        = 0.0
  ft        = -4.0
  fpeak     = 15
  seed      = 200
  shift     = 1 

  strc = Sampling(6*nt,dt,ft)
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,fx)
  #s = Sampling(n1,d1,0.0)

  synth          = zerofloat(nt,nx)
  slope          = zerofloat(nt,nx)
  const_slope    = zerofloat(nt,nx)
  lsf_output     = zerofloat(nt,nx)
  lsf_output_lin = zerofloat(nt,nx)

  trc   = Synthetic.getTrace(strc,fpeak,seed)
  synth = Synthetic.sineData(st,sx,ft,trc,A,slope)
  const_synth = Synthetic.make2D(st,sx,ft,trc,shift,const_slope)
  #const_synth = Synthetic.make2D(trc,nx,1)
  unconformity = False 
  if (unconformity):
      m = 0  
      n1 = 1.0/(4*dt)
      for ix in range(nx):
          n1 += m
          for i1 in range(n1):
              synth[ix][i1] = const_synth[ix][i1]
              slope[ix][i1] = const_slope[ix][i1] 
              m = 0.10*dt/dx

  slope = div(slope,dx)
  print slope[4][4]
  print slope[240][313]
  Synthetic.writeBinary(synth,
          "/Users/earias/Home/git/ea/bench/src/slopes/data/synth.dat");
  Synthetic.writeBinary(slope,
          "/Users/earias/Home/git/ea/bench/src/slopes/data/exact_slope.dat");
  lsf = LocalSlopeFinder(8.0,4.0,3)
  lsf.findSlopes(synth,lsf_output,lsf_output_lin)
  Synthetic.plot(st,sx,synth,"Synthetic Seismic",F,F,F,F)
  #Synthetic.plot(st,sx,const_synth,"Constant Slope",F,F,F,F)
  Synthetic.plot(st,sx,slope,"Exact Slope Values",F,F,T,T)
  #Synthetic.plot(st,sx,synth,"Sythetic Seismic (Slope=2 samples/trace)",
  #  F,F,F,F)
  #Synthetic.plot(st,sx,lsf_output,"Plot of the slope values",F,F,T,T)

######################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
