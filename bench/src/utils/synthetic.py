from imports import *
from utils import *

n1 = 1001
fpeak = 0.1

def main(args):
  test()

def test():
  shift_val = 2
  n_traces  = 101
  fpeak     = 30
  seed      = 200

  s = Sampling(1001,0.002,0.0)

  synth          = zerofloat(1001,n_traces)
  lsf_output     = zerofloat(1001,n_traces)
  lsf_output_lin = zerofloat(1001,n_traces)

  trc   = Synthetic.getTrace(s,fpeak,seed)
  synth = Synthetic.make2D(trc,n_traces,shift_val)
  Synthetic.writeBinary(synth,"synth.dat");

  lsf = LocalSlopeFinder(10.0);
  lsf.findSlopes(synth,lsf_output,lsf_output_lin);
  Synthetic.plot(synth,"Sythetic Seismic (Slope=2 samples/trace)")
  Synthetic.plot(lsf_output,"Plot of the slope values")
  Synthetic.plot(lsf_output_lin,"Plot of the linearities")

######################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
