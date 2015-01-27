package slopes;

import static edu.mines.jtk.util.ArrayMath.*;
import dnp.LocalSlopeFinder;
import util.Plot;

public class StructureTensor extends SetParameters 
                             implements SlopeEstimationMethod {
  public StructureTensor(float sigma1, float sigma2) {
    LocalSlopeFinder _lsf = new LocalSlopeFinder(sigma1,sigma2,_pmax);
  }

  public void estimateSlope() {
    _pe = new float[_n2][_n1];
    _lsf.findSlopes(_f,_pe);
    _pe = mul(_pe,_dg/_d2);
    System.out.println("We did it?");
    System.out.println("slopes: "+_pe[230][225]);    
    // clip, interp, title, paint, colorbar, color
    String hl = "Traces"; //horizontal label
    String vl = "Samples"; //vertical label
    String cbl = "slope (samples/trace)"; //colorbar label
    Plot.plot(_sg,_s2,_pe,"LSF noise= "+_noise,hl,vl,cbl,
        _fw,_fh,-_clipMax,_clipMax,
        _clip,F,_title,_paint,T,T);
  }

  /*
  public static void calculateError() {
    System.out.println("Structure tensor:");
    float error_lsf = Util.rmsError(p_lsf,_pk,_dg,_d2,T);
  }

  public static void plotSlope() {

  }
  
  public static void testOptimalParameters() {

  }

  public static void plotOptimalParameters() {

  }
  */

  /*
  private static Sampling _s1 = new Sampling(1,1,1);  
  private static Sampling _s2 = new Sampling(1,1,1);
  private static Sampling _sg = new Sampling(1,1,1);
  private static Sampling _s1param = new Sampling(1,1,1);
  private static Sampling _s2param = new Sampling(1,1,1);

  private static int _n1,_n2,_ng,_n1param,_n2param;
  private static float _d1,_d2,_dg,_f1,_f2,_fg,_dparam,_fparam,_noise;
  private static float[] _param1 = new float[_n1param];
  private static float[] _param2 = new float[_n2param];
  private static float[][] _f = new float[_n2][_n1];
  private static float[][] _pk = new float[_n2][_n1];
  */
  private static float[][] _pe;
  private LocalSlopeFinder _lsf = new LocalSlopeFinder(1,1,1);
}
