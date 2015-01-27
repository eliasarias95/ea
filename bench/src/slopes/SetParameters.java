package slopes;

import edu.mines.jtk.dsp.Sampling;
import util.Util;
import util.FakeData;

import static edu.mines.jtk.util.ArrayMath.*;

public class SetParameters {

  public static void setChickenTestParameters(float pk) {
    float number = 01.0f;
    _ng = (int)(501*number);
    _n1 = 501;
    _n2 = 501;
    _dg = 1.0f/number;
    _d1 = 1.0f;
    _d2 = 1.0f;
    _fg = 0.0f;
    _f1 = 0.0f;
    _f2 = 0.0f;
    _s1 = new Sampling(_n1,_d1,_f1);
    _s2 = new Sampling(_n2,_d2,_f2);
    _sg = new Sampling(_ng,_dg,_fg);
    float w = 0.1f; //frequency (cycles/sample)
    _f = new float[_n2][_n1];
    _pk = new float[_n2][_n1];
    for (int i2=0; i2<_n2; ++i2){
      for (int i1=0; i1<_n1; ++i1){
        _f[i2][i1] = cos(pk*w*i2-w*i1);
        _pk[i2][i1] = pk;
      }
    }
    for (int i=0; i<_n2; ++i) {
      _f[i] = Util.reSampleSinc(_s1,_f[i],_sg);
      _pk[i] = Util.reSampleSinc(_s1,_pk[i],_sg);
    }
  }

  public static void setSynthParameters(float noise) {
    float number = 1.0f;
    _ng = (int)(501*number);
    _n1 = 501;
    _n2 = 501;
    _dg = 1.0f/number;
    _d1 = 1.0f;
    _d2 = 1.0f;
    _fg = 0.0f;
    _f1 = 0.0f;
    _f2 = 0.0f;
    _s1 = new Sampling(_n1,_d1,_f1);
    _s2 = new Sampling(_n2,_d2,_f2);
    _sg = new Sampling(_ng,_dg,_fg);
    _noise = noise;
    float[][][] fandpk = FakeData.seismicAndSlopes2d2014B(_noise,F);
    _f = fandpk[0];  //synthetic seismic data
    _pk = fandpk[1]; //exact slope values
    for (int i=0; i<_n2; ++i) {
      _f[i] = Util.reSampleSinc(_s1,_f[i],_sg);
      _pk[i] = Util.reSampleSinc(_s1,_pk[i],_sg);
    }
  }

  public static void setGOMParameters() {
    float number = 1.0f;
    _ng = (int)(501*number);
    _n1 = 301;
    _n2 = 920;
    _dg = 1.0f/number;
    _d1 = 0.004f;
    _d2 = .02667f;
    _fg = 0.0f;
    _f1 = 1.6f;
    _f2 = 0;
    _s1 = new Sampling(_n1,_d1,_f1);
    _s2 = new Sampling(_n2,_d2,_f2);
    _sg = new Sampling(_ng,_dg,_fg);
    _f = Util.readImage(_n1,_n2,"data/gom.dat");
    float[][] fs = new float[_n2][_n1];
    mul(_f,.001f,fs);
    Util.writeBinary(fs,"data/gom_scaled.dat");
  }

  public static void setSmoothingParameters(float l_param) {
    _dparam= 1.0f;
    _fparam = 1.0f;
    _n1param = (int)((l_param-_fparam+1)/_dparam);
    _n2param = 10;
    _s1param = new Sampling(_n1param,_dparam,_fparam);
    _s2param = new Sampling(_n2param,_dparam,_fparam);
    _param1 = new float[_n1param];
    _param2 = new float[_n2param];
    for(int i=0; i<_n1param; ++i) {
      _param1[i] = i*_dparam+_fparam;
    }
    for(int i=0; i<_n2param; ++i) {
      _param2[i] = i*_dparam+_fparam;
    }
  }

  public static void setStrainParameters() {
    _dparam = 0.1f; // strain sampling rate
    _fparam = 0.1f;
    _n1param = (int)((1.0f-_fparam+0.1f)/_dparam); //# strain vals to test
    _n2param = _n1param;
    _s1param = new Sampling(_n1param,_dparam,_fparam);
    _s2param = new Sampling(_n2param,_dparam,_fparam);
    _param1 = new float[_n1param];
    _param2 = new float[_n2param];
    for(int i=0; i<_n1param; ++i) {
      _param1[i] = i*_dparam+_fparam;
      _param2[i] = i*_dparam+_fparam;
    }
  }

  protected static final boolean T = true;
  protected static final boolean F = false;

  protected static Sampling _s1 = new Sampling(1,1,1);  
  protected static Sampling _s2 = new Sampling(1,1,1);
  protected static Sampling _sg = new Sampling(1,1,1);
  protected static Sampling _s1param = new Sampling(1,1,1);
  protected static Sampling _s2param = new Sampling(1,1,1);

  protected static int _n1,_n2,_ng,_n1param,_n2param;
  protected static float _d1,_d2,_dg,_f1,_f2,_fg,_dparam,_fparam,_noise;
  protected static float[] _param1 = new float[_n1param];
  protected static float[] _param2 = new float[_n2param];
  protected static float[][] _f = new float[_n2][_n1];
  protected static float[][] _pk = new float[_n2][_n1];
}
