public class Laplacian{
  public static float[][] aniLap(float[][] f  , float[][] d11, 
                                 float[][] d12, float[][] d22) {
    int n2 = f.length;
    int n1 = f[0].length;
    float[][] g = new float[n2][n1];
    float fa, fb, f1, f2, ga, gb, g1, g2;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        fa = 0.0f;
        fb = 0.0f;
        fa += f[i2 ][i1 ];
        fb -= f[i2 ][i1-1]; //gather
        fb += f[i2-1][i1 ];
        fa -= f[i2-1][i1-1];
        f1 = 0.5f*(fa+fb);
        f2 = 0.5f*(fa-fb);
        /***************************/
        g1 = d11[i2][i1]*f1+d12[i2][i1]*f2; //scale
        g2 = d12[i2][i1]*f1+d22[i2][i1]*f2;
        /***************************/
        ga = 0.5f*(g1+g2);
        gb = 0.5f*(g1-g2);
        g[i2 ][i1 ] += ga; //scatter
        g[i2 ][i1-1] -= gb;
        g[i2-1][i1 ] += gb;
        g[i2-1][i1-1] -= ga;
      }
    }
    return g;
  }

  public static float[][] isoLap(float[][] f  , float[][] d11, 
                                 float[][] d12, float[][] d22) {
    int n2 = f.length;
    int n1 = f[0].length;
    float[][] g = new float[n2][n1];
    float fa, fb, f1, f2, ga, gb, g1, g2;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        fa = 0.0f;
        fb = 0.0f;
        fa += f[i2 ][i1 ];
        fb -= f[i2 ][i1-1]; //gather
        fb += f[i2-1][i1 ];
        fa -= f[i2-1][i1-1];
        f1 = 0.5f*(fa+fb);
        f2 = 0.5f*(fa-fb);
        /***************************/
        g1 = f1; // copy
        g2 = f2;
        /***************************/
        ga = 0.5f*(g1+g2);
        gb = 0.5f*(g1-g2);
        g[i2 ][i1 ] += ga; //scatter
        g[i2 ][i1-1] -= gb;
        g[i2-1][i1 ] += gb;
        g[i2-1][i1-1] -= ga;
      }
    }
    return g;
  }
}

