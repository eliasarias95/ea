package pwd;

public class Trianglen {
  public Trianglen(int ndim, int[] nbox, int[] ndat) {
    _dim = ndim;
    _n = new int[_dim];
    _tr = new Triangle[_dim];
    _nd = 1;
    for (int i=0; i<_dim; ++i) {
      _tr[i] = (nbox[i]>1)? new Triangle(nbox[i],ndat[i]): null;
      _s[i] = _nd;
      _n[i] = ndat[i];
      _nd *= ndat[i];
    }
    _tmp = new float[_nd];
  }

  private int _nd,_dim;
  private int[] _n,_s;
  private Triangle[] _tr;
  private float[] _tmp;
}
