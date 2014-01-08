package pwd;

public class Triangle {
  public Triangle(int nbox, int ndat) {
    _nx = ndat;
    _nb = nbox;
    _np = ndat + 2*nbox;
    _tmp = new float[_np];
  }

  private float[] _tmp;
  private int _np,_nb,_nx;
}
