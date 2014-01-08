package pwd;

public class Divn {
  public Divn(int ndim, int nd, int[] ndat, int[] nbox, int niter, 
                                                      boolean verb)  {
    _niter = niter;
    _n = nd;
    Trianglen tr = new Trianglen(ndim,nbox,ndat);

    _p = new float[nd];
  }

  private int _niter,_n;
  private float[] _p;
}
