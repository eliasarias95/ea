package slopes;

public interface SlopeEstimationMethod {
  ///////////////////PRIVATE VARIABLES///////////////////////
  public static final String _path = 
    "/Users/earias/Home/git/ea/bench/src/slopes/data/";
  public static final int _niter = 5;
  public static final float _fw = 0.75f; //fraction width for slide
  public static final float _fh = 0.9f; //fraction height for slide
  public static final float _pmax = 5.0f;
  public static final float _clipMax = 4.0f;
  //public static final boolean T = true;
  //public static final boolean F = false;  
  public static final boolean _title = true;
  public static final boolean _paint = false;  
  public static final boolean _clip = false;
  ///////////////////////METHODS/////////////////////////////
  public void estimateSlope();
  //public static void calculateError();
  //public static void plotSlope();
  //public static void testOptimalParameters();
  //public static void plotOptimalParameters();

}
