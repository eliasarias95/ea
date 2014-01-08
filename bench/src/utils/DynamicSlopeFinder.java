import dnp.*;
import edu.mines.jtk.dsp.*;

public class DynamicSlopeFinder {
  public static void main(String[] args) {
    double sigma1 = 2.0;
    LocalSlopeFinder lsf = new LocalSlopeFinder(sigma1);
  }
}
