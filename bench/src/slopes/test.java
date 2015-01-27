package slopes;

public class test {
  public static void main(String[] args) {
    StructureTensor str = new StructureTensor(24.0f,1.0f);
    SetParameters.setChickenTestParameters(0.7f);
    str.estimateSlope();
  }
}
