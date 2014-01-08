import java.util.Scanner;

public class Grades {
  public static void main(String[] args) {
    Scanner keyboard = new Scanner(System.in);
    System.out.println("This program will tell you what you need on your");
    System.out.println("final exam to get an A, B, or C in your course.");
    System.out.println("This program assumes whatever percent left out");
    System.out.println("of graded values, is the percent your final exam");
    System.out.println("counts for.");
    System.out.println("How many sections is your grade broken down into?");
    System.out.println("For example, 2 midterms, homework, and a final");
    System.out.println("would be 4 sections.");
    int n = keyboard.nextInt();
    float[] g = new float[n]; // grades
    float[] p = new float[n]; // percentages
    System.out.println("Please enter the weighted percentage of your");
    System.out.println("first section of grades (say, your first exam");
    System.out.println("of the course) in decimal form.");
    p[0] = keyboard.nextFloat();
    System.out.println("Please enter the grade you received.");
    g[0] = keyboard.nextFloat();
    for (int i=1; i<n-1; ++i) {
      System.out.println("Please enter the weighted percentage of your");
      System.out.println("next section of grades.");
      p[i] = keyboard.nextFloat();
      System.out.println("Please enter the grade you received.");
      g[i] = keyboard.nextFloat();
    }
    float[] fg = calculate(p,g);
    System.out.println("You need a "+fg[2]+" to get an A in the course.");
    System.out.println("You need a "+fg[1]+" to get a B in the course.");
    System.out.println("You need a "+fg[0]+" to get a C in the course.");
  }

  private static float[] calculate(float[] p, float[] g) {
    int n = p.length;
    float[] fg = new float[3]; // final grade values
    float sum = 0;
    float psum = 0;
    for (int i=0; i<n; ++i) {
      sum += g[i]*p[i];
      psum += p[i];
    }
    fg[0] = (70.0f-sum)/(1.0f-psum);
    fg[1] = (80.0f-sum)/(1.0f-psum);
    fg[2] = (90.0f-sum)/(1.0f-psum);
    return fg;
  }
}
