import javax.swing.JApplet;
import java.awt.Graphics;

public class DrawRings extends JApplet {

  private static final long serialVersionUID = 1L;
  public void paint(Graphics pen) {
    pen.drawOval(200, 200, 150, 150);
    pen.drawOval(375, 200, 150, 150);
    pen.drawOval(550, 200, 150, 150);
    pen.drawOval(287, 300, 150, 150);
    pen.drawOval(463, 300, 150, 150);
  }
  public static void main(String[] args) {
    paint(Graphics);
  }
}
