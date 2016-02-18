package dcolorbar;

// Standard Java packages.
import java.awt.*;
import java.awt.image.*;
import java.io.*;
import java.util.*;
import javax.imageio.*;
import javax.swing.*;
import static java.lang.Math.*;

// Mines Java Toolkit (http://www.mines.edu/~dhale/jtk).
import edu.mines.jtk.awt.*;
import edu.mines.jtk.mosaic.*;

/**
 * Some of what you can do with a well-known product for image processing.
 * This class consists of static methods that transform images. Here, images
 * are represented by 2-D arrays or, more precisely, arrays of arrays. All
 * 2-D arrays are <em>regular</em>, in that each row has the same number of
 * columns and each column has the same number of rows.
 * <p>
 * For simplicity this class assumes that image sample values are in the
 * range 0 (black) to 255 (white). When displayed, values less than 0 are 
 * clipped to 0 and values greater than 255 are clipped to 255.
 *
 * @author your-name-here, Colorado School of Mines
 */
public class PhotoPlay {

	/**
	 * Reads, processes, and plots one or more images.
	 * @param args array of file names specified on the command line.
	 */
	public static void main(String[] args) {

		// If no file names, print a helpful message and exit.
		if (args.length<1) {
			System.out.println("Must specify names of one or more image files!");
			System.exit(1);
		}

		// For each file name in the command line arguments, ...
		for (String fileName:args) {

			// Read and plot red, green, and blue components.
			float[][][] rgb = rgbFromFile(fileName);
			float[][] r = rgb[0];
			float[][] g = rgb[1];
			float[][] b = rgb[2];
			//plot(r); // uncomment to look
			//plot(g); // at each of the
			//plot(b); // RGB components

			// Convert to gray image and transpose.
			float[][] f = transpose(gray(r,g,b));

			// Plot gray image and images obtained after various transforms.
			plot(f);
			plot(transpose(f));
			plot(invert(f));
			plot(stretchContrast(0.5f,f));
			plot(edges(f));
			plot(blur(blur(blur(f)))); // very blurry!
		}
	}

	/**
	 * Returns a copy of the specified array. Image transforms below have 
	 * a similar loop structure, but do something more interesting.
	 */
	public static float[][] copy(float[][] f) {
		int n1 = f[0].length;
		int n2 = f.length;
		float[][] g = new float[n2][n1];
		for (int i2=0; i2<n2; ++i2) {
			for (int i1=0; i1<n1; ++i1) {
				g[i2][i1] = f[i2][i1];
			}
		}
		return g;
	}

	/**
	 * Converts red, green, and blue components to grays. The transform used 
	 * here (and in television) is gray = 0.30*red + 0.59*green + 0.11*blue.
	 * @param r array of red components.
	 * @param g array of green components.
	 * @param b array of blue components.
	 * @return output gray image.
	 */
	public static float[][] gray(float[][] r, float[][] g, float[][] b) {
		int n1 = r[0].length;
		int n2 = r.length;
		float[][] gray = new float[n2][n1];
		for (int i2=0; i2<n2; ++i2) {
			for (int i1=0; i1<n1; ++i1) {
				gray[i2][i1] = 0.30f*r[i2][i1]+ 
						       0.59f*g[i2][i1]+
						       0.11f*b[i2][i1];
			}
		}
		return gray; // TODO: fix this! OK
	}

	/**
	 * Returns the transpose of the specified image. In the transpose,
	 * rows become columns and columns become rows. In other words,
	 * output g[i][j] = f[j][i], for all image sample indices (i,j).
	 * @param f input image.
	 * @return output image.
	 */
	public static float[][] transpose(float[][] f) {
		int n2 = f[0].length;
		int n1 = f.length;
		float[][] g = new float[n2][n1];
		for (int i2=0; i2<n2; ++i2) {
			for (int i1=0; i1<n1; ++i1) {
				g[i2][i1] = f[i1][i2]; 
			}
		}
		return g; // TODO: fix this! OK
	}

	/**
	 * Returns the inverse of the specified image. In the inverse,
	 * black becomes white and white becomes black. In other words,
	 * the inverse is like a negative for a photograph.
	 * @param f input image.
	 * @return output image.
	 */
	public static float[][] invert(float[][] f) {
		int n1 = f[0].length;
		int n2 = f.length;
		float[][] g = new float[n2][n1];
		for (int i2=0; i2<n2; ++i2) {
			for (int i1=0; i1<n1; ++i1) {
				g[i2][i1] = 255f-f[i2][i1]; 
			}
		}
		return g; // TODO: fix this! OK
	}

	/**
	 * Returns a contrast-stretched version of the specified image. 
	 * Contrast stretching is a non-linear transform that does nothing 
	 * to medium gray samples, makes light-gray samples lighter and 
	 * dark-gray samples darker.
	 * <p>
	 * Specifically, the transform is defined by three steps:
	 * <pre>
	 * (1) s = (f-128)/128 
	 * (2) t = sgn(s)*(|s|^p)
	 * (3) g = 128+t*128
	 * </pre>
	 * Here, |s|^p denotes the absolute value of s raised to the power p, 
	 * and sgn(s) = 1 or -1, for s positive or negative, respectively.
	 * (To see how this works, sketch t as a function of s for p = 0.5.)
	 * @param p the power.
	 * @param f input image.
	 * @return output image.
	 */
	public static float[][] stretchContrast(float p, float[][] f) {
		int n1 = f[0].length;
		int n2 = f.length;
		float s, t;
		float[][] g = new float[n2][n1];
		for (int i2=0; i2<n2; ++i2) {
			for (int i1=0; i1<n1; ++i1) {
				s = (f[i2][i1]-128f)/128f;
				if (s>0) {
					t = (float) pow(abs(s),p);
				} else {
					t = (float) pow(abs(s),p)*-1f;
				}
				g[i2][i1] = 128f+t*128f; 
			}
		}
		return g; // TODO: fix this! OK
	}

	/**
	 * Returns an image with edges enhanced.
	 * Output samples are computed in three steps:
	 * <pre>
	 * d1 = f[i2][i1+1]-f[i2][i1-1] (like a derivative in 1st dimension)
	 * d2 = f[i2+1][i1]-f[i2-1][i1] (like a derivative in 2nd dimension)
	 * g[i2][i1] = 2*(|d1|+|d2|) 
	 * </pre>
	 * For samples on image boundaries (where this process would yield an
	 * array-index-out-of-bounds exception), the output is zero.
	 * @param f input image.
	 * @return output image.
	 */
	public static float[][] edges(float[][] f) {
		int n1 = f[0].length;
		int n2 = f.length;
		float[][] g = new float[n2][n1];
		g[0][0] = 0;
		g[n2-1][n1-1] = 0;
		for (int i2=1; i2<n2-1; ++i2) {
			for (int i1=1; i1<n1-1; ++i1) {
				float d1 = f[i2][i1+1]-f[i2][i1-1];
				float d2 = f[i2+1][i1]-f[i2-1][i1];
				g[i2][i1] = 2*(abs(d1)+abs(d2)); 
			}
		}
		return g; // TODO: fix this! OK
	}

	/**
	 * Returns an image with features blurred.
	 * Except at the borders of the image, each output sample is a weighted 
	 * average of the nine nearest input samples. The nine weights used in 
	 * the average are
	 * <pre>
	 *   1/16  2/16  1/16
	 *   2/16  4/16  2/16
	 *   1/16  2/16  1/16
	 * </pre>
	 * At the image borders, the output samples equal the input samples.
	 * The weights sum to one, so that constant regions in an image are 
	 * unchanged by this blurring transform.
	 */
	public static float[][] blur(float[][] f) {
		int n1 = f[0].length;
		int n2 = f.length;
		float[][] g = new float[n2][n1];
		g[0][0] = f[0][0];
		g[n2-1][n1-1] = f[n2-1][n1-1];
		for (int i2=1; i2<n2-1; ++i2) {
			for (int i1=1; i1<n1-1; ++i1) {
				g[i2][i1] = (0.0625f)*f[i2-1][i1-1]+(0.125f)*f[i2][i1-1]+(0.0625f)*f[i2+1][i1-1]
						   +(0.125f)*f[i2-1][i1]   +(0.25f)*f[i2][i1]   +(0.125f)*f[i2+1][i1]
						   +(0.0625f)*f[i2-1][i1+1]+(0.125f)*f[i2][i1+1]+(0.0625f)*f[i2+1][i1+1]; 
			}
		}
		return g; // TODO: fix this!
	}

	///////////////////////////////////////////////////////////////////////////
	// private (You need not change anything below this line.)

	/**
	 * Returns an array of 256 grays from black to white.
	 */
	private static Color[] makeGrayColors() {
		Color[] colors = new Color[256];
		for (int i=0; i<256; ++i) {
			int r = i;
			int g = i;
			int b = i;
			colors[i] = new Color(r,g,b);
		}
		return colors;
	}

	/**
	 * Returns an array of 256 random colors.
	 */
	private static Color[] makeRandomColors() {
		Random random = new Random();
		Color[] colors = new Color[256];
		for (int i=0; i<256; ++i) {
			int r = random.nextInt(256);
			int g = random.nextInt(256);
			int b = random.nextInt(256);
			colors[i] = new Color(r,g,b);
		}
		return colors;
	}

	/**
	 * Plots a 2-D regular array as a gray-scale image.
	 * @param f the array to plot.
	 */
	private static void plot(float[][] f)  {
		plot(f,makeGrayColors());
	}

	/**
	 * Plots a 2-D regular array using specified array of 256 colors.
	 * @param f the array to plot.
	 * @param colors array[256] of colors.
	 */
	private static void plot(float[][] f, Color[] colors)  {
		if (colors.length!=256)
			throw new RuntimeException("Need exactly 256 colors!");

		// Width and height of frame to preserve aspect ratio, approximately.
		int n1 = f[0].length;
		int n2 = f.length;
		int size = 600;
		int w,h;
		if (n1<n2) {
			w = size;
			h = w*n1/n2;
		} else {
			h = size;
			w = h*n2/n1;
		}
		int widthAxis = 100;
		int heightAxis = 50;
		int widthColorBar = 70;
		w += widthAxis+widthColorBar;
		h += heightAxis;

		// Plot panel and frame using the package edu.mines.jtk.mosaic.
		PlotPanel panel = new PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT);
		PixelsView pv = panel.addPixels(f);
		pv.setClips(0.0f,255.0f);
		pv.setInterpolation(PixelsView.Interpolation.NEAREST);
		pv.setColorModel(ColorMap.makeIndexColorModel(colors));
		panel.addColorBar();
		panel.setColorBarWidthMinimum(widthColorBar);
		PlotFrame frame = new PlotFrame(panel);
		frame.setVisible(true);
		frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
		frame.setFontSize(24);
		frame.setSize(w,h);
	}

	/**
	 * Returns RGB (red, green, blue) components of an image read from a file.
	 * @param fileName the name of the file containing a JPEG image.
	 * @return array of RGB components; the components (red,green,blue) are 
	 * stored in {rgb[0],rgb[1],rgb[2]}.
	 */
	private static float[][][] rgbFromFile(String fileName) {
		BufferedImage bi = null;
		try {
			File file = new File(fileName);
			bi = ImageIO.read(file);
		} catch (IOException ioe) {
			throw new RuntimeException(ioe);
		}
		int w = bi.getWidth();
		int h = bi.getHeight();
		int[] pixels = bi.getRGB(0,0,w,h,null,0,w);
		float[][][] rgb = new float[3][h][w];
		for (int i=0,k=0; i<h; ++i) {
			for (int j=0; j<w; ++j,++k) {
				int p = pixels[k];
				int r = (p>>16)&0xff;
				int g = (p>> 8)&0xff;
				int b = (p>> 0)&0xff;
				rgb[0][i][j] = (float)r;
				rgb[1][i][j] = (float)g;
				rgb[2][i][j] = (float)b;
			}
		}
		return rgb;
	}
}
