package dcolorbar;

import java.awt.Color;
import java.io.IOException;
import java.nio.ByteOrder;
import java.util.Random;

import static edu.mines.jtk.util.ArrayMath.cabs;
import static edu.mines.jtk.util.ArrayMath.div;
import static edu.mines.jtk.util.ArrayMath.max;
import static edu.mines.jtk.util.ArrayMath.sqrt;
import static java.lang.Math.*;

import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.SincInterp;
import edu.mines.jtk.io.ArrayFile;
import edu.mines.jtk.io.ArrayInputStream;
import edu.mines.jtk.mosaic.PlotFrame;
import edu.mines.jtk.mosaic.PlotPanel;
import edu.mines.jtk.mosaic.PointsView;
import edu.mines.jtk.mosaic.SimplePlot;

public class SOM2 {
	private int iter, sIter, ic = 0, jc = 0;
	private float alpha = 0, radius, rStart, sigma = 0;
	private static int n_att;
	private float[][][] nodes;

	private float[][][][] distanceTable;
	float[] pickedData;
	private int n1som, n2som;
	private static float[][][] trainData2D;
	private static float[][][][] trainData3D;
	private String filesLocations;
	private static String baseName;
	private static int fileN1, fileN2, fileN3;

	/**
	 * 
	 * @param numIterations
	 * @param hLength
	 * @param vLength
	 * @param numAttributes
	 * @param filesLocations
	 *            For example, if the attribute files location is
	 *            "C:/Users/Chris/workspace/Harmonic_Attributes/TPMorlet0.dat",
	 *            then the filesLocations input will be
	 *            "C:/Users/Chris/workspace/Harmonic_Attributes/". Then the
	 *            baseName input would be TPMorlet, the 0 represents the number
	 *            of the attribute.
	 */
	public SOM2(int numIterations, int hLength, int vLength, int numAttributes,
			String filesLocations, String baseName) {
		this.iter = numIterations;
		this.n_att = numAttributes;
		this.filesLocations = filesLocations;
		this.baseName = baseName;

		n1som = hLength;// hLength;
		n2som = vLength;// vLength;
		nodes = new float[n2som][n1som][numAttributes];
		randomizeAllNodesWeights(nodes);

		if (hLength < vLength) {
			this.radius = vLength * .5f + 1.f;
		} else {
			this.radius = hLength * .5f + 1.f;
		}
		this.rStart = this.radius;

		distanceTable = new float[n2som][n1som][n2som][n1som];
		fillDistanceTable();
	}

	public SOM2(int numIterations, int hLength, int vLength, int numAttributes) {
		this.iter = numIterations;
		this.n_att = numAttributes;
		this.filesLocations = filesLocations;
		this.baseName = baseName;

		n1som = hLength;// hLength;
		n2som = vLength;// vLength;
		nodes = new float[n2som][n1som][numAttributes];
		randomizeAllNodesWeights(nodes);

		if (hLength < vLength) {
			this.radius = vLength * .5f + 1.f;
		} else {
			this.radius = hLength * .5f + 1.f;
		}
		this.rStart = this.radius;

		distanceTable = new float[n2som][n1som][n2som][n1som];
		fillDistanceTable();
	}

	private static void randomizeAllNodesWeights(float[][][] nodes) {

		int n1 = nodes.length;
		int n2 = nodes[0].length;
		for (int i = 0; i < n1; ++i) {
			for (int j = 0; j < n2; ++j) {
				randomizeNodeWeight(nodes[i][j]);
			}
		}
	}

	private static void randomizeNodeWeight(float[] node) {
		Random r = new Random();
		for (int i = 0; i < n_att; ++i) {
			node[i] = 2 * r.nextFloat() - 1;// Initialize floats in square from
			System.out.println(node[i]);
		}
	}

	/**
	 * Normalizes each attribute set of the input data so that the mean is zero
	 * and the standard deviation is 0.
	 * 
	 * @param data
	 *            2D data with a set of attributes for each point
	 * @return
	 */
	public float[][][] norm2DData(float[][][] data) {
		int n3 = data.length;
		int n2 = data[0].length;
		int n1 = data[0][0].length;

		float[] sum = new float[n1];
		float[] mean = new float[n1];

		float[] diffMeanSquarSum = new float[n1];// (x1-xmean)*(x1-xmean) +
													// (x2-xmean)*(x2-xmean) +
													// (xn-xmean)*(xn-xmean)
		float[] stdDev = new float[n1];

		float[][][] normData = new float[n3][n2][n1];

		// Find sum of each attribute
		for (int i1 = 0; i1 < n1; ++i1) {
			// System.out.println(n1);

			for (int i3 = 0; i3 < n3; ++i3) {
				for (int i2 = 0; i2 < n2; ++i2) {
					sum[i1] += data[i3][i2][i1];
				}
			}
		}

		// Find Mean
		mean = div(sum, (float) (n3 * n2));

		// Find (x1-xmean)*(x1-xmean) + (x2-xmean)*(x2-xmean) +
		// (xn-xmean)*(xn-xmean) for each attribute, the inner sum of finding
		// the standard deviation or variance.
		for (int i1 = 0; i1 < n1; ++i1) {
			for (int i3 = 0; i3 < n3; ++i3) {
				for (int i2 = 0; i2 < n2; ++i2) {
					diffMeanSquarSum[i1] += (data[i3][i2][i1] - mean[i1])
							* (data[i3][i2][i1] - mean[i1]);
				}
			}
		}

		// Find Standard Deviation
		stdDev = sqrt(div(diffMeanSquarSum, n3 * n2 - 1));

		// Normalize Data
		for (int i1 = 0; i1 < n1; ++i1) {
			for (int i3 = 0; i3 < n3; ++i3) {
				for (int i2 = 0; i2 < n2; ++i2) {
					normData[i3][i2][i1] = (data[i3][i2][i1] - mean[i1])
							/ stdDev[i1];
				}
			}
		}

		/*
		 * for (int i1=0; i1<n1; ++i1){
		 * System.out.println("mean["+i1+"] = "+mean
		 * [i1]+" sum["+i1+"] = "+sum[i1]); }
		 */

		/****************************** Testing purposes only ************************************/
		/*
		 * float[] sum1 = new float[n1]; float[] mean1 = new float[n1];
		 * 
		 * float[] diffMeanSquarSum1 = new float[n1];//(x1-xmean)*(x1-xmean) +
		 * (x2-xmean)*(x2-xmean) + (xn-xmean)*(xn-xmean) float[] stdDev1 = new
		 * float[n1];
		 * 
		 * 
		 * 
		 * for (int i1=0; i1<n1; ++i1){ for (int i3=0; i3<n3; ++i3){ for (int
		 * i2=0; i2<n2; ++i2){ sum1[i1] += normData[i3][i2][i1]; } } } //Find
		 * Mean mean1 = div(sum1, (float)(n3*n2));
		 * 
		 * //Find (x1-xmean)*(x1-xmean) + (x2-xmean)*(x2-xmean) +
		 * (xn-xmean)*(xn-xmean) for each attribute, the inner sum of finding
		 * the standard deviation or variance. for (int i1=0; i1<n1; ++i1){ for
		 * (int i3=0; i3<n3; ++i3){ for (int i2=0; i2<n2; ++i2){
		 * diffMeanSquarSum1[i1] += (normData[i3][i2][i1] -
		 * mean1[i1])*(normData[i3][i2][i1] - mean1[i1]); } } }
		 * 
		 * //Find Standard Deviation stdDev1 = sqrt(div(diffMeanSquarSum1,
		 * n3*n2-1));
		 * 
		 * for (int i1=0; i1<n1; ++i1){
		 * System.out.println("mean["+i1+"] = "+mean1
		 * [i1]+" stddev["+i1+"] = "+stdDev1[i1]); }
		 */
		return normData;

	}

	/**
	 * Normalizes each attribute set of the input data so that the mean is zero
	 * and the standard deviation is 0.
	 * 
	 * @param data
	 *            3D data with a set of attributes for each point
	 * @return
	 */
	public float[][][][] norm3DData(float[][][][] data) {
		int n4 = data.length;
		int n3 = data[0].length;
		int n2 = data[0][0].length;
		int n1 = data[0][0][0].length;

		float[] sum = new float[n1];
		float[] mean = new float[n1];

		float[] diffMeanSquarSum = new float[n1];// (x1-xmean)*(x1-xmean) +
													// (x2-xmean)*(x2-xmean) +
													// (xn-xmean)*(xn-xmean)
		float[] stdDev = new float[n1];

		// Find sum of each attribute
		for (int i1 = 0; i1 < n1; ++i1) {
			for (int i4 = 0; i4 < n4; ++i4) {
				for (int i3 = 0; i3 < n3; ++i3) {
					for (int i2 = 0; i2 < n2; ++i2) {
						sum[i1] += data[i4][i3][i2][i1];
					}
				}
			}
			// System.out.println(sum[i1]);
		}

		// Find Mean
		mean = div(sum, (float) (n4 * n3 * n2));
		// System.out.println("mean = "+mean);

		// Find (x1-xmean)*(x1-xmean) + (x2-xmean)*(x2-xmean) +
		// (xn-xmean)*(xn-xmean) for each attribute, the inner sum of finding
		// the standard deviation or variance.
		for (int i1 = 0; i1 < n1; ++i1) {
			for (int i4 = 0; i4 < n4; ++i4) {
				for (int i3 = 0; i3 < n3; ++i3) {
					for (int i2 = 0; i2 < n2; ++i2) {
						diffMeanSquarSum[i1] += (data[i4][i3][i2][i1] - mean[i1])
								* (data[i4][i3][i2][i1] - mean[i1]);
					}
				}
			}
		}

		// Find Standard Deviation
		stdDev = sqrt(div(diffMeanSquarSum, n4 * n3 * n2));

		// Normalize Data
		for (int i1 = 0; i1 < n1; ++i1) {
			for (int i4 = 0; i4 < n4; ++i4) {
				for (int i3 = 0; i3 < n3; ++i3) {
					for (int i2 = 0; i2 < n2; ++i2) {
						data[i4][i3][i2][i1] = (data[i4][i3][i2][i1] - mean[i1])
								/ stdDev[i1];
					}
				}
			}
		}

		for (int i1 = 0; i1 < n1; ++i1) {
			System.out.println("mean[" + i1 + "] = " + mean[i1] + " sum[" + i1
					+ "] = " + sum[i1]);
		}

		/****************************** Testing purposes only ************************************/
		float[] sum1 = new float[n1];
		float[] mean1 = new float[n1];

		float[] diffMeanSquarSum1 = new float[n1];// (x1-xmean)*(x1-xmean) +
													// (x2-xmean)*(x2-xmean) +
													// (xn-xmean)*(xn-xmean)
		float[] stdDev1 = new float[n1];

		for (int i1 = 0; i1 < n1; ++i1) {
			for (int i4 = 0; i4 < n4; ++i4) {
				for (int i3 = 0; i3 < n3; ++i3) {
					for (int i2 = 0; i2 < n2; ++i2) {
						sum1[i1] += data[i4][i3][i2][i1];
					}
				}
			}
		}
		// Find Mean
		mean1 = div(sum1, (float) (n4 * n3 * n2));

		// Find (x1-xmean)*(x1-xmean) + (x2-xmean)*(x2-xmean) +
		// (xn-xmean)*(xn-xmean) for each attribute, the inner sum of finding
		// the standard deviation or variance.
		for (int i1 = 0; i1 < n1; ++i1) {
			for (int i4 = 0; i4 < n4; ++i4) {
				for (int i3 = 0; i3 < n3; ++i3) {
					for (int i2 = 0; i2 < n2; ++i2) {
						diffMeanSquarSum1[i1] += (data[i4][i3][i2][i1] - mean1[i1])
								* (data[i4][i3][i2][i1] - mean1[i1]);
					}
				}
			}
		}

		// Find Standard Deviation
		stdDev1 = sqrt(div(diffMeanSquarSum1, n4 * n3 * n2 - 1));

		for (int i1 = 0; i1 < n1; ++i1) {
			System.out.println("mean[" + i1 + "] = " + mean1[i1] + " stddev["
					+ i1 + "] = " + stdDev1[i1]);
		}

		return data;

	}

	private int i1(int node) {
		return node % n1som;
	}

	private int i2(int node) {
		return node / n1som;
	}

	/*
	 * Should this be private? If so how can we do that? Refer to Run.java line
	 * 132
	 */
	public int node(int i1, int i2) {
		return i1 + i2 * n1som;
	}

	private static float[] randomData2D() {
		Random r = new Random();
		int n1 = trainData2D[0].length;
		int n2 = trainData2D.length;
		int i1 = r.nextInt(n1);
		int i2 = r.nextInt(n2);
		return trainData2D[i2][i1];
	}

	private static float[] randomData3D() {
		Random r = new Random();
		int n1 = trainData3D[0][0].length;
		int n2 = trainData3D[0].length;
		int n3 = trainData3D.length;
		int i1 = r.nextInt(n1);
		int i2 = r.nextInt(n2);
		int i3 = r.nextInt(n3);
		return trainData3D[i3][i2][i1];
	}

	/**
	 * Location and name of the data that has been split up into different files
	 * where each file is a different attribute.
	 * 
	 * @return
	 */
	private static float[] randomData3DFromFile() {
		Random r = new Random();
		int n1 = fileN1;
		int n2 = fileN2;
		int n3 = fileN3;
		int i1 = r.nextInt(n1);
		int i2 = r.nextInt(n2);
		int i3 = r.nextInt(n3);
		float[] randomData = new float[n_att];
		ArrayInputStream[] ais = new ArrayInputStream[n_att];
		for (int atr = 0; atr < n_att; ++atr) {
			try {
				ais[atr] = new ArrayInputStream(baseName + atr + ".dat",
						ByteOrder.BIG_ENDIAN);
			} catch (IOException e) {
				System.out.println("IOProblem");
			}
		}
		for (int atr = 0; atr < n_att; ++atr) {

			// Create a writer to write each frequency data set out to a file
			try {

				int skipBytes = (i3 * n2 * n1 + i2 * n1 + i1) * 4;
				ais[atr].skipBytes(skipBytes);
				randomData[atr] = ais[atr].readFloat();
				ais[atr].close();
			} catch (IOException e) {
				System.out.println("IOProblem");
			}

		}

		return randomData;
	}

	public void train2D(float[][][] trainData) {
		this.trainData2D = trainData;// norm2DData(trainData);
		iterationManager2DData();
		System.out.println("rStart = " + rStart);
	}

	public void train3D(float[][][][] trainData) {
		this.trainData3D = norm3DData(trainData);
		iterationManager3DData();
		System.out.println("rStart = " + rStart);
	}

	public void train3DFile(int fileN3, int fileN2, int fileN1) {
		this.fileN3 = fileN3;
		this.fileN2 = fileN2;
		this.fileN1 = fileN1;
		System.out.println("Start Training");
		iterationManager3DDataFile();
		System.out.println("End Training");

		System.out.println("rStart = " + rStart);
	}

	public void iterationManager2DData() {
		for (int t = 0; t < iter; ++t) {
			/*********** Start Plotting Purposes **************/
			if (t == 0 || t == 19 || t == 99 || t == 999 || t == 4999
					|| t == iter - 1) {
				float[][] weightr1 = new float[n1som][n2som];
				float[][] weightr2 = new float[n1som][n2som];
				float[][] weightc1 = new float[n2som][n1som];
				float[][] weightc2 = new float[n2som][n1som];
				for (int i = 0; i < n2som; ++i) {
					for (int j = 0; j < n1som; ++j) {
						weightc1[i][j] = nodes[i][j][0];
						weightc2[i][j] = nodes[i][j][1];
						weightr1[j][i] = nodes[i][j][0];
						weightr2[j][i] = nodes[i][j][1];
					}
				}
				String iterations = Integer.toString(t + 1);
				PointsView pointsr = new PointsView(weightr1, weightr2);
				PointsView pointsc = new PointsView(weightc1, weightc2);

				pointsr.setMarkColor(java.awt.Color.RED);
				pointsr.setMarkStyle(PointsView.Mark.CROSS);
				pointsr.setMarkSize(4);
				pointsc.setMarkColor(java.awt.Color.RED);
				pointsc.setMarkStyle(PointsView.Mark.CROSS);
				pointsc.setMarkSize(4);
				PlotPanel pp = new PlotPanel();
				pp.addTiledView(pointsr);
				pp.addTiledView(pointsc);
				pp.setHLabel("Weight 1");
				pp.setVLabel("Weight 2");
				float max = max(trainData2D) * .75f;
				pp.setLimits(-max, -max, max, max);
				// pp.setLimits(0, 0, 7, 7); // change to 5,5 for seismic data

				pp.addTitle("Iterations = " + iterations);
				// pp.setLimits(0, 0, 1, 1);
				PlotFrame pf = new PlotFrame(pp);
				pf.setVisible(true);
			}
			/*********** End Plotting Purposes **************/

			pickedData = randomData2D();
			findWinner(pickedData);
			if (t < iter * .1f) {
				alpha = (float) (0.9 * (1. - (t / (float) (iter * .1))));
				sigma = (float) (.33355f)
						* (rStart - ((rStart - 1) / (iter * .1f)) * t);// (2-(5.*t)/(float)(3.*iter));//(float)
																		// ((sRadius/3.0-1)/(-(float)iter))*((1-sRadius)/(-sRadius+3))*t+(sRadius/3f);//(2-(5.*t)/(float)(3.*iter));

			} else {
				alpha = .01f;
				sigma = .33355f;
			}

			radius = 3.f * sigma;// (float)
									// (exp(1.5*(1-t/(float)iter)));//(float)3.*sigma;//(exp(1.5*(1-t/(float)iter)));

			// System.out.println("alpha= "+alpha);
			updateNodes(ic, jc);

			// System.out.println("alpha= "+alpha);

		}
	}

	public void iterationManager3DData() {
		for (int t = 0; t < iter; ++t) {
			/*********** Start Plotting Purposes **************/
			/*
			 * if (t==0 || t==19 || t==99 || t==999 || t==4999 || t==iter-1){
			 * float[][] weightr1 = new float[n1som][n2som]; float[][] weightr2
			 * = new float[n1som][n2som]; float[][] weightc1 = new
			 * float[n2som][n1som]; float[][] weightc2 = new
			 * float[n2som][n1som]; for (int i=0; i<n2som; ++i){ for (int j=0;
			 * j<n1som; ++j){ weightc1[i][j] = nodes[i][j][0]; weightc2[i][j] =
			 * nodes[i][j][1]; weightr1[j][i] = nodes[i][j][0]; weightr2[j][i] =
			 * nodes[i][j][1]; } } String iterations = Integer.toString(t+1);
			 * PointsView pointsr = new PointsView(weightr1,weightr2);
			 * PointsView pointsc = new PointsView(weightc1,weightc2);
			 * 
			 * pointsr.setMarkColor(java.awt.Color.RED);
			 * pointsr.setMarkStyle(PointsView.Mark.CROSS);
			 * pointsr.setMarkSize(4); pointsc.setMarkColor(java.awt.Color.RED);
			 * pointsc.setMarkStyle(PointsView.Mark.CROSS);
			 * pointsc.setMarkSize(4); PlotPanel pp = new PlotPanel();
			 * pp.addTiledView(pointsr); pp.addTiledView(pointsc);
			 * pp.setHLabel("Weight 1"); pp.setVLabel("Weight 2"); float pmax =
			 * 0; for (int i4=0; i4<trainData3D.length; ++i4){ for (int i3=0;
			 * i3<trainData3D[0].length; ++i3){ for (int i2=0;
			 * i2<trainData3D[0][0].length; ++i2){ for (int i1=0;
			 * i1<trainData3D[0][0][0].length; ++i1){ if (pmax <
			 * trainData3D[i4][i3][i2][i1]) pmax = trainData3D[i4][i3][i2][i1];
			 * } } } } float max = pmax*.75f; pp.setLimits(-max,-max,max,max);
			 * //pp.setLimits(0, 0, 7, 7); // change to 5,5 for seismic data
			 * 
			 * pp.addTitle("Iterations = "+iterations); //pp.setLimits(0, 0, 1,
			 * 1); PlotFrame pf = new PlotFrame(pp); pf.setVisible(true); }
			 *//*********** End Plotting Purposes **************/

			pickedData = randomData3D();
			findWinner(pickedData);
			if (t < iter * .1f) {
				alpha = (float) (0.9 * (1. - (t / (float) (iter * .1))));
				sigma = (float) (.33355f)
						* (rStart - ((rStart - 1) / (iter * .1f)) * t);// (2-(5.*t)/(float)(3.*iter));//(float)
																		// ((sRadius/3.0-1)/(-(float)iter))*((1-sRadius)/(-sRadius+3))*t+(sRadius/3f);//(2-(5.*t)/(float)(3.*iter));

			} else {
				alpha = .01f;
				sigma = .33355f;
			}

			radius = 3.f * sigma;// (float)
									// (exp(1.5*(1-t/(float)iter)));//(float)3.*sigma;//(exp(1.5*(1-t/(float)iter)));

			// System.out.println("alpha= "+alpha);
			updateNodes(ic, jc);

			// System.out.println("alpha= "+alpha);

		}
	}

	public void iterationManager3DDataFile() {
		for (int t = 0; t < iter; ++t) {
			/*********** Start Plotting Purposes **************/
			if (t == 0 || t == 19 || t == 99 || t == 999 || t == 4999
					|| t == iter - 1) {
				float[][] weightr1 = new float[n1som][n2som];
				float[][] weightr2 = new float[n1som][n2som];
				float[][] weightc1 = new float[n2som][n1som];
				float[][] weightc2 = new float[n2som][n1som];
				for (int i = 0; i < n2som; ++i) {
					for (int j = 0; j < n1som; ++j) {
						weightc1[i][j] = nodes[i][j][0];
						weightc2[i][j] = nodes[i][j][1];
						weightr1[j][i] = nodes[i][j][0];
						weightr2[j][i] = nodes[i][j][1];
					}
				}
				String iterations = Integer.toString(t + 1);
				PointsView pointsr = new PointsView(weightr1, weightr2);
				PointsView pointsc = new PointsView(weightc1, weightc2);

				pointsr.setMarkColor(java.awt.Color.RED);
				pointsr.setMarkStyle(PointsView.Mark.CROSS);
				pointsr.setMarkSize(4);
				pointsc.setMarkColor(java.awt.Color.RED);
				pointsc.setMarkStyle(PointsView.Mark.CROSS);
				pointsc.setMarkSize(4);
				PlotPanel pp = new PlotPanel();
				pp.addTiledView(pointsr);
				pp.addTiledView(pointsc);
				pp.setHLabel("Weight 1");
				pp.setVLabel("Weight 2");

				pp.addTitle("Iterations = " + iterations);
				// pp.setLimits(0, 0, 1, 1);
				PlotFrame pf = new PlotFrame(pp);
				pf.setVisible(true);
			}
			/*********** End Plotting Purposes **************/

			pickedData = randomData3DFromFile();
			findWinner(pickedData);
			if (t < iter * .1f) {
				alpha = (float) (0.9 * (1. - (t / (float) (iter * .1))));
				sigma = (float) (.33355f)
						* (rStart - ((rStart - 1) / (iter * .1f)) * t);// (2-(5.*t)/(float)(3.*iter));//(float)
																		// ((sRadius/3.0-1)/(-(float)iter))*((1-sRadius)/(-sRadius+3))*t+(sRadius/3f);//(2-(5.*t)/(float)(3.*iter));

			} else {
				alpha = .01f;
				sigma = .33355f;
			}

			radius = 3.f * sigma;// (float)
									// (exp(1.5*(1-t/(float)iter)));//(float)3.*sigma;//(exp(1.5*(1-t/(float)iter)));

			// System.out.println("alpha= "+alpha);
			updateNodes(ic, jc);

			// System.out.println("alpha= "+alpha);

		}
	}

	public void updateNodes(int ic, int jc) {
		float h = 0;
		for (int i = 0; i < n2som; ++i) {
			for (int j = 0; j < n1som; ++j) {
				updateNode(ic, jc, i, j);
			}
		}
	}

	public void updateNode(int ic, int jc, int ii, int jj) {
		float h = 0;
		float distance = distanceTable[ic][jc][ii][jj];
		if (distance == 0) {
			for (int i = 0; i < n_att; ++i) {
				// h = (float)
				// (alpha*Math.exp(-Math.pow(distance(ic,jc,ii,jj),2)/Math.pow(alpha,2)));
				h = (float) (alpha * exp(-pow(distanceTable[ic][jc][ii][jj], 2)
						/ (2. * pow(sigma, 2))));

				// System.out.println("h= "+h);
				nodes[ic][jc][i] = nodes[ic][jc][i] + h
						* (pickedData[i] - nodes[ic][jc][i]);

			}
		} else if (distance <= radius) {
			for (int i = 0; i < n_att; ++i) {
				// h = (float)
				// (alpha*Math.exp(-Math.pow(distance(ic,jc,ii,jj),2)/Math.pow(alpha,2)));
				h = (float) (alpha * exp(-pow(distanceTable[ic][jc][ii][jj], 2)
						/ (2. * pow(sigma, 2))));
				// System.out.println("h= "+distance(ic,jc,ii,jj));
				nodes[ii][jj][i] = nodes[ii][jj][i] + h
						* (pickedData[i] - nodes[ii][jj][i]);

			}
		}
	}

	public float[][] classify2DData(float[][][] inputData) {
		int n1 = inputData.length;
		int n2 = inputData[0].length;
		// float[][][] normData = norm2DData(inputData);
		float[][] classifiedData = new float[n1][n2];
		for (int i = 0; i < n1; ++i) {
			for (int j = 0; j < n2; ++j) {
				classifiedData[i][j] = classifyOneData(inputData[i][j]);

			}
		}

		return classifiedData;
	}

	public float[][][] classify3DData(float[][][][] inputData) {
		int n3 = inputData.length;
		int n2 = inputData[0].length;
		int n1 = inputData[0][0].length;
		float[][][][] normData = norm3DData(inputData);
		float[][][] classifiedData = new float[n3][n2][n1];
		for (int i3 = 0; i3 < n3; ++i3) {
			for (int i2 = 0; i2 < n2; ++i2) {
				for (int i1 = 0; i1 < n1; ++i1) {
					classifiedData[i3][i2][i1] = classifyOneData(normData[i3][i2][i1]);
				}
			}
		}

		return classifiedData;
	}

	public float[][][] classify3DDataFile() {
		int n3 = fileN3;
		int n2 = fileN2;
		int n1 = fileN1;
		// float[][][][] normData = norm3DData(inputData);
		float[][][] classifiedData = new float[n3][n2][n1];
		float[] oneSetAtr = new float[n_att];

		ArrayInputStream[] ais = new ArrayInputStream[n_att];
		for (int atr = 0; atr < n_att; ++atr) {
			try {
				ais[atr] = new ArrayInputStream(baseName + atr + ".dat",
						ByteOrder.BIG_ENDIAN);
			} catch (IOException e) {
				System.out.println("IOProblem");
			}
		}
		for (int i3 = 0; i3 < n3; ++i3) {
			for (int i2 = 0; i2 < n2; ++i2) {
				for (int i1 = 0; i1 < n1; ++i1) {
					for (int atr = 0; atr < n_att; ++atr) {

						// Create a writer to read each frequency from a file
						try {
							oneSetAtr[atr] = ais[atr].readFloat();
							// System.out.println("attr"+atr+" = "+oneSetAtr[atr]);
						} catch (IOException e) {
							System.out.println("IOProblem");
						}
					}
					classifiedData[i3][i2][i1] = classifyOneData(oneSetAtr);
					// System.out.println(classifiedData[i3][i2][i1]);
				}
			}
			System.out.println("i3 = " + i3);
		}

		return classifiedData;
	}

	/*
	 * Changed this method to have winCatNum equal the value that is spit out
	 * from the node() method.
	 */
	private int classifyOneData(float[] inputData) {
		int n1 = nodes.length;
		int n2 = nodes[0].length;
		int winCatNum = 0;
		float minEucDis = 999999999;
		float eucDis = 0;
		for (int i = 0; i < n1; ++i) {
			for (int j = 0; j < n2; ++j) {
				eucDis = euclideanDistance(inputData, nodes[i][j]);
				if (eucDis < minEucDis) {
					minEucDis = eucDis;
					winCatNum = node(j, i);
				}

			}
		}
		// System.out.println("Wincatnum= "+winCatNum);
		return winCatNum;
	}

	private void fillDistanceTable() {
		for (int i1 = 0; i1 < n2som; ++i1)
			for (int i2 = 0; i2 < n1som; ++i2)
				for (int i3 = 0; i3 < n2som; ++i3)
					for (int i4 = 0; i4 < n1som; ++i4) {
						distanceTable[i1][i2][i3][i4] = (float) (sqrt(pow(
								(i1 - i3), 2) + pow((i2 - i4), 2)));
					}
	}

	private void findWinner(float[] singleTrainData) {
		float xmin = 99f;

		int n2 = nodes.length;
		int n1 = nodes[0].length;
		for (int i = 0; i < n2; ++i) {
			for (int j = 0; j < n1; ++j) {
				float eucDis = euclideanDistance(singleTrainData, nodes[i][j]);
				if (eucDis < xmin) {
					xmin = eucDis;
					ic = i;
					jc = j;
				}
			}
		}
		// System.out.println(xmin);
	}

	private float euclideanDistance(float[] input, float[] node) {
		float eucDis;
		float diff = 0; // the euclidean distance without the square root
		for (int i = 0; i < n_att; ++i) {
			diff += (input[i] - node[i]) * (input[i] - node[i]);

		}
		eucDis = (float) Math.sqrt((double) diff);
		return eucDis;
	}

	/**
	 * The weights are re-randomized and the category numbers are re-issued.
	 */
	public void resetWeights() {
		nodes = new float[n2som][n1som][n_att];
		randomizeAllNodesWeights(nodes);
	}

	public float[][][] getNodes() {
		return nodes;
	}

}
