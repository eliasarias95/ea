package dcolorbar;

import java.awt.Color;
import java.util.Random;

import static edu.mines.jtk.util.ArrayMath.div;
import static edu.mines.jtk.util.ArrayMath.max;
import static edu.mines.jtk.util.ArrayMath.sqrt;
import static java.lang.Math.*;

import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.SincInterp;
import edu.mines.jtk.mosaic.PlotFrame;
import edu.mines.jtk.mosaic.PlotPanel;
import edu.mines.jtk.mosaic.PointsView;
import edu.mines.jtk.mosaic.SimplePlot;

public class SOM1 {
	private int iter, sIter, ic = 0, jc = 0;
	private float alpha = 0, radius, rStart, sigma = 0, maxH, maxS;
	private static int n_att;
	private float[][] nodes;
	private static float[][][] trainData2D;
	private static float[][][][] trainData3D;

	private float[][] distanceTable;
	float[] pickedData;
	private int n1som;

	public SOM1(int numIterations, int hLength, int numAttributes) {
		this.iter = numIterations;
		this.n_att = numAttributes;

		n1som = hLength;// hLength;
		nodes = new float[n1som][numAttributes];
		randomizeAllNodesWeights(nodes);

		this.radius = hLength * .5f + 1.f;
		this.rStart = this.radius;

		distanceTable = new float[n1som][n1som];
		fillDistanceTable();
	}

	private static void randomizeAllNodesWeights(float[][] nodes) {

		int n1 = nodes.length;
		for (int i = 0; i < n1; ++i) {
			randomizeNodeWeight(nodes[i]);
		}
	}

	private static void randomizeNodeWeight(float[] node) {
		Random r = new Random();
		int sign = 0;
		for (int i = 0; i < n_att; ++i) {
			if (r.nextInt() % 2 == 0)
				sign = 1;
			else {
				sign = -1;
			}
			node[i] = sign * r.nextFloat();// Initialize floats in square from
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
		float[] sum1 = new float[n1];
		float[] mean1 = new float[n1];

		float[] diffMeanSquarSum1 = new float[n1];// (x1-xmean)*(x1-xmean) +
													// (x2-xmean)*(x2-xmean) +
													// (xn-xmean)*(xn-xmean)
		float[] stdDev1 = new float[n1];

		for (int i1 = 0; i1 < n1; ++i1) {
			for (int i3 = 0; i3 < n3; ++i3) {
				for (int i2 = 0; i2 < n2; ++i2) {
					sum1[i1] += normData[i3][i2][i1];
				}
			}
		}
		// Find Mean
		mean1 = div(sum1, (float) (n3 * n2));

		// Find (x1-xmean)*(x1-xmean) + (x2-xmean)*(x2-xmean) +
		// (xn-xmean)*(xn-xmean) for each attribute, the inner sum of finding
		// the standard deviation or variance.
		for (int i1 = 0; i1 < n1; ++i1) {
			for (int i3 = 0; i3 < n3; ++i3) {
				for (int i2 = 0; i2 < n2; ++i2) {
					diffMeanSquarSum1[i1] += (normData[i3][i2][i1] - mean1[i1])
							* (normData[i3][i2][i1] - mean1[i1]);
				}
			}
		}

		// Find Standard Deviation
		stdDev1 = sqrt(div(diffMeanSquarSum1, n3 * n2 - 1));

		for (int i1 = 0; i1 < n1; ++i1) {
			System.out.println("mean[" + i1 + "] = " + mean1[i1] + " stddev["
					+ i1 + "] = " + stdDev1[i1]);
		}

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

		float[][][][] normData = new float[n4][n3][n2][n1];

		// Find sum of each attribute
		for (int i1 = 0; i1 < n1; ++i1) {
			for (int i4 = 0; i4 < n4; ++i4) {
				for (int i3 = 0; i3 < n3; ++i3) {
					for (int i2 = 0; i2 < n2; ++i2) {
						sum[i1] += data[i4][i3][i2][i1];
					}
				}
			}
		}

		// Find Mean
		mean = div(sum, (float) (n4 * n3 * n2));

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
						normData[i4][i3][i2][i1] = (data[i4][i3][i2][i1] - mean[i1])
								/ stdDev[i1];
					}
				}
			}
		}

		/*
		 * for (int i1=0; i1<n1; ++i1){
		 * System.out.println("mean["+i1+"] = "+mean
		 * [i1]+" sum["+i1+"] = "+sum[i1]); }
		 */

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
						sum[i1] += normData[i4][i3][i2][i1];
					}
				}
			}
		}
		// Find Mean
		mean1 = div(sum1, (float) (n3 * n2));

		// Find (x1-xmean)*(x1-xmean) + (x2-xmean)*(x2-xmean) +
		// (xn-xmean)*(xn-xmean) for each attribute, the inner sum of finding
		// the standard deviation or variance.
		for (int i1 = 0; i1 < n1; ++i1) {
			for (int i4 = 0; i4 < n4; ++i4) {
				for (int i3 = 0; i3 < n3; ++i3) {
					for (int i2 = 0; i2 < n2; ++i2) {
						diffMeanSquarSum[i1] += (normData[i4][i3][i2][i1] - mean[i1])
								* (normData[i4][i3][i2][i1] - mean[i1]);
					}
				}
			}
		}

		// Find Standard Deviation
		stdDev1 = sqrt(div(diffMeanSquarSum1, n3 * n2 - 1));

		for (int i1 = 0; i1 < n1; ++i1) {
			System.out.println("mean[" + i1 + "] = " + mean1[i1] + " stddev["
					+ i1 + "] = " + stdDev1[i1]);
		}

		return normData;

	}

	private int i1(int node) {
		return node % n1som;
	}

	private int i2(int node) {
		return node / n1som;
	}

	public int node(int i1) {
		return i1;
	}

	private static float[] randomData() {
		Random r = new Random();
		int n1 = trainData2D[0].length;
		int n2 = trainData2D.length;
		int i1 = r.nextInt(n1);
		int i2 = r.nextInt(n2);
		return trainData2D[i2][i1];
	}

	public void train2D(float[][][] trainData) {
		this.trainData2D = trainData;//norm2DData(trainData);
		iterationManager();
		System.out.println("rStart = " + rStart);
	}

	public void train3D(float[][][][] trainData) {
		this.trainData3D = norm3DData(trainData);
		iterationManager();
		System.out.println("rStart = " + rStart);
	}

	public void iterationManager() {
		for (int t = 0; t < iter; ++t) {
			/*********** Start Plotting Purposes **************/
			if (t == 0 || t == 19 || t == 99 || t == 999 || t == 4999
					|| t == iter - 1) {
				float[] weightr1 = new float[n1som];
				float[] weightr2 = new float[n1som];
				for (int i = 0; i < n1som; ++i) {
					weightr1[i] = nodes[i][0];
					weightr2[i] = nodes[i][1];
				}
				String iterations = Integer.toString(t + 1);
				PointsView pointsr = new PointsView(weightr1, weightr2);

				pointsr.setMarkColor(java.awt.Color.RED);
				pointsr.setMarkStyle(PointsView.Mark.CROSS);
				pointsr.setMarkSize(4);

				PlotPanel pp = new PlotPanel();
				pp.addTiledView(pointsr);
				float max = max(trainData2D) * .75f;
				pp.setLimits(-max, -max, max, max);
				// pp.setLimits(, 0, 7, 7); // change to 5,5 for seismic data

				pp.addTitle("Iterations = " + iterations);
				pp.setHLabel("Weight 1");
				pp.setVLabel("Weight 2");
				// pp.setLimits(0, 0, 1, 1);
				PlotFrame pf = new PlotFrame(pp);
				pf.setVisible(true);
			}
			/*********** End Plotting Purposes **************/

			pickedData = randomData();
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
			updateNodes(ic);

			// System.out.println("alpha= "+alpha);

		}
	}

	public void updateNodes(int ic) {
		float h = 0;
		for (int i = 0; i < n1som; ++i) {
			updateNode(ic, i);
		}
	}

	public void updateNode(int ic, int ii) {
		float h = 0;
		float distance = distanceTable[ic][ii];
		if (distance == 0) {
			for (int i = 0; i < n_att; ++i) {
				// h = (float)
				// (alpha*Math.exp(-Math.pow(distance(ic,jc,ii,jj),2)/Math.pow(alpha,2)));
				h = (float) (alpha * exp(-pow(distanceTable[ic][ii], 2)
						/ (2. * pow(sigma, 2))));

				// System.out.println("h= "+h);
				nodes[ic][i] = nodes[ic][i] + h
						* (pickedData[i] - nodes[ic][i]);

			}
		} else if (distance <= radius) {
			for (int i = 0; i < n_att; ++i) {
				// h = (float)
				// (alpha*Math.exp(-Math.pow(distance(ic,jc,ii,jj),2)/Math.pow(alpha,2)));
				h = (float) (alpha * exp(-pow(distanceTable[ic][ii], 2)
						/ (2. * pow(sigma, 2))));
				// System.out.println("h= "+distance(ic,jc,ii,jj));
				nodes[ii][i] = nodes[ii][i] + h
						* (pickedData[i] - nodes[ii][i]);

			}
		}
	}

	public void train(float[][][][] trainData) {
		// Iteration Manager
		// Update
	}

	/*
	 * public float[][] classifyAllData(float[][][] inputData){ int n1 =
	 * inputData.length; int n2 = inputData[0].length; float[][] classifiedData
	 * = new float[n1][n2];
	 * 
	 * int winner=0;//initialize winner node for (int i=0; i<n1; ++i){ for (int
	 * j=0; j<n2; ++j){ classifiedData[i][j] = classifyOneData(inputData[i][j]);
	 * } }
	 * 
	 * return classifiedData; }
	 */

	public float[][] classify2DData(float[][][] inputData) {
		int n1 = inputData.length;
		int n2 = inputData[0].length;
		float[][][] normData = norm2DData(inputData);
		float[][] classifiedData = new float[n1][n2];
		for (int i = 0; i < n1; ++i) {
			for (int j = 0; j < n2; ++j) {
				classifiedData[i][j] = classifyOneData(normData[i][j]);

			}
		}

		return classifiedData;
	}

	public float[][] classify3DData(float[][][] inputData) {
		int n2 = inputData.length;
		int n1 = inputData[0].length;
		float[][][] normData = norm2DData(inputData);
		float[][] classifiedData = new float[n2][n1];
		int winner = 0;// initialize winner node
		for (int i = 0; i < n2; ++i) {
			for (int j = 0; j < n1; ++j) {
				classifiedData[i][j] = classifyOneData(normData[i][j]);

			}
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
			eucDis = euclideanDistance(inputData, nodes[i]);
			if (eucDis < minEucDis) {
				minEucDis = eucDis;
				winCatNum = node(i);
			}
		}
		return winCatNum;
	}

	private void fillDistanceTable() {
		for (int i1 = 0; i1 < n1som; ++i1)
			for (int i2 = 0; i2 < n1som; ++i2) {
				distanceTable[i1][i2] = (float) (sqrt(pow((i1 - i2), 2)));
			}
	}

	private void findWinner(float[] singleTrainData) {
		float xmin = 99f;

		int n2 = nodes.length;
		for (int i = 0; i < n2; ++i) {
			float eucDis = euclideanDistance(singleTrainData, nodes[i]);
			if (eucDis < xmin) {
				xmin = eucDis;
				ic = i;
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
		nodes = new float[n1som][n_att];
		randomizeAllNodesWeights(nodes);
	}

	public float[][] getNodes() {
		return nodes;
	}

	public ColorMap getColorMap(int min, int max) {
		Color[] color = new Color[256];
		float[][] hsbOrig = new float[3][n1som];
		float[][] hsbInterp = new float[3][256];
		int count = 0;
		float lengthRatio = 256.0f / (float) (n1som);
		float brightness = 8.1f;// scale from 0 through 10, 10 being the
								// brightest
		float saturationAdd = 1f;
		float hueAdd = .3f;
		float saturationStart = 4;

		if (n1som < 256) {
			int ratio1 = (int) ((256.0f / (float) n1som) + .5f);
			// System.out.println("*********:"+ratio2);
			float countC2 = 0.00001f;// small float in the case of divide by
										// zero.
			float countC1 = saturationStart;

			for (int i = 0; i < 256; ++i) {
				hsbInterp[0][i] = countC2;
				if (i % ratio1 == 0 && i > 1) {
					countC2 = countC2 + 1.0f / n1som;
				}
				System.out.println(hsbInterp[0][1]);
				hsbInterp[1][i] = countC1;
			}

		} else {
			System.out.println("****************************");
			for (int i2 = 0; i2 < 1; ++i2) {
				for (int i1 = 0; i1 < n1som; ++i1) {
					hsbOrig[0][count] = (float) i2;
					hsbOrig[1][count] = (float) i1 + saturationAdd;
					// hsbOrig[2][count] = (float) 40;

					// color[(int) (count*lengthRatio)] = getColor(i2, i1, 40);
					++count;
				}
			}

			// extend each array by 50 to get ready for sinc interpolation

			SincInterp si = new SincInterp();
			si.interpolate(n1som, lengthRatio, 0.0, hsbOrig[0], 256, 1.0, 0.0,
					hsbInterp[0]);
			si.interpolate(n1som, lengthRatio, 0.0, hsbOrig[1], 256, 1.0, 0.0,
					hsbInterp[1]);
		}

		maxH = max(hsbInterp[0]);
		maxS = max(hsbInterp[1]) + saturationAdd;

		for (int i = 0; i < 256; ++i) {

			System.out.println("h = " + hsbInterp[0][i] + " s = "
					+ hsbInterp[1][i] + " b = " + brightness);
			color[i] = getColor(hsbInterp[0][i], hsbInterp[1][i], brightness);
		}

		ColorMap cm = new ColorMap(min, max, color);

		return cm;
	}

	/**
	 * Converts a node's coordinates to a specific color using the nodes i2 and
	 * i1 location in the nodes' array. The k value will be set at a constant
	 * value. Process: Rgb coordinates are found from i2, i1, and k. The rgb
	 * coordinates are converted to hue, saturation, and brightness The hue,
	 * saturation, and brightness are used to create a color
	 * 
	 * @param i2
	 * @param i1
	 * @param k
	 * @return
	 */
	public Color getColor(float i2, float i1, float k) {
		// System.out.println("h = "+(i2)/maxH+" s = "+(i1+1)/(maxS)+" b = "+k/10.0f);
		return Color.getHSBColor((i2) / (maxH + .5f), (i1) / (maxS), k / 10.0f);// +1
																				// is
																				// by
																				// experimentation
																				// to
																				// not
																				// get
																				// overlap
																				// in
																				// colors

	}
}
