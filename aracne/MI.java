package aracne;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;

import common.DataVector;
import common.Methods;

/**
 * Computes mutual information between regulators and targets by hybrid adaptive partitioning.
 * Requires ranked data and an array of regulators (optional: array of activators),
 * MI threshold and number of threads to use.
 *
 * @param  rankData HashMap linking gene identifiers with ranked DataVector
 * @param  regulators Array of regulators (e.g. transcription factors or kinases & phosphatases)
 * @param  activators (Optional) array of activators (e.g. kinases)
 * @param  miThreshold MI threshold to use to restrict networks
 * @param  threadCount Number of threads to use
 */
public class MI {
	// Variables
	private String[] genes;
	private String[] regulators;
	private String[] activators;

	// Computed MI threshold
	private double miThreshold;
	// Computed MI between regulators and targets
	private HashMap<String, HashMap<String, Double>> finalNetwork;
	// Computed sign of correlation between regulators and targets (activation: true, deactivation: false)
	private HashMap<String, HashMap<String, Boolean>> finalNetworkSign;

	/**
	 * Constructor for standard ARACNe mode (multi-threaded)
	 *
	 * @param  rankData HashMap linking gene identifiers with ranked DataVector
	 * @param  regulators Array of regulators (e.g. transcription factors or kinases & phosphatases)
	 * @param  activators (Optional) array of activators (e.g. kinases)
	 * @param  miThreshold MI threshold to use to restrict networks
	 * @param  threadCount Number of threads to use
	 */
	public MI(
			HashMap<String, DataVector> rankData,
			String[] regulators,
			String[] activators,
			Double miThreshold,
			Integer threadCount
			) {
		// Set genes
		this.genes = rankData.keySet().toArray(new String[0]);
		Arrays.sort(genes);

		// Set regulators
		this.regulators = regulators;

		// Set activators
		this.activators = activators;

		// Loop to check which edges are kept. It will generate the finalNetwork and finalNetworkSign HashMap
		finalNetwork = new HashMap<String, HashMap<String, Double>>();
		finalNetworkSign = new HashMap<String, HashMap<String, Boolean>>();
		for(int i=0; i<regulators.length; i++){
			HashMap<String, Double> tt = new HashMap<String, Double>();
			HashMap<String, Boolean> ttSign = new HashMap<String, Boolean>();
			finalNetwork.put(regulators[i], tt);
			finalNetworkSign.put(regulators[i], ttSign);
		}

		// Parallelized MI computation
		ExecutorService executor = Executors.newFixedThreadPool(threadCount);
		for(int i=0; i<regulators.length; i++){
			MIThread mt = new MIThread(rankData, genes, regulators, activators, miThreshold, i);
			executor.execute(mt);
		}

		executor.shutdown();
		while (!executor.isTerminated()) {
			// Do not continue before all threads finish
		}

		System.out.println("Regulators processed: "+regulators.length);
	}

	/**
	 * Single thread of MI computation step
	 *
	 * @param  rankData HashMap linking gene identifiers with ranked DataVector
	 * @param  genes Array of genes covered by rankData
	 * @param  regulators Array of regulators (e.g. transcription factors or kinases & phosphatases)
	 * @param  activators (Optional) array of activators (e.g. kinases)
	 * @param  miThreshold MI threshold to use to restrict networks
	 * @param  regulatorIndex Index of processed regulator
	 */
	class MIThread extends Thread {
		// Variables
		private HashMap<String, DataVector> rankData;
		private String[] genes;
		private String[] regulators;
		private String[] activators;
		private double miThreshold;
		private int regulatorIndex;

		// Constructor
		MIThread(
			HashMap<String, DataVector> rankData,
			String[] genes,
			String[] regulators,
			String[] activators,
			double miThreshold,
			int regulatorIndex)
		{
		}

		public void run() {
			for(int j=0; j<genes.length; j++){
				if(!genes[j].equals(regulators[regulatorIndex])){
					// Regulator DataVector
					DataVector vectorX = rankData.get(regulators[regulatorIndex]);
					// Target DataVector
					DataVector vectorY = rankData.get(genes[j]);

					// Compute MI by hybrid adaptive partitioning
					double mi = hapMI(vectorX,vectorY);

					// Only report results if MI is higher than MI threshold
					if(mi >= miThreshold){
						// Compute correlation and sign of interactor (activation or deactivation)
						ArrayList<short[]> splitQuad = vectorX.getQuadrants(vectorY);
						short[] valuesX = splitQuad.get(0);
						short[] valuesY = splitQuad.get(1);
						double[] vX = new double[valuesX.length];		
						vX = castShort2Double(valuesX);
						double[] vY = new double[valuesY.length];
						vY = castShort2Double(valuesY);
						
						// We need at least two data points to assess correlation. If not, we drop the interaction.
						double correlation = 0.0;
						if (vX.length > 1) {
							correlation = new SpearmansCorrelation().correlation(vX,vY);
						}
						
						if (correlation!=0.0){
							// Compute sign from correlation
							boolean sign =( ((int)Math.signum(correlation)) == 1);
							// If activators are specified, ensure that computed mode is correct
							if(activators!=null){
								if(Arrays.asList(activators).contains( regulators[regulatorIndex]) & sign){
										setMI(regulators[regulatorIndex], genes[j], mi);
										setSign(regulators[regulatorIndex], genes[j], sign);
								}else if((!Arrays.asList(activators).contains( regulators[regulatorIndex])) & !sign){
										setMI(regulators[regulatorIndex], genes[j], mi);
										setSign(regulators[regulatorIndex], genes[j], sign);
								}
							// If activators are not specified, just report all results
							}else{
								setMI(regulators[regulatorIndex], genes[j], mi);
								setSign(regulators[regulatorIndex], genes[j], sign);
							}
						}
					}
				}
			}
		}
		private synchronized void setSign(String _tf, String _gene, boolean _sign){
			finalNetworkSign.get(_tf).put(_gene, _sign);
		}
		private synchronized void setMI(String _tf, String _gene, double _mi){
			finalNetwork.get(_tf).put(_gene, _mi);
		}
	}

	/**
	 * Constructor for ARACNe MI threshold mode
	 *
	 * @param  rankData HashMap linking gene identifiers with ranked DataVector
	 * @param  regulators Array of regulators (e.g. transcription factors or kinases & phosphatases)
	 * @param  miPvalue miPvalue for thresholding
	 * @param  seed Seed to use for reproducible results
	 */
	public MI(
			HashMap<String, DataVector> rankData, 
			String[] regulators,
			double miPvalue,
			int seed
			) {
		// Set genes
		this.genes = rankData.keySet().toArray(new String[0]);
		Arrays.sort(genes);

		// Set regulators
		this.regulators = regulators;

		// Compute MI threshold
		this.miThreshold = calibrateMIThreshold(rankData,miPvalue,seed);
	}

	/**
	 * Calibrate MI threshold using permutated matrix considering ranks and NAs.
	 * Regulator - Target relationships are conserved.
	 *
	 * @param  rankData HashMap linking gene identifiers with ranked DataVector
	 * @param  miPvalue miPvalue for thresholding
	 * @param  seed Seed to use for reproducible results
	 */
	public double calibrateMIThreshold(HashMap<String, DataVector> rankData, double miPvalue, int seed){
		int numberOfSamples = rankData.get(genes[0]).values.length;

		System.out.println("Finding threshold for "+genes.length+" genes and "+numberOfSamples+" samples.");

		HashMap<String, DataVector> tempData = new HashMap<String, DataVector>();
		Random r = new Random(seed);

		// Copy data matrix
		for(int i=0; i<genes.length; i++){
			String gene = genes[i];
			tempData.put(gene, rankData.get(gene));
		}

		ArrayList<Double> mit = new ArrayList<Double>();

		// Permutate the values of the data matrix gene-wise
		for(int i=0; i<genes.length; i++){
			for(int j=0; j<numberOfSamples*10; j++){
				int r1 = r.nextInt(numberOfSamples);
				int r2 = r.nextInt(numberOfSamples);

				short temp = tempData.get(genes[i]).values[r1];
				boolean temp_na =  tempData.get(genes[i]).NAs[r1];

				// Flip data points r1 with r2
				tempData.get(genes[i]).values[r1] = tempData.get(genes[i]).values[r2];
				tempData.get(genes[i]).NAs[r1] =  tempData.get(genes[i]).NAs[r2];

				tempData.get(genes[i]).values[r2] = temp;
				tempData.get(genes[i]).NAs[r2] = temp_na;
			}
		}

		// Estimate MI between all regulators and targets in subset
		for(int i=0; i<regulators.length; i++){
			for(int j=0; j<genes.length; j++){
				if(!genes[j].equals(regulators[i])){
					mit.add(hapMI(tempData.get(regulators[i]), tempData.get(genes[j])));
				}
			}
		}

		double[] mis = new double[mit.size()];
		for(int i=0; i<mit.size(); i++){
			mis[i] = mit.get(i);
		}
		double[] rr = fitNull(mis,20);
		double miThreshold = rr[1]*Math.log(miPvalue)+rr[0]; // This looks a lot like the config_kernel.txt parameters
		System.out.println("Parameters for fitted threshold function: "+Arrays.toString(fitNull(mis,100)));
		System.out.println("MI threshold: "+miThreshold);
		return(miThreshold);
	}

	/**
	 * Fit MI null distribution
	 *
	 * @param  _x Array of MI values
	 * @param  _tailLength Length of distribution tail
	 */
	private double[] fitNull(double[] _x, int _tailLength){
		Arrays.sort(_x);
		double[] y = new double[_x.length];
		for(int i=0; i<_x.length; i++){
			y[i] = (_x.length-i)*1.0/_x.length;
		}
		double[] tailx = new double[_tailLength];
		double[] tailMI = new double[_tailLength];
		System.arraycopy(y, y.length-_tailLength, tailx, 0, _tailLength);
		System.arraycopy(_x, _x.length-_tailLength, tailMI, 0, _tailLength);
		for(int i=0; i<tailx.length; i++){
			tailx[i] = Math.log(tailx[i]);
		}
		double[] temp = tailx;
		tailx = tailMI;
		tailMI = temp;
		double sumx = 0.0;
		double sumy = 0.0;
		for(int i=0; i<tailx.length; i++){
			sumx  += tailMI[i];
			sumy  += tailx[i];
		}
		double xbar = sumx / tailx.length;
		double ybar = sumy / tailx.length;
		double xxbar = 0.0;
		double xybar = 0.0;
		for (int i = 0; i < tailx.length; i++) {
			xxbar += (tailMI[i] - xbar) * (tailMI[i] - xbar);
			xybar += (tailMI[i] - xbar) * (y[i] - ybar);
		}
		double beta1 = xybar / xxbar;
		double beta0 = ybar - beta1 * xbar;
		double[] re = new double[2];
		re[0] = beta0;
		re[1] = beta1;
		return re;
	}

	/**
	 * Estimate mutual information between two vectors using hybrid adaptive partitioning.
	 *
	 * Step 1: Split DataVector X and Y into 4 quadrants
	 * 0-1: Ranks of numerical intersection of X and Y
	 * 2: X and Y both are NA
	 * 3. Only X is NA
	 * 4. Only Y is NA
	 * Step 2: Compute MI for each quadrant separately
	 * Step 3: Summarize MI
	 *
	 * @param  vectorX DataVector of regulator
	 * @param  vectorY DataVector of target
	 */
	 public static double hapMI(DataVector vectorX, DataVector vectorY){
		// Preranked matrices for candidate interactors might not overlap perfectly, thus requires reranking
		ArrayList<short[]> splitQuad = vectorX.getQuadrants(vectorY);
		short[] valuesX = splitQuad.get(0);
		short[] valuesY = splitQuad.get(1);

		// Get the counts for NA data points
		int bothNACount = (int)splitQuad.get(2)[0];
		int na1 = (int)splitQuad.get(3)[0];
		int na2 = (int)splitQuad.get(4)[0];

		boolean firstLoop = true;

		double mi = 0;

 		// Perform first split by using the counts and the recursive computeMI call for sample pairs that have both values 
		if(valuesX.length < 8){
			mi = getInformation(valuesX.length, (valuesX.length+na1), (valuesX.length+na2));
		}
		else{
			mi = computeMI(valuesX,valuesY,(short)0,(short)(valuesX.length),(short)0,(short)(valuesY.length),firstLoop);
		}

		int numberOfSamples = vectorX.values.length;
		mi += getInformation(na1, na1+bothNACount, na1+valuesY.length);
		mi += getInformation(na2, na2+bothNACount, na2+valuesY.length);
		mi += getInformation(bothNACount, bothNACount+na1, bothNACount+na2);

		mi = mi/numberOfSamples + Math.log(numberOfSamples); // TODO why this (proof on whiteboard picture by Yishai, late January 2015)

		return mi;
	}

	/**
	 * Estimate mutual information between two numerical vectors using adaptive partitioning.
	 * ToDo: Implement GPU parallelization support
	 *
	 * @param  vectorX Array of ranks of regulator
	 * @param  vectorY Arrau of ranks of target
	 */
	public static double computeMI(short[] vectorX, short[] vectorY){
		boolean firstLoop = true;
		return computeMI(vectorX,vectorY,(short)0,(short)vectorX.length,(short)0,(short)vectorY.length,firstLoop)
				/vectorX.length 
				+ Math.log(vectorX.length); // TODO why this (proof on whiteboard picture by Yishai, late January 2015)
	}

	private static double computeMI(
			short[] vectorX,
			short[] vectorY,
			short limX1,
			short limX2,
			short limY1,
			short limY2,
			boolean firstLoop
			){
		double mi = 0;
		short pointsInQuadrant1 = 0;
		short pointsInQuadrant2 = 0;
		short pointsInQuadrant3 = 0;
		short pointsInQuadrant4 = 0;

		// Pivots on half quadrants
		short splitX = (short)((limX2+limX1)/2);
		short splitY = (short)((limY2+limY1)/2);

		short[] libraryX1 = new short[vectorX.length];
		short[] libraryX2 = new short[vectorX.length];
		short[] libraryX3 = new short[vectorX.length];
		short[] libraryX4 = new short[vectorX.length];
		short[] libraryY1 = new short[vectorX.length];
		short[] libraryY2 = new short[vectorX.length];
		short[] libraryY3 = new short[vectorX.length];
		short[] libraryY4 = new short[vectorX.length];


		for(int i=0; i<vectorX.length; i++){
			if(vectorX[i] <= splitX){
				if(vectorY[i] <= splitY){
					libraryX1[pointsInQuadrant1] = vectorX[i];
					libraryY1[pointsInQuadrant1] = vectorY[i];
					pointsInQuadrant1++;
				}
				else{
					libraryX2[pointsInQuadrant2] = vectorX[i];
					libraryY2[pointsInQuadrant2] = vectorY[i];
					pointsInQuadrant2++;
				}
			}
			else{
				if(vectorY[i] > splitY){
					libraryX3[pointsInQuadrant3] = vectorX[i];
					libraryY3[pointsInQuadrant3] = vectorY[i];
					pointsInQuadrant3++;
				}
				else{
					libraryX4[pointsInQuadrant4] = vectorX[i];
					libraryY4[pointsInQuadrant4] = vectorY[i];
					pointsInQuadrant4++;
				}
			}
		}

		// Chi-square uniformity test
		double expected = vectorX.length/4.0; 
		double tst = (Math.pow(pointsInQuadrant1-expected,2) + Math.pow(pointsInQuadrant2-expected,2) + Math.pow(pointsInQuadrant3-expected,2) + Math.pow(pointsInQuadrant4-expected,2))/expected;

		// 7.815 is the hard-coded threshold, below which you pass the uniformity Chi-Square test with confidence >95% 

		if(tst > 7.815 || firstLoop){
			firstLoop = false;

			// TODO System.arraycopy may be faster for hundreds of samples

			// Operations to put points in different quadrants
			// TODO the splitX-X1 can be precomputed for future optimization
			if(pointsInQuadrant1 > 2){
				short[] x1 = Arrays.copyOfRange(libraryX1, 0, pointsInQuadrant1);
				short[] y1 = Arrays.copyOfRange(libraryY1, 0, pointsInQuadrant1);

				mi += computeMI(x1, y1, limX1, splitX, limY1, splitY, firstLoop);
			}
			else{
				mi += getInformation(pointsInQuadrant1, splitX-limX1, splitY-limY1);
			}

			if(pointsInQuadrant2 > 2){
				short[] x2 = Arrays.copyOfRange(libraryX2, 0, pointsInQuadrant2);
				short[] y2 = Arrays.copyOfRange(libraryY2, 0, pointsInQuadrant2);

				mi += computeMI(x2, y2, limX1, splitX, splitY, limY2, firstLoop);
			}
			else{
				mi += getInformation(pointsInQuadrant2, splitX-limX1, limY2-splitY);
			}

			if(pointsInQuadrant3 > 2){
				short[] x3 = Arrays.copyOfRange(libraryX3, 0, pointsInQuadrant3);
				short[] y3 = Arrays.copyOfRange(libraryY3, 0, pointsInQuadrant3);

				try{
					mi += computeMI(x3, y3, splitX, limX2, splitY, limY2, firstLoop);
				}
				catch(Exception e){
					System.out.println(e);
				}
			}
			else{
				mi += getInformation(pointsInQuadrant3, limX2-splitX, limY2-splitY);
			}

			if(pointsInQuadrant4 > 2){
				short[] x4 = Arrays.copyOfRange(libraryX4, 0, pointsInQuadrant4);
				short[] y4 = Arrays.copyOfRange(libraryY4, 0, pointsInQuadrant4);

				mi += computeMI(x4, y4, splitX, limX2, limY1, splitY, firstLoop);
			}
			else{
				mi += getInformation(pointsInQuadrant4, limX2-splitX, splitY-limY1);
			}
		}
		else{
			mi += getInformation(vectorX.length, limX2-limX1, limY2-limY1);
		}

		return mi;
	}

	private static double getInformation(double pXY, double pX, double pY){
		if(pXY == 0){
			return  0;
		}
		else{
			// This little formula contains the entire MI calculation. Wow
			return pXY*Math.log( pXY/(pX*pY));
		}
	}

	public double[] castShort2Double(short[] vectorX) {
		double[] transformed = new double[vectorX.length];
		short[] buffer = vectorX.clone();
		for (int j=0;j<vectorX.length;j++) {
		    transformed[j] = (double)buffer[j];
		}
		return transformed;
	}

	public Double getThreshold() {
		return miThreshold;
	}

	public HashMap<String, HashMap<String, Double>> getFinalNetwork() {
		return finalNetwork;
	}

	public HashMap<String, HashMap<String, Boolean>> getFinalNetworkSign() {
		return finalNetworkSign;
	}

	public short median(short[] input){
		short median;
		Arrays.sort(input);
		median = (short) input[input.length/2];
		return(median);
	}
}

