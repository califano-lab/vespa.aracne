package aracne;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;

/**
 * Computes mutual information between regulators and targets by hybrid adaptive partitioning.
 * Requires ranked data and an array of regulators (optional: array of activators),
 * MI threshold and number of threads to use.
 *
 * @param  rankData Bootstrapped HashMap linking gene identifiers with ranks
 * @param  rankDataCor HashMap linking gene identifiers with ranks
 * @param  regulators Array of regulators (e.g. transcription factors or kinases & phosphatases)
 * @param  activators (Optional) array of activators (e.g. kinases)
 * @param  targets (Optional) array of targets (e.g. genes)
 * @param  interactionSet (Optional) HashMap of interactions
 * @param  miThreshold MI threshold to use to restrict networks
 * @param  correlationThreshold Correlation threshold to use to trust mode of interaction
 * @param  threadCount Number of threads to use
 * @param  deplete Deplete missing values (no hpMI)
 */
public class MI {
	// Variables
	private String[] genes;
	private String[] regulators;
	private String[] activators;
	private String[] targets;
	private HashMap<String, HashMap<String, Double>> interactionSet;

	// Computed MI threshold
	private double miThreshold;
	// Defined MI threshold
	private double correlationThreshold;
	// Computed MI between regulators and targets
	private HashMap<String, HashMap<String, Double>> finalNetwork;
	// Prior interaction confidence between regulators and targets
	private HashMap<String, HashMap<String, Double>> finalNetworkPrior;
	// Computed sign of correlation between regulators and targets (activation: true, deactivation: false)
	private HashMap<String, HashMap<String, Boolean>> finalNetworkSign;
	// Computed correlation between regulators and targets
	private HashMap<String, HashMap<String, Double>> finalNetworkCorrelation;
	// Computed coverage between regulators and targets
	private HashMap<String, HashMap<String, Double>> finalNetworkCoverage;
	private boolean deplete;

	/**
	 * Constructor for standard ARACNe mode (multi-threaded)
	 *
	 * @param  rankData Bootstrapped HashMap linking gene identifiers with ranks
	 * @param  rankDataCor HashMap linking gene identifiers with ranks
	 * @param  regulators Array of regulators (e.g. transcription factors or kinases & phosphatases)
	 * @param  activators (Optional) array of activators (e.g. kinases)
	 * @param  targets (Optional) array of targets (e.g. genes)
	 * @param  interactionSet (Optional) HashMap of interactions
	 * @param  miThreshold MI threshold to use to restrict networks
	 * @param  correlationThreshold Correlation threshold to use to trust mode of interaction
	 * @param  threadCount Number of threads to use
	 * @param  deplete Deplete missing values (no hpMI)
	 */
	public MI(
			HashMap<String, short[]> rankData,
			HashMap<String, short[]> rankDataCor,
			String[] regulators,
			String[] activators,
			String[] targets,
			HashMap<String, HashMap<String, Double>> interactionSet,
			Double miThreshold,
			Double correlationThreshold,
			Integer threadCount,
			boolean deplete
			) {
		// Set genes
		// this.genes = rankData.keySet().toArray(new String[0]);
		this.genes = targets;
		Arrays.sort(genes);

		// Loop to check which edges are kept. It will generate the finalNetwork and finalNetworkSign HashMap
		finalNetwork = new HashMap<String, HashMap<String, Double>>();
		finalNetworkPrior = new HashMap<String, HashMap<String, Double>>();
		finalNetworkSign = new HashMap<String, HashMap<String, Boolean>>();
		finalNetworkCorrelation = new HashMap<String, HashMap<String, Double>>();
		finalNetworkCoverage = new HashMap<String, HashMap<String, Double>>();
		for(int i=0; i<regulators.length; i++){
			HashMap<String, Double> tt = new HashMap<String, Double>();
			HashMap<String, Double> ttPrior = new HashMap<String, Double>();
			HashMap<String, Boolean> ttSign = new HashMap<String, Boolean>();
			HashMap<String, Double> ttCorrelation = new HashMap<String, Double>();
			HashMap<String, Double> ttCoverage = new HashMap<String, Double>();
			finalNetwork.put(regulators[i], tt);
			finalNetworkPrior.put(regulators[i], ttPrior);
			finalNetworkSign.put(regulators[i], ttSign);
			finalNetworkCorrelation.put(regulators[i], ttCorrelation);
			finalNetworkCoverage.put(regulators[i], ttCoverage);
		}

		// Parallelized MI computation
		ExecutorService executor = Executors.newFixedThreadPool(threadCount);
		for(int i=0; i<regulators.length; i++){
			MIThread mt = new MIThread(rankData, rankDataCor, genes, regulators, activators, interactionSet, miThreshold, correlationThreshold, i, deplete);
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
	 * @param  rankData Bootstrapped HashMap linking gene identifiers with ranks
	 * @param  rankDataCor HashMap linking gene identifiers with ranks
	 * @param  genes Array of genes covered by rankData
	 * @param  regulators Array of regulators (e.g. transcription factors or kinases & phosphatases)
	 * @param  activators (Optional) array of activators (e.g. kinases)
	 * @param  interactionSet (Optional) HashMap of interactions
	 * @param  miThreshold MI threshold to use to restrict networks
	 * @param  correlationThreshold Correlation threshold to use to trust mode of interaction
	 * @param  regulatorIndex Index of processed regulator
	 * @param  deplete Deplete missing values (no hpMI)
	 */
	class MIThread extends Thread {
		// Variables
		private HashMap<String, short[]> rankData;
		private HashMap<String, short[]> rankDataCor;
		private String[] genes;
		private String[] regulators;
		private String[] activators;
		private HashMap<String, HashMap<String, Double>> interactionSet;
		private double miThreshold;
		private double correlationThreshold;
		private int regulatorIndex;
		private boolean deplete;

		// Constructor
		MIThread(
			HashMap<String, short[]> rankData,
			HashMap<String, short[]> rankDataCor,
			String[] genes,
			String[] regulators,
			String[] activators,
			HashMap<String, HashMap<String, Double>> interactionSet,
			double miThreshold,
			double correlationThreshold,
			int regulatorIndex,
			boolean deplete)
		{
			this.rankData = rankData;
			this.rankDataCor = rankDataCor;
			this.genes = genes;
			this.regulators = regulators;
			this.activators = activators;
			this.interactionSet = interactionSet;
			this.miThreshold = miThreshold;
			this.correlationThreshold = correlationThreshold;
			this.regulatorIndex = regulatorIndex;
			this.deplete = deplete;
		}

		public void run() {
			for(int j=0; j<genes.length; j++){
				if (interactionSet.containsKey(regulators[regulatorIndex])) {
					if (interactionSet.get(regulators[regulatorIndex]).containsKey(genes[j])) {
						// Regulator ranks
						short[] vectorX = rankData.get(regulators[regulatorIndex]);
						short[] vectorX_Cor = rankDataCor.get(regulators[regulatorIndex]);
						// Target ranks
						short[] vectorY = rankData.get(genes[j]);
						short[] vectorY_Cor = rankDataCor.get(genes[j]);

						// Compute MI by hybrid adaptive partitioning
						double mi = hapMI(vectorX,vectorY,deplete);

						// Only report results if MI is higher than MI threshold
						if(mi >= miThreshold){
							// Compute correlation and sign of interactor (activation or deactivation)
							ArrayList<short[]> splitQuad = splitQuadrants(vectorX_Cor,vectorY_Cor);
							short[] valuesX = splitQuad.get(0);
							short[] valuesY = splitQuad.get(1);
							double[] vX = new double[valuesX.length];		
							vX = castShort2Double(valuesX);
							double[] vY = new double[valuesY.length];
							vY = castShort2Double(valuesY);
							double coverage = ((double)valuesX.length) / vectorX_Cor.length;
							
							// We need at least two data points to assess correlation. If not, we drop the interaction.
							double correlation = Double.NaN;
							if (vX.length > 1) {
								correlation = new SpearmansCorrelation().correlation(vX,vY);
							}
							
							if (!Double.isNaN(correlation)){
								// Compute sign from correlation
								boolean sign =( ((int)Math.signum(correlation)) == 1);
								// If activators are specified, ensure that computed mode is correct
								if(activators!=null){
									if(Arrays.asList(activators).contains( regulators[regulatorIndex]) & (sign | Math.abs(correlation) < correlationThreshold)){
											setMI(regulators[regulatorIndex], genes[j], mi);
											setSign(regulators[regulatorIndex], genes[j], true);
											setCorrelation(regulators[regulatorIndex], genes[j], (1.0*Math.abs(correlation)));
											setCoverage(regulators[regulatorIndex], genes[j], coverage);
									}else if(!Arrays.asList(activators).contains( regulators[regulatorIndex]) & (!sign | Math.abs(correlation) < correlationThreshold)){
											setMI(regulators[regulatorIndex], genes[j], mi);
											setSign(regulators[regulatorIndex], genes[j], false);
											setCorrelation(regulators[regulatorIndex], genes[j], (-1.0*(Math.abs(correlation))));
											setCoverage(regulators[regulatorIndex], genes[j], coverage);
									}
								// If activators are not specified, just report all results
								}else{
									setMI(regulators[regulatorIndex], genes[j], mi);
									setSign(regulators[regulatorIndex], genes[j], sign);
									setCorrelation(regulators[regulatorIndex], genes[j], correlation);
									setCoverage(regulators[regulatorIndex], genes[j], coverage);
								}
							}
							setPrior(regulators[regulatorIndex], genes[j], interactionSet.get(regulators[regulatorIndex]).get(genes[j]));
						}
					}
				}
			}
		}
		private synchronized void setPrior(String _tf, String _gene, double _prior){
			finalNetworkPrior.get(_tf).put(_gene, _prior);
		}
		private synchronized void setCorrelation(String _tf, String _gene, Double _correlation){
			finalNetworkCorrelation.get(_tf).put(_gene, _correlation);
		}
		private synchronized void setCoverage(String _tf, String _gene, Double _coverage){
			finalNetworkCoverage.get(_tf).put(_gene, _coverage);
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
	 * @param  rankData HashMap linking gene identifiers with ranks
	 * @param  regulators Array of regulators (e.g. transcription factors or kinases & phosphatases)
	 * @param  targets (Optional) array of targets (e.g. genes)
	 * @param  interactionSet (Optional) HashMap of interactions
	 * @param  miPvalue miPvalue for thresholding
	 * @param  maximumInteractions Maximum number of interactions to assess
	 * @param  seed Seed to use for reproducible results
	 * @param  deplete Deplete missing values (no hpMI)
	 */
	public MI(
			HashMap<String, short[]> rankData, 
			String[] regulators,
			String[] targets,
			HashMap<String, HashMap<String, Double>> interactionSet,
			double miPvalue,
			int maximumInteractions,
			int seed,
			boolean deplete
			) {
		// Set genes
		this.genes = rankData.keySet().toArray(new String[0]);
		Arrays.sort(genes);

		// Set regulators
		this.regulators = regulators;

		// Set targets
		this.targets = targets;

		// Set interactions
		this.interactionSet = interactionSet;

		// Compute MI threshold
		this.miThreshold = calibrateMIThreshold(rankData,miPvalue,maximumInteractions,seed,deplete);
	}

	/**
	 * Calibrate MI threshold using permutated matrix considering ranks and NAs.
	 * Regulator - Target relationships are conserved.
	 *
	 * @param  rankData HashMap linking gene identifiers with ranks
	 * @param  miPvalue miPvalue for thresholding
	 * @param  maximumInteractions Maximum number of interactions to assess
	 * @param  seed Seed to use for reproducible results
	 * @param  deplete Deplete missing values (no hpMI)
	 */
	public double calibrateMIThreshold(HashMap<String, short[]> rankData, double miPvalue, int maximumInteractions, int seed, boolean deplete){
		int numberOfSamples = rankData.get(genes[0]).length;

		// Compute number of interactions
		int interactionslength = 0;
		for (String regulator : interactionSet.keySet()) {
			for (String target : interactionSet.get(regulator).keySet()) {
				interactionslength += 1;
			}
		}

		HashMap<String, short[]> tempData = new HashMap<String, short[]>();
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

				short temp = tempData.get(genes[i])[r1];

				// Flip data points r1 with r2
				tempData.get(genes[i])[r1] = tempData.get(genes[i])[r2];

				tempData.get(genes[i])[r2] = temp;
			}
		}

		// Estimate MI between all regulators and targets in subset

		if (interactionslength > maximumInteractions) {
			System.out.println("Reducing "+interactionslength+" interactions to be below threshold of "+maximumInteractions+" maximum interactions.");
		}

		int interactionCounter = 0;
		for(int i=0; i<regulators.length; i++){
			for(int j=0; j<targets.length; j++){
				if (interactionCounter < maximumInteractions) {
					if (interactionSet.containsKey(regulators[i])) {
						if (interactionSet.get(regulators[i]).containsKey(targets[j])) {
							mit.add(hapMI(tempData.get(regulators[i]), tempData.get(targets[j]), deplete));
						}
					}
				}
				interactionCounter += 1;
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
	 * Returns quadrants for MI estimation of x against y
	 *
	 * @param  x short[] ranks with 0 designated as NaN
	 * @param  y short[] ranks with 0 designated as NaN
	 * @return Arraylist with the following vectors:
	 * 0-1: Intersecting reranked ranks
	 * 2: Nr. of samples with NAs in both vectors
	 * 3: Nr. of samples with NAs in vector of object
	 * 4: Nr. of samples with NAs in vector x
	 */
	public static ArrayList<short[]> splitQuadrants(short[] x, short[] y){
		ArrayList<short[]> output = new ArrayList<short[]>();

		// Quadrant counts
		short[] bothNA = new short[1];
		short[] xNA = new short[1];
		short[] yNA = new short[1];
		bothNA[0] = xNA[0] = yNA[0] = 0;
		
		// First loop does the counting
		int lengthOfNotNAsInBoth = 0;
		for(int i=0; i<x.length; i++){
			if(x[i]!=0){
				if (y[i]!=0){
					lengthOfNotNAsInBoth++;
				} else {
					yNA[0]++;
				}
			}  else {
				if (y[i]!=0){
					xNA[0]++;
				} else {
					bothNA[0]++;
				}
			}
		}
	
		// Second loop generates the not NAs vectors
		short [] xNotNA = new short[lengthOfNotNAsInBoth];
		short [] yNotNA = new short[lengthOfNotNAsInBoth];
		int j = 0;
		for(int i=0; i<x.length; i++){
			if(x[i]!=0 & y[i]!=0){
				xNotNA[j] = x[i];
				yNotNA[j] = y[i];
				j++;		
			}
		}

		// Rerank if necessary
		if (lengthOfNotNAsInBoth != x.length){
			xNotNA = reRankVector(xNotNA);
			yNotNA = reRankVector(yNotNA);
		}

		output.add(xNotNA);
		output.add(yNotNA);
		output.add(bothNA);
		output.add(xNA);
		output.add(yNA);

		return output;
	}

	public static short[] reRankVector(short[] inputVector){
		short[] rankVector = new short[inputVector.length];
		for(int i=0; i<inputVector.length; i++){
			int counter = 1;
			for(int j=0; j<inputVector.length; j++){
				if(inputVector[i] > inputVector[j]){
					counter++;
				}
			}
			rankVector[i] = (short)counter;
		}
		return rankVector;
	}	

	/**
	 * Estimate mutual information between two vectors using hybrid adaptive partitioning.
	 *
	 * Step 1: Split ranks X and Y into 4 quadrants
	 * 0-1: Ranks of numerical intersection of X and Y
	 * 2: X and Y both are NA
	 * 3. Only X is NA
	 * 4. Only Y is NA
	 * Step 2: Compute MI for each quadrant separately
	 * Step 3: Summarize MI
	 *
	 * @param  vectorX short[] of regulator
	 * @param  vectorY short[] of target
	 * @param  deplete Deplete missing values (no hpMI)
	 */
	 public static double hapMI(short[] vectorX, short[] vectorY, boolean deplete){
		// Preranked matrices for candidate interactors might not overlap perfectly, thus requires reranking
		ArrayList<short[]> splitQuad = splitQuadrants(vectorX, vectorY);
		short[] valuesX = splitQuad.get(0);
		short[] valuesY = splitQuad.get(1);

		// Get the counts for NA data points
		int bothNACount = 0;
		int na1 = 0;
		int na2 = 0;

		if (!deplete) {
			bothNACount = (int)splitQuad.get(2)[0];
			na1 = (int)splitQuad.get(3)[0];
			na2 = (int)splitQuad.get(4)[0];
		}

		boolean firstLoop = true;

		double mi = 0;

 		// Perform first split by using the counts and the recursive computeMI call for sample pairs that have both values 
		if(valuesX.length < 8){
			mi = getInformation(valuesX.length, (valuesX.length+na1), (valuesX.length+na2));
		}
		else{
			mi = computeMI(valuesX,valuesY,(short)0,(short)(valuesX.length),(short)0,(short)(valuesY.length),firstLoop);
		}

		int numberOfSamples = vectorX.length;
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

	public HashMap<String, HashMap<String, Double>> getFinalNetworkPrior() {
		return finalNetworkPrior;
	}

	public HashMap<String, HashMap<String, Boolean>> getFinalNetworkSign() {
		return finalNetworkSign;
	}

	public HashMap<String, HashMap<String, Double>> getFinalNetworkCorrelation() {
		return finalNetworkCorrelation;
	}

	public HashMap<String, HashMap<String, Double>> getFinalNetworkCoverage() {
		return finalNetworkCoverage;
	}

	public short median(short[] input){
		short median;
		Arrays.sort(input);
		median = (short) input[input.length/2];
		return(median);
	}
}

