package common;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

import common.DataVector;

/**
 * Generates ExpressionMatrix object from expression matrix file accounting for missing values
 *
 * @param  file Input expression matrix file
 * @return ExpressionMatrix
 */
public class ExpressionMatrix {
	// Variable declaration
	double[][] data;
	boolean[][] naBoolean;
	ArrayList<String> genes;
	ArrayList<String> samples;
	private int[] bootsamples;

	// Constructor (inherited objects) with the addition of boolean information about NA status of values
	public ExpressionMatrix(double[][] data, boolean[][] naBoolean, ArrayList<String> genes, ArrayList<String> samples){
		this.data=data;
		this.naBoolean=naBoolean;
		this.genes=genes;
		this.samples=samples;
	}

	// Constructor (parses a file) with missing values
	public ExpressionMatrix(File file) throws NumberFormatException, IOException, Exception{
		BufferedReader br0 = new BufferedReader(new FileReader(file));
		String strLine0;
		// Count number of lines (features, variables) and columns (arrays, samples)
		int rows = -1; // We won't count the first line
		int columns = 0;
		while ((strLine0 = br0.readLine()) != null) {
			if (strLine0.length()>0) { // This statement will skip empty lines.
				rows++;
				if (rows == 0) {
					String[] splitter = strLine0.split("\t");
					for (int j = 1; j<splitter.length; j++) {
						if (splitter[j].trim().length() > 0) {
							columns++;
						}
					}
				}
			}
		}
		br0.close();

		// Now reloop on the file to load it
	
		BufferedReader br = new BufferedReader(new FileReader(file));
		String strLine;
		// Read File Line By Line
		data = new double[rows][columns];
		naBoolean = new boolean[rows][columns];
		genes = new ArrayList<String>();
		samples = new ArrayList<String>();

		if(samples.size()>32767){
			throw new IOException("The number of samples is higher than the maximum supported by the short data type.");
		}

		int i = 0;
		while ((strLine = br.readLine()) != null) {
			String[] splitter = strLine.split("\t");
			if (i == 0) { // Fill the column names
				for (int j = 1; j<splitter.length; j++) {
					samples.add(splitter[j]);
				}
			} else {
				String gene = null;
				for (int j = 0; j<splitter.length; j++) {
					if (j == 0) { // If this is the first column
						gene=splitter[j];
						genes.add(gene);
					} else {
						if(!splitter[j].toUpperCase().equals("NA") && !splitter[j].toUpperCase().equals("NAN") && !splitter[j].toUpperCase().equals("")){
							data[i-1][j-1] = Double.parseDouble(splitter[j]);
							naBoolean[i-1][j-1] = false;
						} else {
							data[i-1][j-1] = 0;
							naBoolean[i-1][j-1] = true;
						}
					}
				}
			}
			i++;
		}
		br.close();
	}
		
	// Rank with white noise
	public HashMap<String,short[]> rank(Random random){
		HashMap<String,short[]> rankData = new HashMap<String, short[]>();

		int i = 0;
		for(String gene : genes){
			double[] inputVector = data[i];
			boolean[] naVector = naBoolean[i];
			// Add white noise to break ties
			for(int ii = 0; ii<inputVector.length; ii++){
				 inputVector[ii] = inputVector[ii] + random.nextDouble() / 1e5;
			}
			
			short[] values = Methods.rankVector(inputVector,naVector);
			rankData.put(gene, values);
			i++;
		}

		return(rankData);
	}

	// Rank with white noise
	public HashMap<String,DataVector> rankDV(Random random){
		HashMap<String,DataVector> rankData = new HashMap<String, DataVector>();

		int i = 0;
		for(String gene : genes){
			double[] inputVector = data[i];
			boolean[] naVector = naBoolean[i];
			// Add white noise to break ties
			for(int ii = 0; ii<inputVector.length; ii++){
				 inputVector[ii] = inputVector[ii] + random.nextDouble() / 1e5;
			}
			
			short[] values = Methods.rankVector(inputVector,naVector);
			rankData.put(gene, new DataVector(values));
			i++;
		}
		return(rankData);
	}
	
	// Rank without white noise
	public HashMap<String,short[]> rank(){
		HashMap<String,short[]> rankData = new HashMap<String, short[]>();

		int i = 0;
		for(String gene : genes){
			double[] inputVector = data[i];
			boolean[] naVector = naBoolean[i];
			short[] values = Methods.rankVector(inputVector,naVector);
			rankData.put(gene, values);
			i++;
		}
		return(rankData);
	}

	// Rank without white noise
	public HashMap<String,DataVector> rankDV(){
		HashMap<String,DataVector> rankData = new HashMap<String, DataVector>();

		int i = 0;
		for(String gene : genes){
			double[] inputVector = data[i];
			boolean[] naVector = naBoolean[i];
			short[] values = Methods.rankVector(inputVector,naVector);
			rankData.put(gene, new DataVector(values));
			i++;
		}
		return(rankData);
	}

	public ExpressionMatrix bootstrap(Random random){

		// Which samples to get
		int[] bootsamples = new int[samples.size()];
		for(int i = 0; i<bootsamples.length; i++){
			bootsamples[i] = random.nextInt(samples.size()-1);
		}
		this.bootsamples=bootsamples;
		// Generate a bootstrapped data
		double[][] newdata = new double[genes.size()][samples.size()];
		boolean[][] newdataNA = new boolean[genes.size()][samples.size()];
		for(int i = 0; i<genes.size(); i++){
			for(int j = 0; j<samples.size(); j++){
				// Add white noise to break ties
				newdata[i][j] = data[i][bootsamples[j]] + random.nextDouble() / 1e5;
				newdataNA[i][j] = naBoolean[i][bootsamples[j]];
			}
		}
		ExpressionMatrix em = new ExpressionMatrix(newdata, newdataNA, genes, samples);
		return(em);
	}

	public ExpressionMatrix bootstrap(int[] bootsamples){
		Random noiseGenerator = new Random();
		double[][] newdata = new double[genes.size()][samples.size()];
		boolean[][] newdataNA = new boolean[genes.size()][samples.size()];
		for(int i = 0; i<genes.size(); i++){
			for(int j = 0; j<samples.size(); j++){
				// Add white noise to break ties
				newdata[i][j] = data [i][bootsamples[j]] + noiseGenerator.nextDouble() / 1e5;
				newdataNA[i][j] = naBoolean[i][bootsamples[j]];
			}		
		}
		ExpressionMatrix em = new ExpressionMatrix(newdata, newdataNA, genes, samples);
		return(em);
	}

	public double[][] getData() {
		return data;
	}

	public void setData(int i, int j, double value) {
		data[i][j] = value;
	}

	public boolean[][] getNA() {
		return naBoolean;
	}

	public void setNA(int i, int j, boolean value) {
		naBoolean[i][j] = value;
	}

	public ArrayList<String> getSamples() {
		return samples;
	}


	public ArrayList<String> getGenes() {
		return genes;
	}

	// This method generates a new expressionMatrix which is a subset of the current object supporting NA values
	public ExpressionMatrix selectSamples(ArrayList<String> common) {
		double[][] newdata = new double [genes.size()][];
		boolean[][] newdataNA = new boolean [genes.size()][];
		int igene = 0;
		for (String gene : genes){
			double vector[] = new double[common.size()];
			int i = 0;
			for (String sample : common){
				int j = samples.indexOf(sample);
				vector[i] = data[igene][j];
				i++;
			}
			newdata[igene]=vector;
			igene++;
		}
		ExpressionMatrix newEm = new ExpressionMatrix(newdata, newdataNA, genes, common);
		return(newEm);
	}

}
