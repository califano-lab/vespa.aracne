package aracne;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

import org.apache.commons.math3.distribution.PoissonDistribution;

import java.util.Arrays;
import java.util.Comparator;
import java.util.stream.IntStream;

import org.apache.commons.lang3.ArrayUtils;

public class BootstrapConsolidator {

	private Double[] poissonLibrary = new Double[100];
	
	private HashMap<String, Integer> edgesOccurrences = new HashMap<String, Integer>();
	private HashMap<String, Double> mi = new HashMap<String, Double>();
	private HashMap<String, Double> pvalue = new HashMap<String, Double>();
	private HashMap<String, Double> pvalue_adjusted = new HashMap<String, Double>();
	private HashMap<String, Double> correlation = new HashMap<String, Double>();
	private HashMap<String, Double> prior = new HashMap<String, Double>();

	private HashSet<String> tfs = new HashSet<String>();
	private HashSet<String> targets = new HashSet<String>();
	
	private int maxCount = 0;

	private String multipletesting;

	
	public BootstrapConsolidator(String multipletesting) {
		this.multipletesting=multipletesting;
	}

	public void mergeFiles(File folder) throws IOException{
		FilenameFilter filter = new FilenameFilter() {
			@Override
			public boolean accept(File dir, String name) {
				return name.matches("bootstrapNetwork_.*");
			}
		};
		File[] listOfFiles = folder.listFiles(filter);
	
		// Calculate the occurrence of each edge (count object)
		// And sum the global mi for each edge
		System.out.println("Integrating "+listOfFiles.length+" bootstraps...");
		for (File file : listOfFiles) {
			if (file.isFile()) {
				System.out.println(file.getName());
				try {
					BufferedReader br = new BufferedReader(new FileReader(file));
					String line = "";
					boolean firstline = true;
					while((line = br.readLine()) != null){
						if(firstline){
							firstline=false;
						} else {
							String[] sp = line.split("\t");
							String key = sp[0]+"#"+sp[1]; // edge key
							tfs.add(sp[0]);
							targets.add(sp[1]);
							if(edgesOccurrences.containsKey(key)){
								edgesOccurrences.put(key, edgesOccurrences.get(key)+1);
								mi.put(key, mi.get(key)+Double.parseDouble(sp[2]));
								correlation.put(key, correlation.get(key)+Double.parseDouble(sp[3]));
								prior.put(key, prior.get(key)+Double.parseDouble(sp[4]));
							}
							else{
								edgesOccurrences.put(key, 1);
								mi.put(key, Double.parseDouble(sp[2]));
								correlation.put(key, Double.parseDouble(sp[3]));
								prior.put(key, Double.parseDouble(sp[4]));
							}
						}
					}
					br.close();
				}
				catch(Exception e){
					e.printStackTrace();
				}
			}
		}
		
		// Initialize the poisson distribution based on the average number of edge appearance
		long allCount = 0;
		//FileWriter fw = new FileWriter(new File("file.txt"));
		for(String key : edgesOccurrences.keySet()){
			// Calculate the average mi of each edge
			int edgeOccurrence=edgesOccurrences.get(key);
			mi.put(key, mi.get(key)/edgeOccurrence);
			correlation.put(key, correlation.get(key)/edgeOccurrence);
			prior.put(key, prior.get(key)/edgeOccurrence);
			allCount += edgeOccurrence; // Total number of edges in all the bootstrap files
			if(edgeOccurrence > maxCount){
				maxCount = edgeOccurrence; // Edge with the most count
			}
			//fw.write(edgeOccurrence+"\n");
		}
		//fw.close();
		// Lambda
		//long allObservedEdges = edgesOccurrences.size();
		//Double meanEdgeNumber = allCount/allObservedEdges; // Wrong way: observed edges in the denominator
		Double meanEdgeNumber = (allCount*1.0)/(targets.size()*tfs.size()); // Right way: all possible edges
		generatePoissonPvalues(meanEdgeNumber);

		// Compute p-values
		for(String key : edgesOccurrences.keySet()){
			int occurrence = edgesOccurrences.get(key);
			Double thisPvalue;
			thisPvalue = poissonLibrary[occurrence];
			if(thisPvalue>1){
				thisPvalue=1.0;
			}
			pvalue.put(key, thisPvalue);
		}

		// Multiple testing correction
		String [] tmp_pvalues_keys = pvalue.keySet().toArray(new String[0]);
		Double [] tmp_pvalues = pvalue.values().toArray(new Double[0]);
		Double [] tmp_pvalues_adjusted = new Double[tmp_pvalues.length];

		if (multipletesting.equals("BH")) {
			tmp_pvalues_adjusted = pAdjustBH(tmp_pvalues);
		}
		else if (multipletesting.equals("none")) {
			tmp_pvalues_adjusted = tmp_pvalues;
		}
		else if (multipletesting.equals("Bonferroni")) {
			// Bonferroni correction
			// In the C++ version, the correction was done using the observed number of edges
			// Here, the correction is done based on ALL the possible combinations
			int totalCount = tfs.size()*targets.size();
			
			for(int i = 0; i < tmp_pvalues.length; i++)
			{
				tmp_pvalues_adjusted[i] = tmp_pvalues[i]*totalCount;
				if(tmp_pvalues_adjusted[i]>1){
					tmp_pvalues_adjusted[i]=1.0;
				}
			}
		}

		for(int j = 0; j < tmp_pvalues.length; j++)
		{
			pvalue_adjusted.put(tmp_pvalues_keys[j],tmp_pvalues_adjusted[j]);
		}

	}
	
	public void writeSignificant(String outFile, Double poissonPvalue){
		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outFile)));

			// Header
			bw.write("Regulator\tTarget\tMI\tCorrelation\tPrior\tPvalue\n");

			String[] keys = edgesOccurrences.keySet().toArray(new String[0]);
			for(String key : keys){
				int occurrence = edgesOccurrences.get(key);

				if(pvalue_adjusted.get(key) < poissonPvalue){
					String[] sp = key.split("#");
					bw.write(sp[0]+"\t"+sp[1]+"\t"+mi.get(key)+"\t"+correlation.get(key)+"\t"+prior.get(key)+"\t"+pvalue_adjusted.get(key)+"\n");
				}
			}
			bw.close();
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public void generatePoissonPvalues(Double mean){
		poissonLibrary = new Double[maxCount+1];
		PoissonDistribution pdist = new PoissonDistribution(mean);
		for(int i=0; i<maxCount+1; i++){
			poissonLibrary[i] = 1.0-pdist.cumulativeProbability(i);
		}
	}

	private static Double[] pAdjustBH(Double[] pvalues) {
		// Generate decreasing index
		int[] ix = IntStream.rangeClosed(1, pvalues.length).toArray();
		ArrayUtils.reverse(ix);

		// Decreasing order
		int[] dorder = IntStream.range(0, pvalues.length)
		.boxed().sorted((i, j) -> Double.compare(pvalues[j], pvalues[i]))
		.mapToInt(i -> i).toArray();

		// Increasing order
		int[] iorder = IntStream.range(0, dorder.length)
		.boxed().sorted(Comparator.comparingInt(i -> dorder[i]))
		.mapToInt(i -> i).toArray();

		// Main equation
		Double[] v = new Double[pvalues.length];
		for (int i = 0; i < pvalues.length; ++i) {
			v[i] = ((double) pvalues.length / ix[i]) * pvalues[dorder[i]];
		}

		// Compute cumulative minimum
		Double[] vcm = new Double[v.length];
		Double min = v[0];

		for (int i = 0; i < v.length; ++i) {
			if (v[i] < min) {
				min = v[i];
			}
			vcm[i] = min;
		}

		// Threshold to maximum values of 1
		for (int i = 0; i < vcm.length; ++i) {
			if (vcm[i] > 1) {
				vcm[i] = 1.0;
			}
		}

		// Restore order
		Double[] pvalues_adj = new Double[pvalues.length];
		for (int i = 0; i < pvalues.length; ++i) {
			pvalues_adj[i] = vcm[iorder[i]];
		}

		return pvalues_adj;
	}
	
}








