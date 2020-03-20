package common;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.HashMap;
import java.util.Arrays;


/**
 * Class defining static methods to load specific data formats
 */
public class DataParser {
	static ArrayList<String> genes;
	private int[] bootRandomSamples = null;
	private int sampleNumber = 0;

	public String[] getGenes() {
		String[] genearray = new String[genes.size()];
		genearray = genes.toArray(genearray);
		return genearray;
	}
	/**
	 * Method to read a file with a single number
	 * @throws IOException 
	 */
	public static Double readValue(File inFile) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(inFile));
		String string = br.readLine();
		Double value = Double.valueOf(string);
		br.close();
		return(value);
	}
	/**
	 * Method to write a single value into a file
	 */
	public static void writeValue(Double value, File outFile) throws IOException{
		BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));
		bw.write(String.valueOf(value));
		bw.close();
	}

	/**
	 * Method to read a gene list
	 */
	public static String[] readGeneSet(File regulatorsFile){
		try{
			BufferedReader br = new BufferedReader(new FileReader(regulatorsFile));
			String l = "";
			HashSet<String> regulators = new HashSet<String>();
			while((l = br.readLine()) != null){
				regulators.add(l);
			}
			br.close();
			String[] regulatorsList = regulators.toArray(new String[0]);
			Arrays.sort(regulatorsList);

			return(regulatorsList);
		}
		catch(Exception e){
			e.printStackTrace();
		}
		return null;
	}

	/**
	 * Method to read a gene list and ensure that all values are covered by the expression matrix
	 */
	public static String[] readGeneSet(File regulatorsFile, ArrayList<String> genes){
		try{
			BufferedReader br = new BufferedReader(new FileReader(regulatorsFile));
			String l = "";
			HashSet<String> regulators = new HashSet<String>();
			while((l = br.readLine()) != null){
				regulators.add(l);
			}
			br.close();

			// Only keep regulators present in expression matrix
			regulators.retainAll(new HashSet<String>(genes));

			String[] regulatorsList = regulators.toArray(new String[0]);
			Arrays.sort(regulatorsList);

			if(regulatorsList.length==0){
				throw new Exception("The regulators file is not valid.");
			}

			return(regulatorsList);
		}
		catch(Exception e){
			e.printStackTrace();
		}
		return null;
	}

	/**
	 * Method to read an interaction list and ensure that all values are covered by the expression matrix
	 */
	public static HashMap<String, HashMap<String, Double>> readInteractionSet(File interactionsFile){
		try{
			BufferedReader br = new BufferedReader(new FileReader(interactionsFile));
			String l = "";
			HashMap<String, HashMap<String, Double>> interactionNetwork = new HashMap<String, HashMap<String, Double>>();
			boolean firstline = true;
			while((l = br.readLine()) != null){
				if(firstline){
					firstline=false;
				} else {
					String[] sp = l.split("\t");
						// protein A: regulator, protein B: target
						HashMap<String, Double> interactionMap = new HashMap<String, Double>();
						if (interactionNetwork.containsKey(sp[0])){
							interactionMap = interactionNetwork.get(sp[0]);
						}
						interactionMap.put(sp[1], Double.parseDouble(sp[2]));
						interactionNetwork.put(sp[0], interactionMap);

						// protein B: regulator, protein A: target
						HashMap<String, Double> interactionMap2 = new HashMap<String, Double>();
						if (interactionNetwork.containsKey(sp[1])){
							interactionMap2 = interactionNetwork.get(sp[1]);
						}
						interactionMap.put(sp[0], Double.parseDouble(sp[2]));
						interactionNetwork.put(sp[1], interactionMap);
				}
			}
			br.close();

			return(interactionNetwork);
		}
		catch(Exception e){
			e.printStackTrace();
		}
		return null;
	}
}