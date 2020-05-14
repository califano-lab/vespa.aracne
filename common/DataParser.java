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
			HashMap<String, HashMap<String, Double>> interactionSet = new HashMap<String, HashMap<String, Double>>();
			boolean firstline = true;
			int interactionCounter = 0;
			int regulatorCounter = 0;
			while((l = br.readLine()) != null){
				if(firstline){
					firstline=false;
				} else {
					String[] sp = l.split("\t");

					// skip potential self interactions
					if (!sp[0].equals(sp[1])) {
						HashMap<String, Double> interaction = new HashMap<String, Double>();
						if (interactionSet.containsKey(sp[0])){
							interaction = interactionSet.get(sp[0]);
						}
						else {
							regulatorCounter += 1;
						}
						interaction.put(sp[1], Double.parseDouble(sp[2]));
						interactionSet.put(sp[0], interaction);

						interactionCounter += 1;
					}
				}
			}
			br.close();

			System.out.println("Info: Parsed "+regulatorCounter+" regulators.");
			System.out.println("Info: Parsed "+interactionCounter+" interactions.");

			return(interactionSet);
		}
		catch(Exception e){
			e.printStackTrace();
		}
		return null;
	}
}