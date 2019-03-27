package aracne;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.ParseException;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigInteger;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;

import aracne.BootstrapConsolidator;
import common.DataParser;
import common.ExpressionMatrix;

import common.DataVector;

public class Aracne {
	// Variable definition
	static NumberFormat formatter = new DecimalFormat("0.###E0");
	static float simCut = 0; // In the original paper this was a parameter. E.g. if two TFs have a very high MI, DPI is not calculated
	static Random random = new Random();
	private static boolean singlemode = false;

	// Main Method
	public static void main(String[] args) throws Exception {
		
		Logger.getGlobal().setLevel(Level.OFF);
		
		Locale.setDefault(Locale.US);

		System.out.println("ARACNe (revision: " + Aracne.class.getPackage().getImplementationVersion() + ")");

		//// Parse arguments
		CommandLineParser parser = new DefaultParser();
		Options options = new Options();

		// Flag arguments
		options.addOption("c", "consolidate", false, "");
		options.addOption("j", "calculateThreshold", false, "");
		options.addOption("n", "nodpi", false, "");
		options.addOption("b", "nobootstrap", false, "");
		options.addOption("r", "nobonferroni", false, "");

		// Arguments with values
		options.addOption("e", "expfile", true, "");
		options.addOption("o", "output", true, "");
		options.addOption("k", "kinases", true, "");
		options.addOption("t", "tfs", true, "");
		options.addOption("p", "pvalue", true, "");
		options.addOption("g", "geneNumber", true, "");
		options.addOption("s", "seed", true, "");
		options.addOption("m", "threads", true, "");
		options.addOption("v", "consolidatepvalue", true, "");

		// Default arguments
		boolean isConsolidate = false;
		boolean isThreshold = false;
		boolean noDPI = false;
		boolean nobootstrap = false;
		boolean nobonferroni = false;

		String expPath = null;
		String outputPath = null;
		String tfsPath = null;
		String kinasesPath = null;
		Double miPvalue = 1E-8;
		Integer geneNumber = 3000;
		Integer seed = null;
		Integer threadCount = 1;
		Double consolidatePvalue = 0.05; 

		// Parse arguments
		try {
			CommandLine cmd = parser.parse( options, args);

			if (cmd.hasOption("consolidate")) {
				isConsolidate = true;
			}
			if (cmd.hasOption("calculateThreshold")) {
				isThreshold = true;
			}
			if (cmd.hasOption("nodpi")) {
				noDPI = true;
			}
			if (cmd.hasOption("nobootstrap")) {
				nobootstrap = true;
			}
			if (cmd.hasOption("nobonferroni")) {
				nobonferroni = true;
			}
			if (cmd.hasOption("consolidatepvalue")) {
				consolidatePvalue = Double.parseDouble(cmd.getOptionValue("consolidatepvalue"));
			}
			if (cmd.hasOption("pvalue")) {
				miPvalue = Double.parseDouble(cmd.getOptionValue("pvalue"));
			}
			if (cmd.hasOption("geneNumber")) {
				geneNumber = Integer.parseInt(cmd.getOptionValue("geneNumber"));
			}
			if (cmd.hasOption("seed")) {
				seed = Integer.parseInt(cmd.getOptionValue("seed"));
			}
			if (cmd.hasOption("threads")) {
				threadCount = Integer.parseInt(cmd.getOptionValue("threads"));
			}
			if (cmd.hasOption("expfile")) {
				expPath = (String)cmd.getOptionValue("expfile_upstream");
			}
			if (cmd.hasOption("output")) {
				outputPath = (String)cmd.getOptionValue("output");
			}
			if (cmd.hasOption("tfs")) {
				tfsPath = (String)cmd.getOptionValue("tfs");
			}
			if (cmd.hasOption("kinases")) {
				kinasesPath = (String)cmd.getOptionValue("kinases");
			}

			if (isConsolidate && outputPath==null) {
				throw new ParseException("Missing required option for consolidation: output");
			}
			else if (isThreshold && (outputPath==null || expPath==null)) {
				throw new ParseException("Missing required options for finding MI threshold: expfile_upstream or output");
			}
			else if (!isConsolidate && !isThreshold && (outputPath==null || expPath==null || tfsPath==null)) {
				throw new ParseException("Missing required options for ARACNe: expfile_upstream, tfs or output");
			}

		} catch( ParseException exp ) {
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("aracne.jar", options);
		    System.out.println( exp.getMessage() );
			System.exit(1);
		}

		File outputFolder = new File(outputPath);
		
		if(seed!=null){
			random = new Random(seed);
		} else{
			random = new Random();
		}

		// Here the program forks
		// You can calculate the MI threshold
		// You can run a single bootstrap (or non bootstrap)
		// You can consolidate bootstraps

		if(isThreshold){
			File expressionFile = new File(expPath);
			outputFolder.mkdir();

			runThreshold(expressionFile,outputFolder,geneNumber,miPvalue,seed);
		}
		else if(!isConsolidate){
			File expressionFile = new File(expPath);
			File tfFile = new File(tfsPath);

			// if kinases file is present, proteomic-mode is enabled.
			File kinasesFile = null;
			if(kinasesPath != null){
				kinasesFile = new File(kinasesPath);
			}

			String processId = new BigInteger(130, random).toString(32);
			runAracne(
					expressionFile,
					tfFile,
					kinasesFile,
					outputFolder,
					processId,
					miPvalue,
					threadCount,
					singlemode,
					noDPI,
					nobootstrap
					);
		} else {
			runConsolidate(outputFolder,nobonferroni,consolidatePvalue);
		}
	}

	// Calculate Threshold mode
	private static void runThreshold(File expressionFile, File outputFolder, int geneNumber, double miPvalue, int seed) throws NumberFormatException, Exception{
		// Read expression matrix and transcription factor lists
		ExpressionMatrix em = new ExpressionMatrix(expressionFile);

		// Generate ranked data
		HashMap<String, DataVector> rankData = em.rankDV(random);
		
		// Check if the sample size
		int sampleNumber = em.getSamples().size();
		if(sampleNumber>32767){
			System.err.println("Warning: sample number is higher than the short data limit");
			System.exit(1);
		}

		//// Calculate threshold for the required p-value
		// Don't if a threshold file already exists
		File miThresholdFile = new File(outputFolder+"/miThreshold_p"+formatter.format(miPvalue)+"_samples"+sampleNumber+".txt");
		double miThreshold;

		if(miThresholdFile.exists()){
			System.out.println("MI threshold file was already there, but I am recalculating it.");
		}
		MI miCPU = new MI(rankData);
		miThreshold = miCPU.calibrateMIThresholdNA(rankData,geneNumber,miPvalue,seed);
		DataParser.writeValue(miThreshold, miThresholdFile);
	}

	// Single run mode
	private static void runAracne(
			File expressionFile,
			File transcriptionFactorsFile,
			File kinasesFile,
			File outputFolder, 
			String processId,
			Double miPvalue,
			Integer threadCount,
			boolean singlemode, // Only one expression matrix was provided
			boolean noDPI, // Do not use DPI
			boolean nobootstrap // Do not use bootstrap
			) throws NumberFormatException, Exception {
		long initialTime = System.currentTimeMillis();

		// Read expression matrix 
		ExpressionMatrix em = new ExpressionMatrix(expressionFile);

		// Generate ranked data
		HashMap<String, short[]> rankData;

		// Bootstrap matrix
		if(!nobootstrap){
			System.out.println("Bootstrapping input matrix with "+em.getGenes().size()+" genes and "+em.getSamples().size()+" samples");
			ExpressionMatrix bootstrapped = em.bootstrap(random);
			rankData = bootstrapped.rank(random);
		} else {
			rankData = em.rank(random);
		}

		// Check if the sample size is less than the short limit
		int sampleNumber = rankData.get(rankData.keySet().toArray()[0]).length;
		if(sampleNumber>32767){
			System.err.println("Warning: sample number is higher than the short data limit");
			System.exit(1);
		}

		// TF list
		HashSet<String> tfSet = DataParser.readGeneSet(transcriptionFactorsFile);
		String[] tfList = tfSet.toArray(new String[0]);
		Arrays.sort(tfList);

		if(tfList.length==0){
			System.err.println("The regulator file is badly formatted or empty");
			System.exit(1);
		}

		// kinase HashSet;
		HashSet<String> kinaseSet = null;
		if(kinasesFile!=null){
			kinaseSet = DataParser.readGeneSet(kinasesFile);
			if(kinaseSet.size()==0){
				System.err.println("The kinase regulator file is badly formatted or empty");
				System.exit(1);
			}
		}

		// Check if the threshold file exists
		File miThresholdFile = new File(outputFolder+"/miThreshold_p"+formatter.format(miPvalue)+"_samples"+sampleNumber+".txt");
		double miThreshold;
		if(!miThresholdFile.exists()){
			System.err.println("MI threshold file is not present.");
			System.err.println("Please run ARACNE in --calculateThreshold mode first");
			System.exit(1);
		}
		System.out.println("MI threshold file is present");
		miThreshold = DataParser.readValue(miThresholdFile);

		// Calculate a single ARACNE (using the APAC4 implementation by Alex)
		long time1 = System.currentTimeMillis();
		System.out.println("Calculate network from: "+expressionFile);
		HashMap<String, HashMap<String, Double>> finalNetwork = new HashMap<String, HashMap<String, Double>>();
		HashMap<String, HashMap<String, Boolean>> finalNetworkSign = new HashMap<String, HashMap<String, Boolean>>();

		MI miCPU = new MI(rankData,tfList,kinaseSet,miThreshold,threadCount);

		finalNetwork = miCPU.getFinalNetwork();
		finalNetworkSign = miCPU.getFinalNetworkSign();

		System.out.println("Time elapsed for calculating MI: "+(System.currentTimeMillis() - time1)/1000+" sec\n");

		// And then calculate DPI
		long time2 = System.currentTimeMillis();
		HashMap<String, HashSet<String>> removedEdges;
		if(noDPI){
			removedEdges = new HashMap<String, HashSet<String>>();
		}else {
			removedEdges = dpi(finalNetwork, finalNetworkSign, threadCount);
			System.out.println("DPI time elapsed: "+(System.currentTimeMillis() - time2)/1000+" sec");
		}

		// And write out the single bootstrap network
		File outputFile;
		if(nobootstrap){
			outputFile = new File(outputFolder.getAbsolutePath()+"/nobootstrap_network.txt");
		} else {
			outputFile = new File(outputFolder.getAbsolutePath()+"/bootstrapNetwork_"+processId+".txt");
		}
		writeFinal(outputFile,removedEdges,finalNetwork,finalNetworkSign);

		long finalTime = System.currentTimeMillis();
		System.out.println("Total time elapsed: "+(finalTime - initialTime)/1000+" sec");
	}

	// This method consolidates the several bootstraps
	private static void runConsolidate(File outputFolder, boolean nobonferroni, Double consolidatePvalue) throws IOException {
		BootstrapConsolidator c = new BootstrapConsolidator(nobonferroni);
		c.mergeFiles(outputFolder);
		String outputFile = outputFolder+"/network.txt";

		
		// consolidatePvalue is the P-value for the Poisson distribution. Aka how many times an edge has to appear in the bootstraps to be kept.
		// Hard-coded to 0.3 in the original ARACNe
		c.writeSignificant(outputFile, consolidatePvalue);

		System.out.println("\n        :");
		System.out.println("       :");
		System.out.println("        :       Cool! All done!");
		System.out.println("    /\\('')/\\");
		System.out.println("    \\      /");
	}


	// Method to read the expression file
	// Method to load the TF list
	// DPI method
	private static HashMap<String, HashSet<String>> dpi(
			HashMap<String,	HashMap<String, Double>> finalNet, 
			HashMap<String, HashMap<String, Boolean>> finalNetSign, 
			int _threadNumber){
		DPI dpi = new DPI();
		return dpi.dpi(finalNet,finalNetSign,_threadNumber);
	}
	

	public static void writeFinal(
	File finalDPIfile, 
		HashMap<String, HashSet<String>> removedEdges,
		HashMap<String, HashMap<String, Double>> finalNet,
		HashMap<String, HashMap<String, Boolean>> finalNetSign){
		try{
			int left = 0;
			int removed = 0;

			BufferedWriter bw = new BufferedWriter(new FileWriter(finalDPIfile));

			// Header
			bw.write("Regulator\tTarget\tMI\tSign\n");

			for(String k : finalNet.keySet()){
				HashSet<String> tr = null;
				if(removedEdges.containsKey(k)){
					tr = removedEdges.get(k);
				}
				else{
					tr = new HashSet<String>();
				}

				for(String kk : finalNet.get(k).keySet()){
					if(!tr.contains(kk)){
						bw.write(k+"\t"+kk+"\t"+finalNet.get(k).get(kk)+"\t" + finalNetSign.get(k).get(kk) +"\n");
						left++;
					} else {
						removed++;
					}
				}
			}
			bw.close();
			System.out.println("Edges removed by DPI:\t"+removed);
			System.out.println("Final Network size:\t"+left);
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
}
