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

import common.DataParser;
import common.ExpressionMatrix;


public class Aracne {
	// Variable definition
	static NumberFormat formatter = new DecimalFormat("0.###E0");
	static float simCut = 0; // In the original paper this was a parameter. E.g. if two TFs have a very high MI, DPI is not calculated
	final static int geneNumberInMIthresholding = 3000;
	static Random random = new Random();
	private static boolean singlemode = false;

	// Main Method
	public static void main(String[] args) throws Exception {
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
		options.addOption("e", "expfile_upstream", true, "");
		options.addOption("d", "expfile_downstream", true, "");
		options.addOption("o", "output", true, "");
		options.addOption("t", "tfs", true, "");
		options.addOption("p", "pvalue", true, "");
		options.addOption("s", "seed", true, "");
		options.addOption("m", "threads", true, "");
		options.addOption("v", "consolidatepvalue", true, "");



		// Default arguments
		boolean isConsolidate = false;
		boolean isThreshold = false;
		boolean noDPI = false;
		boolean nobootstrap = false;
		boolean nobonferroni = false;

		String exp1File = null;
		String exp2File = null;
		String outputFolderPath = null;
		String tfsPath = null;
		Double miPvalue = 1E-8;
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
			if (cmd.hasOption("seed")) {
				seed = Integer.parseInt(cmd.getOptionValue("seed"));
			}
			if (cmd.hasOption("threads")) {
				threadCount = Integer.parseInt(cmd.getOptionValue("threads"));
			}
			if (cmd.hasOption("expfile_upstream")) {
				exp1File = (String)cmd.getOptionValue("expfile_upstream");
			}
			if (cmd.hasOption("expfile_downstream")) {
				exp2File = (String)cmd.getOptionValue("expfile_downstream");
			}
			if (cmd.hasOption("output")) {
				outputFolderPath = (String)cmd.getOptionValue("output");
			}
			if (cmd.hasOption("tfs")) {
				tfsPath = (String)cmd.getOptionValue("tfs");
			}

			if (isConsolidate && outputFolderPath==null) {
				throw new ParseException("Missing required option for consolidation: output");
			}
			else if (isThreshold && (outputFolderPath==null || exp1File==null)) {
				throw new ParseException("Missing required options for finding MI threshold: expfile_upstream or output");
			}
			else if (!isConsolidate && !isThreshold && (outputFolderPath==null || exp1File==null || tfsPath==null)) {
				throw new ParseException("Missing required options for ARACNe: expfile_upstream, tfs or output");
			}

		} catch( ParseException exp ) {
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("aracne.jar", options);
		    System.out.println( exp.getMessage() );
			System.exit(1);
		}

		File outputFolder = new File(outputFolderPath);
		
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
			File expressionFile = new File(exp1File);
			outputFolder.mkdir();

			runThreshold(expressionFile,outputFolder,miPvalue,seed);
		}
		else if(!isConsolidate){
			File expressionFile1 = new File(exp1File);
			File transcriptionFactorsFile = new File(tfsPath);

			if(exp2File == null){
				exp2File = exp1File;
				singlemode  = true;
			}

			File expressionFile2 = new File(exp2File);

			String processId = new BigInteger(130, random).toString(32);
			runAracne(
					expressionFile1,
					expressionFile2,
					transcriptionFactorsFile,
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
	private static void runThreshold(File _expressionFile, File _outputFolder, double miPvalue, int seed) throws NumberFormatException, Exception{
		// Read expression matrix and transcription factor lists
		ExpressionMatrix em = new ExpressionMatrix(_expressionFile);

		// Rank transformation through an external Normalization class
		HashMap<String, short[]> rankData = em.rank(random);

		// Check if the sample size
		int sampleNumber = rankData.get(rankData.keySet().toArray()[0]).length;
		if(sampleNumber>32767){
			System.err.println("Warning: sample number is higher than the short data limit");
			System.exit(1);
		}

		//// Calculate threshold for the required p-value
		// Don't if a threshold file already exists
		File miThresholdFile = new File(_outputFolder+"/miThreshold_p"+formatter.format(miPvalue)+"_samples"+sampleNumber+".txt");
		double miThreshold;

		if(miThresholdFile.exists()){
			System.out.println("MI threshold file was already there, but I am recalculating it.");
		}
		MI miCPU = new MI(rankData);
		miThreshold = miCPU.calibrateMIThreshold(sampleNumber,geneNumberInMIthresholding,miPvalue,seed);
		DataParser.writeValue(miThreshold, miThresholdFile);
	}

	// Single run mode
	private static void runAracne(
			File expressionFile1,
			File expressionFile2,
			File transcriptionFactorsFile,
			File outputFolder, 
			String processId,
			Double miPvalue,
			Integer threadCount,
			boolean singlemode, // Only one expression matrix was provided
			boolean noDPI, // Do not use DPI
			boolean nobootstrap // Do not use bootstrap
			) throws NumberFormatException, Exception {
		long initialTime = System.currentTimeMillis();
		// Read expression matrices 
		ExpressionMatrix em1 = new ExpressionMatrix(expressionFile1);
		ExpressionMatrix em2;
		if(singlemode){
			em2 = em1;
		} else {
			em2 = new ExpressionMatrix(expressionFile2);
		}
		
		HashMap<String, short[]> rankData1;
		HashMap<String, short[]> rankData2;

		
		
		
		// Don't bootstrap them
		if(nobootstrap){
			rankData1 = em1.rank(random);
			rankData2 = em2.rank(random);
		}
			
		// Or Bootstrap them
		else {
			System.out.println("Bootstrapping input matrix 1 with "+em1.getGenes().size()+" genes and "+em1.getSamples().size()+" samples");
			ExpressionMatrix bootstrapped1 = em1.bootstrap(random);
			ExpressionMatrix bootstrapped2;
			if(!singlemode){
				System.out.println("Bootstrapping input matrix 2 with "+em2.getGenes().size()+" genes and "+em2.getSamples().size()+" samples");
				if(em2.getSamples().size()!=em1.getSamples().size()){
					System.err.println("The two input files have different lengths!");
					System.exit(1);
				}
				bootstrapped2 = em2.bootstrap(em1.getBootsamples());

			} else {
				bootstrapped2 = bootstrapped1;
			}
			rankData1 = bootstrapped1.rank(random);
			rankData2 = bootstrapped2.rank(random);
		}

		// Check if the sample size is less than the short limit
		int sampleNumber = rankData1.get(rankData1.keySet().toArray()[0]).length;
		if(sampleNumber>32767){
			System.err.println("Warning: sample number is higher than the short data limit");
			System.exit(1);
		}

		// TF list
		HashSet<String> tfSet = DataParser.readGeneSet(transcriptionFactorsFile);
		String[] tfList = tfSet.toArray(new String[0]);
		Arrays.sort(tfList);

		if(tfList.length==0){
			System.err.println("The transcription factor file is badly formatted or empty");
			System.exit(1);
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
		if(singlemode){
			System.out.println("Calculate network from: "+expressionFile1);
		} else {
			System.out.println("Calculate network from: "+expressionFile1+ " and " +expressionFile2);
		}
		HashMap<String, HashMap<String, Double>> finalNetwork = new HashMap<String, HashMap<String, Double>>();
		MI miCPU = new MI(rankData1,rankData2,tfList,miThreshold,threadCount);
		finalNetwork = miCPU.getFinalNetwork();
		System.out.println("Time elapsed for calculating MI: "+(System.currentTimeMillis() - time1)/1000+" sec\n");

		// And then calculate DPI
		long time2 = System.currentTimeMillis();
		HashMap<String, HashSet<String>> removedEdges;
		if(noDPI){
			removedEdges = new HashMap<String, HashSet<String>>();
		}else {
			removedEdges = dpi(finalNetwork, threadCount);
			System.out.println("DPI time elapsed: "+(System.currentTimeMillis() - time2)/1000+" sec");
		}

		// And write out the single bootstrap network
		File outputFile;
		if(nobootstrap){
			outputFile = new File(outputFolder.getAbsolutePath()+"/nobootstrap_network.txt");
		} else {
			outputFile = new File(outputFolder.getAbsolutePath()+"/bootstrapNetwork_"+processId+".txt");
		}
		writeFinal(outputFile,removedEdges,finalNetwork);

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
			HashMap<String,
			HashMap<String, Double>> finalNet,
			int threadNumber
		){
		DPI dpi = new DPI();
		return dpi.dpi(finalNet, threadNumber);
	}

	public static void writeFinal(File finalDPIfile, HashMap<String, HashSet<String>> removedEdges, HashMap<String, HashMap<String, Double>> finalNet){
		try{
			int left = 0;
			int removed = 0;

			BufferedWriter bw = new BufferedWriter(new FileWriter(finalDPIfile));
			
			// Header
			bw.write("Regulator\tTarget\tMI\n");

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
						bw.write(k+"\t"+kk+"\t"+finalNet.get(k).get(kk)+"\n");
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






