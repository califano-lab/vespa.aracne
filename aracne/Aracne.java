package aracne;

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

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomAdaptor;

import aracne.BootstrapConsolidator;

import common.DataParser;
import common.ExpressionMatrix;
import common.DataVector;

public class Aracne {
	// Variable definition
	static NumberFormat formatter = new DecimalFormat("0.###E0");
	static MersenneTwister random = new MersenneTwister();

	// Main Method
	public static void main(String[] args) throws Exception {

		Locale.setDefault(Locale.US);

		System.out.println("ARACNe (" + Aracne.class.getPackage().getImplementationVersion() + ")");

		//// Parse arguments
		CommandLineParser parser = new DefaultParser();
		Options options = new Options();

		// Flag arguments
		options.addOption("c", "consolidate", false, "Run ARACNe in consolidation mode");
		options.addOption("t", "threshold", false, "Run ARACNe in MI threshold estimation mode");
		options.addOption("nd", "nodpi", false, "Run ARACNe without DPI");
		options.addOption("nb", "nobootstrap", false, "Run ARACNe without bootstrapping");
		options.addOption("nm", "nobonferroni", false, "Run ARACNe without Bonferroni correction");

		// Arguments with values
		options.addOption("e", "expfile", true, "Expression Matrix (M x N); M=genes, N=samples; Designate missing values with NA");
		options.addOption("o", "output", true, "Output directory");
		options.addOption("r", "regulators", true, "Regulator identifier file (e.g. transcription factors or kinases & phosphatases)");
		options.addOption("a", "activators", true, "Activator identifier file (e.g. kinases)");
		options.addOption("f", "fwer", true, "Threshold estimation mode: family-wise error-rate [default: 0.05]");
		options.addOption("s", "seed", true, "Optional seed for reproducible results [default: random]");
		options.addOption("j", "threads", true, "Number of threads to use [default: 1]");
		options.addOption("p", "pvalue", true, "P-value threshold for the Poisson test of edge significance [default: 0.05]");

		// Default arguments
		boolean isConsolidate = false;
		boolean isThreshold = false;
		boolean noDPI = false;
		boolean nobootstrap = false;
		boolean nobonferroni = false;

		String expPath = null;
		String outputPath = null;
		String regulatorsPath = null;
		String activatorsPath = null;
		Double fwer = 0.05;
		Integer seed = null;
		Integer threadCount = 1;
		Double pvalue = 0.05; 

		// Parse arguments
		try {
			CommandLine cmd = parser.parse( options, args);

			if (cmd.hasOption("consolidate")) {
				isConsolidate = true;
			}
			if (cmd.hasOption("threshold")) {
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
			if (cmd.hasOption("pvalue")) {
				pvalue = Double.parseDouble(cmd.getOptionValue("pvalue"));
			}
			if (cmd.hasOption("fwer")) {
				fwer = Double.parseDouble(cmd.getOptionValue("fwer"));
			}
			if (cmd.hasOption("seed")) {
				seed = Integer.parseInt(cmd.getOptionValue("seed"));
			}
			if (cmd.hasOption("threads")) {
				threadCount = Integer.parseInt(cmd.getOptionValue("threads"));
			}
			if (cmd.hasOption("expfile")) {
				expPath = (String)cmd.getOptionValue("expfile");
			}
			if (cmd.hasOption("output")) {
				outputPath = (String)cmd.getOptionValue("output");
			}
			if (cmd.hasOption("regulators")) {
				regulatorsPath = (String)cmd.getOptionValue("regulators");
			}
			if (cmd.hasOption("activators")) {
				activatorsPath = (String)cmd.getOptionValue("activators");
			}

			if (isConsolidate && outputPath==null) {
				throw new ParseException("Required option: output");
			}
			else if (!isConsolidate && (outputPath==null || expPath==null || regulatorsPath==null)) {
				throw new ParseException("Required options: expfile, regulators and output");
			}

		} catch( ParseException exp ) {
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("aracne.jar", options);
		    System.out.println( exp.getMessage() );
			System.exit(1);
		}

		File outputFolder = new File(outputPath);
		
		// Set seed if specified
		if(seed!=null){
			random = new MersenneTwister(seed);
		} else{
			random = new MersenneTwister();
		}

		// Read expression matrix if present
		File expressionFile = null;
		ExpressionMatrix em = null;
		if(expPath != null){
			expressionFile = new File(expPath);
			em = new ExpressionMatrix(expressionFile);
		}

		// Read regulators file if present
		File regulatorsFile = null;
		String[] regulators = null;
		if(regulatorsPath != null){
			regulatorsFile = new File(regulatorsPath);
			regulators = DataParser.readGeneSet(regulatorsFile, em.getGenes());
		}

		// Read activators file if present
		File activatorsFile = null;
		String[] activators = null;
		if(activatorsPath != null){
			activatorsFile = new File(activatorsPath);
			activators = DataParser.readGeneSet(activatorsFile, em.getGenes());
		}

		// Main ARACNe routines
		// ARACNe threshold estimation mode
		if(isThreshold){
			outputFolder.mkdir();

			runThreshold(em,regulators,outputFolder,fwer,seed);
		}
		// ARACNe bootstrapping / standard mode
		else if(!isConsolidate){
			String processId = new BigInteger(130, new RandomAdaptor(random)).toString(32);

			runAracne(
					em,
					regulators,
					activators,
					outputFolder,
					processId,
					fwer,
					threadCount,
					noDPI,
					nobootstrap
					);
		// ARACNe consolidation mode
		} else {
			runConsolidate(outputFolder,nobonferroni,pvalue);
		}
	}

	// ARACNe MI Threshold mode
	private static void runThreshold(ExpressionMatrix em, String[] regulators, File outputFolder, double fwer, int seed) throws NumberFormatException, Exception{
		// Generate ranked data
		HashMap<String, DataVector> rankData = em.rank(random);

		// Get number of samples and genes
		int sampleNumber = em.getSamples().size();
		int geneNumber = em.getGenes().size();

		// Initialize MI threshold file
		File miThresholdFile = new File(outputFolder+"/fwer_"+formatter.format(fwer)+"_samples"+sampleNumber+".txt");
		double miThreshold;

		if(miThresholdFile.exists()){
			System.out.println("Warning: Overwriting existing MI threshold file.");
		}

		// Compute miPvalue for thresholding
		double miPvalue = fwer / (regulators.length * geneNumber-1);
		System.out.println("Estimating p-value threshold for "+regulators.length+" regulators, "+geneNumber+" genes and FWER="+fwer+": "+miPvalue);

		// Compute miThreshold and write to file
		MI miCPU = new MI(rankData, regulators, miPvalue, seed);
		DataParser.writeValue(miCPU.getThreshold(), miThresholdFile);
	}

	// ARACNe bootstrapping / main mode
	private static void runAracne(
			ExpressionMatrix em,
			String[] regulators,
			String[] activators,
			File outputFolder, 
			String processId,
			Double fwer,
			Integer threadCount,
			boolean noDPI, // Do not use DPI
			boolean nobootstrap // Do not use bootstrap
			) throws NumberFormatException, Exception {
		long initialTime = System.currentTimeMillis();

		HashMap<String, DataVector> rankData;
		// Bootstrap matrix
		if(!nobootstrap){
			System.out.println("Bootstrapping input matrix with "+em.getGenes().size()+" genes and "+em.getSamples().size()+" samples.");
			ExpressionMatrix bootstrapped = em.bootstrap(random);
			rankData = bootstrapped.rank(random);
		} else {
			rankData = em.rank(random);
		}

		// Apply directional DPI if activators were specified
		boolean dDPI = false;
		if(activators!=null){
			dDPI = true;
		}

		// Check if the threshold file exists
		int sampleNumber = em.getSamples().size();
		File miThresholdFile = new File(outputFolder+"/fwer_"+formatter.format(fwer)+"_samples"+sampleNumber+".txt");
		double miThreshold;
		if(!miThresholdFile.exists()){
			throw new IOException("MI threshold file is not present. Run ARACNe in threshold mode first.");
		}
		miThreshold = DataParser.readValue(miThresholdFile);

		// Compute a single ARACNE (using the APAC4 implementation by Alex)
		long time1 = System.currentTimeMillis();
		System.out.println("Compute network.");

		MI miCPU = new MI(rankData,regulators,activators,miThreshold,threadCount);

		// MI between two interactors
		HashMap<String, HashMap<String, Double>> finalNetwork = miCPU.getFinalNetwork();
		// Correlation sign (directionality) between two interactors
		HashMap<String, HashMap<String, Boolean>> finalNetworkSign = miCPU.getFinalNetworkSign();

		System.out.println("Time elapsed for calculating MI: "+(System.currentTimeMillis() - time1)/1000+" sec\n");

		// Compute DPI
		long time2 = System.currentTimeMillis();
		HashMap<String, HashSet<String>> removedEdges;
		if(noDPI){
			removedEdges = new HashMap<String, HashSet<String>>();
		}else {
			removedEdges = dpi(finalNetwork, finalNetworkSign, threadCount, dDPI);
			System.out.println("DPI time elapsed: "+(System.currentTimeMillis() - time2)/1000+" sec");
		}

		// Write bootstrapped network to file
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

	// ARACNe consolidation mode
	private static void runConsolidate(File outputFolder, boolean nobonferroni, Double pvalue) throws IOException {
		BootstrapConsolidator c = new BootstrapConsolidator(nobonferroni);
		c.mergeFiles(outputFolder);
		String outputFile = outputFolder.getAbsolutePath()+"/network.txt";

		
		// P-value for the Poisson distribution, also describes how many times an edge has to appear in the bootstraps to be kept.
		// Hard-coded to 0.3 in the original ARACNe
		c.writeSignificant(outputFile, pvalue);

		System.out.println("\n        :");
		System.out.println("       :");
		System.out.println("        :       Cool! All done!");
		System.out.println("    /\\('')/\\");
		System.out.println("    \\      /");
	}

	// Compute DPI
	private static HashMap<String, HashSet<String>> dpi(
			HashMap<String,	HashMap<String, Double>> finalNet, 
			HashMap<String, HashMap<String, Boolean>> finalNetSign, 
			int threadNumber,
			boolean dDPI){
		DPI dpi = new DPI();
		return dpi.dpi(finalNet, finalNetSign, threadNumber, dDPI);
	}
	
	// Write results
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
