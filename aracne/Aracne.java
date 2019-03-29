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

import aracne.BootstrapConsolidator;
import common.DataParser;
import common.ExpressionMatrix;

import common.DataVector;

public class Aracne {
	// Variable definition
	static NumberFormat formatter = new DecimalFormat("0.###E0");
	static Random random = new Random();

	// Main Method
	public static void main(String[] args) throws Exception {

		Locale.setDefault(Locale.US);

		System.out.println("ARACNe (" + Aracne.class.getPackage().getImplementationVersion() + ")");

		//// Parse arguments
		CommandLineParser parser = new DefaultParser();
		Options options = new Options();

		// Flag arguments
		options.addOption("c", "consolidate", false, "Run ARACNe in consolidation mode");
		options.addOption("n", "nodpi", false, "Run ARACNe without DPI");
		options.addOption("b", "nobootstrap", false, "Run ARACNe without bootstrapping");
		options.addOption("r", "nobonferroni", false, "Run ARACNe without Bonferroni correction");

		// Arguments with values
		options.addOption("e", "expfile", true, "Expression Matrix (M x N); M=genes, N=samples; Designate missing values with NA");
		options.addOption("o", "output", true, "Output directory");
		options.addOption("k", "kinases", true, "Kinase identifier file; enables phosphoproteomics dDPI");
		options.addOption("t", "tfs", true, "Regulator identifier file (transcription factors or kinases/phosphatases)");
		options.addOption("f", "fwer", true, "Threshold mode: family-wise error-rate");
		options.addOption("g", "numberOfGenes", true, "Threshold mode: Number of randomly sampled genes [default: 3000]");
		options.addOption("s", "seed", true, "Optional seed for reproducible results [default: random]");
		options.addOption("m", "threads", true, "Number of threads to use [default: 1]");
		options.addOption("v", "consolidatepvalue", true, "Bootstrapping: p-value threshold for the Poisson test of edge significance [default: 0.05]");

		// Default arguments
		boolean isConsolidate = false;
		boolean noDPI = false;
		boolean nobootstrap = false;
		boolean nobonferroni = false;

		String expPath = null;
		String outputPath = null;
		String tfsPath = null;
		String kinasesPath = null;
		Double fwer = 0.05;
		Integer numberOfGenes = 3000;
		Integer seed = null;
		Integer threadCount = 1;
		Double consolidatePvalue = 0.05; 

		// Parse arguments
		try {
			CommandLine cmd = parser.parse( options, args);

			if (cmd.hasOption("consolidate")) {
				isConsolidate = true;
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
			if (cmd.hasOption("fwer")) {
				fwer = Double.parseDouble(cmd.getOptionValue("fwer"));
			}
			if (cmd.hasOption("numberOfGenes")) {
				numberOfGenes = Integer.parseInt(cmd.getOptionValue("numberOfGenes"));
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
			if (cmd.hasOption("tfs")) {
				tfsPath = (String)cmd.getOptionValue("tfs");
			}
			if (cmd.hasOption("kinases")) {
				kinasesPath = (String)cmd.getOptionValue("kinases");
			}

			if (isConsolidate && outputPath==null) {
				throw new ParseException("Required option: output");
			}
			else if (!isConsolidate && (outputPath==null || expPath==null || tfsPath==null)) {
				throw new ParseException("Required options: expfile, tfs and output");
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
		// You can run a single bootstrap (or non bootstrap)
		// You can consolidate bootstraps
		if(!isConsolidate){
			File expressionFile = new File(expPath);
			File tfFile = new File(tfsPath);
			outputFolder.mkdir();

			// If kinases file is present, proteomic-mode is enabled.
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
					fwer,
					threadCount,
					noDPI,
					nobootstrap,
					seed
					);
		} else {
			runConsolidate(outputFolder,nobonferroni,consolidatePvalue);
		}
	}

	// ARACNe mode
	private static void runAracne(
			File expressionFile,
			File transcriptionFactorsFile,
			File kinasesFile,
			File outputFolder, 
			String processId,
			Double fwer,
			Integer threadCount,
			boolean noDPI, // Do not use DPI
			boolean nobootstrap, // Do not use bootstrap
			Integer seed
			) throws NumberFormatException, Exception {
		long initialTime = System.currentTimeMillis();

		// Read expression matrix 
		ExpressionMatrix em = new ExpressionMatrix(expressionFile);

		// Generate ranked data
		HashMap<String, DataVector> rankData;

		// Bootstrap matrix
		if(!nobootstrap){
			System.out.println("Bootstrapping input matrix with "+em.getGenes().size()+" genes and "+em.getSamples().size()+" samples");
			ExpressionMatrix bootstrapped = em.bootstrap(random);
			rankData = bootstrapped.rankDV(random);
		} else {
			rankData = em.rankDV(random);
		}

		// Check if the sample size is less than the short limit
		int sampleNumber = rankData.get(rankData.keySet().toArray()[0]).values.length;
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
		String[] kinList = null;
		boolean dDPI = false;
		if(kinasesFile!=null){
			HashSet<String> kinSet = DataParser.readGeneSet(kinasesFile);
			if(kinSet.size()==0){
				System.err.println("The kinase regulator file is badly formatted or empty");
				System.exit(1);
			}
			kinList = kinSet.toArray(new String[0]);
			dDPI = true;
		}

		// Compute MI threshold
		int geneNumber = em.getGenes().size();

		// Compute miPvalue
		double miPvalue = fwer / (tfSet.size() * geneNumber-1);
		System.out.println("Estimating MI threshold p-value for "+tfSet.size()+" regulators, "+geneNumber+" genes and FWER="+fwer+": "+miPvalue);

		MI miTCPU = new MI(rankData,tfList);
		double miThreshold = miTCPU.calibrateMIThreshold(rankData,miPvalue,seed);
		File miThresholdFile = new File(outputFolder.getAbsolutePath()+"/bootstrapMIThreshold_"+processId+".txt");
		DataParser.writeValue(miThreshold, miThresholdFile);

		// Calculate a single ARACNE (using the APAC4 implementation by Alex)
		long time1 = System.currentTimeMillis();
		System.out.println("Calculate network from: "+expressionFile);
		HashMap<String, HashMap<String, Double>> finalNetwork = new HashMap<String, HashMap<String, Double>>();
		HashMap<String, HashMap<String, Boolean>> finalNetworkSign = new HashMap<String, HashMap<String, Boolean>>();

		MI miCPU = new MI(rankData,tfList,kinList,threadCount);

		finalNetwork = miCPU.getFinalNetwork();
		finalNetworkSign = miCPU.getFinalNetworkSign();

		System.out.println("Time elapsed for calculating MI: "+(System.currentTimeMillis() - time1)/1000+" sec\n");

		// And then calculate DPI
		long time2 = System.currentTimeMillis();
		HashMap<String, HashSet<String>> removedEdges;
		if(noDPI){
			removedEdges = new HashMap<String, HashSet<String>>();
		}else {
			removedEdges = dpi(finalNetwork, finalNetworkSign, threadCount, dDPI);
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

	// Consolidate mode
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

	// Compute DPI
	private static HashMap<String, HashSet<String>> dpi(
			HashMap<String,	HashMap<String, Double>> finalNet, 
			HashMap<String, HashMap<String, Boolean>> finalNetSign, 
			int _threadNumber,
			boolean dDPI){
		DPI dpi = new DPI();
		return dpi.dpi(finalNet,finalNetSign,_threadNumber,dDPI);
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
