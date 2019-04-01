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

import common.DataParser;

public class BootstrapConsolidator {
	private Double mit = 0.0;

	private HashMap<String,	HashMap<String, Double>> net = new HashMap<String,	HashMap<String, Double>>();
	private HashMap<String, HashMap<String, Boolean>> netSign = new HashMap<String, HashMap<String, Boolean>>();
	private HashMap<String, HashMap<String, Double>> netSignMI = new HashMap<String, HashMap<String, Double>>();

	public BootstrapConsolidator() {
	}

	public void mergeMIT(File folder) throws IOException{
		FilenameFilter filter = new FilenameFilter() {
			@Override
			public boolean accept(File dir, String name) {
				return name.matches("mithreshold_.*");
			}
		};
		File[] listOfFiles = folder.listFiles(filter);
	
		System.out.println("Integrating "+listOfFiles.length+" MI threshold bootstraps.");
		for (File file : listOfFiles) {
			if (file.isFile()) {
				System.out.println(file.getName());
				try {
					mit = mit + DataParser.readValue(file);
				}
				catch(Exception e){
					e.printStackTrace();
				}
			}
		}
		mit = mit/listOfFiles.length;
		System.out.println("Mean MI threshold: "+mit);
	}

	public void mergeNetworks(File folder) throws IOException{
		FilenameFilter filter = new FilenameFilter() {
			@Override
			public boolean accept(File dir, String name) {
				return name.matches("bootstrapNetwork_.*");
			}
		};
		File[] listOfFiles = folder.listFiles(filter);
	
		System.out.println("Integrating "+listOfFiles.length+" network bootstraps.");
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
							if(net.containsKey(sp[0])){
								HashMap<String, Double> targetNet = net.get(sp[0]);
								HashMap<String, Boolean> targetNetSign = netSign.get(sp[0]);
								HashMap<String, Double> targetNetSignMI = netSignMI.get(sp[0]);

								if(targetNet.containsKey(sp[1])){
									targetNet.put(sp[1], targetNet.get(sp[1]) + Double.parseDouble(sp[2]));
									// Replace sign if MI of current bootstrap is higher then others.
									if (Double.parseDouble(sp[2]) > targetNetSignMI.get(sp[1])){
										targetNetSign.put(sp[1], Boolean.parseBoolean(sp[3]));
										targetNetSignMI.put(sp[1], Double.parseDouble(sp[2]));
									}
								}
								else{
									targetNet.put(sp[1], Double.parseDouble(sp[2]));
									targetNetSign.put(sp[1], Boolean.parseBoolean(sp[3]));
									targetNetSignMI.put(sp[1], Double.parseDouble(sp[2]));
								}

								net.put(sp[0], targetNet);
								netSign.put(sp[0], targetNetSign);
								netSignMI.put(sp[0], targetNetSignMI);
							}
							else{
								HashMap<String, Double> targetNet = new HashMap<String, Double>();
								HashMap<String, Boolean> targetNetSign = new HashMap<String, Boolean>();
								HashMap<String, Double> targetNetSignMI = new HashMap<String, Double>();

								targetNet.put(sp[1], Double.parseDouble(sp[2]));
								targetNetSign.put(sp[1], Boolean.parseBoolean(sp[3]));
								targetNetSignMI.put(sp[1], Double.parseDouble(sp[2]));

								net.put(sp[0], targetNet);
								netSign.put(sp[0], targetNetSign);
								netSignMI.put(sp[0], targetNetSignMI);
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
		
		for(String regulatorKey : net.keySet()){
			HashMap<String, Double> targetNet = net.get(regulatorKey);
			for(String targetKey : targetNet.keySet()){
				targetNet.put(targetKey, targetNet.get(targetKey) / listOfFiles.length);
			}
		}
	}

	public double getMIT() {
		return mit;
	}

	public HashMap<String, HashMap<String, Double>> getFinalNetwork() {
		return net;
	}

	public HashMap<String, HashMap<String, Boolean>> getFinalNetworkSign() {
		return netSign;
	}
}
