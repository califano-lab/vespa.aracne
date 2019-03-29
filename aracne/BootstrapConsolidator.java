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

	private HashMap<String, Integer> edgesOccurrences = new HashMap<String, Integer>();
	private HashMap<String, Double> mi = new HashMap<String, Double>();
	
	private HashSet<String> tfs = new HashSet<String>();
	private HashSet<String> targets = new HashSet<String>();
	
	public BootstrapConsolidator() {
	}

	public void mergeMIT(File folder) throws IOException{
		FilenameFilter filter = new FilenameFilter() {
		    @Override
		    public boolean accept(File dir, String name) {
		        return name.matches("bootstrapMIThreshold_.*");
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
	
		// Calculate the occurrence of each edge (count object)
		// And sum the global mi for each edge
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
			        		String key = sp[0]+"#"+sp[1]; // edge key
			        		tfs.add(sp[0]);
			        		targets.add(sp[1]);
			        		if(edgesOccurrences.containsKey(key)){
			        			edgesOccurrences.put(key, edgesOccurrences.get(key)+1);
			        			mi.put(key, mi.get(key)+Double.parseDouble(sp[2]));
			        		}
			        		else{
			        			edgesOccurrences.put(key, 1);
			        			mi.put(key, Double.parseDouble(sp[2]));
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
		
		for(String key : edgesOccurrences.keySet()){
			// Calculate the average mi of each edge
			int edgeOccurrence=edgesOccurrences.get(key);
			mi.put(key, mi.get(key)/edgeOccurrence);
		}
	}
	
	public void writeNetwork(String outFile){
		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outFile)));

			// Header
			bw.write("Regulator\tTarget\tMI\n");

			String[] keys = edgesOccurrences.keySet().toArray(new String[0]);
			for(String key : keys){
				String[] sp = key.split("#");
				bw.write(sp[0]+"\t"+sp[1]+"\t"+mi.get(key)+"\n");
			}
			bw.close();
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
}








