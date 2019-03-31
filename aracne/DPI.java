package aracne;

import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * Apply data processing inequality to a network
 *
 * @param  finalNet HashMap linking regulators with targets and their MI
 * @param  finalNetSign HashMap linking regulators with targets and their mode of interaction
 * @param  threadCount Number of threads to use
 * @param  dDPI Whether directional DPI, accouting for mode of interaction should be used
 * @return DPI-filtered HashMap linking regulators with targets
 */
public class DPI {

	public DPI(){
		
	}

	/**
	 * Apply data processing inequality to a network
	 *
	 * @param  finalNet HashMap linking regulators with targets and their MI
	 * @param  finalNetSign HashMap linking regulators with targets and their mode of interaction
	 * @param  threadCount Number of threads to use
	 * @param  dDPI Whether directional DPI, accouting for mode of interaction should be used
	 * @return DPI-filtered HashMap linking regulators with targets
	 */
	public HashMap<String, HashSet<String>> dpi(
			HashMap<String, HashMap<String, Double>> finalNet,
			HashMap<String, HashMap<String, Boolean>> finalNetSign,
			int threadCount,
			boolean dDPI
			){
			
			// Keys in finalNet are regulators
			String[] regulators = finalNet.keySet().toArray(new String[0]);

			// First level keys are regulators, second level keys are regulators, doubles are MI values
			HashMap<String, HashMap<String, Double>> regregNetwork = new HashMap<String, HashMap<String, Double>>();
			HashMap<String, HashMap<String, Boolean>> regregNetworkSign = new HashMap<String, HashMap<String, Boolean>>();
			HashMap<String, HashSet<String>> removedEdges = new HashMap<String,HashSet<String>>();

			for(String regulator : regulators){
				HashSet<String> temp = new HashSet<String>();
				removedEdges.put(regulator, temp);
			}
			
			// This loop populates the regulator-regulator MI HashMap (regregNetwork)
			for (int i = 0; i < regulators.length; i++) {
				for (int j = 0; j < regulators.length; j++) {
					if (finalNet.get(regulators[i]).containsKey(regulators[j])) {
						
						if (regregNetwork.containsKey(regulators[i])) {
							regregNetwork.get(regulators[i]).put(regulators[j], finalNet.get(regulators[i]).get(regulators[j]));
							regregNetworkSign.get(regulators[i]).put(regulators[j], finalNetSign.get(regulators[i]).get(regulators[j]));

						} else {
							HashMap<String, Double> dtemp = new HashMap<String, Double>();
							dtemp.put(regulators[j], finalNet.get(regulators[i]).get(regulators[j]));
							regregNetwork.put(regulators[i], dtemp);
							HashMap<String, Boolean> dtempSign = new HashMap<String, Boolean>();
							dtempSign.put(regulators[j], finalNetSign.get(regulators[i]).get(regulators[j]));
							regregNetworkSign.put(regulators[i], dtempSign);
						}
					}
				}
			}

			// Check common targets between regulators and apply (d)DPI on all regulator-target-regulator clique triplets
			try{
				ExecutorService executor = Executors.newFixedThreadPool(threadCount);
				
				for (int i = 0; i < regulators.length; i++) {
					// Apply DPI if the first regulator is in our regregMI
					if (regregNetwork.containsKey(regulators[i])) {
						DPIThread mt = new DPIThread(regulators, finalNet, finalNetSign, regregNetwork, regregNetworkSign, dDPI, i, removedEdges);
						executor.execute(mt);
					}
				}
				
				executor.shutdown();
		        while (!executor.isTerminated()) {
		        	// do not continue before all threads finish
		        }
				
				return(removedEdges);
			}
			catch(Exception e){
				e.printStackTrace();
			}
			return(null);
	}

	/**
	 * Single thread of DPI computation step
	 *
	 * @param  regulators Array of regulators (e.g. transcription factors or kinases & phosphatases)
	 * @param  finalNet HashMap linking regulators with targets and their MI
	 * @param  finalNetSign HashMap linking regulators with targets and their mode of interaction
	 * @param  regregNetwork HashMap linking regulators with regulators and their mode of interaction
	 * @param  regregNetworkSign HashMap linking regulators with regulators and their mode of interaction
	 * @param  dDPI Whether directional DPI, accouting for mode of interaction should be used
	 * @param  regulatorIndex Index of processed regulator
	 * @param  removedEdges Edges to remove
	 */
	public class DPIThread extends Thread{
		// Variables
		private String[] regulators;
		private HashMap<String, HashMap<String, Double>> finalNet;
		private HashMap<String, HashMap<String, Boolean>> finalNetSign;
		HashMap<String, HashMap<String, Double>> regregNetwork;
		HashMap<String, HashMap<String, Boolean>> regregNetworkSign;
		private boolean dDPI = false;
		private int regulatorIndex = 0;
		HashMap<String, HashSet<String>> removedEdges;
		
		// Constructor
		public DPIThread(
				String[] regulators,
				HashMap<String, HashMap<String, Double>> finalNet, 
				HashMap<String, HashMap<String, Boolean>> finalNetSign, 
				HashMap<String, HashMap<String, Double>> regregNetwork, 
				HashMap<String, HashMap<String, Boolean>> regregNetworkSign, 
				boolean dDPI,
				int regulatorIndex,
				HashMap<String, HashSet<String>> removedEdges
				){
			this.regulators = regulators;
			this.finalNet = finalNet;
			this.finalNetSign = finalNetSign;
			this.regregNetwork = regregNetwork;
			this.regregNetworkSign = regregNetworkSign;
			this.dDPI = dDPI;
			this.regulatorIndex = regulatorIndex;
			this.removedEdges = removedEdges;
		}
		
		public void run(){	
			HashSet<String> targetsOfI = new HashSet<String>(finalNet.get(regulators[regulatorIndex]).keySet());
			HashMap<String, Double> fin1 = finalNet.get(regulators[regulatorIndex]);
			HashMap<String, Double> tft1 = regregNetwork.get(regulators[regulatorIndex]);
			HashMap<String, Boolean> fin1Sign = finalNetSign.get(regulators[regulatorIndex]);
			HashMap<String, Boolean> tft1Sign = regregNetworkSign.get(regulators[regulatorIndex]);
			HashSet<String> rem1 = removedEdges.get(regulators[regulatorIndex]);
			
			for (int j = regulatorIndex + 1; j < regulators.length; j++) {
				// And if the second regulator has an edge with the first regulator...
				if (tft1.containsKey(regulators[j])) {
					
					HashMap<String, Double> fin2 = finalNet.get(regulators[j]);
					HashMap<String, Boolean> fin2Sign = finalNetSign.get(regulators[j]);

					HashSet<String> rem2 = removedEdges.get(regulators[j]);
					
					// targetsOfJ is all genes having a significant edge with the regulator j
					HashSet<String> targetsOfJ = new HashSet<String>(fin2.keySet());
					
					// Intersection: targetsOfJ now contains only targets in common between the two regulators
					targetsOfJ.retainAll(targetsOfI);		// this eats most of the time needed in the DPI calculation
					
					// Obtain the regulator-regulator MI
					double regregMI = tft1.get(regulators[j]);
					boolean regregMISign = tft1Sign.get(regulators[j]);
					
					// Loop over the common targets
					for (String target : targetsOfJ){
						double v1 = fin1.get(target);
						double v2 = fin2.get(target);
						
						boolean s1 = fin1Sign.get(target);
						boolean s2 = fin2Sign.get(target);

						// apply directional DPI for regulators with known mode of interaction
						if (dDPI)
						{
							// regulator is positive correlation
							if (regregMISign & (s1 == s2)){
								if (v1 < regregMI && v1 < v2) {
									synchronized(rem1){
										rem1.add(target);
									}
								} else if (v2 < regregMI && v2 < v1) {
									synchronized(rem2){
										rem2.add(target);
									}
								}
							// regulator is negative correlation
							}else if( (!regregMISign) & (s1 != s2)){
								if (v1 < regregMI && v1 < v2) {
									synchronized(rem1){
										rem1.add(target);
									}
								} else if (v2 < regregMI && v2 < v1) {
									synchronized(rem2){
										rem2.add(target);
									}
								}
							}else{
								synchronized(rem1){
									rem1.add(target);
								}
								synchronized(rem2){
									rem2.add(target);
								}
							}
						} else {
							if (v1 < regregMI && v1 < v2) {
								synchronized(rem1){
									rem1.add(target);
								}
							} else if (v2 < regregMI && v2 < v1) {
								synchronized(rem2){
									rem2.add(target);
								}
							}
						}
					}
				}
			}
		}
	}
}
