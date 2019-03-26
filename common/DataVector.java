package common;

import java.util.ArrayList;

public class DataVector {
	// Variable Declaration
	public boolean [] NAs;
	public short [] values;
	
	// Constructor
	public DataVector(short[] values){
		boolean[] na = new boolean[values.length];
		for (int i = 0; i<values.length; i++){
			if (values[i]==0){
				na[i] = true;
			} else {
				na[i] = false;
			}
		}
		this.NAs=na;
		this.values=values;
	}
	
	// Methods

	/*
	 * This Method provides two vectors (for not NAs samples) and three counts:
	 * Nr. of samples with NAs in both vectors
	 * Nr. of samples with NAs in vector 1
	 * Nr. of samples with NAs in vector 2
	 * 
	*/
	public ArrayList<short[]> getQuadrants(DataVector x){
		ArrayList<short[]> output = new ArrayList<short[]>();

		// Quadrant counts
		short[] bothNAs = new short[1];
		short[] NAs1 = new short[1];
		short[] NAs2 = new short[1];
		bothNAs[0] = NAs1[0] = NAs2[0] = 0;
		
		// First loop does the counting
		int lengthOfNotNAsInBoth = 0;
		for(int i=0; i<NAs.length; i++){
			if(NAs[i]==false){
				if (x.NAs[i]==false){
					lengthOfNotNAsInBoth++;
				} else {
					NAs2[0]++;
				}
			}  else {
				if (x.NAs[i]==false){
					NAs1[0]++;
				} else {
					bothNAs[0]++;
				}
			}
		}
	
		// Second loop generates the not NAs vectors
		short [] notNA1 = new short[lengthOfNotNAsInBoth];
		short [] notNA2 = new short[lengthOfNotNAsInBoth];
		int j = 0;
		for(int i=0; i<NAs.length; i++){
			if(NAs[i]==false & x.NAs[i]==false){
				notNA1[j] = values[i];
				notNA2[j] = x.values[i];
				j++;		
			}
		}
		output.add(notNA1);
		output.add(notNA2);
		output.add(bothNAs);
		output.add(NAs1);
		output.add(NAs2);

		return output;
	}
	
	
	public short[] getNotNAvalues() {
		int length = 0;
		for(int i=0; i<NAs.length; i++){
			if(NAs[i]==false){
				length++;
			}
		}
		
		short [] notNAvalues = new short[length];
		int j = 0;
		for(int i=0; i<NAs.length; i++){
			if(NAs[i]==false){
				notNAvalues[j] = values[i];
				j++;		
			}
		}
		return notNAvalues;
	}
	
	public ArrayList<short[]> getNotNAsInBoth(DataVector x){
		int length = 0;
		for(int i=0; i<NAs.length; i++){
			if(NAs[i]==false & x.NAs[i]==false){
				length++;
			}
		}
	
		short [] notNA1 = new short[length];
		short [] notNA2 = new short[length];
		int j = 0;
		for(int i=0; i<NAs.length; i++){
			if(NAs[i]==false & x.NAs[i]==false){
				notNA1[j] = values[i];
				notNA2[j] = x.values[i];
				j++;		
			}
		}
		ArrayList<short[]> output = new ArrayList<short[]>();
		output.add(notNA1);
		output.add(notNA2);
		return output;
	}

	
	public ArrayList<short[]> getNotNAsInFirst(DataVector x){
		int length = 0;
		for(int i=0; i<NAs.length; i++){
			if(NAs[i]==false & x.NAs[i]==true){
				length++;
			}
		}
	
		short [] notNA1 = new short[length];
		short [] NA2 = new short[length];
		int j = 0;
		for(int i=0; i<NAs.length; i++){
			if(NAs[i]==false & x.NAs[i]==true){
				notNA1[j] = values[i];
				NA2[j] = x.values[i];
				j++;		
			}
		}
		ArrayList<short[]> output = new ArrayList<short[]>();
		output.add(notNA1);
		output.add(NA2);
		return output;
	}

	public ArrayList<short[]> getNotNAsInSecond(DataVector x){
		int length = 0;
		for(int i=0; i<NAs.length; i++){
			if(NAs[i]==true & x.NAs[i]==false){
				length++;
			}
		}
	
		short [] NA1 = new short[length];
		short [] notNA2 = new short[length];
		int j = 0;
		for(int i=0; i<NAs.length; i++){
			if(NAs[i]==true & x.NAs[i]==false){
				NA1[j] = values[i];
				notNA2[j] = x.values[i];
				j++;		
			}
		}
		ArrayList<short[]> output = new ArrayList<short[]>();
		output.add(NA1);
		output.add(notNA2);
		return output;
	}
	
}
