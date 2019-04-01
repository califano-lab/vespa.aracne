package common;

import java.util.ArrayList;

/**
 * DataVector is a generic structure to store ranked vectors and NA values
 *
 * @param  values A ranked vector, where NaN == NA and all other ranks are integers.
 * @return DataVector
 */
public class DataVector {
	// Variable Declaration
	public boolean [] NAs;
	public short [] values;
	
	// Constructor
	public DataVector(short[] values){
		boolean[] na = new boolean[values.length];
		for (int i = 0; i<values.length; i++){
			if (values[i]==Double.NaN){
				na[i] = true;
			} else {
				na[i] = false;
			}
		}
		this.NAs=na;
		this.values=values;
	}
	
	/**
	 * Returns quadrants for MI estimation of object against DataVector x
	 *
	 * @param  x DataVector
	 * @return Arraylist with the following vectors:
	 * 0-1: Intersecting reranked ranks
	 * 2: Nr. of samples with NAs in both vectors
	 * 3: Nr. of samples with NAs in vector of object
	 * 4: Nr. of samples with NAs in vector x
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

		// Rerank if necessary
		if (lengthOfNotNAsInBoth != values.length){
			notNA1 = reRankVector(notNA1);
			notNA2 = reRankVector(notNA2);
		}

		output.add(notNA1);
		output.add(notNA2);
		output.add(bothNAs);
		output.add(NAs1);
		output.add(NAs2);

		return output;
	}

	public static short[] reRankVector(short[] inputVector){
		short[] rankVector = new short[inputVector.length];
		for(int i=0; i<inputVector.length; i++){
			int counter = 1;
			for(int j=0; j<inputVector.length; j++){
				if(inputVector[i] > inputVector[j]){
					counter++;
				}
			}
			rankVector[i] = (short)counter;
		}
		return rankVector;
	}	
}

