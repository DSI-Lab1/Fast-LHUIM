package algo;

import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.net.URL;
import java.util.HashMap;



public class MainTestFast_LHUIM {

	public static void main(String [] arg) throws IOException{

		// the input and output file paths
		String input = fileToPath("T10I4D100K_50.txt");
		String output = ".//output.txt";
		
		// the minutil threshold
		int lminUtil = 10000; 
		int minLen = 5000;

//		int lminUtil = 200000; 
//		int minLen = 360;


		AlgoFast_LHUIM algo = new AlgoFast_LHUIM();
		algo.runAlgorithm(lminUtil, minLen, input, output, true, Integer.MAX_VALUE, true);
		// Print statistics
		algo.printStats();

		
	}
	
	public static String fileToPath(String filename) throws UnsupportedEncodingException{
		URL url = MainTestFast_LHUIM.class.getResource(filename);
		 return java.net.URLDecoder.decode(url.getPath(),"UTF-8");
	}
}
