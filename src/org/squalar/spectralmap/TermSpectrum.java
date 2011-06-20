/**
 * Spectral mapping of index terms
 *  Copyright (C) 2011 Peter Wittek
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

package org.squalar.spectralmap;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;

import org.squalar.util.DenseVector;

/**
 * Command line utility for generating the spectrum of a given term
 * in the visible range.
 */


public class TermSpectrum {

	/** The eigenvectors. */
	private double[][] eigenVectors;
	
	/** The term matrix. */
	private double[][] termMatrix;
	
	/** The constant cutOff value of similarity over which term vectors
	 *  and eigenvectors are considered similar. */
	final static double cutOff = 0.05;

	/**
	 * The main method. The output will be written in figures/spectrum_term.png.
	 *
	 * @param args[0] the co-occurrence file name
	 * @param args[1] the eigenvectors file name
	 * @param args[2] the eigenvalues file name
	 * @param args[3] the index terms file name
	 * @param args[4] the term for which the spectrum is to be drawn
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
	public static void main(String[] args) throws IOException {
	    if (args.length != 5) {
		      throw (new IllegalArgumentException("There were " + args.length + " arguments, instead of the expected 5."));
		    }

		String cooccurFilename=args[0];
		String eigenVectorsFilename=args[1];
		String eigenValuesFilename=args[2];
		String indexTermFilename=args[3];
		String term=args[4];
		TermSpectrum ts = new TermSpectrum(eigenVectorsFilename,
				cooccurFilename);
		String[] keywords = readWordList(indexTermFilename);
			int termIndex = Arrays.binarySearch(keywords,term);
			double[] cosines = ts.findCosineValuesForTerm(termIndex);
			int[] maxCosines = DenseVector.argMax(cosines, 20);
			for (int j = 0; j < maxCosines.length; j++) {
				System.out.print(maxCosines[j] + ":" + cosines[maxCosines[j]]
						+ " ");
				if (cosines[maxCosines[j]] < cutOff) {
					maxCosines[j] = 0;
				}
			System.out.println();
			double[] singularValues=DenseVector.readDoubleArray(eigenValuesFilename);
			GenerateVisibleSpectrum.drawTermSpectrum(singularValues, maxCosines, term);
		}
	}

	/**
	 * Instantiates a new term spectrum.
	 *
	 * @param eigenVectorsFilename the eigen ectors filename
	 * @param cooccurFilename the cooccurrence filename
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
	public TermSpectrum(String eigenVectorsFilename, String cooccurFilename)
			throws IOException {
		this.eigenVectors = DenseVector.readMatrix(eigenVectorsFilename);
		this.termMatrix = DenseVector.readMatrix(cooccurFilename);
	}

	/**
	 * Finds cosine values with every eigenvector for a given term.
	 *
	 * @param termIndex the index of the term
	 * @return the cosine values 
	 */
	private double[] findCosineValuesForTerm(int termIndex) {
		double[] result = new double[eigenVectors.length];
		for (int i = 0; i < eigenVectors.length; i++) {
			result[i] = DenseVector.cos(termMatrix[termIndex], eigenVectors[i]);
		}
		return result;
	}

	/**
	 * Reads a term list.
	 *
	 * @param filename the filename
	 * @return the sorted term list
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
	private static String[] readWordList(String filename) throws IOException {
		if (filename == null) {
			return null;
		}
		ArrayList<String> words = new ArrayList<String>();
		Scanner s = new Scanner(new BufferedReader(new FileReader(filename)))
				.useDelimiter("[\n\r]");
		while (s.hasNext()) {
			String tmp = s.next();
			tmp = tmp.toLowerCase();
			if (!tmp.equals("")) {
				words.add(tmp);
			}
		}
		s.close();

		String[] result = new String[words.size()];
		for (int i=0;i<words.size();i++){
			result[i]=words.get(i);
		}
		java.util.Arrays.sort(result);
		return result;
	}

}
