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

package org.squalar.util;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.StringTokenizer;

/**
 * The utility class DenseVector deals with double[][] matrices.
 */
public class DenseVector {

	/**
	 * The k arg max of an int array.
	 *
	 * @param values the values
	 * @param k the number of arg max values to find
	 * @return the arg max values
	 */
	public static int[] argMax(double[] values, int k) {
		int[] keys = new int[values.length];
		for (int i = 0; i < keys.length; i++)
			keys[i] = i;
		int[] result = new int[k];
		double[] copyValues = new double[values.length];
		for (int i = 0; i < values.length; i++)
			copyValues[i] = values[i];
		for (int i = 0; i < k; i++) {
			double max = Double.NEGATIVE_INFINITY;
			int argmax = 0;
			for (int j = 0; j < copyValues.length; j++)
				if (copyValues[j] > max) {
					max = copyValues[j];
					argmax = j;
				}
			copyValues[argmax] = Double.NEGATIVE_INFINITY;
			result[i] = keys[argmax];
		}
		return result;
	}

	/**
	 * Convert from sparse matrix to a dense double[][] matrix.
	 *
	 * @param mx the sparse matrix
	 * @return the double[][] matrix
	 */
	public static double[][] convertFromSparseMatrix(VectorNode[][] mx) {
		if (SparseVector.findMinColumnIndex(mx)==0){
			mx=SparseVector.shiftColumns(mx,1);
		}
		int maxIndex = SparseVector.findMaxColumnIndex(mx);
		double[][] result = new double[mx.length][];
		for (int i = 0; i < mx.length; i++) {
			result[i] = new double[maxIndex];
			if (mx[i] != null) {
				for (int j = 0; j < mx[i].length; j++) {
					result[i][mx[i][j].index - 1] = mx[i][j].value;
				}
			}
		}
		return result;
	}
	
	/**
	 * Cosine value of two double arrays.
	 *
	 * @param x the x vector 
	 * @param y the y vector
	 * @return the cosine
	 */
	public static double cos(double[] x, double[] y){
		return dotProduct(x,y)/(norm(x)*norm(y));
	}
	
	/**
	 * Dot product of two double arrays.
	 *
	 * @param x the x vector
	 * @param y the y vector
	 * @return the dot product
	 */
	public static double dotProduct(double[] x, double[] y){
		double result=0;
		if (x.length!=y.length){
			return Double.NaN;
		}
		for (int i=0;i<x.length;i++){
			result+=x[i]*y[i];
		}
		return result;
	}

	/**
	 * Maximum in a double array.
	 *
	 * @param x the x vector
	 * @return the max value
	 */
	public static double max(double[] x) {
		double result = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < x.length; i++) {
			if (result < x[i]) {
				result = x[i];
			}
		}
		return result;
	}

	/**
	 * Norm of a double array.
	 *
	 * @param x the x vector
	 * @return the norm
	 */
	public static double norm(double[] x){
		double result=0;
		for (int i=0;i<x.length;i++){
			result+=x[i]*x[i];
		}
		return Math.sqrt(result);
	}

	/**
	 * Reads a dense double[][] matrix.
	 *
	 * @param fileName the file name
	 * @return the matrix
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
	public static double[][] readMatrix(String fileName) throws IOException {
		Scanner scn = new Scanner(new BufferedReader(new FileReader(fileName)))
				.useDelimiter("\n");
		int nRow = 0;
		int nCol = 0;
		while (scn.hasNext()) {
			if (nRow == 0) {
				StringTokenizer st = new StringTokenizer(scn.next(),", ");
				while (st.hasMoreTokens()) {
					st.nextToken();
					nCol++;
				}
			} else {
				scn.next();
			}
			nRow++;
		}
		scn.close();
		double[][] result = new double[nRow][nCol];
		scn = new Scanner(new BufferedReader(new FileReader(fileName)))
				.useDelimiter("\n");
		int i = 0, j = 0;
		while (scn.hasNext()) {
			StringTokenizer st = new StringTokenizer(scn.next(),", ");
			while (st.hasMoreTokens()) {
				result[i][j] = Double.valueOf(st.nextToken()).doubleValue();
				j++;
			}
			j = 0;
			i++;
		}
		scn.close();
		return result;
	}

	/**
	 * Reads a double array.
	 *
	 * @param filename the file name
	 * @return the double array
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
	public static double[] readDoubleArray(String filename) throws IOException {
		ArrayList<Double> tmp = new ArrayList<Double>();
		Scanner s = new Scanner(new BufferedReader(new FileReader(filename)))
		.useDelimiter("\n");
		while (s.hasNext()) {
			tmp.add(Double.valueOf(s.next()));
		}
		s.close();
		double[] result = new double[tmp.size()];
		for (int i = 0; i < result.length; i++)
			result[i] = (double) tmp.get(i);
		return result;

	}

	/**
	 * Writes a dense matrix.
	 *
	 * @param mx the dense matrix
	 * @param filename the file name
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
	public static void writeDenseMatrix(double[][] mx, String filename)
	throws IOException {
		writeDenseMatrix(mx, filename," ");
	}

	
	/**
	 * Write dense matrix.
	 *
	 * @param mx the dense matrix
	 * @param filename the file name
	 * @param delimiter the delimiter between values
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
	public static void writeDenseMatrix(double[][] mx, String filename, String delimiter)
			throws IOException {
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(
				filename)));
		for (int i = 0; i < mx.length; i++) {
			for (int j = 0; j < mx[i].length; j++) {
				out.write(mx[i][j] + "");
				if (j != mx[i].length - 1) {
					out.write(delimiter);
				}
			}
			out.write("\n");
		}
		out.close();
	}

}
