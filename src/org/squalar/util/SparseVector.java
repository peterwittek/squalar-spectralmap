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
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.StringTokenizer;

/**
 * The utility class SparseVector deals with libsvm-formatted sparse matrices.
 */
public class SparseVector {

	private static VectorNode[] addNode(VectorNode[] v, int index, double value) {
		if (v == null || v[0] == null) {
			VectorNode[] w = new VectorNode[1];
			w[0] = new VectorNode();
			w[0].index = index;
			w[0].value = value;
			return w;
		}
		VectorNode[] w = new VectorNode[v.length + 1];
		boolean inserted = false;
		int j = 0;
		for (int i = 0; i < v.length; i++) {
			w[j] = new VectorNode();
			w[j].index = v[i].index;
			w[j].value = v[i].value;
			if (v[i].index == index) {
				w[j].value = w[j].value + value;
				inserted = true;
			}
			if (!inserted && index < v[i].index) {
				w[j].index = index;
				w[j].value = value;
				inserted = true;
				i--;
			}
			j++;
		}
		if (!inserted) {
			w[j] = new VectorNode(index, value);
			j++;
		}
		if (j == v.length + 1)
			return w;
		VectorNode[] result = new VectorNode[v.length];
		for (int i = 0; i < v.length; i++) {
			result[i] = new VectorNode();
			result[i].index = w[i].index;
			result[i].value = w[i].value;
		}
		return result;
	}

	/**
	 * Extracts a column vector.
	 *
	 * @param mx the sparse matrix
	 * @param columnIndex the column index
	 * @return sparse vector of the requested column
	 */
	public static VectorNode[] columnVector(VectorNode[][] mx, int columnIndex) {
		ArrayList<VectorNode> v = new ArrayList<VectorNode>();
		for (int i = 0; i < mx.length; i++) {
			if (mx[i] != null) {
				for (int j = 0; j < mx[i].length; j++) {
					if (mx[i][j].index == columnIndex) {
						v.add(new VectorNode(i, mx[i][j].value));
						break;
					}
				}
			}
		}
		VectorNode[] result = new VectorNode[v.size()];
		v.toArray(result);
		return result;
	}

	/**
	 * Dot product of two sparse vectors.
	 *
	 * @param x the x vector
	 * @param y the y vector
	 * @return the dot product
	 */
	public static double dotProduct(VectorNode[] x, VectorNode[] y) {
		double sum = 0;
		if (x == null || y == null) {
			return 0;
		}
		int xlen = x.length;
		int ylen = y.length;
		int i = 0;
		int j = 0;
		while (i < xlen && j < ylen) {
			if (x[i].index == y[j].index)
				sum += x[i++].value * y[j++].value;
			else {
				if (x[i].index > y[j].index)
					++j;
				else
					++i;
			}
		}
		return sum;
	}

	
	/**
	 * Sparse matrix multiply with a transpose of the second sparse matrix.
	 *
	 * @param mx1 the first sparse matrix
	 * @param mx2 the second sparse matrix
	 * @return the product sparse matrix
	 */
	public static VectorNode[][] matrixMultiplyWithTranspose(VectorNode[][] mx1,
			VectorNode[][] mx2) {
		int m = mx1.length;
		if (findMaxColumnIndex(mx1) != findMaxColumnIndex(mx2)) {
			return null;
		}
		VectorNode[][] result = new VectorNode[m][];
		for (int i = 0; i < m; i++) {
			result[i] = null;
			for (int j = 0; j < m; j++) {
				double tmp = 0;
				tmp = dotProduct(mx1[i], mx2[j]);
				if (tmp != 0) {
					result[i] = addNode(result[i], j, tmp);
				}
			}
		}
		return result;
	}

	/**
	 * Reads a sparse matrix.
	 *
	 * @param filename the file name
	 * @return the sparse matrix
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
	public static VectorNode[][] readSparseMatrix(String filename)
			throws IOException {
		ArrayList<VectorNode[]> result = new ArrayList<VectorNode[]>();
		Scanner scn = new Scanner(new BufferedReader(new FileReader(filename))).useDelimiter("[\n\r]");
		int m = 0;
		while (scn.hasNext()) {
			String tmp = scn.next();
			if (tmp.length() > 0) {
				result.add(parseSparseVectorString(tmp));
				m++;
			}
		}
		scn.close();
		VectorNode[][] resultArray=new VectorNode[result.size()][];
		for (int i=0;i<result.size();i++){
			resultArray[i]=result.get(i);
		}
		return resultArray;
	}

	private static VectorNode[] parseSparseVectorString(String s){
		StringTokenizer st = new StringTokenizer(s, "[ :]");
		int nTokens=st.countTokens();
		if (nTokens%2!=0){
			st.nextToken();
			nTokens--;
		}
		VectorNode[] result = new VectorNode[nTokens / 2];
		int n = 0;
		while (st.hasMoreTokens()) {
			result[n] = new VectorNode();
			result[n].index = Integer.valueOf(st.nextToken());
			result[n].value = Double.valueOf(st.nextToken());
			n++;
		}
		return result;
	}
	
	/**
	 * Finds the maximum column index (often the dimension of the space).
	 *
	 * @param mx the sparse matrix
	 * @return the maximum column index
	 */
	public static int findMaxColumnIndex(VectorNode[][] mx) {
		int result = 0;
		for (int i = 0; i < mx.length; i++) {
			if (mx[i] != null) {
				for (int j = 0; j < mx[i].length; j++) {
					if (mx[i][j].index > result)
						result = mx[i][j].index;
				}
			}
		}
		return result;
	}

	
	/**
	 * Finds the minimum column index. This is important to determine
	 * whether the matrix elements are zero-indexed.
	 *
	 * @param mx the sparse matrix
	 * @return the minimum index
	 */
	public static int findMinColumnIndex(VectorNode[][] mx) {
		int result = Integer.MAX_VALUE;
		for (int i = 0; i < mx.length; i++) {
			if (mx[i] != null) {
				if (mx[i][0].index < result)
					result = mx[i][0].index;
			}
		}
		return result;
	}

	/**
	 * Shift columns to the right. This method is useful when
	 * converting from a zero-indexed matrix to a one-indexed matrix. 
	 *
	 * @param mx the sparse matrix
	 * @param k the number of shifts to the right 
	 * @return the shifted matrix
	 */
	public static VectorNode[][] shiftColumns(VectorNode[][] mx, int k) {
		VectorNode[][] result = new VectorNode[mx.length][];
		for (int i = 0; i < mx.length; i++) {
			result[i] = null;
			if (mx[i] != null) {
				result[i] = new VectorNode[mx[i].length];
				for (int j = 0; j < mx[i].length; j++) {
					result[i][j] = new VectorNode(mx[i][j].index + k,
							mx[i][j].value);
				}
			}
		}
		return result;
	}

	/**
	 * Transposes a sparse matrix.
	 *
	 * @param mx the sparse matrix
	 * @return the transposed matrix
	 */
	public static VectorNode[][] transpose(VectorNode[][] mx) {
		int minCol = findMinColumnIndex(mx);
		int adjust=0;
		if (minCol>0){
			adjust=1;
		}
		int mci = findMaxColumnIndex(mx)+1-adjust;
		VectorNode[][] result = new VectorNode[mci][];
		for (int i = 0; i < mx.length; i++) {
			if (mx[i] != null) {
				for (int j = 0; j < mx[i].length; j++) {
					result[mx[i][j].index-adjust] = addNode(result[mx[i][j].index-adjust], i,
							mx[i][j].value);
				}
			}
		}
		return result;
	}

}
