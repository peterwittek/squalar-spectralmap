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

import java.io.IOException;
import java.util.logging.Logger;

import org.squalar.util.DenseVector;
import org.squalar.util.SparseVector;
import org.squalar.util.VectorNode;

import ch.akuhn.edu.mit.tedlab.SMat;
import ch.akuhn.edu.mit.tedlab.SVDRec;
import ch.akuhn.edu.mit.tedlab.Svdlib;

/**
 * Command line utility for creating a co-occurrence matrix and its
 * singular value decomposition.
 */

public class Decompose {
	
	/** The logger. */
	public static Logger logger = Logger.getLogger("org.squalar.spectralmap");
	
	/** The number of singular values to calculate in the decomposition */
	private static int nSingularValues = 1000;
	
	
	/**
	 * The main method. 
	 *
	 * @param args[0] the libsvm-formatted term-document matrix file name
	 * @param args[1] the co-occurrence file name
	 * @param args[2] the left singular vectors file name
	 * @param args[3] the right singular vector file name
	 * @param args[4] the singular values file name
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
	public static void main(String[] args) throws IOException {
	    if (args.length != 5) {
	      throw (new IllegalArgumentException("There were " + args.length + " arguments, instead of the expected 5."));
	    }
		String tdMatrixFilename=args[0];
		String cooccurFilename=args[1];
		String leftSingularVectorsFilename=args[2];
		String rightSingularVectorsFilename=args[3];
		String singularValuesFilename=args[4];

		VectorNode[][] mx = SparseVector.transpose(SparseVector.readSparseMatrix(tdMatrixFilename));
		logger.info("Calculating co-occurrence matrix");
		mx=SparseVector.matrixMultiplyWithTranspose(mx, mx);
		DenseVector.writeDenseMatrix(DenseVector.convertFromSparseMatrix(mx), cooccurFilename);
		SVDRec svdResult=decompose(mx);
		DenseVector.writeDenseMatrix(svdResult.Ut.value,leftSingularVectorsFilename);
		DenseVector.writeDenseMatrix(svdResult.Vt.value,rightSingularVectorsFilename);
		DenseVector.writeDenseMatrix(new double[][]{ svdResult.S} ,singularValuesFilename);
	}

	/**
	 * Decompose.
	 *
	 * @param mx the rectangular matrix to be decomposed 
	 * @return the SVD decompositon of the matrix
	 */
	private static SVDRec decompose(VectorNode[][] mx) {
		logger.info("Converting matrix");
		SMat smx=convertLibsvmToSvdlibj(mx);
		logger.info("Starting SVD using algorithm LAS2");
		Svdlib svd = new Svdlib();
		return svd.svdLAS2A(smx, nSingularValues);
	}

	/**
	 * Converts a libsvm sparse matrix to a svdlibj sparse one
	 *
	 * @param mx the libsvm sparse matrix
	 * @return the svdlibj matrix
	 */
	private static SMat convertLibsvmToSvdlibj(VectorNode[][] mx) {
		SMat S;
		int i, j, n;
		int maxColIndex = 0;
		// n = number of non-zero elements
		for (i = 0, n = 0; i < mx.length; i++) {
			for (j = 0; j < mx[i].length; j++) {
				n++;
				if (mx[i][j].index > maxColIndex) {
					maxColIndex = mx[i][j].index;
				}
			}
		}
		S = new SMat(mx.length, maxColIndex, n);
		for (j = 0, n = 0; j < maxColIndex+1; j++) {
			VectorNode[] column = SparseVector.columnVector(mx, j);
			S.pointr[j] = n;
			for (i = 0; i < column.length; i++) {
				S.rowind[n] = column[i].index;
				S.value[n] = column[i].value;
				n++;
			}
		}
		S.pointr[S.cols] = S.vals;
		return S;
	}

}
