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

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

import org.squalar.util.DenseVector;

/**
 * Utility class to map a term's eigenvalue to
 * the visible range and draw the term's spectrum.
 */

public class GenerateVisibleSpectrum{
	
	private final static String path = "figures/linear-spectrum.png";
	private final static int minVisible=400;
	private final static int maxVisible=700;
		
	/**
	 * Draws a term spectrum.
	 *
	 * @param singularValues the singular values
	 * @param singularVectorIndices the indices of singular vectors to which 
	 *        the given term is similar
	 * @param term the term for which the spectrum is to be drawn
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
	public static void drawTermSpectrum(double[] singularValues, int[] singularVectorIndices, String term) throws IOException {
		double[] spectrum=scale(singularValues);
		double[] subSpectrum=new double[singularVectorIndices.length];
		for (int i=0;i<subSpectrum.length;i++){
			subSpectrum[i]=spectrum[singularVectorIndices[i]];
		}
		draw(subSpectrum,"figures/spectrum_"+term+".png");
	}
	
	private static double[] scale(double[] singularValues){
		double[] result=new double[singularValues.length];
		double max=DenseVector.max(singularValues);
		for (int i=0;i<singularValues.length;i++){
			result[i]=minVisible+singularValues[i]*(maxVisible-minVisible)/max;
		}
		return result;
	}
	
	private static void draw(double[] frequencies, String outputFilename) throws IOException{
		BufferedImage image = ImageIO.read(new File(path));
		int width=image.getWidth();
		int height=image.getHeight();
		int stripeWidth = width / 370;
		BufferedImage bufferedImage = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);

		Graphics2D g=bufferedImage.createGraphics();
		g.setColor(Color.black);
		g.fillRect(0, 0, width, height);
		for (int i = 0; i < frequencies.length; i++) {
			int x = stripeWidth*((int)frequencies[i]-minVisible);
			g.drawImage(image, 
					x, 0, x + stripeWidth, height, 
					x, 0, x + stripeWidth, height,
					null);
		}
		g.dispose();
		ImageIO.write(bufferedImage, "png", new File(outputFilename));		
	}
}
