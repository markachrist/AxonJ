/** 
	Developer: Mark Christopher & Kasra Zarei, University of Iowa 
	Copyright - University of Iowa, 2014
**/ 

import java.awt.Color;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.measure.Measurements;
import ij.measure.ResultsTable;
import ij.measure.Calibration; 
import ij.plugin.ImageCalculator;
import ij.plugin.filter.ParticleAnalyzer;
import ij.plugin.frame.RoiManager;
import ij.process.BinaryProcessor;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.plugin.filter.GaussianBlur;

public class AxonJ_ implements PlugIn{
	
	/** Default Resolution/Scale from Anderson Scope
		Known Distance in Pixels */
	public double pixels = 373;
	
	/** Known Distance in Micrometers */
	public double microns = 20; 

	/** Hessian scales */
	public double hess_scales[] = {1.6}; 
	
	/** Min area of axon, in pixels^2 */
	public double min_pixels = 50;
	
	/** Min area of axon, in pixels^2 */
	public double max_pixels = 5000;
	
	/** Line width used to draw overlay on original image */
	public int overlayLineWidth = 1;
	
	/** Line width used to draw overlay on original image */
	public Color overlayColor = new Color(0, 255, 255);
	
	/** Half-width of filter used to compute heat map */
	public int h_width = 50;
	
	/** Heat map scale factor */
	public double h_scale = 0.5;
	
	/** Block size used for CLAHE local contrast adjustment */
	int contrastBlockSize = 127;
	
	/** Boolean indicating whether CLAHE local contrast adjustment should be performed */
	boolean runContrast = true;
	
	/** Option to generate overlay */
	boolean createOverlay = true;
	
	/** Option to generate heat map */
	boolean createHeatMap = false; 
	
	/** Dead and Sick Axons filtered out of total axon count */
	public int nonAxons = 0;
	
	/** The commands that can be called to run CLAHE contrast adjustment */
	String fijiContrast = "Enhance Local Contrast (CLAHE)";
	String ijContrast = "CLAHE ";
	String contrastCommand = ijContrast;
	
	/** The args that can be called to run CLAHE contrast adjustment */
	String fijiArgs = "blocksize=BLOCKSIZE histogram=256 maximum=3 mask=*None*";
	String ijArgs = "blocksize=BLOCKSIZE histogram=256 maximum=3";
	String contrastArgs = ijArgs;
	
	/**
	 * Counts axons of the currently selected image. Counts axons by averaging Hessian features
	 * over several scales to find cell boundaries. Result is then thresholded to produce binary
	 * cell boundary-background image. Resulting connected regions are then counted to estimate
	 * the number of axons in the image.
	 * 
	 * @param imp
	 *   Should be an optic nerve slice microscopy image.
	 */
	public void run(String args){
		
		ImagePlus p = WindowManager.getCurrentImage(); 
		ImageProcessor imp = p.getProcessor(); 
	
		int heat_w = (int)(h_scale * imp.getWidth() + 0.5);
		int heat_h = (int)(h_scale * imp.getHeight() + 0.5);
		FloatProcessor outputImage = new FloatProcessor(heat_w, heat_h);
		float axonIntensity = 0; 
		boolean contrast = this.runContrast;
		boolean continueProcessing = false; 
		
		//Convert to 32-bit and create accumulator
		ImagePlus working = new ImagePlus("Running AxonJ...", imp);
		ImageProcessor gray = imp.convertToFloat();
		ImagePlus sum = new ImagePlus("Sum", new FloatProcessor(imp.getWidth(), imp.getHeight()));
		working.setProcessor(gray); 
		working.show();
		
		Calibration cal = p.getCalibration(); 
		double pixWidth = 1.0 / cal.pixelWidth;
//		double pixHeight = 1.0 / cal.pixelHeight;
		
		if (pixWidth != 1.0){
			IJ.showMessage("pixels/micron: " + pixWidth); 
			
			/**
			 * Determine min and max axon sizes, hessian scales using entered resolution scale 
			 */
			min_pixels = 2.681 * pixWidth;
			max_pixels = 268.0965 * pixWidth; 
			
			for (int i = 0; i < hess_scales.length; i++)
			{
				hess_scales[i] = 0.08579 * pixWidth; 
			}
			
			continueProcessing = true; 
		}
		else if(this.createDialog()){						
			/**
			 * Determine min and max axon sizes, hessian scales using entered resolution scale 
			 */
			min_pixels = 2.681 * pixels / microns;
			max_pixels = 268.0965 * pixels / microns; 
			
			for (int i = 0; i < hess_scales.length; i++)
			{
				hess_scales[i] = 0.08579 * pixels / microns; 
			}
			
			continueProcessing = true; 
		}
			
		if (continueProcessing == true){
			WindowManager.setCurrentWindow(working.getWindow());
			
			//Run contrast adjustment
			if(contrast){
				contrast = this.adjustContrast(contrastBlockSize);
			}
				
			ImageCalculator calc = new ImageCalculator();
				
			//Compute Hessians at each indicated scale and sum them
			for(int i = 0; i < hess_scales.length; ++i){

				String hess_args = "largest smoothing=" + hess_scales[i];
				IJ.run("FeatureJ Hessian", hess_args);
					
				ImagePlus cur_hess = WindowManager.getCurrentImage();
				cur_hess.setTitle(hess_args);
					
				sum = calc.run("Add create 32-bit", cur_hess, sum);
					
				cur_hess.close();
				WindowManager.setCurrentWindow(working.getWindow());
			}
				
			//Convert sum to 8-bit and threshold it			
			working.close();
			sum.setProcessor(sum.getProcessor().convertToByte(true));
			sum.setTitle("Axons");
			sum.show();
						
			IJ.run("Auto Threshold", "method=Mean white");
			RoiManager rm = new RoiManager();
			
			//Count the particles (axons) in the resulting thresholed image
			ResultsTable resultsTable = new ResultsTable(); 
			BinaryProcessor bp = new BinaryProcessor((ByteProcessor)sum.getProcessor());
			sum.setProcessor(bp);
			int measure = Measurements.AREA | Measurements.CENTROID | Measurements.MODE | Measurements.STD_DEV | Measurements.FERET;
//			int partOpts = ParticleAnalyzer.SHOW_OVERLAY_OUTLINES | ParticleAnalyzer.EXCLUDE_EDGE_PARTICLES;
			int partOpts = ParticleAnalyzer.SHOW_OVERLAY_OUTLINES;
//			int partOpts = ParticleAnalyzer.SHOW_OVERLAY_MASKS;
			ParticleAnalyzer analyze = new ParticleAnalyzer(partOpts, measure, null, min_pixels, max_pixels, 0, 1.0);
			ParticleAnalyzer.setRoiManager(rm);
			ParticleAnalyzer.setResultsTable(resultsTable);
			analyze.analyze(sum);
			
			if(createOverlay){
				ImagePlus overlay = drawROIsOnImage(rm, p);
				overlay.show();
			}
			rm.close();
			
			resultsTable.show("Results");
				
			int xIndex = resultsTable.getColumnIndex("X");
			int yIndex = resultsTable.getColumnIndex("Y");
			int fIndex = resultsTable.getColumnIndex("MinFeret");
				
			double[] xCentroids = resultsTable.getColumnAsDoubles(xIndex);
			double[] yCentroids = resultsTable.getColumnAsDoubles(yIndex);
			double[] feretD = resultsTable.getColumnAsDoubles(fIndex);
//			double[] cellArea = resultsTable.getColumnAsDoubles(aIndex);
				
			int axonCount = resultsTable.getCounter(); 
					
			// Traverse through counted axons and sample from pixels to determine average intensity
			// note: this allows the darker, dying/sick axons to be excluded 
			for (int j = 0; j < axonCount; j++){
				double xCen = xCentroids[j];
				double yCen = yCentroids[j]; 
				double feret = feretD[j]; 
				int blackPixels = 0; 
				axonIntensity = 0; 
					
				// sample from pixels of axons
				for (double r = yCen - feret/2; r <= yCen + feret/2; r++){
					for (double c = xCen - feret/2; c <= xCen + feret/2; c++){
						if((Math.pow((double) r - yCen,2.0) + Math.pow((double) c-xCen,2.0) <= Math.pow(feret/2,2.0))){
							axonIntensity += imp.getPixelValue((int) c, (int) r);
							blackPixels++; 
						}
					}
				}
				axonIntensity /= blackPixels; 
				if (axonIntensity < 100){
					nonAxons++; 
				}
				else{
					/** If necessary, display "+" marks on counted axons **/ 
					//imp.putPixelValue((int) xCen, (int) yCen, 125);
					//imp.putPixelValue((int) xCen+1, (int) yCen, 125);
					//imp.putPixelValue((int) xCen-1, (int) yCen, 125);
					//imp.putPixelValue((int) xCen, (int) yCen+1, 125);
					//imp.putPixelValue((int) xCen, (int) yCen-1, 125);
				}
			}
			IJ.showMessage("Total axons: " + axonCount);
			
			// Create Heat Map if option is selected	
			if (createHeatMap){	
				for (int j = 0; j < axonCount; j++){
					double xCen = h_scale * xCentroids[j] + 0.5;
					double yCen = h_scale * yCentroids[j] + 0.5; 
					
					outputImage.putPixelValue((int) xCen, (int) yCen, 1.0); 
				}
				
				//Display cell centroids
				new ImagePlus("Cells map" , outputImage).show();
				
				//Create mask to count cells
				FloatProcessor mask = this.createCircularMask((int)(h_scale * h_width + 0.5));
				new ImagePlus("Count mask", mask).show();
				
				//Now do the counting via convolution
				outputImage.convolve((float[])mask.getPixels(), mask.getWidth(), mask.getHeight());
				//Compute min/max values so it displays correctly
				outputImage.resetMinAndMax();
				ImagePlus heatMap = new ImagePlus("Heat map", outputImage);
				
				// Apply Gaussian blur to generated axon heat map 
				GaussianBlur gaussian = new GaussianBlur(); 
				gaussian.blur(outputImage, 50); 
				/* user has to apply Fire look-up table manually */ 
				heatMap.show();			
			}
		}
	}
	
	/**
	 * Performs local contrast enhancement on the current image. Relies on the CLAHE method
	 * This is built into Fiji, but requires a plugin for vanilla ImageJ. This method
	 * attempts to run this method, but fails if it can't find it in the expected places.
	 * 
	 * Find info on CLAHE here:
	 * http://rsbweb.nih.gov/ij/plugins/clahe/index.html
	 * http://fiji.sc/Enhance_Local_Contrast_(CLAHE)
	 * 
	 * @return
	 *   True if contrast adjustment was performed, false otherwise 
	 */
	public boolean adjustContrast(int blockSize){
		
		boolean ran = true;
		String cmd, args;
		
		//Try using CLAHE
		try{
			cmd = this.contrastCommand;
			args = this.contrastArgs.replace("BLOCKSIZE", blockSize + "");
			IJ.run(cmd, args);
		}
		catch(Exception e1){
			//It didn't work, use alternate command
			try{
				cmd = this.fijiContrast;
				args = this.fijiArgs.replace("BLOCKSIZE", blockSize + "");
				IJ.run(cmd, args);
				
				//Remember that this was the correct command
				this.contrastCommand = fijiContrast;
				this.contrastArgs = fijiArgs;
			}
			catch(Exception e2){
				ran = false;
			}
		}
		
		return ran;
	}
	
	/**
	 * Creates a dialog used to confirm parameter values with the user. When/if the "OK"
	 * button is clicked, the analysis begans.
	 * 
	 * @return
	 *   True if the user entered valid params and clicked OK, false otherwise
	 */
	public boolean createDialog(){	
		GenericDialog d = new GenericDialog("AxonJ Setup");
		
		d.addMessage("Set the image scale values used for segmenting axons.");
		d.addNumericField("Known Distance in Pixels:", this.pixels, 2);
		d.addNumericField("Known Distance in Micrometers:", this.microns, 2);
		d.addCheckbox("Create Axon Overlay", createOverlay);
		d.addCheckbox("Create Axon Density Heat Map", createHeatMap);
		d.showDialog();
		
		if(!d.wasOKed()){
			return false;
		}
		
		double params[] = new double[2];
		
		if(!getParams(d, params)){
			IJ.showMessage("Error", "Missing parameter values!");
			return false;
		}
		
		this.pixels = params[0];
		this.microns = params[1];
		this.createOverlay = d.getNextBoolean();
		this.createHeatMap = d.getNextBoolean();
		
		return true;
	}

	/**
	 * Gets the parameter values input by the user.
	 * 
	 * @param d
	 *   Dialog from which parameter values are retrieved
	 * @param params
	 *   Output destination
	 * @return
	 *   True if all values were available, false otherwise 
	 */
	public boolean getParams(GenericDialog d, double params[]){	
		for(int i = 0; i < params.length; ++i){
			double n = d.getNextNumber();
			if(Double.isNaN(n)){
				return false;
			}
			params[i] = n;
		}
		return true;
	}
	
	/**
	 * Draws the set of ROIs controlled by RoiManager onto the given 
	 */
	public ImagePlus drawROIsOnImage(RoiManager rm, ImagePlus imp){
		Roi rois[] = rm.getRoisAsArray();
		
		ImageProcessor rgb = imp.getProcessor().convertToRGB();
		
		rgb.setLineWidth(overlayLineWidth);
		rgb.setColor(overlayColor);
		
		for(Roi r : rois){
			rgb.draw(r);
		}
		
		return new ImagePlus(imp.getTitle() + " axon outlines", rgb);
	}
	
	/**
	 * Converts comma-separated string of numbers to double array.
	 * 
	 * @param list
	 *   Comma separated list of numeric values (eg. 1,0.4,1e6)
	 * @return
	 *   Primitive double array containing values indicated by list
	 */
	public double[] listToDoubleArray(String list){	
		String data[] = list.split(",");
		double result[] = new double[data.length];
		
		for(int i = 0; i < data.length; ++i){
			result[i] = Double.parseDouble(data[i]);
		}
		
		return result;
	}
	
	/**
	 * Converts double array into comma-separated string of numbers
	 * 
	 * @param array
	 *   Array of double values to be converted
	 * @return
	 *   String containing comma-separated list of values in 
	 */
	public String doubleArrayToList(double array[]){	
		String list = "";
		
		for(int i = 0; i < array.length; ++i){
			list += array[i] + ",";
		}
		
		return list.substring(0, list.length() - 1);
	}
	
	/**
	 * Creates circular mask with the given radius. Pixels within mask have value of 1.0
	 * and all others have value of 0.0.
	 * 
	 * @param radius
	 *   Radius of circular mask
	 * @return
	 *   Float-valued, circular mask image with size of (2*radius + 1) x (2*radius + 1)
	 */
	public FloatProcessor createCircularMask(int radius){
		
		int r2 = radius*radius;
		FloatProcessor mask = new FloatProcessor(2*radius + 1, 2*radius + 1);
		
		for(int x = 0; x < mask.getWidth(); ++x){
			for(int y = 0; y < mask.getHeight(); ++y){
				if((x - radius)*(x - radius) + (y - radius)*(y - radius) <= r2){
					mask.putPixelValue(x, y, 1.0);
				}
			}
		}
		return mask;
	}
}
