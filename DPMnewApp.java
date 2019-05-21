package org.opensourcephysics.sip.wyatt;
import java.awt.*;
import java.awt.geom.Ellipse2D;
import java.awt.geom.AffineTransform;
import org.opensourcephysics.display.*;
import org.opensourcephysics.frames.*;
import org.opensourcephysics.numerics.*;
import java.util.stream.*;
import org.opensourcephysics.controls.*;
import org.opensourcephysics.display.GUIUtils;
import java.io.*;
import java.io.FileWriter; 
import java.util.List;
import java.util.ArrayList;
import java.util.Calendar; 

/**
 * DPMnewApp does a Monte Carlo simulation for a 2D mixture of fluctuating penetrable 
 * ellipses and hard disks (ellipses also translate and rotate)
 *
 * @author Alan Denton & Wyatt Davis based on code ch08/hd/HardDisksApp.java 
 * by Jan Tobochnik, Wolfgang Christian, Harvey Gould
 * @version 1.1 revised 21/05/19
 */

 
public class DPMnewApp extends AbstractSimulation {
	DPMnew hd = new DPMnew();
  	DisplayFrame display = new DisplayFrame("x", "y", "Hard Disk Particles");
  	int i = 0; // tracks number of steps taken
  	int j = 0; // tracks number of data points written 
  	int k = 0; // tracks number of independent runs 
  	Calendar cal = Calendar.getInstance();
  	ArrayList<Double> L1dat = new ArrayList<Double>(10000);
  	ArrayList<Double> L2dat = new ArrayList<Double>(10000);
  
  	/**
   	* Initializes the simulation 
   	*/
  	public void initialize() {
    		hd.Nd = control.getInt("Nd");
    		hd.Np = control.getInt("Np");
    		hd.phi = control.getDouble("phi");
    		hd.Lx = 0.5*Math.sqrt(hd.Nd*Math.PI/hd.phi);
    		hd.Ly = 0.5*Math.sqrt(hd.Nd*Math.PI/hd.phi);
    		hd.condition = false;
    		hd.tolerance = control.getDouble("tolerance");
    		hd.thetaTolerance = control.getDouble("theta tolerance");
    		hd.shapeTolerance = control.getDouble("shape tolerance");
    		String configuration = control.getString("initial configuration");
    		hd.sizeRatio = control.getDouble("size ratio");
                hd.penetrationCost= control.getDouble("penetration cost");
                hd.solventQuality = control.getString("solvent quality");
    		hd.initialize(configuration);
    		display.addDrawable(hd);
    		display.setPreferredMinMax(0, hd.Lx, 0, hd.Ly);
    		display.setSquareAspect(true);
    
  	}

  	//Dictates how data are taken over duration of simulation
  	public void doStep() {
		//if(k<6){ // 5 independent runs
		if(k<2){ // a single run
	  	if(k==0 && i==0){
			File dir1 = new File("data/DPM/phiN="+hd.phi);
		  	boolean successful1 = dir1.mkdir();
    			if (successful1){
      				System.out.println("directory was created successfully");
    			}
    			else{
      				System.out.println("failed trying to create the directory");
    			} 
		  	//for(int i=1; i<6; i++){ // 5 independent runs
		  	for(int i=1; i<2; i++){ // a single run
		  		File dir2 = new File("data/DPM/phiN="+hd.phi+"/"+i);
				boolean successful2 = dir2.mkdir();
    				if (successful2){
      					System.out.println("directory was created successfully");
    				}
    				else{
      					System.out.println("failed trying to create the directory");
    				} 
    				try{    
           				FileWriter fw=new FileWriter("data/DPM/phiN="+hd.phi+"/"+i+"/"+"phiN="+hd.phi+"_eX.dat");      
          			}catch(Exception e){System.out.println(e);}    
          			System.out.println("Success..."); 
    				try{    
           				FileWriter fw=new FileWriter("data/DPM/phiN="+hd.phi+"/"+i+"/"+"phiN="+hd.phi+"_eY.dat");       
          			}catch(Exception e){System.out.println(e);}    
          			System.out.println("Success...");    
			}
			k++;
		}
	  //check if enough data points have been taken
		if(j<1e6){ // number of data points (NOT number of steps)
			if(i<1e5){ // number of equilibration steps
				hd.step();
				i++;
			}
			else{ // data collection 
				hd.step(); 
				i++;
				//if (i%100 == 0){ // interval between samples
				if(i%10 == 0){
					if(hd.l1[0]>=hd.l2[0]){
                                        	L1dat.add(hd.l1[0]);
                                        	L2dat.add(hd.l2[0]);
                                	}
                                	else{
                                        	L1dat.add(hd.l2[0]);
                                        	L2dat.add(hd.l1[0]);
                                	}
					j++;
					//System.out.print(j);
				}
			}
		}
     		//WRITE DATA 
		else{
			String strFilePathX = "data/DPM/phiN="+hd.phi+"/"+k+"/"+"phiN="+hd.phi+"_eX.dat";
			String strFilePathY = "data/DPM/phiN="+hd.phi+"/"+k+"/"+"phiN="+hd.phi+"_eY.dat";
			System.out.println(strFilePathX); 
			try{    
          			FileWriter fwX = new FileWriter(strFilePathX);
		  		BufferedWriter bw = new BufferedWriter(fwX);
				for (double dat: L1dat){
					bw.write(dat+"");
					bw.newLine();
				}
				bw.close();  
          		}catch(Exception e){System.out.println(e);} 
          		try{    
          			FileWriter fwY = new FileWriter(strFilePathY);
				BufferedWriter bw = new BufferedWriter(fwY);
				for (double dat: L2dat){
					bw.write(dat+""); 
					bw.newLine();
				}
				bw.close();  
          		}catch(Exception e){System.out.println(e);}   
			k++;
        		L1dat.clear();
        		L2dat.clear();
        		i=0;
        		j=0;
        		initialize();
        	}
		}
	}

  /**
   * Resets the hard disks model to its default state.
   */
	public void reset() {
    		enableStepsPerDisplay(true);
    		control.setValue("Nd", 100);
    		control.setValue("Np", 1);
    		control.setValue("phi", 0.3);
   	        control.setValue("tolerance", 0.1);
   	        control.setValue("theta tolerance", 0.05);
   	        control.setValue("shape tolerance", 0.001); // about 10% of most probable principal radius
    		control.setValue("initial configuration", "square");
    		control.setValue("size ratio", 1);
                control.setValue("penetration cost", 0.5);
                control.setValue("solvent quality", "good");
                //control.setValue("solvent quality", "theta");
    		initialize();
  	}

  /**
   * Resets the hard disks model and the data graphs.
   *
   * This method is invoked using a custom button.
   */
	public void resetData() {
    		hd.resetAverages();
    		GUIUtils.clearDrawingFrameData(false); // clears old data from the plot frames
  	}

  /**
   * Starts the Java application.
   * @param args  command line parameters
   */
	public static void main(String[] args) { // set up animation control structure using this class
    		SimulationControl control = SimulationControl.createApp(new DPMnewApp());
    		control.addButton("resetData", "Reset Data");
	}
}
