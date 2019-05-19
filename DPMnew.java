package org.opensourcephysics.sip.wyatt;
import java.awt.*;
import java.awt.geom.Ellipse2D;
import java.awt.geom.AffineTransform;
import org.opensourcephysics.display.*;
import org.opensourcephysics.frames.*;
import org.opensourcephysics.numerics.*;
import java.util.stream.*;
import java.util.Arrays;

/**
 * DPMnew does a Monte Carlo simulation for a 2D mixture of fluctuating penetrable ellipses 
 * and hard disks (ellipses do not translate or rotate)
 *
 * @author Alan Denton & Wyatt Davis based on OSP code by Jan Tobochnik, Wolfgang Christian, Harvey Gould
 * @version 1.1 revised 19/05/19
 */
 

public class DPMnew implements Drawable {
	public double x[], y[];
        public double xp[], yp[];
  	public double a[], b[];
  	public double l1[], l2[];
        public double theta[];
  	public double sizeRatio;
        public double penetrationCost;
        public double al1 = 1.3690*10001/10000/Math.PI/Math.PI;
        public double al2 = 1.0972*10001/10000/2/2/Math.PI/Math.PI;
        public double v1 = 1.8769;
        public double v2 = 3.1012;
        public double b1 = 3.62051;
        public double b2 = 2.75939;
        public double c1 = 49.3780;
        public double c2 = 256.564;
        public int over[][];
        public int overtest[][];
  	public double prob; 
  	public int Nd;
  	public int Np;
        String solventQuality;
  	public double phi;
  	public double Lx, Ly;
  	public double steps = 0;
  	public double tolerance;
  	public double shapeTolerance;
  	public boolean condition;

  	public void initialize(String configuration) {
        	resetAverages();
    		x = new double[Nd];
    		y = new double[Nd];
    		xp = new double[Np];
    		yp = new double[Np];
    		a = new double[Np];
    		b = new double[Np];
    		l1 = new double[Np];
    		l2 = new double[Np];
                over = new int[Np][Nd];
                overtest = new int[Np][Nd];
		for(int i=0; i<Np; ++i){
			for(int j=0; j<Nd; ++j){
				over[i][j] = 0;
				overtest[i][j] = 0;
			}
		}
                theta = new double[Np];
                theta[0] = 0;
    		if(configuration.equals("square")) {
      			setSquarePositions();
    		} 
		else {
      			setRandomPositions();
    		}
   	}

	public void resetAverages() {
        	steps = 0;
  	}
  
	public void setRandomPositions() {
        	boolean overlap;
    		for(int i=0; i<Nd; ++i) {
      			do {
        			overlap = false;
        			x[i] = Lx*Math.random();
        			y[i] = Ly*Math.random();
        			int j = 0;
        			while(j<i && !overlap) {
          				double dx = PBC.separation(x[i]-x[j], Lx);
          				double dy = PBC.separation(y[i]-y[j], Ly);
          				if(dx*dx+dy*dy < 1) {
            					overlap = true;
          				}
          				j++;
        			}
      			} while(overlap);
    		}
  	}
	
	// set initial positions of disks and ellipses
	public void setSquarePositions() {
    		double dnx = Math.sqrt(Nd);
    		int nx = (int) dnx;
    		if(dnx-nx > 0.00001) {
      			nx++; // N is not a perfect square
    		}
    		double ax = Lx/nx; // distance between columns of disks
    		double ay = Ly/nx; // distance between rows
    		int i = 0;
    		int iy = 0;
    		while(i < Nd) { // place disks on square lattice
      			for(int ix = 0; ix<nx; ++ix) { // loops through disks in a row
        			if(i < Nd) {
            				x[i] = ax*ix;
          				y[i] = ay*iy;
          				i++;
        			}
      			}
    			iy++;
    		}
    
    		boolean overlap;
		for(int k=0; k<Np; ++k){ // place ellipses at interstitial sites
			do{
				overlap = false;
				xp[k] = (nx/2+0.5)*ax;
				yp[k] = (nx/2+0.5)*ay;
				a[k] = 0.2;
				b[k] = 0.2;

				if (solventQuality.equals("theta")){
                                        l1[k] = Math.pow(a[k]/sizeRatio/.5,2)/12;
                                        l2[k] = Math.pow(b[k]/sizeRatio/.5,2)/12;
                                }
                                else{
                                        l1[k] = Math.pow(a[k]/sizeRatio/2.15016,2);
                                        l2[k] = Math.pow(b[k]/sizeRatio/2.15016,2);
                                }

				i = 0;
				while(!overlap && i < Nd){
					double dx=PBC.separation(xp[k]-x[i], Lx);
					double dy=PBC.separation(yp[k]-y[i], Ly);
					if (dx*dx+dy*dy < (a[k]+0.5)*(a[k]+0.5)){
						overlap = true;
					}
					i++;
				}
			} while(overlap); 
		}
	}

        //Probability distribution for polymer shape in a theta solvent
        public void shapeProb_theta(double L1O, double L2O, double L1N, double L2N){

                if(L1N < L2N){ // swap eigenvalues
			double temp = L1O;
			L1O = L2O;
			L2O = temp;
			temp = L1N;
			L1N = L2N;
			L2N = temp;
		}

                prob = Math.pow(L1N/L1O,v1-1)*Math.pow(L2N/L2O,v2-1)*Math.exp(-v1*(L1N-L1O)/al1-v2*(L2N-L2O)/al2);
        }

        //Probability distribution for polymer shape in a good solvent
        public void shapeProb_good(double L1O, double L2O, double L1N, double L2N){

                if(L1N < L2N){ // swap eigenvalues
			double temp = L1O;
			L1O = L2O;
			L2O = temp;
			temp = L1N;
			L1N = L2N;
			L2N = temp;
		}

		prob = Math.pow(L1N/L1O,b1)*Math.pow(L2N/L2O,b2)*Math.exp(-c1*(L1N-L1O)-c2*(L2N-L2O));
        }

         //Make B a copy of array A (A and B must be of the same dimensions) 
        public void arrayCopy(int A[][], int B[][]){
                for(int i=0; i<A.length; i++)
                        for(int j=0; j<A[i].length; j++)
                                B[i][j]=A[i][j];
        }

	//Determines whether overlap occurs between disk and ellipse 
  	public void Overcond(int d, int p){ // disk (d) and polymer (p) indices
	  	// Step 1: Coordinate transform
		condition = false; 
		double Tx = PBC.separation(x[d]-xp[p], Lx); // no rotations
		double Ty = PBC.separation(y[d]-yp[p], Ly); 
		// Step 2: Calculate coefficients of polynomial
		double A = Tx*Tx/a[p]/a[p];
		double B = b[p]*b[p]*Ty*Ty/a[p]/a[p]/a[p]/a[p];
		double C = b[p]*b[p]/a[p]/a[p] - 1;
		// Step 3: Find roots of polynomial
	 	double coef[] = new double[5];
                coef[0] = -1.;
                coef[1] = -2.*C;
                coef[2] = A+B-C*C;
                coef[3] = 2.*A*C;
                coef[4] = A*C*C;
	  	Polynomial polynomial = new Polynomial(coef);
		double[] realroots;
                realroots = polynomial.rootsReal();
                int nRoots = realroots.length;
                for(int i=0; i<nRoots; ++i){ // check each root for overlap
			if(!condition){
                        	double xroot = Tx*realroots[i];
                        	double yroot = b[p]*b[p]*Ty*xroot/a[p]/a[p]/(Tx+(b[p]*b[p]/a[p]/a[p]-1)*xroot);
                        	if((xroot-Tx)*(xroot-Tx)+(yroot-Ty)*(yroot-Ty) < 0.25){
                                	condition = true;
                        	}
                        }
                }
	}
 
	// Monte Carlo step 
	public void step() {
    		boolean overlap;
    		double dxtrial, dytrial;
    		double dl1trial, dl2trial;
    		double L;
		
                // make trial moves and check for overlaps
                for(int i=0; i<Nd; ++i) {

                        if(Math.random() > 0.1){ // try to move a disk 

      				overlap = false;
      				dxtrial = tolerance*2.*(Math.random()-0.5);
      				dytrial = tolerance*2.*(Math.random()-0.5);
      				x[i] = PBC.position(x[i]+dxtrial, Lx); 
      				y[i] = PBC.position(y[i]+dytrial, Ly); 

				for(int j=0; j<Nd; ++j) { // check for overlaps with other disks
        				if(j!=i && !overlap) {
          					double dx = PBC.separation(x[i]-x[j], Lx);
          					double dy = PBC.separation(y[i]-y[j], Ly);
          					if(dx*dx+dy*dy < 1) { // reject disk move
            						overlap = true; 
            						x[i] = PBC.position(x[i]-dxtrial, Lx); 
            						y[i] = PBC.position(y[i]-dytrial, Ly); 
          					}
        				}
      				} // end of loop over other disks

				if(overlap){
					continue; // next move
				}

      				for(int k=0; k<Np; k++){ // check for overlaps with polymers
      					condition = false;
                                        L = Math.max(a[k], b[k]);
                                        double dxp = PBC.separation(x[i]-xp[k], Lx); // no rotations
                                        double dyp = PBC.separation(y[i]-yp[k], Ly);
			  		if(dxp*dxp+dyp*dyp < (L+0.5)*(L+0.5)){ // overlap candidate
                                                if(dxp*dxp/(a[k]+0.5)/(a[k]+0.5)+dyp*dyp/(b[k]+0.5)/(b[k]+0.5) < 1){ // disk inside fattened ellipse
                                                        condition = true; 
							overtest[k][i] = 1;
							continue; // go to next polymer
                                                }					
						Overcond(i,k); // check for disk-ellipse overlap 
			  		}
					if(condition){ 
						overtest[k][i] = 1;
					}
					else{
						overtest[k][i] = 0;
					}
		   		} // end of loop over polymers

				double rand = Math.random();
				if(rand < Math.exp(-penetrationCost*(IntStream.of(overtest[0]).sum()-IntStream.of(over[0]).sum()))){ // accept disk move
                               		arrayCopy(overtest,over);
                        	}
                        	else{ // reject disk move
                               		arrayCopy(over,overtest); 
                               		x[i] = PBC.position(x[i]-dxtrial, Lx); 
                               		y[i] = PBC.position(y[i]-dytrial, Ly);
                        	}
			} // end of disk move

			else{ // make a polymer move

			// move polymers and check for overlaps with disks
      				for(int k=0; k<Np; k++){ // loop over polymers

					dl1trial = 5.*shapeTolerance*2.*(Math.random()-0.5);
					dl2trial = shapeTolerance*2.*(Math.random()-0.5);
					l1[k] += dl1trial;
					l2[k] += dl2trial;

					if(l1[k] < 0 || l2[k] < 0){ // reject polymer move
						l1[k] -= dl1trial;
						l2[k] -= dl2trial;
						continue; // next polymer
					}

 					if(solventQuality.equals("theta")){
						shapeProb_theta(l1[k]-dl1trial, l2[k]-dl2trial, l1[k], l2[k]);
                                		a[k] = .5*sizeRatio*Math.pow(l1[k]*12,.5);
                                		b[k] = .5*sizeRatio*Math.pow(l2[k]*12,.5);
					}
					else{
						shapeProb_good(l1[k]-dl1trial, l2[k]-dl2trial, l1[k], l2[k]);
                                		a[k] = 2.15016*sizeRatio*Math.pow(l1[k],.5);
                                		b[k] = 2.15016*sizeRatio*Math.pow(l2[k],.5);
					}

                        		L = Math.max(a[k], b[k]);

		        		for(int j=0; j<Nd; ++j){ // check for overlaps with disks
						condition = false;
                                        	double dxp = PBC.separation(x[j]-xp[k], Lx); // no rotations
                                        	double dyp = PBC.separation(y[j]-yp[k], Ly);
						if(dxp*dxp+dyp*dyp < (L+0.5)*(L+0.5)){ // overlap candidate
                                               		if(dxp*dxp/(a[k]+0.5)/(a[k]+0.5)+dyp*dyp/(b[k]+0.5)/(b[k]+0.5) < 1){ // overlap: disk inside fattened ellipse
                                                       		condition = true; 
								overtest[k][j] = 1;
								continue; // go to next disk
                                               		}
							Overcond(j,k);  
						} 
						if(condition){ 
							overtest[k][j] = 1;
						}
						else{
							overtest[k][j] = 0;
						}
					} // end of loop over disks

					double rand = Math.random();
					if(rand < prob*Math.exp(-penetrationCost*(IntStream.of(overtest[0]).sum()-IntStream.of(over[0]).sum()))){ // accept polymer move
						arrayCopy(overtest,over);
					}
					else{ // reject polymer move 
						arrayCopy(over,overtest); 
						l1[k] -= dl1trial;
						l2[k] -= dl2trial;
						if(solventQuality.equals("theta")){
                                        		a[k] = .5*sizeRatio*Math.pow(l1[k]*12,.5);
                                        		b[k] = .5*sizeRatio*Math.pow(l2[k]*12,.5);
                                		}
                                		else{
                                        		a[k] = 2.15016*sizeRatio*Math.pow(l1[k],.5);
                                        		b[k] = 2.15016*sizeRatio*Math.pow(l2[k],.5);
                                		}
					}
				}
	    		}
		}
	}
	
	// Draw the mixture
	public void draw(DrawingPanel drawingPanel, Graphics g) {
    		double radius = 0.5;
    		int pxRadiusP[] = new int[Np];
    		int pyRadiusP[] = new int[Np];
    		Graphics2D g2 = (Graphics2D) g; // New
    		if(x==null) {
      			return;
    		}
    		int pxRadius = Math.abs(drawingPanel.xToPix(radius)-drawingPanel.xToPix(0));
    		int pyRadius = Math.abs(drawingPanel.yToPix(radius)-drawingPanel.yToPix(0));
	        g.setColor(Color.blue);
    		for(int i=0; i<Nd; i++) {
      			int xpix = drawingPanel.xToPix(x[i])-pxRadius;
      			int ypix = drawingPanel.yToPix(y[i])-pyRadius;
      			g.fillOval(xpix, ypix, 2*pxRadius, 2*pyRadius);
    		} 
    		g.setColor(Color.red);
    		for(int i=0; i<Np; i++) {
      			pxRadiusP[i] = Math.abs(drawingPanel.xToPix(a[i])-drawingPanel.xToPix(0)); // Alan
      			pyRadiusP[i] = Math.abs(drawingPanel.yToPix(b[i])-drawingPanel.yToPix(0)); // Alan
      			int xpixP = drawingPanel.xToPix(xp[i])-pxRadiusP[i];
      			int ypixP = drawingPanel.yToPix(yp[i])-pyRadiusP[i];
      			Ellipse2D e = new Ellipse2D.Double(xpixP, ypixP, 2*pxRadiusP[i], 2*pyRadiusP[i]);
      			AffineTransform at = AffineTransform.getRotateInstance(theta[i], xpixP+pxRadiusP[i], ypixP+pyRadiusP[i]);
      			g2.fill(at.createTransformedShape(e));
  		} 
    		g.setColor(Color.black);
    		int xpix = drawingPanel.xToPix(0);
    		int ypix = drawingPanel.yToPix(Ly);
    		int lx = drawingPanel.xToPix(Lx)-drawingPanel.xToPix(0);
    		int ly = drawingPanel.yToPix(0)-drawingPanel.yToPix(Ly);
    		g.drawRect(xpix, ypix, lx, ly);
  	}
}
