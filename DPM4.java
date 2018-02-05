package org.opensourcephysics.sip.ch15;
import java.awt.*;
import java.awt.geom.Ellipse2D;
import java.awt.geom.AffineTransform;
import org.opensourcephysics.display.*;
import org.opensourcephysics.frames.*;
import org.opensourcephysics.numerics.*;
import java.util.stream.*;
import java.util.Arrays;

/**
 * Does a Monte Carlo simulation for a 2D system of hard disks and polymer
 *
 *@author Alan Denton & Wyatt Davis based on MD code ch08/hd/HardDisksApp.java by Jan Tobochnik, Wolfgang Christian, Harvey Gould
 *@version 4.0 revised 02/5/18 
 */
 
 
 ////////GOAL/////////
 /**
  * Incorporate ordered eigenvalues. 
*/
public class DPM4 implements Drawable {
  public double x[], y[];
  public double xp[],yp[];
  public double a[],b[];
  public double l1[],l2[];
  public double theta[];
  public double sizeRatio;
  public double penetrationCost;
  public double al1 = 1.3690*10001/10000/Math.PI/Math.PI;
  public double al2 = 1.0972*10001/10000/2/2/Math.PI/Math.PI;
  public double v1 = 1.8769;
  public double v2 = 3.1012;
  public int over[][];
  public int overtest[][];
  public double prob; 
  public int Nd;
  public int Np;
  public double phi;
  public double Lx,Ly;
  public double keSum = 0, virialSum = 0;
  public double steps = 0;
  public double tolerance;
  public double thetolerance = .05;
  public double l1tolerance = .01;
  public double l2tolerance = .001;
  public boolean condition = false;
  // NEW, when true the nearest point distance calculation is valid. When false (more than two real roots), the calculation is invalid and the move needs 
  //to be skipped.
  public boolean tworoots = true;
  //public int d = 0;
  //public int p = 0;
  // end break
  // start break
  // initialize
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
    theta = new double[Np];
    theta[0] = 0.; 
    if(configuration.equals("regular")) {
      setRegularPositions();
    } else {
      setRandomPositions();
    }
  }
  // end break

  // start break
  // reset averages
  public void resetAverages() {
    steps = 0;
    virialSum = 0;
  }
  // end break

  // start break
  // set particle coordinates
  public void setRandomPositions() {
    boolean overlap;
    for(int i = 0;i<Nd;++i) {
      do {
        overlap = false;
        x[i] = Lx*Math.random();
        y[i] = Ly*Math.random();
        int j = 0;
        while((j<i)&&!overlap) {
          double dx = PBC.separation(x[i]-x[j], Lx);
          double dy = PBC.separation(y[i]-y[j], Ly);
          if(dx*dx+dy*dy<1.0) {
            overlap = true;
          }
          j++;
        }
      } while(overlap);
    }
  }

	
  public void setRegularPositions() {
    double dnx = Math.sqrt(Nd);
    int nx = (int) dnx;
    if(dnx-nx>0.00001) {
      nx++; // N is not a perfect square
    }
    double ax = Lx/nx; // distance between columns of disks
    double ay = Ly/nx; // distance between rows
    int i = 0;
    int iy = 0;
    while(i<Nd) {
      for(int ix = 0;ix<nx;++ix) { // loops through disks in a row
        if(i<Nd) {
          y[i] = ay*(iy+0.5);
          if(iy%2==0) {            // even rows displaced from odd rows
            x[i] = ax*(ix+0.25);
          } else {
            x[i] = ax*(ix+0.75);
          }
          i++;
        }
      }
      iy++;
    }
    
    
    boolean overlap;
	for(int j = 0; j<Np; ++j){
		do{
			overlap = false;
			xp[j] = Lx*Math.random();
			yp[j] = Ly*Math.random();
			a[j] = .5;
			b[j] = .5;
			l1[j] = Math.pow(a[j]/sizeRatio/.5,2)/12;
			l2[j] = Math.pow(b[j]/sizeRatio/.5,2)/12;
			int k = 0;
			while(overlap==false&&k<Nd){
				double dx=PBC.separation(xp[j]-x[k], Lx);
				double dy=PBC.separation(yp[j]-y[k], Ly);
				if (dx*dx+dy*dy<(Math.pow(.5+.5,2))){
					overlap = true;
				}
				k++;
			}
		} while(overlap); 
	}
					
					//Continue with random placement 
				
  }
  // end break
  
  // Incorperate ordered eigenvalues here for ideal polymer
  public void sciutto(int p, double dl1, double dl2){
	  double L1O = 0; double L2O = 0;
	  double L1N = 0; double L2N = 0;
	  
	  if(l1[p]>=l2[p]){
		  L1O = l1[p]; L1N = L1O + dl1;
		  L2O = l2[p]; L2N = L2O + dl2;
		  
	  }
	  else if(l1[p]<l2[p]){
		  L1O = l2[p]; L1N = L1O + dl2;
		  L2O = l1[p]; L2N = L2O + dl1;
	  }
	  
	  prob = Math.pow((L1N)/L1O,v1-1)*Math.pow((L2N)/L2O,v2-1)*Math.exp(-(v1*(L1N-L1O)/al1)-(v2*(L2N-L2O)/al2));	  
  }
  
  //Make B a copy of array A (A and B must be of the same dimensions) 
  public void arrayCopy(int A[][], int B[][]){
	  for(int i=0; i<A.length; i++)
        for(int j=0; j<A[i].length; j++)
         B[i][j]=A[i][j];
      }   
  
  // Determines whether a crowder and polymer are overlapping
  public void Overcond(int d, int p){
	  // Step 1: Coordinate Transform
	  tworoots = true;
	  condition = false; 
	  double Txa = PBC.separation(x[d]-xp[p], Lx);
	  double Tya = PBC.separation(y[d]-yp[p], Ly); 
	  double Tx = Math.cos(theta[p])*Txa-Math.sin(theta[p])*Tya;
	  double Ty = Math.sin(theta[p])*Txa+Math.cos(theta[p])*Tya;
	  
	  // Step 2: Calcualte Coefficients of Polynomial
	  double A = Tx*Tx/a[p]/a[p];
	  double B = b[p]*b[p]*Ty*Ty/a[p]/a[p]/a[p]/a[p];
	  double C = (b[p]*b[p]/a[p]/a[p]) - 1;
	  
	  // Step 3: Find Roots of Polynomial  
	  double coef[] = new double[5]; 
	  coef[0] = A*C*C;
	  coef[1] = 2*A*C;
	  coef[2] = A+B-(C*C);
	  coef[3] = -2*C;
	  coef[4] = -1;
	  
	  double roots[][] = Root.polynomial(coef);
	  double realroot1 = 0; 
	  double realroot2 = 0;
	  
	// start of Alan's edits
	  double realroots[] = new double[4];
          int nReal = 0;
          for(int i = 0; i < 4; i++){
	     if(roots[1][i]*roots[1][i] < 1.e-40){ 
                realroots[nReal] = roots[0][i];
                nReal++;
             }
          }
          if(nReal != 2){
             tworoots = false;
          }                        
          else{
             realroot1 = realroots[0];
             realroot2 = realroots[1];
	     double xone = (Tx/realroot1);
	     double xtwo = (Tx/realroot2);
	     double yone = b[p]*b[p]*Ty*xone/a[p]/a[p]/(Tx+((b[p]*b[p]/a[p]/a[p])-1)*xone);
	     double ytwo = b[p]*b[p]*Ty*xtwo/a[p]/a[p]/(Tx+((b[p]*b[p]/a[p]/a[p])-1)*xtwo);
	 
	     if((xone-Tx)*(xone-Tx)+(yone-Ty)*(yone-Ty) < .25){
		condition = true;
	     }
	     else if((xtwo-Tx)*(xtwo-Tx)+(ytwo-Ty)*(ytwo-Ty) < .25){
		condition = true;
	     }
          } 
	// end of Alan's edits
}


	  
	  
  // PROBLEMS //
/**
 * 
 * Seems to ignore overlap condition. 
 * 
*/
     ///// TO DO IN THIS METHOD /////
  // Incorporate ordered eigenvalues // 
  public void step() {
// make trial displacements and check for overlap
	tworoots = true;
    boolean overlap;
    double dxtrial, dytrial;
    double l1trial, l2trial;
    double L;
    double dtheta; 
    
    for(int i = 0;i<Nd;++i) {
      overlap = false;
      tworoots = true;
      // Wherever overlap = true, make sure to change overtest.
      // The 
      dxtrial = tolerance*2.*(Math.random()-0.5);
      dytrial = tolerance*2.*(Math.random()-0.5);
      x[i] = PBC.position(x[i]+dxtrial, Lx); 
      y[i] = PBC.position(y[i]+dytrial, Ly); 
      // Leave this loop unchanged 
      for(int j = 0;j<Nd;++j) {
        if((j!=i)&&!overlap) {
          double dx = PBC.separation(x[i]-x[j], Lx);
          double dy = PBC.separation(y[i]-y[j], Ly);
          if(dx*dx+dy*dy<1.0) {
            overlap = true; 
            
            x[i] = PBC.position(x[i]-dxtrial, Lx); // reject displacement 
            y[i] = PBC.position(y[i]-dytrial, Ly); 
          }
        }
      }
     
      condition = false;
      tworoots = true;
      for(int k = 0; k<Np; k++){
		  
		 if(overlap != true && tworoots!=false){
			 L = a[k];
			  if(b[k]>=a[k]){
				L = b[k];
			}
			
			if(over[k][i]==1){
				Overcond(i,k);
				if(tworoots==true && condition==false){
					overtest[k][i] = 0;
				}
			}
			    double Txa = PBC.separation(x[i]-xp[k],Lx);
				double Tya = PBC.separation(y[i]-yp[k],Ly);
				double dxp = Math.cos(theta[k])*Txa-Math.sin(theta[k])*Tya;
				double dyp = Math.sin(theta[k])*Txa+Math.cos(theta[k])*Tya; 
			  if(dxp*dxp + dyp*dyp<(L+.5)*(L+.5)){
				  Overcond(i,k);
				  if(tworoots==false){
					  x[i] = PBC.position(x[i]-dxtrial, Lx); // reject displacement 
					  y[i] = PBC.position(y[i]-dytrial, Ly);
				  }
				  if(condition == true && tworoots==true){
					   overtest[k][i] = 1;
					   
			  }
				  if(condition == false && tworoots==true){
					  overtest[k][i]=0; 
				  }
				  if((dxp*dxp/a[k]/a[k])+(dyp*dyp/b[k]/b[k])<1){
					  overtest[k][i]=1;
					  
					  }
			  
			  }
		   }
      }
      // rng here. 
	 
	  double rand = Math.random();
	  if(tworoots){
	  if(rand<Math.exp(-penetrationCost*(IntStream.of(overtest[0]).sum()-IntStream.of(over[0]).sum())) && tworoots==true){
						  arrayCopy(overtest,over);
					  }
					  else if(rand>Math.exp(-penetrationCost*(IntStream.of(overtest[0]).sum()-IntStream.of(over[0]).sum())) && tworoots==true){
						  arrayCopy(over,overtest);
						  x[i] = PBC.position(x[i]-dxtrial, Lx); // reject displacement 
						  y[i] = PBC.position(y[i]-dytrial, Ly);
					  }
					  
  }
}
 
  
 
      condition = false;
      
      
      // The condition for accepting a move change 
      for(int k = 0; k<Np; k++){
			  condition = false;
			  tworoots = true;
			  dxtrial = tolerance*2.*(Math.random()-0.5);
			  dytrial = tolerance*2.*(Math.random()-0.5);
			  l1trial = l1tolerance*2.*(Math.random()-0.5);
			  l2trial = l2tolerance*2.*(Math.random()-0.5);
			  dtheta = thetolerance *2.*(Math.random()-0.5);
			  //dtheta = 0.;
			  sciutto(k,l1trial,l2trial);
			  
			  l1[k] = l1[k]+l1trial;
			  l2[k] = l2[k]+l2trial;
			  //System.out.print(Math.sin(Math.PI/2));
			  if (l1[k]<0 || l2[k]<0){
				  condition = true;
			  }
			  a[k] = .5*sizeRatio*Math.pow(l1[k]*12,.5);
			  b[k] = .5*sizeRatio*Math.pow(l2[k]*12,.5);
			  xp[k] = PBC.position(xp[k]+dxtrial,Lx); 
			  yp[k] = PBC.position(yp[k]+dytrial,Ly);
			  theta[k] += dtheta;
			  if (condition == true){
					 l1[k] = l1[k]-l1trial;
					 l2[k] = l2[k]-l2trial;
					 a[k] = .5*sizeRatio*Math.pow(l1[k]*12,.5);
					 b[k] = .5*sizeRatio*Math.pow(l2[k]*12,.5);
					 xp[k] = PBC.position(xp[k]-dxtrial, Lx); // reject displacement 
					 yp[k] = PBC.position(yp[k]-dytrial, Ly);
					 theta[k] -= dtheta;
					 continue;
					 
			  }
			  L = a[k];
				if(b[k]>a[k]){
			  L = b[k];
		      }
		      
			  for(int j = 0;j<Nd;++j){
				 
				  if(tworoots==true){
					   
				double Txa = PBC.separation(x[j]-xp[k],Lx);
				double Tya = PBC.separation(y[j]-yp[k],Ly);
				double dxp = Math.cos(theta[k])*Txa-Math.sin(theta[k])*Tya;
				double dyp = Math.sin(theta[k])*Txa+Math.cos(theta[k])*Tya;
				
				if(over[k][j]==1){ 
					Overcond(j,k);
					if(!tworoots){
						continue;
					}
					if(tworoots==true){
						if(condition==false){
							overtest[k][j]=0;
						}
						if (condition==true){
							overtest[k][j]=1;
							}
				    }
				}
				if(dxp*dxp+dyp*dyp<(L+.5)*(L+.5)){
					Overcond(j,k);
					if(tworoots==false){
						continue;
					}
					if(condition==true){ 
						overtest[k][j]=1;
					}
					
					if(condition == false){
					  overtest[k][j]=0; 
				    }
					
					 
					if ((dxp*dxp/a[k]/a[k])+(dyp*dyp/b[k]/b[k])<1){
						overtest[k][j]=1;
					}
					
					
					
				  

		
}

  

					  
				  }
		

				  
}

if (tworoots){
 double rand = Math.random();


 if(rand<prob*Math.exp(-penetrationCost*(IntStream.of(overtest[0]).sum()-IntStream.of(over[0]).sum())) && tworoots==true){
						  arrayCopy(overtest,over);
						  
					  }
 else if(rand>prob*Math.exp(-penetrationCost*(IntStream.of(overtest[0]).sum()-IntStream.of(over[0]).sum())) && tworoots==true){
					  l1[k] = l1[k]-l1trial;
					  l2[k] = l2[k]-l2trial;
					  arrayCopy(over,overtest);
					  a[k] = .5*sizeRatio*Math.pow(l1[k]*12,.5);
					  b[k] = .5*sizeRatio*Math.pow(l2[k]*12,.5);
					  xp[k] = PBC.position(xp[k]-dxtrial, Lx); //reject displacement 
					  yp[k] = PBC.position(yp[k]-dytrial, Ly);
					  theta[k] -= dtheta;
					  
					  }
}
					  
  
  
 if (!tworoots){
					  
					  l1[k] = l1[k]-l1trial;
					  l2[k] = l2[k]-l2trial;
					  arrayCopy(over,overtest);
					  a[k] = .5*sizeRatio*Math.pow(l1[k]*12,.5);
					  b[k] = .5*sizeRatio*Math.pow(l2[k]*12,.5);
					  xp[k] = PBC.position(xp[k]-dxtrial, Lx); //reject displacement 
					  yp[k] = PBC.position(yp[k]-dytrial, Ly);
					  theta[k] -= dtheta;
				  } 
  
}


}



  

  /**
   * Draws the hard disks by painting circles at each hard disk location.
   *
   * @param drawingPanel DrawingPanel
   * @param g Graphics
   * 
   */
  public void draw(DrawingPanel drawingPanel, Graphics g) {
    double radius = 0.5;
    int pxRadiusP[] = new int[Np];
    int pyRadiusP[] = new int[Np];
    Graphics2D g2 = (Graphics2D) g; //New
    if(x==null) {
      return;
    }
    int pxRadius = Math.abs(drawingPanel.xToPix(radius)-drawingPanel.xToPix(0));
    int pyRadius = Math.abs(drawingPanel.yToPix(radius)-drawingPanel.yToPix(0));
    g.setColor(Color.blue);
    for(int i = 0;i<Nd;i++) {
      int xpix = drawingPanel.xToPix(x[i])-pxRadius;
      int ypix = drawingPanel.yToPix(y[i])-pyRadius;
      g.fillOval(xpix, ypix, 2*pxRadius, 2*pyRadius);
    } // draw cell boundaries
    g.setColor(Color.red);
    for(int i = 0;i<Np;i++) {
      pxRadiusP[i] = Math.abs(drawingPanel.xToPix(a[i])-drawingPanel.xToPix(0)); // Alan
      pyRadiusP[i] = Math.abs(drawingPanel.yToPix(b[i])-drawingPanel.yToPix(0)); // Alan
      int xpixP = drawingPanel.xToPix(xp[i])-pxRadiusP[i];
      int ypixP = drawingPanel.yToPix(yp[i])-pyRadiusP[i];
      Ellipse2D e = new Ellipse2D.Double(xpixP, ypixP, 2*pxRadiusP[i], 2*pyRadiusP[i]);
      AffineTransform at = AffineTransform.getRotateInstance(theta[i], xpixP+pxRadiusP[i], ypixP+pyRadiusP[i]);
      g2.fill(at.createTransformedShape(e));
      
  } // draw cell boundaries
    g.setColor(Color.black);
    int xpix = drawingPanel.xToPix(0);
    int ypix = drawingPanel.yToPix(Ly);
    int lx = drawingPanel.xToPix(Lx)-drawingPanel.xToPix(0);
    int ly = drawingPanel.yToPix(0)-drawingPanel.yToPix(Ly);
    g.drawRect(xpix, ypix, lx, ly);
  }
}

