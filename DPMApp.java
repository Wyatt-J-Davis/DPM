package org.opensourcephysics.sip.DPM;

import java.awt.*;
import java.awt.geom.Ellipse2D;
import java.awt.geom.AffineTransform;
import java.text.SimpleDateFormat;

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
import java.util.Date;

/**
 * HardDisksApp does a Monte Carlo simulation for a 2D system of uncharged hard disks crowding an ellipsoidal polymer
 *
 * @author Alan Denton & Wyatt Davis based on MD code ch08/hd/HardDisksApp.java by Jan Tobochnik, Wolfgang Christian, Harvey Gould
 * @version 1.0 revised 12/28/17
 * Init.s a second species of disk randomly amongst the first species of disk. Avoids interspecific overlap. Updates positions of both species,  species two behaves ideally in isolation.
 */

public class DPMApp extends AbstractSimulation {
    DPM hd = new DPM();
    DisplayFrame display = new DisplayFrame("x", "y", "Hard-Disk Polymer Mixture");
    int i = 0; // Will track number of steps taken
    int j = 0; // Will track number of data points written
    int k = 0; // Will track number of independant runs
    ArrayList<Double> L1dat = new ArrayList<Double>(10000);
    ArrayList<Double> L2dat = new ArrayList<Double>(10000);
    Calendar cal = Calendar.getInstance();
    SimpleDateFormat dateFormat = new SimpleDateFormat("MM-dd-yyyy hh.mm.ss");
    Date date = cal.getTime();
    String strDate = dateFormat.format(date);


    /**
     * Initializes the simulation
     */
    public void initialize() {
        hd.Nd = control.getInt("Nd");
        hd.Np = control.getInt("Np");
        hd.phi = control.getDouble("phi");
        hd.Lx = .5 * Math.sqrt(hd.Nd * Math.PI / hd.phi);
        hd.Ly = .5 * Math.sqrt(hd.Nd * Math.PI / hd.phi);
        hd.condition = false;
        hd.tolerance = control.getDouble("tolerance");
        String configuration = control.getString("initial configuration");
        hd.sizeRatio = control.getDouble("sizeRatio");
        hd.penetrationCost = control.getDouble("penetrationCost");
        hd.solventQuality = control.getString("solventQuality");
        hd.initialize(configuration);
        display.addDrawable(hd);
        display.setPreferredMinMax(0, hd.Lx, 0, hd.Ly);
        display.setSquareAspect(true);

    }

    //Dictates how data is taken over duration of simulation
    public void doStep() {
        if (k < 6) {
            if (k == 0 && i == 0) {
                String topLevelPath = "data/phiN=" + hd.phi + ",q=" + hd.sizeRatio + "_" + strDate;
                File dir1 = new File(topLevelPath);
                boolean successful1 = dir1.mkdir();
                if (successful1) {
                    System.out.println("directory was created successfully");
                } else {
                    System.out.println("failed trying to create the directory");
                }
                for (int i = 1; i < 6; i++) {
                    File dir2 = new File(topLevelPath + "/" + i);
                    boolean successful2 = dir2.mkdir();
                    if (successful2) {
                        System.out.println("directory was created successfully");
                    } else {
                        System.out.println("failed trying to create the directory");
                    }
                    try {
                        FileWriter fw = new FileWriter(topLevelPath + "/" + i + "/" + "phiN=" + hd.phi + ",q=" + hd.sizeRatio + "_eX.dat");
                    } catch (Exception e) {
                        System.out.println(e);
                    }
                    System.out.println("Success...");
                    try {
                        FileWriter fw = new FileWriter(topLevelPath + "/" + i + "/" + "phiN=" + hd.phi + ",q=" + hd.sizeRatio + "_eY.dat");
                    } catch (Exception e) {
                        System.out.println(e);
                    }
                    System.out.println("Success...");
                }
                k++;
            }
            //check if enough datapoints have been taken
            if (j < 10000) {
                //allow steps for equilibriation
                if (i < 10000) {
                    hd.step();
                    i++;
                }
                //after steps, start data collection procedure
                else {
                    hd.step();
                    i++;
                    if (i % 1000 == 0) {
                        if (hd.l1[0] >= hd.l2[0]) {
                            L1dat.add(hd.l1[0]);
                            L2dat.add(hd.l2[0]);
                        } else {
                            L1dat.add(hd.l2[0]);
                            L2dat.add(hd.l1[0]);
                        }
                        j++;
                        //System.out.print(j);
                    }
                }
            }
            //WRITE DATA
            else {
                String topLevelPath = "data/phiN=" + hd.phi + ",q=" + hd.sizeRatio + "_" + strDate;
                String strFilePathX = topLevelPath + "/" + k + "/" + "phiN=" + hd.phi + ",q=" + hd.sizeRatio + "_eX.dat";
                String strFilePathY = topLevelPath + "/" + k + "/" + "phiN=" + hd.phi + ",q=" + hd.sizeRatio + "_eY.dat";
                System.out.println(strFilePathX);
                try {
                    FileWriter fwX = new FileWriter(strFilePathX);
                    BufferedWriter bw = new BufferedWriter(fwX);
                    for (double dat : L1dat) {
                        bw.write(dat + "");
                        bw.newLine();
                    }
                    bw.close();
                } catch (Exception e) {
                    System.out.println(e);
                }
                try {
                    FileWriter fwY = new FileWriter(strFilePathY);
                    BufferedWriter bw = new BufferedWriter(fwY);
                    for (double dat : L2dat) {
                        bw.write(dat + "");
                        bw.newLine();
                    }
                    bw.close();
                } catch (Exception e) {
                    System.out.println(e);
                }
                k++;
                L1dat.clear();
                L2dat.clear();
                i = 0;
                j = 0;
                initialize();
            }
        }
    }

    /**
     * Resets the hard disks model to its default state.
     */
    public void reset() {
        enableStepsPerDisplay(true);
        control.setValue("Nd", 36);
        control.setValue("Np", 1);
        control.setValue("phi", .2);
        control.setValue("tolerance", 0.1);
        control.setValue("initial configuration", "regular");
        control.setValue("sizeRatio", 1);
        control.setValue("penetrationCost", .5);
        control.setValue("solventQuality", "theta");
        initialize();
    }

    /**
     * Resets the hard disks model and the data graphs.
     * <p>
     * This method is invoked using a custom button.
     */
    public void resetData() {
        hd.resetAverages();
        GUIUtils.clearDrawingFrameData(false); // clears old data from the plot frames
    }

    /**
     * Starts the Java application.
     *
     * @param args command line parameters
     */
    public static void main(String[] args) { // set up animation control structure using this class
        SimulationControl control = SimulationControl.createApp(new DPMApp());
        control.addButton("resetData", "Reset Data");

    }
}

