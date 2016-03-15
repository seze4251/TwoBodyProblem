//
//  functions.c
//  functions
//
//  Created by Seth on 3/14/16.
//  Copyright © 2016 Seth. All rights reserved.
//


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>




int fowardEuler(double * t, double IC []) {
    extern double mM, mE, mS, rM, rE, G, thetaS, thetaM, clear, xE, yE, vEx, vEy;
    // Break down IC vector
    double vS, vM, dES, dEM, xS, yS, vSx, vSy, xM, yM, vMx, vMy;
    vS = IC[0];
    vM = IC[1];
    dES = IC[2];
    dEM = IC[3];
    xS = IC[4];
    yS = IC[5];
    vSx = IC[6];
    vSy = IC[7];
    xM = IC[8];
    yM = IC[9];
    vMx = IC[10];
    vMy = IC[11];
    
    double dMS = sqrt( pow(xS - xM, 2) + pow(yS - yM, 2) );
    // printf("vS = %5.5f vM = %5.5f dES = %5.5f dEM = %5.5f xS = %5.5f \n ", vS, vM, dES, dEM, xS);
    // printf("yS = %5.5f vSx = %5.5f vSy = %5.5f xM = %5.5f yM = %5.5f  vMx = %5.5f vMy = %5.5f \n", yS, vSx, vSy, xM, yM, vMx, vMy);
    // printf("dMS = %5.5f\n", dMS);
    // Distance between Spaceship and Moon
    // Availble Stuff from [IC] vS, vM, dES, dEM, dSM, xS, yS, vSx, vSy, xM, yM, vMx, vMy
    // Globals (constants): mM, mE, mS, rM, rE, G, thetaS, thetaM, clear, xE, yE, vEx, vEy
    
    // Declare time;
    double dt = 1; // Second
    *t = 0;
    
    // Declare Derivatives:
    double aSx, aSy, aMx, aMy;
    
    // Declare Forces: F = Ma, a = F/m
    double fSMx, fSMy, fEMx, fEMy, fESx, fESy;
    //  printf("power: %5.5f\n",powf(dMS,3));
    // Start Integration Loop
    //int i = 1;
    while (1) {
        // If spacecraft is heading home, Increase timestep!
        if (vSx < 0) {
        }
        
        // Change Time
        *t = *t + dt;
        
        // Find Forces
        // fMSx and fMSy are forces on Spacecraft by moon
        // Draw this out and check my Signs!!!
        // Mase sure fMSx is the one I want it to be
        fSMx = G * mS * mM * (xM - xS) / pow(dMS,3);
        fSMy = G * mS * mM * (yM - yS) / pow(dMS,3);
        fEMx = G * mM * mE * (-xM) / pow(dEM,3);
        fEMy = G * mM * mE * (-yM) / pow(dEM,3);
        fESx = G * mS * mE * (-xS) / pow(dES,3);
        fESy = G * mS * mE * (-yS) / pow(dES,3);
        
        /*
         if (i == 1) {
         printf("fMSx = %5.5f  fMSy = %5.5f FEMx = %5.5f FEMy = %5.5f  fMSx = %5.5f  fMSx = %5.5f\n",fSMx, fSMy, fEMx, fEMy, fESx, fESy);
         
         }
         
         */
        // Find acelerations
        aSx = (fSMx + fESx)/mS;
        aSy = (fSMy + fESx)/mS;
        aMx = (fEMx - fSMx)/mM;
        aMy = (fEMy - fSMy)/mM;
        /*
         if (i == 1) {
         printf("aSx = %5.5f  aSy = %5.5f aMx = %5.5f aMy = %5.5f\n",aSx, aSy, aMx, aMy);
         
         }
         i += 1;
         */
        // New Velocities with Foward Euler
        vSx = vSx + dt * aSx;
        vSy = vSy + dt * aSy;
        vMx = vMx + dt * aMx;
        vMy = vMy + dt * aMy;
        
        // New Positions with Foward Euler
        xS = xS + dt * vSx;
        yS = yS + dt * vSy;
        xM = xM + dt * vMx;
        yM = yM + dt * vMy;
        
        // Solve for Ending Conditions
        dMS =  sqrt(pow(xS - xM, 2) + pow(yS - yM, 2));
        dEM =  sqrt(pow(xM,2) + pow(yM,2));
        dES =  sqrt(pow(xS,2) + pow(yS,2));
        
        //   printf("vS = %5.5f vM = %5.5f dES = %5.5f dEM = %5.5f xS = %5.5f \n ", vS, vM, dES, dEM, xS);
        //     printf("yS = %5.5f vSx = %5.5f vSy = %5.5f xM = %5.5f yM = %5.5f  vMx = %5.5f vMy = %5.5f \n", yS, vSx, vSy, xM, yM, vMx, vMy);
        // printf("dMS = %5.5f\n", dMS);
        
        if (dMS <= rM + clear) {
            //printf("Exiting Simulation because Spacecraft Crashed with Moon\n");
            return -1;
            
        } else if (dES < rE) {
            //printf("Exiting Simulation because Spacecraft Returned to Earth\n");
            return 0;
            
        } else if(dES > 2 * dEM) {
            // printf("Exiting Simulation because Spacecraft Lost to Space\n");
            return - 2;
        }
        //printf("Distance Between SC and Moon %5.3f\n",powf(xS - xM, 2) + powf(yS - yM, 2));
    }
}

int rk4(double * t, double IC []) {
    extern double mM, mE, mS, rM, rE, G, thetaS, thetaM, clear, xE, yE, vEx, vEy;
    // Break down IC vector
    double vS, vM, dES, dEM, xS, yS, vSx, vSy, xM, yM, vMx, vMy;
    vS = IC[0];
    vM = IC[1];
    dES = IC[2];
    dEM = IC[3];
    xS = IC[4];
    yS = IC[5];
    vSx = IC[6];
    vSy = IC[7];
    xM = IC[8];
    yM = IC[9];
    vMx = IC[10];
    vMy = IC[11];
    
    double dMS = sqrt( pow(xS - xM, 2) + pow(yS - yM, 2) );
    // printf("vS = %5.5f vM = %5.5f dES = %5.5f dEM = %5.5f xS = %5.5f \n ", vS, vM, dES, dEM, xS);
    // printf("yS = %5.5f vSx = %5.5f vSy = %5.5f xM = %5.5f yM = %5.5f  vMx = %5.5f vMy = %5.5f \n", yS, vSx, vSy, xM, yM, vMx, vMy);
    // printf("dMS = %5.5f\n", dMS);
    // Distance between Spaceship and Moon
    // Availble Stuff from [IC] vS, vM, dES, dEM, dSM, xS, yS, vSx, vSy, xM, yM, vMx, vMy
    // Globals (constants): mM, mE, mS, rM, rE, G, thetaS, thetaM, clear, xE, yE, vEx, vEy
    
    // Declare time;
    double dt = 0.5; // Second
    *t = 0;
    
    // Declare Derivatives:
    double aSx, aSy, aMx, aMy;
    
    // Declare Forces: F = Ma, a = F/m
    double fSMx, fSMy, fEMx, fEMy, fESx, fESy;
    //  printf("power: %5.5f\n",powf(dMS,3));
    // Start Integration Loop
    //int i = 1;
    
    // Run Kutta Parameters
    double k1, k2, k3, k4;
    
    while (1) {
        // If spacecraft is heading home, Increase timestep!
        if (vSx < 0) {
            dt = 3;
        }
        
        // Change Time
        *t = *t + dt;
        
        // Find Forces
        // fMSx and fMSy are forces on Spacecraft by moon
        // Draw this out and check my Signs!!!
        // Mase sure fMSx is the one I want it to be
        fSMx = G * mS * mM * (xM - xS) / pow(dMS,3);
        fSMy = G * mS * mM * (yM - yS) / pow(dMS,3);
        fEMx = G * mM * mE * (-xM) / pow(dEM,3);
        fEMy = G * mM * mE * (-yM) / pow(dEM,3);
        fESx = G * mS * mE * (-xS) / pow(dES,3);
        fESy = G * mS * mE * (-yS) / pow(dES,3);
        
        /*
         if (i == 1) {
         printf("fMSx = %5.5f  fMSy = %5.5f FEMx = %5.5f FEMy = %5.5f  fMSx = %5.5f  fMSx = %5.5f\n",fSMx, fSMy, fEMx, fEMy, fESx, fESy);
         
         }
         
         */
        // Find acelerations
        aSx = (fSMx + fESx)/mS;
        aSy = (fSMy + fESx)/mS;
        aMx = (fEMx - fSMx)/mM;
        aMy = (fEMy - fSMy)/mM;
        /*
         if (i == 1) {
         printf("aSx = %5.5f  aSy = %5.5f aMx = %5.5f aMy = %5.5f\n",aSx, aSy, aMx, aMy);
         
         }
         i += 1;
         */
        
        k1 = vSx + dt * aSx;
        k2 = vSx + dt/2 * k1;
        k3 = vSx + dt/2  * k2;
        k4 = vSx + dt * k3;
        vSx = vSx + dt / 6 * (k1 + 2* k2 + 2* k3 + k4);  // aSx;
        
        k1 = vSy + dt * aSy;
        k2 = vSy + dt/2 * k1;
        k3 = vSy + dt/2  * k2;
        k4 = vSy + dt * k3;
        vSy = vSy + dt  / 6 * (k1 + 2* k2 + 2* k3 + k4); //* aSy;
        
        k1 = vMx + dt * aMx;
        k2 = vMx + dt/2 * k1;
        k3 = vMx + dt/2  * k2;
        k4 = vMx + dt * k3;
        vMx = vMx + dt  / 6 * (k1 + 2* k2 + 2* k3 + k4); //* aMx;
        
        k1 = vMy + dt * aMy;
        k2 = vMy + dt/2 * k1;
        k3 = vMy + dt/2  * k2;
        k4 = vMy + dt * k3;
        vMy = vMy + dt / 6 * (k1 + 2* k2 + 2* k3 + k4);  // * aMy;
        
        // New Positions with MidPoint
        
        k1 = xS + dt * vSx;
        k2 = xS + dt/2 * k1;
        k3 = xS + dt/2  * k2;
        k4 = xS + dt * k3;
        xS = xS + dt  / 6 * (k1 + 2* k2 + 2* k3 + k4);  // * vSx;
        
        k1 = yS + dt * vSy;
        k2 = yS + dt/2 * k1;
        k3 = yS + dt/2  * k2;
        k4 = yS + dt * k3;
        yS = yS + dt  / 6 * (k1 + 2* k2 + 2* k3 + k4);  // * vSy;
        
        k1 = xM + dt * vMx;
        k2 = xM + dt/2 * k1;
        k3 = xM + dt/2  * k2;
        k4 = xM + dt * k3;
        xM = xM + dt  / 6 * (k1 + 2* k2 + 2* k3 + k4);  //* vMx;
        
        k1 = yM + dt * vMy;
        k2 = yM + dt/2 * k1;
        k3 = yM + dt/2  * k2;
        k4 = yM + dt * k3;
        yM = yM + dt / 6 * (k1 + 2* k2 + 2* k3 + k4); // * vMy;
        
        // Solve for Ending Conditions
        dMS =  sqrt(pow(xS - xM, 2) + pow(yS - yM, 2));
        dEM =  sqrt(pow(xM,2) + pow(yM,2));
        dES =  sqrt(pow(xS,2) + pow(yS,2));
        
        //   printf("vS = %5.5f vM = %5.5f dES = %5.5f dEM = %5.5f xS = %5.5f \n ", vS, vM, dES, dEM, xS);
        //     printf("yS = %5.5f vSx = %5.5f vSy = %5.5f xM = %5.5f yM = %5.5f  vMx = %5.5f vMy = %5.5f \n", yS, vSx, vSy, xM, yM, vMx, vMy);
        // printf("dMS = %5.5f\n", dMS);
        
        if (dMS <= rM + clear) {
            //printf("Exiting Simulation because Spacecraft Crashed with Moon\n");
            return -1;
            
        } else if (dES < rE) {
            //printf("Exiting Simulation because Spacecraft Returned to Earth\n");
            return 0;
            
        } else if(dES > 2 * dEM) {
            // printf("Exiting Simulation because Spacecraft Lost to Space\n");
            return - 2;
        }
        //printf("Distance Between SC and Moon %5.3f\n",powf(xS - xM, 2) + powf(yS - yM, 2));
    }
}
