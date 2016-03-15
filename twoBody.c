//
//  twoBodyMain.h
//
//
//  Created by Seth on 3/8/16.
//  Copyright © 2016 Seth. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrixFix.h"

int fowardEuler(double * t, double IC []);
int rk4(double * t, double IC []);

// Declare Global Constants
const double mM = 7.34767309E22, mE = 5.97219E24, mS = 28833, rM = 1737100, rE = 6371000, G = 6.675E-11;
const double  thetaS = 50 * M_PI / 180, thetaM = 42.5 * M_PI/180, clear = 10000, xE = 0, yE = 0, vEx = 0, vEy = 0;


int main(int argc, const char * argv[]) {
    // Declare Initial Conditions
    double dES = 340000000, dEM = 384403000;
    double vM =  sqrt(G * pow(mE,2) / ((mE + mM) *dEM) ), vS = 1000;
    double IC [12];
    IC[0] = vS; // Velocity of Spaceship vS
    IC[1] = vM; // Velocity of Moon vM
    IC[2] = dES; // Distance Between Earth Spaceship dES
    IC[3] = dEM; // Distance Between Earth Moon dEM
    IC[4] = dES * cos(thetaS); // X position of Spaceship xS
    IC[5] = dES * sin(thetaS); // Y position of Spaceship yS
    IC[6] = vS * cos(thetaS);  // Velocity in x direction of Spaceship vSx
    IC[7] = vS * sin(thetaS);  // Velcity in y direction of Spaceship vSy
    IC[8] = dEM * cos(thetaM); // X position of Moon xM
    IC[9] = dEM * sin(thetaM); // Y position of Moon yM
    IC[10] = -vM * sin(thetaM); // Velocity in x direction of Moon vMx
    IC[11] = vM * cos(thetaM);  // Velocity in y direction of Moon vMy
    
    // Test Inital Conditions
    double time;
    int ret = fowardEuler(& time, IC);
    printf("time = %7.7f, return value = %d\n", time, ret);
    
    
    // Solve for smallest Velocity
    /*    int i, j, steps = 1, minV = 50, maxV = 100, angleStep = 10;
     double angleStepSize = 2*M_PI/angleStep, velocityStepSize = (double) maxV/steps;
     printf("angleStep= %5.5f, max vel = %5.5f, velocityStep = %5.5f", angleStepSize * steps, steps * velocityStepSize, velocityStepSize);
     
     // Rows represent Velocities and Cols Represent Theta
     matrix * mtxVx = newMatrix(steps, angleStep);
     matrix * mtxVy = newMatrix(steps, angleStep);
     matrix * mtxResult = newMatrix(steps, angleStep);
     printf("Filling Matrix\n");
     // Fill mtxVx and mtxVy
     for (i = 0; i < steps; i++) { // coresponds to rows / Velocities
     for (j = 0; j < angleStep; j++) { // cols / theta
     ELEM(mtxVx, i, j) = vS * cos(thetaS) +  velocityStepSize * cos(angleStepSize * j);
     ELEM(mtxVy, i, j) = vS * sin(thetaS) +  velocityStepSize * sin(angleStepSize * j);
     printf("cos(anglestep * j) = %5.5f  sin(anglestep * j) = %5.5f\n", i*velocityStepSize * cos(angleStepSize * j),i*  velocityStepSize * sin(angleStepSize * j));
     }
     }
     
     /*  printMatrix(mtxVx);
     
     printf("Starting Integration\n");
     // Compute integration and save results
     for (i = 0; i < steps; i++) {
     for (j = 0; j < angleStep; j++) {
     IC[6] = ELEM(mtxVx, i, j);
     IC[7] = ELEM(mtxVy, i, j);
     ret = fowardEuler(& time, IC);
     printf("ret = %d, Velocity %5.5f  deltaVx = %5.5f  delta Vy = %5.5f \n",ret,  sqrt(pow(ELEM(mtxVx, i, j),2) + pow(ELEM(mtxVy, i, j) ,2)), ELEM(mtxVx, i, j) , ELEM(mtxVy, i, j));
     if (ret == 0) {
     ELEM(mtxResult,i,j) = 1;
     printf("Magnitude Velocity = %5.5f, deltaVx = %5.5f  delta Vy = %5.5f\n", sqrt(pow(ELEM(mtxVx, i, j),2) + pow(ELEM(mtxVy, i, j) ,2)),ELEM(mtxVx, i, j) , ELEM(mtxVy, i, j));
     break;
     }
     }
     if (ret == 0) {
     break;
     }
     }
     
     */
    
    // Find the fastest way back to Earth
    double fastestTime = 0;
    int i, j, steps = 50, minV = 50, maxV = 100, angleStep = 10;
    
    double angleStepSize = 2*M_PI/angleStep, velocityStepSize = (double) (maxV - minV)/steps;
    printf("angleStep= %5.5f, max vel = %5.5f, velocityStep = %5.5f", angleStepSize * steps, steps * velocityStepSize, velocityStepSize);
    
    // Rows represent Velocities and Cols Represent Theta
    matrix * mtxVx = newMatrix(steps, angleStep);
    matrix * mtxVy = newMatrix(steps, angleStep);
    matrix * mtxResult = newMatrix(steps, angleStep);
    printf("Filling Matrix\n");
    // Fill mtxVx and mtxVy
    for (i = 0; i < steps; i++) { // coresponds to rows / Velocities
        for (j = 0; j < angleStep; j++) { // cols / theta
            ELEM(mtxVx, i, j) = vS * cos(thetaS) + 50 + ( i + 1) * velocityStepSize * cos(angleStepSize * j);
            ELEM(mtxVy, i, j) = vS * sin(thetaS) +  50 +(i + 1) * velocityStepSize * sin(angleStepSize * j);
            printf("cos(anglestep * j) = %5.5f  sin(anglestep * j) = %5.5f\n", i*velocityStepSize * cos(angleStepSize * j),i*  velocityStepSize * sin(angleStepSize * j));
        }
    }
    
    
    for (i = 0; i < steps; i++) {
        for (j = 0; j < angleStep; j++) {
            IC[6] = ELEM(mtxVx, i, j);
            IC[7] = ELEM(mtxVy, i, j);
            ret = fowardEuler(& time, IC);
            if (ret == 0) {
                if (fastestTime == 0) {
                    fastestTime = time;
                }else if (time < fastestTime) {
                    fastestTime = time;
                    printf("New Fastest Time is %5.5f", fastestTime);
                }
                
            }
        }
    }
    
    
    
    
    return 0;
}

int fowardEuler(double * t, double IC []) {
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


