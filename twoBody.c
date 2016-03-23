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
#include "functions.h"
#include <time.h>

// Declare Global Constants
const double mM = 7.34767309E22, mE = 5.97219E24, mS = 28833, rM = 1737100, rE = 6371000, G = 6.674E-11;
const double  thetaS = 50 * M_PI / 180, thetaM = 42.5 * M_PI/180, xE = 0, yE = 0, vEx = 0, vEy = 0;
double clear;

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
    
    double time, x, y;
    int ret = fowardEuler(& time, IC);
    printf("time = %7.7f, return value = %d\n", time, ret);
    
    
    // INPUT BLOCK
    double numA, numV, vStart, vMax = 100, angMax = 2 * M_PI;
    double * velX, * velY;
    int i, j, count, size, err;
    
    while (1) {
        printf("# Angles, # Velocities, #Velocity Start, #clearence \n");
        err = scanf("%lf %lf %lf %lf",&numA, &numV, &vStart, &clear);
        printf("err = %d  clear %f\n",err, clear);
        if (err == 4) {
            if (numA > 0 && numV > 0 && vStart > 0 && clear >= 0) {
                break;
            } else {
                printf("Please Enter Three Numbers that Are greater than 0 \n");
            }
        }
        fseek(stdin,0,SEEK_END);
    }
    
    count = 0;
    size = (int) (numA * numV);
    velX = (double *) malloc(size * sizeof(double));
    velY = (double *) malloc(size * sizeof(double));
    
    double aSS = angMax / numA , vSS = (vMax - vStart) / numV;
    
    for (i = 0; i < numV ; i++) {
        for (j = 0; j < numA; j++) {
            velX[count] = (vStart + vSS * (i + 1)) * cos(j * aSS) + vS * cos(thetaS);
            velY[count] = (vStart + vSS * (i + 1)) * sin(j * aSS) + vS * sin(thetaS);
            count++;
        }
    }
    
    clock_t start_t, end_t;
    double total_t;
    start_t = clock();
    
    // Compute Block
    double currentTime, bestTime = 0, topVX = 0, topVY = 0;
    
    for (i = 0; i < numA*numV; i++) {
        IC[6] = velX[i];
        IC[7] = velY[i];
        ret = fowardEuler(& currentTime, IC);
        //  printf("ret: %d  i: %d \n",ret, i);
        if (ret == 0) {
            if (bestTime == 0){
                bestTime = currentTime;
                topVX = velX[i];
                topVY = velY[i];
            } else if (currentTime < bestTime) {
                bestTime = currentTime;
                topVX = velX[i];
                topVY = velY[i];
            }
        }
    }
    
    x = topVX - vS * cos(thetaS);
    y = topVY - vS * sin(thetaS);
    
    printf("Best Time: %5.5f, vx %5.5f,  vy %5.5f, magnitude %5.5f theta %5.5f \n", bestTime, x, y, sqrt(pow(x,2) + pow(y,2)), atan2(y,x) * 180/M_PI);
    end_t = clock();
    total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
    printf("Total Time: %5.5f\n",total_t);
    
    // Print to File
    FILE * file = fopen("SingleCore","a");
    fprintf(file, " %d \t\t %f \t\t %f \n",1,total_t, numA * numV);
    fclose(file);
    
    return 0;
}


