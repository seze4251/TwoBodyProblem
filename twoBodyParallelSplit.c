//
//  distParaAsync.c
//  MatrixFix
//
//  Created by Seth on 2/2/16.
//  Copyright © 2016 Seth. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "functions.h"

const double mM = 7.34767309E22, mE = 5.97219E24, mS = 28833, rM = 1737100, rE = 6371000, G = 6.675E-11;
const double  thetaS = 50 * M_PI / 180, thetaM = 42.5 * M_PI/180, clear = 10000, xE = 0, yE = 0, vEx = 0, vEy = 0;


int main(int argc, char * argv[]) {
    // Initialize MPI
    MPI_Init(&argc, &argv);
    
    //Initialize Base Vars on all processes
    int nprocs, myrank, tagX = 1, tagY = 2, tagResult = 3, mpi_error;
    const int serverRank = 0;
    double starttime, endtime, totaltime;
    
    
    // Determine # of procs and my rank
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    
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
    
    
    // Read in Data
    double numA, numV, vStart, vMax = 100, angMax = 2 * M_PI;
    double * velX, * velY;
    int i, count, loopCount;
    
    if (myrank == serverRank) {
        
        while (1) {
            printf("Please enter 3 #'s larger than 0,  #angles, #Velocities and start velocity \n");
            int check = scanf("%lf %lf %lf",&numA, &numV, &vStart);
            if (numA > 0 && numV > 0 && vStart > 0) {
                break;
            } else {
                printf("Please Enter Three Numbers that Are greater than 0 \n");
            }
        }
        int i, j, count = 0, size;
        size = (int) (numA * numV);
        double holderX[size], holderY[size], aSS = angMax / numA , vSS = (vMax - vStart) / numV;
        
        for (i = 0; i < numV ; i++) {
            for (j = 0; j < numA; j++) {
                holderX[count] = (vStart + vSS * i) * cos(j * aSS) + vS * cos(thetaS);
                holderY[count] = (vStart + vSS * i) * sin(j * aSS) + vS * sin(thetaS);
                count++;
            }
        }
        
        velX = holderX;
        velY = holderY;
        starttime = MPI_Wtime();
    }
    
    
    // Determine necissary variables for splitting by # of processes
    int rem, scatterSize;
    
    MPI_Request reqX [nprocs - 1];
    MPI_Request reqY [nprocs - 1];
    MPI_Status status;
    
    // Sends data
    if (myrank == serverRank) {
        loopCount = count + rem;
        rem = (int)(numA * numV) % nprocs;
        scatterSize = (numA * numV) / nprocs;
        int count;
        for (i = 1; i < nprocs; i++) {
            mpi_error = MPI_Isend(velX + rem +scatterSize * (i + 1), scatterSize, MPI_DOUBLE, i, tagX, MPI_COMM_WORLD, (reqX+i));
            mpi_error = MPI_Isend(velY + rem + scatterSize * (i + 1), scatterSize, MPI_DOUBLE, i, tagY, MPI_COMM_WORLD, (reqY+i));
        }
    } else {
        int i;
        for (i = 0; i < 2; i++) {
            mpi_error = MPI_Probe(serverRank, MPI_ANY_TAG, MPI_COMM_WORLD, & status);
            mpi_error = MPI_Get_count(& status, MPI_DOUBLE, & count);
            double data [count];
            if (status.MPI_TAG == tagX) {
                MPI_Recv(data, count, MPI_DOUBLE, serverRank, tagX, MPI_COMM_WORLD, & status);
                velX = data;
            } else {
                MPI_Recv(data, count, MPI_DOUBLE, serverRank, tagY, MPI_COMM_WORLD, & status);
                velY = data;
            }
        }
        loopCount = count;
    }
    
    // Perform Calculation
    double currentTime, bestTime = 0, topVX, topVY;
    int ret;
    
    for (i = 0; i < loopCount; i++) {
        IC[6] = velX[i];
        IC[7] = velY[i];
        ret = fowardEuler(& currentTime, IC);
        if (ret == 0) {
            if (bestTime == 0) {
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
    
    // Wait on server for Everything to have gone through
    if (myrank == serverRank) {
        MPI_Waitall(nprocs, reqX, MPI_STATUS_IGNORE);
        MPI_Waitall(nprocs, reqY, MPI_STATUS_IGNORE);
    }
    
    
    // Collect Results
    if (myrank == serverRank) {
        double data[3];
        for (i = 1; i < nprocs; i++) {
            mpi_error = MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, & status);
            mpi_error = MPI_Recv(data, 3, MPI_DOUBLE, status.MPI_SOURCE, tagResult, MPI_COMM_WORLD, &status);
            if (data[1] < bestTime) {
                bestTime = data[0];
                topVX = data[1];
                topVY = data[2];
            }
        }
        
    } else {
        double send [] = {bestTime, topVX, topVY};
        mpi_error =  MPI_Send(send, 3, MPI_DOUBLE, serverRank, tagResult, MPI_COMM_WORLD);
    }
    
    MPI_Finalize();
    if (myrank == serverRank) {
            printf("Best Time: %5.3f  MagV: %5.3f,  theta: %5.3f degrees", bestTime, sqrt(pow(topVX,2) + pow(topVY,2)), 180 / M_PI *atan2(topVY, topVX));
    }

    return 0;
}



