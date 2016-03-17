//
//  distParaAsync.c
//  MatrixFix
//
//  Created by Seth on 2/2/16.
//  Copyright © 2016 Seth. All rights reserved.
//

// Takes in scanF inputs: Number of Angles, Number of Velocities, Starting Velocity and Clearence Height

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "functions.h"

// Declare Globals
const double mM = 7.34767309E22, mE = 5.97219E24, mS = 28833, rM = 1737100, rE = 6371000, G = 6.675E-11;
const double  thetaS = 50 * M_PI / 180, thetaM = 42.5 * M_PI/180, xE = 0, yE = 0, vEx = 0, vEy = 0;
double clear;

int main(int argc, char * argv[]) {
    // Initialize MPI
    MPI_Init(&argc, &argv);
    
    //Initialize Base Vars on all processes
    int nprocs, myrank, tagX = 5, tagY = 8, tagResult = 3, mpi_error;
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
    int i, j, count, size, loopCount;
    
    if (myrank == serverRank) {
        // Read in data
        while (1) {
            printf("# Angles, # Velocities, #Velocity Start, #clearence \n");
            scanf("%lf %lf %lf %lf",&numA, &numV, &vStart, &clear);
            if (numA > 0 && numV > 0 && vStart > 0 && clear >= 0) {
                break;
            } else {
                printf("Please Enter Three Numbers that Are greater than 0 \n");
            }
        }
        // Progress Tracker
        printf("Out of Scanloop\n");
        
        count = 0;
        size = (int) (numA * numV);
        velX = (double *) malloc(size * sizeof(double));
        velY = (double *) malloc(size * sizeof(double));

        double aSS = angMax / numA , vSS = (vMax - vStart) / numV;
        
        for (i = 0; i < numV ; i++) {
            for (j = 0; j < numA; j++) {
                velX[count] = (vStart + vSS * (i + 1)) * cos(j * aSS) + vS * cos(thetaS);
                velY[count] = (vStart + vSS * (i + 1)) * sin(j * aSS) + vS * sin(thetaS);
               // printf(" velX = %5.5f, velY = %5.5f \n", velX[count], velY[count]); //(vStart + vSS * (i + 1)) * cos(j * aSS), (vStart + vSS * (i + 1)) * sin(j * aSS) );
                count++;
            }
        }
        
        /* Test Inputs
        for (i = 0; i < size; i ++) {
            printf(" velX = %5.5f, velY = %5.5f  i  %d\n", velX[i], velY[i], i);
        }
         */
        
        starttime = MPI_Wtime();
    }
    
    // Broadcast Clearence
    MPI_Request req;
    mpi_error = MPI_Ibcast(& clear, 1, MPI_DOUBLE, serverRank, MPI_COMM_WORLD, & req);
    
    // Determine necissary variables for splitting by # of processes
    int rem, scatterSize;
    
    MPI_Request reqX [nprocs - 1];
    MPI_Request reqY [nprocs - 1];
    MPI_Status status;
    
    // Sends data
    if (myrank == serverRank) {
        rem = size % nprocs;
        scatterSize = size / nprocs;
        loopCount = scatterSize + rem;
       // printf("SERVER: scatersize: %d rem: %d  loopCount %d\n", scatterSize, rem, loopCount);
       /*
        for (i = 0; i < size ; i++) {
            printf("SERVER: velX = %5.5f, velY = %5.5f \n", velX[i], velY[i]);
        } */
        
        for (i = 1; i < nprocs; i++) {
            mpi_error = MPI_Isend(velX + rem +scatterSize * i, scatterSize, MPI_DOUBLE, i, tagX, MPI_COMM_WORLD, (reqX+i));
            mpi_error = MPI_Isend(velY + rem + scatterSize * i , scatterSize, MPI_DOUBLE, i, tagY, MPI_COMM_WORLD, (reqY+i));
        }
    } else {
        double * data;
        
        for (i = 0; i < 2; i++) {
            mpi_error = MPI_Probe(serverRank, MPI_ANY_TAG, MPI_COMM_WORLD, & status);
            mpi_error = MPI_Get_count(& status, MPI_DOUBLE, & count);
           // printf("count:%d  tag: %d  tagX = %d,  tagY = %d \n",count, status.MPI_TAG, tagX, tagY);
            data = (double *) malloc(count * sizeof(double));
            
            if (status.MPI_TAG == tagX) {
                MPI_Recv(data, count, MPI_DOUBLE, serverRank, tagX, MPI_COMM_WORLD, & status);
                velX = data;
                
            } else if (status.MPI_TAG == tagY) {
                MPI_Recv(data, count, MPI_DOUBLE, serverRank, tagY, MPI_COMM_WORLD, & status);
                velY = data;
            }
        }
        
        loopCount = count;
    }
    
    // Wait for CLearence
    MPI_Wait(&req, MPI_STATUS_IGNORE);
    /*
   if (myrank != serverRank) {
        for (i = 0; i < loopCount; i++) {
            printf(" velX = %5.5f, velY = %5.5f \n", velX[i], velY[i]);
        }
    }
     */
     // Perform Calculation
    printf("loopcount: %d, myrank %d",loopCount, myrank);
    printf("Begining Calculation: myrank = %d\n",myrank);
    double currentTime, bestTime = 0, topVX, topVY;
    int ret;
    
    for (i = 0; i < loopCount; i++) {
        IC[6] = velX[i];
        IC[7] = velY[i];
        ret = fowardEuler(& currentTime, IC);
     //  printf("ret: %d, myrank: %d  i = %d\n", ret, myrank, i);
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
    MPI_Barrier(MPI_COMM_WORLD);
    
    printf("Finished Calculation: myrank = %d \n", myrank);
    
    // Wait on server for Everything to have gone through
    if (myrank == serverRank) {
        MPI_Waitall(nprocs - 1, reqX + 1, MPI_STATUS_IGNORE);
        MPI_Waitall(nprocs - 1, reqY + 1, MPI_STATUS_IGNORE);
        printf("Server made it past waitall\n");
    }
    
    
    // Collect Results
    if (myrank == serverRank) {
        printf("Server Begining to Collect Results \n");
        double data[3];
        for (i = 1; i < nprocs; i++) {
            mpi_error = MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, & status);
            mpi_error = MPI_Recv(data, 3, MPI_DOUBLE, status.MPI_SOURCE, tagResult, MPI_COMM_WORLD, &status);
            printf("SLAVE: t: %5.3f Vx: %5.3f  Vy: %5.3f \n", data[0], data[1], data[2]);
            if (data[1] < bestTime) {
                bestTime = data[0];
                topVX = data[1];
                topVY = data[2];
            }
        }
        
    } else {
        printf("Worker Process Begining to send results \n");
        double send [] = {bestTime, topVX, topVY};
        mpi_error =  MPI_Send(send, 3, MPI_DOUBLE, serverRank, tagResult, MPI_COMM_WORLD);
    }
    
    // Free Data
    free(velX);
    free(velY);
    
    // Print Results
    if (myrank == serverRank) {
      double  x = topVX - vS * cos(thetaS);
       double y = topVY - vS * sin(thetaS);
        printf("Best Time: %5.5f, vx %5.5f,  vy %5.5f, magnitude %5.5f theta %5.5f \n", bestTime, x, y, sqrt(pow(x,2) + pow(y,2)), atan2(y,x) * 180/M_PI);
        
        endtime = MPI_Wtime();
        totaltime = endtime - starttime;
        printf("Total Time = %5.5f\n", totaltime);
        
        // Print to File
        /*FILE * file = fopen("OutputParallel","w");
         fprintf(file, " %d \t\t %lu /n",m,total_t);
         fclose(file);
         */
        
        
    }
    
    MPI_Finalize();

    
    return 0;
}



