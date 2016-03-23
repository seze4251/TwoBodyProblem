//
//  twoBodySlaveMaster.c
//  SlaveMaster
//
//  Created by Seth on 03/15/16.
//  Copyright © 2016 Seth. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "functions.h"

// Master Logic
void handleMasterInit(double * velX, double * velY, int * trash, int place [][2], int load, int tagX, int tagY, MPI_Status status, MPI_Request req []);
void handleMasterBody(double * bestTime, double * topVX, double * topVY, double * transfer, int place[][2], MPI_Status status);
void handleMasterRequestMore(double * velX, double * velY, int size, int trash[], int place[][2], int load, int tagX, int tagY, MPI_Status status, MPI_Request req []);

void handleMasterCompute(double * IC, double * velX, double * velY, double * currentTime, double * bestTime, double * topVX, double * topVY, int place[][2]);
void handleMasterFinishShort(int trash [], int tagFinilize, int nprocs, MPI_Request req []);
void handleMasterFinishLong(double * bestTime, double * topVX, double * topVY, double * transfer, int nprocs, int trash [], int place[][2], int hasData, int tagResult, int tagFinilize, MPI_Status status, MPI_Request req []);
void finish(int trash [], int tagFinilize, int nprocs, MPI_Request req []);

// Slave Logic
int handleSlaveInit (int * slaveCount, double ** velX, double ** velY, double ** velXtwo, double ** velYtwo, int * trash, int serverRank, int tagInit, int tagX, int tagY, int tagFinilize, int myrank, MPI_Status status);
void handleSlaveBody(int * slaveCount, double * IC, double * velX, double * velY, double * velXtwo, double * velYtwo, double * currentTime, double * bestTime, double * topVX, double * topVY, int serverRank, int tagX, int tagY, int tagResult, int tagFinilize, int tagMoreData, int myrank);
void handleServerFinish(double starttime, double endtime, double totaltime, double * bestTime, double * topVX, double * topVY, double * vS, int nprocs, int load, int size);

// Private Functions
int * setPlace( int place[][2], int status);
int * findPlace(int place[][2], int status);
void printPlace(int place[][2], int nprocs);

// Globals
const double mM = 7.34767309E22, mE = 5.97219E24, mS = 28833, rM = 1737100, rE = 6371000, G = 6.675E-11;
const double  thetaS = 50 * M_PI / 180, thetaM = 42.5 * M_PI/180, xE = 0, yE = 0, vEx = 0, vEy = 0;
double clear;

int main(int argc, char * argv[]) {
    // Initialize MPI
    MPI_Init(&argc, &argv);
    
    //Initialize Base Vars on all processes
    int nprocs, myrank, mpi_error;
    const int serverRank = 0;
    double starttime = 0, endtime = 0, totaltime = 0;
    
    // Determine # of procs and my rank
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    
    // Declare Constants
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
    double * velX = NULL, * velY = NULL;
    int i, j, count, load, size;
    
    if (myrank == serverRank) {
        // Read in data
        while (1) {
        //    printf(" #Load Distribution # Angles, # Velocities, #Velocity Start, #clearence \n");
            scanf("%d %lf %lf %lf %lf", &load, &numA, &numV, &vStart, &clear);
            if (numA > 0 && numV > 0 && vStart > 0 && clear >= 0 && load > 1) {
                break;
            } else {
                printf("Please Enter Three Numbers that Are greater than 0 \n");
            }
        }
        // Progress Tracker
      //  printf("Out of Scanloop\n");
        
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
        starttime = MPI_Wtime();
    }
    
    
    // Status and Request Variables
    MPI_Status status;
    MPI_Request req [nprocs];
    
    // Broadcast Clearence
    mpi_error = MPI_Ibcast(& clear, 1, MPI_DOUBLE, serverRank, MPI_COMM_WORLD, req);
    
    // Slave Master Variables
    int flag,  tagX = 7, tagY = 8, tagResult =9, tagInit = 3, tagMoreData = 4, tagFinilize = 5;
    int place [nprocs][2], trash[] = {-1};
    
    // Slave Specific Variables
    double * velXtwo = NULL, * velYtwo = NULL;
    
    //  Init Place
    if (myrank == serverRank) {
      //  printf("tagX = 7, tagY = 8, tagResult =9, tagInit = 3, tagMoreData = 4, tagFinilize = 5\n");
        for (j = 0; j < 2; j++) {
            for (i = 0; i < nprocs; i++) {
                if (i == 0 && j == 0) {
                    place[i][j] = 0;
                } else {
                    place[i][j] = -1;
                }
            }
        }
    }
    
    // Declare Variables Needed for Parallel Optimization
    double currentTime, bestTime = 0, topVX, topVY, transfer[3];
    
    // Main Loop
    i = 0;
    //printf("Entering Main loop rank: %d \n",myrank);
    while (1) {
        if (myrank == serverRank) {
            // Master
            if (place[0][0] != size) {
                MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, & flag, & status);
                //   printf("Master: flag: %d  Tag: %d  From: %d\n", flag, status.MPI_TAG, status.MPI_SOURCE);
                
                if (flag == 1) {
                    if (status.MPI_TAG == tagInit ) {
                        // Master Handle Initialize Request
                       // printf("Master: Enter Init\n");
                        handleMasterInit(velX, velY, trash, place, load, tagX, tagY, status, req);
                        
                    } else if ( status.MPI_TAG == tagResult ) {
                        // Master Data Recive Request
                     //   printf("Master: Enter Data Recive\n");
                        handleMasterBody(& bestTime, & topVX, &topVY, transfer, place, status);
                        
                    } else if(status.MPI_TAG == tagMoreData) {
                        // Master Handle Request for more data
                       // printf("MASTER: Enter Request More Data \n");
                        handleMasterRequestMore(velX, velY, size, trash, place, load, tagX, tagY, status, req);
                        
                    } else {
                        // Master Log Unexpected Tag
                        printf("Error: Log not implemented yet, Message Tag inncorrect \n");
                    }
                    
                } else if( flag == 0) {
                    // Master Compute
                      //printf("MASTER: Enter Compute\n");
                      handleMasterCompute(IC, velX, velY, & currentTime, & bestTime, & topVX, & topVY, place);
                }
                
            } else {
                // Master Finish
              //  printf("MASTER: Start FINISH Section\n");
                int hasData = 0, j;
                
                for (i = 1; i < nprocs; i++) {
                    for (j = 0; j < 2; j++) {
                        if (place[i][j] != -1) {
                            hasData += 1;
                        }
                    }
                }
                
                if (!hasData) {
                    // Case 1: Master did all the work (small matrix)
                 //   printf("MASTER: FINISH SHORT\n");
                    handleMasterFinishShort(trash, tagFinilize, nprocs, req);
                    break;
                } else {
                    // Case 2: Others Have data
                   // printf("MASTER: FINISH LONG\n");
                    handleMasterFinishLong(& bestTime, & topVX, & topVY, transfer, nprocs, trash, place, hasData, tagResult, tagFinilize, status, req);
                    break;
                }
            }
        } else {
            // All Slaves
            int slaveCount;
            if (i == 0) {
                //Initialize Slave
                i = 1;
               // printf("Slave: %d enter Init\n", myrank);
                int checkExit = handleSlaveInit(& slaveCount, & velX, & velY, & velXtwo, & velYtwo, trash, serverRank, tagInit, tagX, tagY, tagFinilize, myrank, status);
                if (checkExit) {
                    break;
                }
            }
           // printf("Slave: %d waiting on Clearance\n",myrank);
            MPI_Wait(req, MPI_STATUS_IGNORE); // Wait on Clearence
           // printf("Slave: %d has Clearence of %f -> Entering Body Slave Loop\n",myrank, clear);
            
            handleSlaveBody(& slaveCount, IC, velX, velY, velXtwo, velYtwo, & currentTime, & bestTime, & topVX, & topVY, serverRank, tagX, tagY, tagResult, tagFinilize, tagMoreData, myrank);
            break;
        }
    }
    
   // printf("I rank %d made it out of the infinite loop\n", myrank);
    if (myrank == serverRank) {
        handleServerFinish(starttime, endtime, totaltime, & bestTime, & topVX, & topVY, & vS, nprocs, load, size);
        
    }
    
    MPI_Finalize();
    return 0;
}


/*
 * Start Sub Function Definitions
 *
 */

void handleMasterInit(double * velX, double * velY, int * trash, int place [][2], int load, int tagX, int tagY, MPI_Status status, MPI_Request req []) {
    int mpi_error;
   // printf("SERVER: Initilize Message from Process %d \n",status.MPI_SOURCE);
    mpi_error = MPI_Irecv(trash, 1, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, req);
    mpi_error = MPI_Isend(velX + place[0][0], load, MPI_DOUBLE, status.MPI_SOURCE, tagX, MPI_COMM_WORLD, req);
    mpi_error = MPI_Isend(velY + place[0][0], load, MPI_DOUBLE, status.MPI_SOURCE, tagY, MPI_COMM_WORLD, req);
    MPI_Request_free(req);
    
    place[status.MPI_SOURCE][0] = place[0][0];
    place[0][0] += load;
}

// Returns a pointer to the Lowest NON -1 entry of place
// If both are 0, print an error
int * findPlace(int place[][2], int status) {
    int * holder = NULL;
    if (place[status][0] != -1 && place[status][0] < place[status][1]) {
        holder = & place[status][0];
    } else if( place[status][1] != -1 && place[status][0] > place[status][1]) {
        holder = & place[status][1];
    } else {
        printf("Function Find lowest Place requires both entries be non 0\n");
    }
    return holder;
}

void handleMasterBody(double * bestTime, double * topVX, double * topVY, double * transfer, int place[][2], MPI_Status status) {
    int count, mpi_error, * holder;
    mpi_error = MPI_Get_count(&status, MPI_DOUBLE, &count);
    holder = findPlace(place, status.MPI_SOURCE);
    
    // Blocking Recv
    mpi_error = MPI_Recv(transfer, 3, MPI_DOUBLE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, & status);
   // printf("MASTER: RECV-- BestTime %f   TopVX  %f  TopVY %f\n",transfer[0], transfer[1], transfer[2]);
    if (*bestTime == 0) {
        * bestTime = transfer[0];
        * topVX = transfer[1];
        * topVY = transfer[2];
    } else if (transfer [0] < * bestTime && transfer[0] > 0) {
        * bestTime = transfer[0];
        * topVX = transfer[1];
        * topVY = transfer[2];
    }
    
    *holder = -1;
}

// Return either place with a 0 or the lowest place
int * setPlace( int place[][2], int status) {
    int * holder = NULL;
    
    if (place[status][0] < place[status][1]) {
        holder = & place[status][0];
        
    } else if (place[status][0] > place[status][1]) {
        holder = & place[status][1];
        
    } else {
        printf("ERROR: in setPlace\n");
    }
    return holder;
}


void handleMasterRequestMore(double * velX, double * velY, int size, int trash[], int place[][2], int load, int tagX, int tagY, MPI_Status status, MPI_Request req []) {
    int mpi_error;
    
    if (place[0][0] + load >= size) {
        load = size - place[0][0];
    }
    
    // Recive message and throw it away
    mpi_error = MPI_Irecv(trash, 1, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, req);
    MPI_Request_free(req);
    
    // Send the new Message
    mpi_error = MPI_Isend(velX + place[0][0], load, MPI_DOUBLE, status.MPI_SOURCE, tagX, MPI_COMM_WORLD, req);
    mpi_error = MPI_Isend(velY + place[0][0], load, MPI_DOUBLE, status.MPI_SOURCE, tagY, MPI_COMM_WORLD, req);
    
    int * holder = setPlace(place, status.MPI_SOURCE);
    * holder = place[0][0];
    place[0][0] += load;
    MPI_Request_free(req);
}

void handleMasterCompute(double * IC, double * velX, double * velY, double * currentTime, double * bestTime, double * topVX, double * topVY, int place[][2]) {
    
    IC[6] = velX[place[0][0]];
    IC[7] = velY[place[0][0]];
    int ret = fowardEuler(currentTime, IC);
    
    if (ret == 0) {
       // printf("MASTER COMPUTE ret = 0, currentT = %f, topVX %f\n", *currentTime, velX[place[0][0]]);
        if (*bestTime == 0) {
            * bestTime = * currentTime;
            * topVX = velX[place[0][0]];
            * topVY = velY[place[0][0]];
        } else if (* currentTime < * bestTime && * currentTime > 0) {
            * bestTime = * currentTime;
            * topVX = velX[place[0][0]];
            * topVY = velY[place[0][0]];
        }
    }
    
    place[0][0] += 1;
}

void finish(int trash [], int tagFinilize, int nprocs, MPI_Request req []) {
    // Send termination message to all processes
    int i, mpi_error;
    for (i = 1; i < nprocs; i++) {
        mpi_error = MPI_Isend(trash, 1, MPI_INT, i, tagFinilize, MPI_COMM_WORLD, req + i);
    }
    
    // Wait for all messages to go through to avoid seg fault
    MPI_Waitall(nprocs -1, req + 1, MPI_STATUS_IGNORE);
 //   printf("SERVER: Completed Finish Function\n");
}

void handleMasterFinishShort(int trash [], int tagFinilize, int nprocs, MPI_Request req []) {
    finish(trash, tagFinilize, nprocs, req);
   // printf("SERVER: Finish Short\n");
}

//Cases: Slaves have 1 data or 2 sets of data
void handleMasterFinishLong(double * bestTime, double * topVX, double * topVY, double * transfer, int nprocs, int trash [], int place[][2], int hasData, int tagResult, int tagFinilize, MPI_Status status, MPI_Request req []) {
    
    int  mpi_error, count, i;
    
    for (i = 0; i < hasData; i++) {
        MPI_Probe(MPI_ANY_SOURCE, tagResult, MPI_COMM_WORLD, & status);
        mpi_error = MPI_Get_count(&status, MPI_DOUBLE, &count);
        
        //Blocking Recv
        mpi_error = MPI_Recv(transfer, count, MPI_DOUBLE, status.MPI_SOURCE, tagResult, MPI_COMM_WORLD, & status);
        if (transfer [0] < * bestTime && transfer[0] > 0) {
            * bestTime = transfer[0];
            * topVX = transfer[1];
            * topVY = transfer[2];
        }
    }
    
    finish(trash, tagFinilize, nprocs, req);
}


int handleSlaveInit(int * slaveCount, double ** velX, double ** velY, double ** velXtwo, double ** velYtwo, int * trash, int serverRank, int tagInit, int tagX, int tagY, int tagFinilize, int myrank, MPI_Status status) {
    //printf("SLAVE: %d \n",myrank);
    
    // Blocking Send
    int mpi_error = MPI_Send(trash, 1, MPI_INT, serverRank, tagInit, MPI_COMM_WORLD);
    
    MPI_Probe(serverRank, MPI_ANY_TAG, MPI_COMM_WORLD, & status);
    
    //printf("SLAVE: %d Finished Probe tag: %d \n",myrank, status.MPI_TAG);
    
    if (status.MPI_TAG == tagX || status.MPI_TAG == tagY) {
        mpi_error = MPI_Get_count(&status, MPI_DOUBLE, slaveCount);
        * velX = (double *) malloc(*slaveCount * sizeof(double));
        * velY = (double *) malloc(*slaveCount * sizeof(double));
        * velXtwo = (double *) malloc(*slaveCount * sizeof(double));
        * velYtwo = (double *) malloc(*slaveCount * sizeof(double));
        
        //Blocking Recv
        mpi_error = MPI_Recv(* velX, *slaveCount, MPI_DOUBLE, serverRank, tagX, MPI_COMM_WORLD, & status);
        mpi_error = MPI_Recv(* velY, *slaveCount, MPI_DOUBLE, serverRank, tagY, MPI_COMM_WORLD, & status);
        return 0;
        
    } else {
        mpi_error = MPI_Recv(trash, 1, MPI_INT, serverRank, tagFinilize, MPI_COMM_WORLD, & status);
        return 1;
    }
}


void handleSlaveBody(int * slaveCount, double * IC, double * velX, double * velY, double * velXtwo, double * velYtwo, double * currentTime, double * bestTime, double * topVX, double * topVY, int serverRank, int tagX, int tagY, int tagResult, int tagFinilize, int tagMoreData, int myrank) {
    
    MPI_Request req[4];
    MPI_Status status[1];
    int calc, rem, mpi_error, i = 0, j = 0, flag = 0, ret, trash[] = {-1};
    double transfer [3], * holder;
    
    // Need to start loop with mtxA_one being full and mtxC_one initialized to 0
    // j is loop counter
    
    while (1) {
        // Debug Stuff
        //  printf("SLAVE: %d Top of Loop j = %d\n",myrank, j);
        for (i = 0; i < *slaveCount; i++) {
            //    printf("Vx: %5.5f  Vy: %5.5f\n", velX[i], velY[i]);
        }
        
        // Determine sizes of data
        calc = *slaveCount / 2;
        rem = *slaveCount % 2;
        
        // Ask For More Data
        mpi_error = MPI_Isend(trash, 1, MPI_INT, serverRank, tagMoreData, MPI_COMM_WORLD, req);  // Ask for more data req
        
        // Reset BestTime
        * bestTime = 0;
        
        // Compute First Data
        for (i = 0; i < calc; i++) {
            IC[6] = velX[i];
            IC[7] = velY[i];
            ret = fowardEuler(currentTime, IC);
            if (ret == 0) {
                if (* bestTime == 0) {
                    * bestTime = * currentTime;
                    * topVX = velX[i];
                    * topVY = velY[i];
                } else if (* currentTime < *bestTime) {
                    * bestTime = * currentTime;
                    * topVX = velX[i];
                    * topVY = velY[i];
                }
            }
        }
        
        // printf("Slave %d: Completed first compute j = %d\n",myrank, j);
        // Attempt to Recieve Information
        MPI_Iprobe(serverRank, tagX, MPI_COMM_WORLD, & flag, status);
        
        if (flag == 1) {
            //    printf("SLAVE: Recv X Short\n");
            mpi_error = MPI_Get_count( status, MPI_DOUBLE, slaveCount);
            mpi_error = MPI_Irecv(velXtwo, * slaveCount, MPI_DOUBLE, serverRank, tagX, MPI_COMM_WORLD, req + 1); // req + 1 Recv VelX
            mpi_error = MPI_Irecv(velYtwo, * slaveCount, MPI_DOUBLE, serverRank, tagY, MPI_COMM_WORLD, req + 2); // req + 2 Recv VelY
            
        }
        
        
        // Calculate Rest of Data
        for (i = calc; i < 2 * calc + rem; i++) {
            IC[6] = velX[i];
            IC[7] = velY[i];
            ret = fowardEuler(currentTime, IC);
            if (ret == 0) {
                if (* bestTime == 0) {
                    * bestTime = * currentTime;
                    * topVX = velX[i];
                    * topVY = velY[i];
                } else if (* currentTime < *bestTime) {
                    * bestTime = * currentTime;
                    * topVX = velX[i];
                    * topVY = velY[i];
                }
            }
        }
        
        //printf("Slave %d: Completed Second Compute j = %d\n",myrank, j);
        
        
        // Wait for previous send mtxC request to go through
        if (j != 0) {
            MPI_Wait(req+3, MPI_STATUS_IGNORE);
        }
        //   printf("Slave %d: Completed req + 3 Wait j = %d\n",myrank, j);
        
        // Set Up Result
        transfer[0] = *bestTime;
        transfer[1] = *topVX;
        transfer[2] = *topVY;
        
        
        // Send Answer Back
        mpi_error = MPI_Isend(transfer, 3, MPI_DOUBLE, serverRank, tagResult, MPI_COMM_WORLD, req + 3); // req + 3
        
        MPI_Wait(req, MPI_STATUS_IGNORE); // Wait for the more data request to be received
        
        // Finish Recive or block Recv
        if (flag == 1) {
            // Finish Recv
            // printf("Short Recv\n");
            MPI_Wait(req + 1, MPI_STATUS_IGNORE);
            MPI_Wait(req + 2, MPI_STATUS_IGNORE);
            
        } else {
            //printf("Long Recv\n");
            //Blocking Probe -- Need not to touch (status + 0)
            MPI_Probe(serverRank, MPI_ANY_TAG, MPI_COMM_WORLD, status);
            
            if (status[0].MPI_TAG == tagX || status[0].MPI_TAG == tagY) {
                //  printf("BAD BAD Here\n");
                mpi_error = MPI_Get_count(status, MPI_DOUBLE, slaveCount);
                // Blocking Recv
                mpi_error = MPI_Recv(velXtwo, * slaveCount, MPI_DOUBLE, serverRank, tagX, MPI_COMM_WORLD, status);
                mpi_error = MPI_Recv(velYtwo, * slaveCount, MPI_DOUBLE, serverRank, tagY, MPI_COMM_WORLD, status);
                //printf("SLAVE: %d Long Recv!!!!!!! j = %d\n",myrank, j);
                for (i = 0; i < *slaveCount; i++) {
                    //  printf("Vx: %5.5f  Vy: %5.5f\n", velXtwo[i], velYtwo[i]);
                }
                
            } else if(status[0].MPI_TAG == tagFinilize) {
                //printf("Slave: Recv Tag Finialize\n");
                
                MPI_Wait(req + 3, MPI_STATUS_IGNORE); //
                MPI_Iprobe(serverRank, tagX, MPI_COMM_WORLD, & flag, status + 1);
                
                if (flag == 1) {
                    // Case not supported
                    printf("Unsupported Cased occured Slave myrank: %d\n", myrank);
                }
                //  printf("SLAVE: %d About to enter Final Recv\n",myrank);
                mpi_error = MPI_Recv(trash, 1, MPI_INT, serverRank, tagFinilize, MPI_COMM_WORLD, status);
                break;
            }
        }
        
        //    printf("Slave %d: Completed Recv i = %d\n",myrank, i);
        //Switch Pointers
        holder = velY;
        velY = velYtwo;
        velYtwo = holder;
        
        holder = velX;
        velX = velXtwo;
        velXtwo = holder;
        j++;
    }
    
    free(velX);
    free(velY);
    free(velXtwo);
    free(velYtwo);
}


void handleServerFinish(double starttime, double endtime, double totaltime, double * bestTime, double * topVX, double * topVY, double * vS, int nprocs, int load, int size) {
    //Time
    double  x = * topVX - *vS * cos(thetaS);
    double y = * topVY - *vS * sin(thetaS);
    printf("Best Time: %5.5f, vx %5.5f,  vy %5.5f, magnitude %5.5f theta %5.5f \n", * bestTime, x, y, sqrt(pow(x,2) + pow(y,2)), atan2(y,x) * 180/M_PI);
    endtime = MPI_Wtime();
    totaltime = endtime - starttime;
    printf("Total Running Time: %5.3f\n",totaltime);
    // Print to File
     FILE * file = fopen("slaveSpeed.txt","a");
     fprintf(file, " %d \t\t %f  \t %d \t %d \n",nprocs, totaltime, size, load);
     fclose(file);
}

void printPlace(int place[][2], int nprocs) {
    int i, j;
    for (i = 0; i < nprocs; i ++) {
        for (j = 0 ; j < 2; j++) {
            printf("%d\t",place[i][j]);
        }
        printf("\n");
    }
}




