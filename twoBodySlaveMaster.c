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
#include <unistd.h>

// Master Logic
void handleMasterInit(double * velX, double * velY, int * trash, int place [][2], int load, int tagX, int tagY, Status status, Request req []);
void handleMasterBody(double * bestTime, double * topVX, double * topVY, double * transfer, int place[][2], MPI_Status status);
handleMasterRequestMore(double * velX, double * velY, int size, int trash[], int place[][2], int load, int tagX, int tagY, MPI_Status status, MPI_Request req []);

void handleMasterCompute(double * IC, double * velX, double * velY, double * currentTime, double * bestTime, double * topVX, double * topVY, int place[][2]);
void handleMasterFinishShort(int trash [], int tagFinilize, int nprocs, MPI_Request req []);
void handleMasterFinishLong(double * bestTime, double * topVX, double * topVY, double * transfer, int nprocs, int trash [], int place[][2], int hasData, int tagResult, int tagFinilize, MPI_Status status, MPI_Request req []);
void finish(int trash [], int tagFinilize, int nprocs, MPI_Request req []);

// Slave Logic
int handleSlaveInit(double ** velX, double ** velY, double ** velXtwo, double ** velYtwo, int * slaveCount int * trash, int serverRank, int tagInit, int tagX, int tagY, int tagFinilize, int myrank, MPI_Status status)
void handleSlaveBody(double * IC, double * velX, double * velY, double * velXtwo, double * velYtwo, double * currentTime, double * bestTime, double * topVX, double * topVY, int * slaveCount, int serverRank, int tagX, int tagY, int tagResult, int tagFinilize, int tagMoreData, int myrank);
void handleServerFinish(matrix * mtxA, matrix * mtxB, matrix * mtxC, int n, int p, double starttime, double endtime, double totaltime);

// Private Functions
int * setPlace( int place[][2], int status);
int * findPlace(int place[][2], int status);
void printPlace(int place[][2], int nprocs);

// Globals
const double mM = 7.34767309E22, mE = 5.97219E24, mS = 28833, rM = 1737100, rE = 6371000, G = 6.675E-11;
const double  thetaS = 50 * M_PI / 180, thetaM = 42.5 * M_PI/180, clear = 10000, xE = 0, yE = 0, vEx = 0, vEy = 0;

int main(int argc, char * argv[]) {
    // Initialize MPI
    MPI_Init(&argc, &argv);
    
    //Initialize Base Vars on all processes
    int nprocs, myrank, mpi_error;
    const int serverRank = 0;
    double starttime, endtime, totaltime;
    
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
    double * velX, * velY, * velXtwo, * velYtwo;
    int i, count, load, size = (int) numA * numV;
    
    if (myrank == serverRank) {
        
        while (1) {
            printf("Please enter 4 #'s larger than 0,  #angles, #Velocities and start velocity and load size \n");
            int check = scanf("%lf %lf %lf %d",&numA, &numV, &vStart, & load);
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
    
    // Status and Request Variables
    MPI_Status status;
    MPI_Request req [nprocs];
    
    // Slave Master Variables
    int i = 0, * slaveCount;
    int flag,  tagX = 7, tagY = 8, tagResult =9, tagInit = 3, tagMoreData = 4, tagFinilize = 5;
    int place [nprocs][2], trash [1];
    trash[0] = -1;
    
    //  Init Place
    if (myrank == serverRank) {
        int j;
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
    while (1) {
        if (myrank == serverRank) {
            // Master
            if (place[0][0] != size) {
                MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, & flag, & status);
                
                if (flag == 1) {
                    if (status.MPI_TAG == tagInit ) {
                        // Master Handle Initialize Request
                        handleMasterInit(velX, velY, trash, place, load, tagX, tagY, status, req);
                        
                    } else if ( status.MPI_TAG == tagResult ) {
                        // Master Data Recive Request
                        handleMasterBody(& bestTime, & topVX, &topVY, transfer, place, status);
                        
                    } else if(status.MPI_TAG == tagMoreData) {
                        // Master Handle Request for more data
                        handleMasterRequestMore(velX, velY, size, trash, place, load, tagX, tagY, status, req);
                        
                    } else {
                        // Master Log Unexpected Tag
                        printf("Error: Log not implemented yet, Message Tag inncorrect \n");
                    }
                    
                } else if( flag == 0) {
                    // Master Compute
                    handleMasterCompute(IC, velX, velY, & currentTime, & bestTime, & topVX, & topVY, place);
                }
                
            } else {
                // Master Finish
                printf("MASTER: Start FINISH Section\n");
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
                    handleMasterFinishShort(trash, tagFinilize, nprocs, req);
                    break;
                } else {
                    // Case 2: Others Have data
                    handleMasterFinishLong(& bestTime, & topVX, & topVY, transfer, nprocs, trash, place, hasData, tagResult, tagFinilize, status, req);
                    break;
                }
            }
        } else {
            // All Slaves
            if (i == 0) {
                //Initialize Slave
                i = 1;
                int checkExit = handleSlaveInit(& velX, & velY, & velXtwo, % velYtwo, slaveCount, trash, serverRank, tagInit, tagX, tagY, tagFinilize, myrank, status);
                if (checkExit) {
                    break;
                }
            }
            
            handleSlaveBody(IC, velX, velY, velXtwo, velYtwo, & currentTime, & bestTime, & topVX, & topVY, slaveCount, serverRank, tagX, tagY, tagResult, tagFinilize, tagMoreData, myrank);
            break;
        }
    }
    
    printf("I rank %d made it out of the infinite loop\n", myrank);
    if (myrank == serverRank) {
        handleServerFinish(mtxA, mtxB, mtxC, n, p, starttime, endtime, totaltime);
    }
    
    MPI_Finalize();
    if (myrank == serverRank) {
        printf("Best Time: %5.3f  MagV: %5.3f,  theta: %5.3f degrees", bestTime, sqrt(pow(topVX,2) + pow(topVY,2)), 180 / M_PI *atan2(topVY, topVX));
    }
    return 0;
}


/*
 * Start Sub Function Definitions
 *
 */

void handleMasterInit(double * velX, double * velY, int * trash, int place [][2], int load, int tagX, int tagY, Status status, Request req []) {
    int mpi_error;
    printf("SERVER: Initilize Message from Process %d \n",status.MPI_SOURCE);
    mpi_error = MPI_Irecv(trash, 1, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, req);
    mpi_error = MPI_Isend(velX[place[0][0]], load, MPI_DOUBLE, status.MPI_SOURCE, tagX, MPI_COMM_WORLD, req);
    mpi_error = MPI_Isend(velY[place[0][0]], load, MPI_DOUBLE, status.MPI_SOURCE, tagY, MPI_COMM_WORLD, req);
    MPI_Request_free(req);
    
    place[status.MPI_SOURCE][0] = place[0][0];
    place[0][0] += load;
}

// Returns a pointer to the Lowest NON -1 entry of place
// If both are 0, print an error
int * findPlace(int place[][2], int status) {
    int * holder;
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
    
    if (transfer [0] < * bestTime && transfer[0] > 0) {
        * bestTime = data[0];
        * topVX = data[1];
        * topVY = data[2];
    }

    *holder = -1;
}

// Return either place with a 0 or the lowest place
int * setPlace( int place[][2], int status) {
    int * holder;
    
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
    mpi_error = MPI_Isend(velX[place[0][0]], load, MPI_DOUBLE, status.MPI_SOURCE, tagX, MPI_COMM_WORLD, req);
    mpi_error = MPI_Isend(velY[place[0][0]], load, MPI_DOUBLE, status.MPI_SOURCE, tagY, MPI_COMM_WORLD, req);
    
    int * holder = setPlace(place, status.MPI_SOURCE);
    * holder = place[0][0];
    place[0][0] += load;
    MPI_Request_free(req);
}

void handleMasterCompute(double * IC, double * velX, double * velY, double * currentTime, double * bestTime, double * topVX, double * topVY, int place[][2]) {
    IC[6] = velX[i];
    IC[7] = velY[i];
    int ret = fowardEuler(& currentTime, IC);
    if (ret == 0) {
        if (bestTime == 0) {
            * bestTime = currentTime;
            * topVX = velX[i];
            * topVY = velY[i];
        } else if (currentTime < bestTime) {
            * bestTime = currentTime;
            * topVX = velX[i];
            * topVY = velY[i];
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
    printf("SERVER: Completed Finish Function\n");
}

void handleMasterFinishShort(int trash [], int tagFinilize, int nprocs, MPI_Request req []) {
    finish(trash, tagFinilize, nprocs, req);
    printf("SERVER: Finish Short\n");
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
            * bestTime = data[0];
            * topVX = data[1];
            * topVY = data[2];
        }
    }
    
    finish(trash, tagFinilize, nprocs, req);
}


int handleSlaveInit(double ** velX, double ** velY, double ** velXtwo, double ** velYtwo, int * slaveCount int * trash, int serverRank, int tagInit, int tagX, int tagY, int tagFinilize, int myrank, MPI_Status status) {
    printf("SLAVE: %d Send Initialize\n",myrank);
    
    // Blocking Send
    int mpi_error = MPI_Send(trash, 1, MPI_INT, serverRank, tagInit, MPI_COMM_WORLD);
    
    MPI_Probe(serverRank, MPI_ANY_TAG, MPI_COMM_WORLD, & status);
    
    if (status.MPI_TAG == tagX || status.MPI_TAG == tagY) {
        mpi_error = MPI_Get_count(&status, MPI_DOUBLE, slaveCount);
        * velX = (double) malloc(*slaveCount * sizeof(double));
        * velY = (double) malloc(*slaveCount * sizeof(double));
        * velXtwo = (double) malloc(*slaveCount * sizeof(double));
        * velYtwo = (double) malloc(*slaveCount * sizeof(double));
        
        //Blocking Recv
        mpi_error = MPI_Recv(** velX, *slaveCount, MPI_DOUBLE, serverRank, tagX, MPI_COMM_WORLD, & status);
        mpi_error = MPI_Recv(** velY, *slaveCount, MPI_DOUBLE, serverRank, tagY, MPI_COMM_WORLD, & status);
        return 0;
        
    } else {
        mpi_error = MPI_Recv(trash, 1, MPI_INT, serverRank, tagFinilize, MPI_COMM_WORLD, & status);
        return 1;
    }
}


void handleSlaveBody(double * IC, double * velX, double * velY, double * velXtwo, double * velYtwo, double * currentTime, double * bestTime, double * topVX, double * topVY, int * slaveCount, int serverRank, int tagX, int tagY, int tagResult, int tagFinilize, int tagMoreData, int myrank) {
    
    MPI_Request req[4];
    MPI_Status status [2];
    
    int calc, rem, err, count, mpi_error, i = 0, flagX = 0, flagY = 0, j, ret;
    int trash[1];
    trash[0] = -1;
    double transfer [3];
    double * holder;
    
    // Need to start loop with mtxA_one being full and mtxC_one initialized to 0
    printf("SLAVE: %d entering main body loop\n",myrank);
    while (1) {
        // Determine sizes of data
        calc = * slaveCount / 2;
        rem = * slaveCount % 2;
        
        // Ask For More Data
        mpi_error = MPI_Isend(trash, 1, MPI_INT, serverRank, tagMoreData, MPI_COMM_WORLD, req);
        
        // Reset BestTime
        * bestTime = 0;
        
        // Wait for data to send if
        if (i != 0) {
            MPI_Wait(req + 1, MPI_STATUS_IGNORE);
            MPI_Wait(req + 2, MPI_STATUS_IGNORE);
        }
        // Compute First Data
        for (i = 0; j < calc; j++) {
            IC[6] = velX[i];
            IC[7] = velY[i];
            ret = fowardEuler(currentTime, IC);
            if (ret == 0) {
                if (* bestTime == 0) {
                   * bestTime = * currentTime;
                   * topVX = velX[i];
                   * topVY = velY[i];
                } else if (* currentTime < *bestTime) {
                   * bestTime = currentTime;
                   * topVX = velX[i];
                   * topVY = velY[i];
                }
            }
        }
        
        // Attempt Recv
        MPI_Iprobe(serverRank, tagX, MPI_COMM_WORLD, & flagX, status);
        
        if (flagX == 1) {
            mpi_error = MPI_Get_count( status, MPI_DOUBLE, slaveCount);
            mpi_error = MPI_Irecv(velXtwo, count, MPI_DOUBLE, serverRank, tagX, MPI_COMM_WORLD, req + 1); // req + 1
            
        }
        
        MPI_Iprobe(serverRank, tagY, MPI_COMM_WORLD, & flagY, status + 1);
        
        if (flagY == 1) {
            mpi_error = MPI_Get_count( status, MPI_DOUBLE, slaveCount);
            mpi_error = MPI_Irecv(velYtwo, count, MPI_DOUBLE, serverRank, tagY, MPI_COMM_WORLD, req + 2); // req + 2
            
        }
        
        // Calculate Rest of Data
        for (i = calc; j < 2 * calc + rem; j++) {
            IC[6] = velX[i];
            IC[7] = velY[i];
            ret = fowardEuler(currentTime, IC);
            if (ret == 0) {
                if (* bestTime == 0) {
                    * bestTime = * currentTime;
                    * topVX = velX[i];
                    * topVY = velY[i];
                } else if (* currentTime < *bestTime) {
                    * bestTime = currentTime;
                    * topVX = velX[i];
                    * topVY = velY[i];
                }
            }
        }

        // Wait for previous send mtxC request to go through
        if (i != 0) {
            MPI_Wait(req+3, MPI_STATUS_IGNORE);
        }
        
        
        // Set Up Result
        transfer[0] = *bestTime;
        transfer[1] = *topVX;
        transfer[2] = *topVY;
        
        // Send mtxC_two
        mpi_error = MPI_Isend(transfer, 3, MPI_DOUBLE, serverRank, tagResult, MPI_COMM_WORLD, req + 3); // req + 3
        
        MPI_Wait(req, MPI_STATUS_IGNORE); // Wait for the more data request to be received
        
        // Finish Recive or block Recv
        if (flagY == 1 && flagX == 1 ) {
            // Finish Recv
            MPI_Wait(req + 1, MPI_STATUS_IGNORE);
            MPI_Wait(req + 2, MPI_STATUS_IGNORE);
            
        } else if (flagY == 1 && flagX == 0) {
            MPI_Wait(req + 2, MPI_STATUS_IGNORE);
            
        } else if (flagY == 0 && flagX == 1) {
            MPI_Wait(req + 1, MPI_STATUS_IGNORE);
            
        } else {
            //Blocking Probe -- Need not to touch (status + 0)
            MPI_Probe(serverRank, MPI_ANY_TAG, MPI_COMM_WORLD, status);
            
            if (status[0].MPI_TAG == tagMoreData) {
                mpi_error = MPI_Get_count(status, MPI_DOUBLE, slaveCount);
                // Blocking Recv
                
                mpi_error = MPI_Irecv(velXtwo, *slaveCount, MPI_DOUBLE, serverRank, tagX, MPI_COMM_WORLD, req + 1);
                mpi_error = MPI_Irecv(velYtwo, *slaveCount, MPI_DOUBLE, serverRank, tagY, MPI_COMM_WORLD, req + 2);
                
            } else if(status[0].MPI_TAG == tagFinilize) {
                //Finishing Cases for slaves:
                // Normal - have a finishing message with no A waiting
                // Bad - Have a finishing message + have mtxA data waiting
                
                MPI_Wait(req + 2, MPI_STATUS_IGNORE); //
                MPI_Iprobe(serverRank, tagA, MPI_COMM_WORLD, & flagX, status + 1);
                
                if (flagX == 1) {
                    // Case not supported
                    printf("Unsupported Cased occured Slave myrank: %d\n", myrank);
                }
                    printf("SLAVE: %d About to enter Final Recv\n",myrank);
                mpi_error = MPI_Recv(trash, 1, MPI_INT, serverRank, tagFinilize, MPI_COMM_WORLD, status);
                break;
            }
        }
        
        //Switch Pointers
        holder = velY;
        veylY = velYtwo;
        velYtwo = holder;
        
        holder = velX;
        velX = velXtwo;
        velXtwo = holder;
        i++;
    }
    free(velX);
    free(velY);
    free(velXtwo);
    free(velYtwo);
}




void handleServerFinish(matrix * mtxA, matrix * mtxB, matrix * mtxC, int n, int p, double starttime, double endtime, double totaltime) {
    //Time
    endtime = MPI_Wtime();
    totaltime = endtime - starttime;
    printf("Total Running Time: %5.3f\n",totaltime);
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




