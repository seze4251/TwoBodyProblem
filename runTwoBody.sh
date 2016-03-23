#!/bin/bash
# Runs Full Project

#Initializes Files
gcc -O3 -c functions.c
make
make -f makeParaSplit
make -f makeSlaveMast

rm *.txt

# Run script for Number of Processes vs Compute Time for fixed load
./numberProcesses.sh
mv slaveSpeed.txt a.txt
mv distParalellSpeed.txt b.txt
mv SingleCore c.txt


# Run script for Number of Tasks vs Computetime for Fixed amount of 2 cores
./numberTasks.sh
mv slaveSpeed.txt d.txt
mv distParalellSpeed.txt e.txt
mv SingleCore f.txt

# Runs script for Amount of Load vs Compute time for 2 cores and large load
./numberLoad.sh

mv slaveSpeed.txt g.txt
