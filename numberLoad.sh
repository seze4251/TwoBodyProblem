#!/bin/bash
# Purpose: Creates Num Processes Vs Performance

load="3"
nAngles="720"
nVelocity="10"
vStart="99"
clearence="100"
counter="2"

while [ $counter -lt 20 ]
do
#&load, &numA, &numV, &vStart, &clear
mpirun -np 2 ./twoBodySlave  << EOF
$load $nAngles $nVelocity $vStart $clearence
EOF


((load+=3))
((counter+=1))
done
