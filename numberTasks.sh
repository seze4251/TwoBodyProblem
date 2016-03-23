#!/bin/bash
# Purpose: Creates Num Processes Vs Performance

load="50"
nAngles="4"
nVelocity="10"
vStart="99"
clearence="100"
counter="1"



while [ $counter -lt 7 ]
do

((nAngles=nAngles*counter))

#&load, &numA, &numV, &vStart, &clear
mpirun -np 2 ./twoBodySlave  << EOF
$load $nAngles $nVelocity $vStart $clearence
EOF

#&numA, &numV, &vStart, &clear
mpirun -np 2 ./paraSplit << EOF
$nAngles $nVelocity $vStart $clearence
EOF


./twoBody << EOF
$nAngles $nVelocity $vStart $clearence
EOF


((counter+=1))
done
