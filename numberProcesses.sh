#!/bin/bash
# Purpose: Creates Num Processes Vs Performance

load="50"
nAngles="720"
nVelocity="10"
vStart="99"
clearence="100"
counter="2"

./twoBody << EOF
$nAngles $nVelocity $vStart $clearence
EOF

while [ $counter -lt 4 ]
do
#&load, &numA, &numV, &vStart, &clear
mpirun -np $counter ./twoBodySlave  << EOF
$load $nAngles $nVelocity $vStart $clearence
EOF

#&numA, &numV, &vStart, &clear
mpirun -np $counter ./paraSplit << EOF
$nAngles $nVelocity $vStart $clearence
EOF

((counter+=1))
done
