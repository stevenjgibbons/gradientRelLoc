#!/bin/sh
event1=DPRK3
event2=DPRK4
statphasefile=DPRK_ak135_slovecs.txt
CCtimesfile=DPRK_CC_times.txt
NSP=`wc ${statphasefile} | awk '{print $1}'`
cp ${statphasefile} gradientRelLoc.input
cat ${CCtimesfile} >> gradientRelLoc.input
./bin/gradientRelLoc ${NSP} ${event1} ${event2} < gradientRelLoc.input > out.txt
