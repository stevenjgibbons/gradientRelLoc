#!/bin/sh
event1=DPRK3
event2=DPRK4
statphasefile=DPRK_ak135_slovecs.txt
# CCtimesfile=DPRK3.DPRK4.CC_times_multiple.txt
CCtimesfile=dummy.txt
NSP=`wc ${statphasefile} | awk '{print $1}'`
NREAL=50
cp ${statphasefile} mc_gradientRelLoc.input
cat ${CCtimesfile} >> mc_gradientRelLoc.input
./bin/mc_gradientRelLoc ${NSP} ${NREAL} ${event1} ${event2} < mc_gradientRelLoc.input 
