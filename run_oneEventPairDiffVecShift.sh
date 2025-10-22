#!/bin/sh
# -0.94237    1.55647
DX0KM=-0.94237
DY0KM=1.55647
event1=DPRK2
event2=DPRK6
statphasefile=DPRK_ak135_slovecs.txt
# CCtimesfile=DPRK_CC_times.txt
CCtimesfile=manual_DPRK2_DPRK6.txt
NSP=`wc ${statphasefile} | awk '{print $1}'`
cp ${statphasefile} shift_gradientRelLoc.input
cat ${CCtimesfile} >> shift_gradientRelLoc.input
./bin/shift_gradientRelLoc ${NSP} ${event1} ${event2} ${DX0KM} ${DY0KM} < shift_gradientRelLoc.input > out.txt
