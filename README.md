# gradientRelLoc
A set of programs for calculating the relative location of seismic events using differential time measurements.

Compile by running the *compile_all.sh* script.  

```
sh compile_all.sh
```

and then, for a couple of events, e.g. DPRK3 and DPRK4, run using a script of the form *run_oneEventPairDiffVec.sh*  

```
#!/bin/sh
event1=DPRK3
event2=DPRK4
statphasefile=DPRK_ak135_slovecs.txt
CCtimesfile=DPRK_CC_times.txt
NSP=`wc ${statphasefile} | awk '{print $1}'`
cp ${statphasefile} gradientRelLoc.input
cat ${CCtimesfile} >> gradientRelLoc.input
./bin/gradientRelLoc ${NSP} ${event1} ${event2} < gradientRelLoc.input > out.txt
``` 
