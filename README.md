# gradientRelLoc
A set of programs for calculating the relative location of seismic events using differential time measurements.  

This code is available on Zenodo with the DOI https://zenodo.org/doi/10.5281/zenodo.12680277  

<a href="https://zenodo.org/doi/10.5281/zenodo.12680277"><img src="https://zenodo.org/badge/825336703.svg" alt="DOI"></a>  

[![SQAaaS badge shields.io](https://img.shields.io/badge/sqaaas%20software-silver-lightgrey)](https://api.eu.badgr.io/public/assertions/WYpfzBaHTOCew0H1QWAdpA "SQAaaS silver badge achieved")  
[![SQAaaS badge](https://github.com/EOSC-synergy/SQAaaS/raw/master/badges/badges_150x116/badge_software_silver.png)](https://api.eu.badgr.io/public/assertions/WYpfzBaHTOCew0H1QWAdpA "SQAaaS silver badge achieved")  

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
