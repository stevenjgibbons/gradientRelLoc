#!/bin/sh
# set -x
scriptname=./solve_all_pairs.sh
if [ $# != 1 ]
then
  echo
  echo "USAGE: "
  echo "bash $scriptname   CCtimesfile "
  echo "bash $scriptname   DPRK_CC_times.txt "
  echo
  exit 1
fi
#
filename=$1
if test ! -r ${filename}
then
  echo No file ${filename} found ...
  exit 1
fi
statphasefile=DPRK_ak135_slovecs.txt
if test ! -r ${statphasefile}
then
  echo No file ${statphasefile} found ...
  exit 1
fi
NSP=`wc ${statphasefile} | awk '{print $1}'`
cp ${statphasefile} gradientRelLoc.input
cat ${filename} >> gradientRelLoc.input
stem=`basename $filename .txt`_solutions
outfile=${stem}.txt
if test -r ${outfile}
then
  rm ${outfile}
fi
touch ${outfile}
#
for event1 in DPRK1 DPRK2 DPRK3 DPRK4 DPRK5 DPRK6
do
  for event2 in DPRK1 DPRK2 DPRK3 DPRK4 DPRK5 DPRK6
  do
    if [ "$event1" != "$event2" ]
    then
      # echo $event1 $event2
      ./bin/gradientRelLoc ${NSP} ${event1} ${event2} < gradientRelLoc.input > out.txt
      grep Result out.txt | tail -n 1 | cut -c9-70 >> ${outfile}
    fi
  done
done
cat ${outfile}
