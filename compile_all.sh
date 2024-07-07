#!/bin/sh
cd SUBS
make
cd ..
if test ! -r subslib.a
then
  echo "Failed to create library subslib.a"
fi
make gradientRelLoc
for file in \
   bin/gradientRelLoc 
do
  if test ! -r $file
  then
    echo Failed to create executable $file
  fi
  chmod 755 $file
done
chmod 755 run_oneEventPairDiffVec.sh
