#!/bin/bash
for i in 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50
do
   val=`cat AM${i}/lambdaValue.txt`
   val2=4.375
   result=$(expr $val*$val2 | bc)
   touch AM${i}/lambdaValue2.txt
   echo "0$result" > AM${i}/lambdaValue2.txt
done
