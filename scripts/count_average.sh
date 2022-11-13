#!/bin/bash

count=0;
total=0;

for i in $(awk ' { print $2; }' $1 )
   do
     i=$(echo $i | sed 's/\([0-9.]*\)[eE]\([-+]*[0-9]*\)/(\1*10^\2)/g' )

     total=$(echo $total+$i | bc -l )

     ((count++))
   done
echo "scale=5; $total / $count" | bc -l
