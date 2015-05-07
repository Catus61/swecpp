#!/bin/bash 
for i in {0..8000}; do 
if [ $[$i%500] -eq 0 ]; then 
python plotSWE.py 16 $i 
fi
done 

