#!/bin/bash

for MODE in CUT SHAPE COMB ALLC
do
  for MASS in {110,115,120,125,130,135,140,150,160,170,180,190,200,250,300,350,400,450,500,550,600}
    do
    cat logs/${MASS}/log_${MODE}_125_${MASS}_* >& logs/${MASS}/log_${MODE}_125_${MASS}.log
  done
done
