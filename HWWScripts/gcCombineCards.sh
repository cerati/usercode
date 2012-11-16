#!/bin/bash

#this is cut, shape, combine, ...
MODE=$1

echo 'mode:'$MODE
#for MASS in {110,115,120,125,130,135,140,145,150,160,170,180,190,200,250,300,350,400,450,500,550,600}
for MASS in {110,125,160,250,300,600}
  do
  cd cards/${MASS}/
  if [ "${MODE}" == "cut" ] ; then  
      CARDS="hwwof_0j_cut.txt hwwsf_0j_cut.txt hwwof_1j_cut.txt hwwsf_1j_cut.txt hwwof_2j_cut.txt hwwsf_2j_cut.txt"
      COMBC="hww_012j_cut.txt"
  fi
  if [ "${MODE}" == "shape" ] ; then  
      CARDS="hwwof_0j_shape.txt hwwof_1j_shape.txt"
      COMBC="hwwof_01j_shape.txt"
  fi
  if [ "${MODE}" == "combine" ] ; then  
      CARDS="hwwof_0j_shape.txt hwwsf_0j_cut.txt hwwof_1j_shape.txt hwwsf_1j_cut.txt hwwof_2j_cut.txt hwwsf_2j_cut.txt"
      COMBC="hww_012j_combine.txt"
  fi
  if [ "${MODE}" == "allcomb" ] ; then  
      CARDS="hwwof_0j_shape.txt hwwsf_0j_cut.txt hwwof_1j_shape.txt hwwsf_1j_cut.txt hwwof_2j_cut.txt hwwsf_2j_cut.txt hwwof_0j_shape_7TeV.txt hwwsf_0j_shape_7TeV.txt hwwof_1j_shape_7TeV.txt hwwsf_1j_shape_7TeV.txt hww_2j_cut_7TeV.txt"
      COMBC="hww_012j_allcomb_7p8TeV.txt"
  fi
  combineCards.py -S ${CARDS} > ${COMBC}
  cd -
done

#./gcCombineCards.sh combine
#for mode in cut shape combine allcomb; do ./gcCombineCards.sh $mode; done
