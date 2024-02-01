#!/bin/bash

rm -r plots/*

paraname=""
rwstart=0
rwstep=0
rwnstep=0

while read -a line
do
paraname=${line[0]}
rwstart=${line[1]}
rwstep=${line[2]}
rwnstep=${line[3]}

#if [ $paraname != "T0" ]
#then
#continue
#fi

if [ $paraname == "M0" ] ; then continue ; echo $paraname ; fi
if [ $paraname == "M1" ] ; then continue ; echo $paraname ; fi
if [ $paraname == "M2" ] ; then continue ; echo $paraname ; fi
if [ $paraname == "M3" ] ; then continue ; echo $paraname ; fi
if [ $paraname == "M4" ] ; then continue ; echo $paraname ; fi
if [ $paraname == "M5" ] ; then continue ; echo $paraname ; fi
if [ $paraname == "M6" ] ; then continue ; echo $paraname ; fi
if [ $paraname == "S0" ] ; then continue ; echo $paraname ; fi
if [ $paraname == "S1" ] ; then continue ; echo $paraname ; fi
if [ $paraname == "S2" ] ; then continue ; echo $paraname ; fi


root -l -q -b 'muoncolliderdelphes_selection_draw.C("'${paraname}'", '${rwstart}', '${rwstep}', '${rwnstep}')' > muoncolliderdelphes_selection_draw_${paraname}.log
done < ../QCKM_5_reweightscan_setups.list
