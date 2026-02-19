#!/bin/bash


#FIELD=("00" "01" "02" "03" "04" "05" "06" "07" "08" "09" "11" "13" "15" "17" "19" "22" "25" "30" "35" "40" "45")
FIELD=("00" "01" "02" "03" "04" "05" "06" "07" "08" "09" "11" "13" "15" "17" "19" "22" "25" "30" "35" "40" "45" "50" "55" "60" "65" "70" "75" "80" "85" "90" "95" "100" "105" "110" "115")
len=${#FIELD[@]}

#name="LiPF6inACN"
#name="LiIinH2O"
#name="TLiPF6inACN"
name="LiPF6inH2O"

len=${#FIELD[@]}
#name="LiIinH2O"
#name="StockD"
#name="KPF6inACN"
density="05"
rm diel${name}D$density
for(( i=0; i<len; i++));
do
  field=${FIELD[i]}
  cat ${name}D${density}E${field}.dat >> diel${name}D$density
  mv ${name}D${density}E${field}.dat ./RAW
done
