#!/bin/bash
#cd ./bin && bsub -I -q q_x86_cn_cess -node 60842-60843,60848,60860,60862,60909,60922,60944,60947 -N $1 -np $2 ./nh_gomo 2>&1 |tee ../nonhydro_$1_$2.log  
bsub -I -q $1 -N $2 -np $3 ./nh_gomo 2>&1 |tee ./log.txt
#cd ./bin && ./nh_gomo_v2  2>&1 |tee ../tiny.log
