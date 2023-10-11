#!/bin/sh
timeTopk="/home/dataFiles/"
for i in '16' # '16' # '04'
do
    for j in  '05' # '05' '25' '5' '001' '0025' '005' '01' '025' # '05' # '01' '025' '05'
    do  
        for k in '1024' # '16' '32' '64' '128' '256' '512' '1024' '2048' '4096' '8192' # '16'
        do
            for l in 'emd' 'euclidean' 'ks' 'cd' # 'emd' 'euclidean'
            do
                stdbuf -oL ./test topk_query --metric ${l} --prefix /data/ --dataset pus0920_10000 --abs_error 0.${j} --topk ${i} --query_num 10 --initial_size ${k} > ${timeTopk}pus0920Top${i}absError0${j}initSize${k}${l}.txt     
                stdbuf -oL ./test topk_query --metric ${l} --prefix /data/ --dataset hus0920_10000 --abs_error 0.${j} --topk ${i} --query_num 10 --initial_size ${k} > ${timeTopk}hus0920Top${i}absError0${j}initSize${k}${l}.txt 
                stdbuf -oL ./test topk_query --metric ${l} --prefix /data/ --dataset enem1121_10000 --abs_error 0.${j} --topk ${i} --query_num 10 --initial_size ${k} > ${timeTopk}enem1121Top${i}absError0${j}initSize${k}${l}.txt 
                stdbuf -oL ./test topk_query --metric ${l} --prefix /data/ --dataset airline_10000 --abs_error 0.${j} --topk ${i} --query_num 10 --initial_size ${k} > ${timeTopk}airlineTop${i}absError0${j}initSize${k}${l}.txt 
            done
        done
    done
done