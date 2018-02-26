#!/user/bin/env bash

for file in ./*.bed; do
    remove_pre=${file##*/}
    remove_suf=${remove_pre%%.bed*}
    #echo ${remove_suf}
    cut -f1,2,3 $file > ${remove_suf}.bed
done
