#!/user/bin/env bash


for line in `cat ./cnv_vcfs.txt`; do
    remove_pre=${line##*/}
    #echo ${remove_pre}
    remove_suf=${remove_pre%%_*}
    #echo ${remove_suf}
    dir=${line%%/f*}
    echo ${dir}
    new=${dir}"/"${remove_suf}"_cnv.vcf"  
    #echo ${new}
    mv ${line} ${new} 
done 
