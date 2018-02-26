#!/user/bin/env bash

for line in `cat ./cnv_beds.txt`; do
    remove_pre=${line##*/}
    remove_suf=${remove_pre%%.bed}
    #echo ${remove_suf}
    for text in `cat ./indel_beds.txt`; do
        reduce=${text##*/}
        reduce2=${reduce%%.bed}
        #echo ${reduce2}
        if [[ ${remove_suf} == ${reduce2} ]]; then
            #bcftools merge ${line} ${text} > /mnt/data/jeremy/projects/smri/mosaic_hunter/filter_files/merged/${remove_suf}_merged.vcf.gz
            #echo ${remove_suf} ${reduce2}
            cat ${line} ${text} > /mnt/data/jeremy/projects/smri/mosaic_hunter/filter_files/100X/merged/${remove_pre}_cnv_indel.bed
        else
            continue
        fi
    done
done
