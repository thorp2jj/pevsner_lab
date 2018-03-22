#!/user/bin/env bash


for sample in ./*.cns; do
    remove_pre=${sample##*/}
    #echo ${remove_pre}
    remove_suf=${remove_pre%%_*}
    #echo ${remove_suf}
    mkdir ${remove_suf}
    sleep 2
    srun cnvkit.py call ${remove_suf}_recal.cns -y -m threshold -t=-1.1,-0.4,0.3,0.7 -o ./${remove_suf}/${remove_suf}_recal.call.cns
    wait
    srun cnvkit.py call ./${remove_suf}/${remove_suf}_recal.call.cns -y -m clonal --purity 1 -o ./${remove_suf}/${remove_suf}_recal_purity.call.cns
    wait
    srun cnvkit.py export bed ./${remove_suf}/${remove_suf}_recal_purity.call.cns --show all -y -o ./${remove_suf}/${remove_suf}_recal_purity_call.bed
    wait
    mv ${sample} ./${remove_suf}
done 
