#!/usr/bin/bash

#======================================
# https://github.com/yexianingyue
#======================================


set -e
shopt -s expand_aliases

if [ $# -lt 3 ] || [ $# -gt 4 ];then
    echo -e "\n\tUsage:\t$0 <out_prefix> <megahit_parameters> <fq1>[,fq1.2] [<fq2>[,fq2.2]]"
    echo -e "\n\t\t\t<megahit_parameters> 必须使用一对\e[1;31m双/单引号\e[0m引起来. 示例如下"
    echo -e "\n\teg.1. for paired end seq"
    echo -e "\t\t$0 output_dir/output_file \e[1;31m'\e[0m--k-list 21,41,61,81,101,121,141\e[1;31m'\e[0m input.1.fq.gz input.2.fq.gz"
    echo -e "\n\teg.2. for single seq"
    echo -e "\t\t$0 output_dir/output_file \e[1;31m\"\e[0m--k-list 21,41,61,81,101,121,141\e[1;31m\"\e[0m input.1.fq.gz"
    echo -e "\n\teg.3. for other megahit parameters"
    echo -e "\t\t$0 output_dir/output_file \e[1;31m\"\e[0m--k-list 21,41,61,81,101,121,141 -t 20\e[1;31m\"\e[0m input.1.fq.gz input.2.fq.gz"
    echo -e "\n\n\t如果使用parallel直接生成命令时，注意将\e[1;31m双/单引号\e[0m使用\e[1;31m反斜杠\e[0m转义。否则出错\n"
    exit 127
fi


#--------------------------------------------
# version：2
# 组装完成后，会重命名contigs
#



## software
alias megahit='/usr/local/bin/megahit'
alias gzip='/usr/local/bin/pigz'

## defaulte parameters
mode="single"

out=`realpath -s $1`
ps=$2
fq1=$3

if [ -f ${out}.tmplock ];then
    echo "${out} is running."
    echo "you can delete ${out}.tmplock to run it again."
    exit 127
elif [ -f ${out}.ok ];then
    echo "${out} success."
    exit 0
fi

touch ${out}.tmplock

if [ $# -eq 4 ];then
    mode="paired"
    fq2=$4
fi

cleanup() {
    rm -f "${out}.tmplock"
    exit 1
}

# 捕获终止信号
trap cleanup SIGINT SIGTERM EXIT


outp=${out##*/}
## 单双端测序数据公用的命令
comm_cmd=" ${ps} -o ${out} --out-prefix ${outp} 2> ${out}.run.err > ${out}.run.log \
    && mv ${out}/${outp}.contigs.fa ${out}.contigs.fa && gzip ${out}.contigs.fa \
    && mv ${out}/${outp}.log ${out}.megahit.log \
    && mv ${out}/options.json ${out}.megahit.opt.json \
    && seqkit replace -p \"^\" -r \"${outp}_ctg_{nr} \" -w 0 -o ${out}.contigs.fa.rename.gz ${out}.contigs.fa.gz && mv ${out}.contigs.fa.rename.gz ${out}.contigs.fa.gz \
    && chmod 444 ${out}.megahit.opt.json ${out}.megahit.log ${out}.contigs.fa.gz \
    && rm -rf ${out} \
    && mv ${out}.tmplock ${out}.ok \
    && mv ${out}.run.err ${out}.run.log \
    && chmod 444 ${out}.run.log"

if [ $mode == "single" ];then
    cmd="megahit -r $fq1 ${comm_cmd}"
    echo $cmd | sh
else
    # megahit -1 ${fq1} -2 ${fq2} ${ps} -o ${out} --out-prefix ${outp} 2> ${out}.run.err > ${out}.run.log
    cmd="megahit -1 ${fq1} -2 ${fq2} ${comm_cmd}"
    echo $cmd | sh
fi

