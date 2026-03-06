#!/usr/bin/bash

#======================================
# https://github.com/yexianingyue
#======================================


##############################################
#       焊死的参数，不然每次手动指定太麻烦了
##############################################
env_ckv=/share/data1/zhangy2/conda/envs/ckv.1.0.3
db_ckv=/share/data1/Database/checkv/v1.5/checkv-db-v1.5

env_vib="/share/data1/software/miniconda3/envs/vibrant"

env_dvf="/share/data1/software/miniconda3/envs/dvf"
bin_dvf="/share/data1/lvqb/software/DeepVirFinder/dvf.py"

env_genomad="/share/data2/guorc/Software/conda/genomad.v1.11.1"
db_genomad="/share/data1/Database/geNomad/genomad_db"

db_busco="/share/data1/lish/virome/2020_kideny/27.scaffolding/busco/bacteria_odb10/tot.hmm"
cutoff_busco="/share/data1/lish/virome/2020_kideny/27.scaffolding/busco/bacteria_odb10/scores_cutoff"

bin_hmm=/usr/local/bin/hmmsearch

bin_o1=/path/to/print_id_len.py
bin_o2=prodigal-gv-parallel.py
export PATH=$PATH:$bin_o1:$bin_o2



# 初始化默认值
minlen=1000 # 最小序列
maxlen=0 # 最长序列
threads=80
force=0 # 是否强制重跑
clean=0 # 是删除临时文件
skip_deepvirfinder=0 # 是否跳过deepvirfinder
skip_virbrant=0 ## 是否跳过virbrant

base_path=$(realpath -s $0)
base_path=${base_path%/*}

# 定义颜色
GREEN='\033[0;32m'  # 绿色
LIGHT_GRAY='\033[1;30m'   # 浅灰色
YELLOW='\033[0;33m'  # 黄色
RED='\033[0;31m'  # 红色
NC='\033[0m'        # 无颜色

help=$(cat << EOF

    Usage:  ${GREEN}$0 ${YELLOW}<input.fasta> <output_prefix> [Parameters] ${NC}

    Parameters:

    -m|--minlen         Min size.[${GREEN}${minlen}${NC}]
    -M|--maxlen         Max size.[${GREEN}${maxlen}${NC}]
    -p|--threads        Number of threads.[${GREEN}${threads}${NC}]

    -skip_dvf|--skip_deepvirfinder
                        skip deepvirfinder [Flag]
    -skip_vib|--skip_virbrant
                        skip virbrabt [Flag]

    -bp|--base_path     dependes dir.(includes, flow.virus.{ckv.dvf,genomad,vib,busco}.sh)
                        [${GREEN}${base_path}${NC}]
    -f|--force          overwrite [Flag]
    -c|--clean          remove tmp directory [Flag]. (includes: <output_dir>/tmp)
    \n
EOF
)

if [[ $# -lt 2 ]];then
    echo -e "$help"
    exit 1
fi

# 初始化位置参数
positional=()

# 解析参数
while [[ $# -gt 0 ]]; do
    case $1 in
        -m|--minlen)
            minlen=$2
            shift 2
            ;;
        -M|--maxlen)
            maxlen=$2
            shift 2
            ;;
        -p|--threads)
            threads=$2
            shift 2
            ;;
        -skip_dvf|--skip_deepvirfinder)
            skip_deepvirfinder=1
            shift 1
            ;;
        -skip_vib|--skip_virbrant)
            skip_virbrant=1
            shift 1
            ;;
        -f|--force)
            force=1
            shift 1
            ;;
        -c|--clean)
            clean=1
            shift 1
            ;;
        --) # 终止选项参数解析，后续全部为位置参数
            shift
            positional+=("$@")
            break
            ;;
        -*|--*) # 处理未知选项
            echo "Unknown option $1"
            exit 1
            ;;
        *) # 处理位置参数
            positional+=("$1")
            shift
            ;;
    esac
done

set -- "${positional[@]}"

set -euo pipefail

input=$(realpath -s $1)
output=$(realpath -s $2)

cleanup() {
    rm -f "${output}.lock"
    exit 1
}

# 捕获终止信号
trap cleanup SIGINT SIGTERM

export PATH=$PATH:${base_path}

exec > >(tee ${output}.run.log) 2> >(tee ${output}.run.err >&2)

if [ $force -eq 1 ];then
    [ -f ${output} ] && rm -rf ${output}
    [ -d ${output}.tmp ] && rm -rf ${output}.tmp
    [ -f ${output}.lock ] && rm -rf ${output}.lock
    [ -f ${output}.ok ] && rm -rf ${output}.ok
else
    [ -f ${output}.lock ] && echo -e "${RED}running${NC}: ${output}" && exit 127
    [ -f ${output}.ok ] && echo -e "${GREEN}success${NC}: ${output}" && exit 0
fi

# 创建临时目录
[ ! -d ${output}.tmp ] && mkdir -p ${output}.tmp
touch ${output}.lock

## 过滤长度
# 根据长度进行过滤
if [ $minlen -eq 0 ];then
        ln -s ${input} ${output}.tmp/input_raw
else
        seqkit seq -m $minlen -w 0 -o  ${output}.tmp/input_raw ${input}
fi


####################################################
#               find virus
## 跑ckv on raw
flow.virus.ckv.sh ${output}.tmp/input_raw ${output}.tmp/st1.ckv --threads ${threads} --conda_env $env_ckv --db_path $db_ckv

## 跑ckv on provieus
flow.virus.ckv.sh ${output}.tmp/st1.ckv/proviruses.fna ${output}.tmp/st2.ckv --threads  ${threads} --conda_env $env_ckv --db_path $db_ckv

## merge ckv results
cat ${output}.tmp/st1.ckv/virus.select.fna  ${output}.tmp/st2.ckv/virus.select.fna > ${output}.tmp/ckv.virus.fna
## 因为seqkit head会强行终端输出，所以无论怎么样，都会exit 141，因此这边强行跳过就可以
nseq=$(seqkit seq -m $minlen ${output}.tmp/ckv.virus.fna | seqkit head -n 10 | seqkit seq -n | wc -l) || true

if [[ $nseq -eq 0 ]];then
    echo -e "[WARN] can't find any virus from:\n\n${input}\n\n"
    rm ${output}.lock
    touch ${output}.warning ${output}.ok
    exit 0
fi

## vibrant
if [[ $skip_virbrant -ne 1 ]];then
    flow.virus.vib.sh ${output}.tmp/ckv.virus.fna  ${output}.tmp/st3.vib \
        --threads ${threads} --conda_env ${env_vib}
else
    echo "skip VIRBRANT info"
fi

## deepvirfinder
if [[ $skip_deepvirfinder -ne 1 ]];then
    echo "extract DeepVirFinder inof"
    flow.virus.dvf.sh ${output}.tmp/ckv.virus.fna ${output}.tmp/st4.dvf --threads 1 --conda_env ${env_dvf} --bin_path ${bin_dvf}
else
    echo "skip DeepVirFinder info"
fi

####################################################
#               remove contaminations
## genomad
flow.virus.genomad.sh ${output}.tmp/ckv.virus.fna ${output}.tmp/st5.genomad \
    --threads ${threads}  --conda_env ${env_genomad} --db_path ${db_genomad}

## busco
prodigal-gv-parallel.py -t ${threads} -m -i ${output}.tmp/ckv.virus.fna -o /dev/null -q -a ${output}.tmp/ckv.virus.faa
flow.virus.busco.sh ${output}.tmp/ckv.virus.faa ${output}.tmp/st6.busco \
    --db_path ${db_busco} -cutoff ${cutoff_busco} --bin_path ${bin_hmm}
rm  ${output}.tmp/ckv.virus.faa

## intergate
echo "remove plasmid or virus with a high proportion of bacterial genes"
cat $(find ${output}.tmp -name "need_rm") |  seqkit grep -v -f - -o ${output}.fna.gz ${output}.tmp/ckv.virus.fna

echo "predict genes use prodigal-gv"
prodigal-gv-parallel.py -t ${threads} -m -q -i ${output}.fna.gz -o ${output}.gff -d ${output}.ffn -a ${output}.faa
print_id_len.py  ${output}.fna.gz > ${output}.fna.len

echo "match checkv quality"
head -n 1 ${output}.tmp/st1.ckv/quality_summary.tsv > ${output}.ckv
cat $(find ${output}.tmp -name "virus.select.ckv") | grep -v miuvig_quality >> ${output}.ckv
awk -F "\t" 'ARGIND==1{a[$1]=1;next};ARGIND==2{if($1 in a)print}' ${output}.fna.len ${output}.ckv > ${output}.ckv.tmp && mv ${output}.ckv.tmp ${output}.ckv
grep "^>" ${output}.faa | perl -ne 'chomp; $_=~s/^>//; @l=split/ # /; $l[0]=~/(\S+)_\d+/; print "$1\t$l[0]\t$l[1]\t$l[2]\n"' > ${output}.faa.feature

echo "match VIBRANT inof"
if [[ $skip_virbrant -ne 1 ]];then
    cat $(find ${output}.tmp -name "VIBRANT_genome_quality_*.tsv") | perl -ne 'chomp; @l=split/\t/; @n=split(/\s+/, $l[0]); print "$n[0]\t$l[1]\n"' >> ${output}.vib.lifestyle
    awk -F "\t" 'ARGIND==1{a[$1]=1;next};ARGIND==2{if($1 in a)print}' ${output}.fna.len ${output}.vib.lifestyle > ${output}.vib.lifestyle.tmp && mv ${output}.vib.lifestyle.tmp ${output}.vib.lifestyle
fi

if [[ $skip_deepvirfinder -ne 1 ]];then
    cat $(find ${output}.tmp -name "*bp_dvfpred.txt")  | perl -ne 'chomp; @l=split/\t/; @n=split(/\s+/, $l[0]); print "$n[0]\t$l[2]\t$l[3]\t\n"' > ${output}.dvf.info
    awk -F "[[:space:]]" 'ARGIND==1{a[$1]=1;next};ARGIND==2{if($1 in a)print}' ${output}.fna.len ${output}.dvf.info > ${output}.dvf.info.tmp && mv ${output}.dvf.info.tmp ${output}.dvf.info
fi

rm ${output}.run.err
rm ${output}.lock ${output}.fna.len  ${output}.gff  ${output}.tmp/input_raw
rm ${output}.tmp/ckv.virus.fna
touch ${output}.ok
pigz -f ${output}.{faa,ffn}
echo "Finished!"
exit 0
