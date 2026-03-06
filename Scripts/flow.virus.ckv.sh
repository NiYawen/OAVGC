#!/usr/bin/bash


# 初始化默认值
threads=8
conda_env="/share/data1/zhangy2/conda/envs/ckv.1.0.3"
db_path="/share/data1/Database/checkv/v1.5/checkv-db-v1.5"
force=0 # 是否强制重跑
clean=0 # 是删除临时文件

# 定义颜色
GREEN='\033[0;32m'  # 绿色
LIGHT_GRAY='\033[1;30m'   # 浅灰色
YELLOW='\033[0;33m'  # 黄色
RED='\033[0;31m'  # 红色
NC='\033[0m'        # 无颜色

help=$(cat << EOF

    Usage:  ${GREEN}$0 ${YELLOW}<input.fasta> <output_dir> [Parameters] ${NC}


    Parameters:

    -p|--threads        Number of threads.[${GREEN}${threads}${NC}]
    -env|--conda_env    /path/to/btn.index.csp.[default: ${GREEN}${conda_env}${NC}]
    -db|--db_path       /path/to/btn.index.genome.[default: ${GREEN}${db_path}${NC}]
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
        -p|--threads)
            threads=$2
            shift 2
            ;;
        -env|--conda_env)
            conda_env=$2
            shift 2
            ;;
        -db|--db_path)
            db_path=$2
            shift 2
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

input=$(realpath -s $1)
output=$(realpath -s $2)

cleanup() {
    rm -f "${output}.lock"
    exit 1
}

# 捕获终止信号
trap cleanup SIGINT SIGTERM

if [ $force -eq 1 ];then
    [ -d ${output} ] && rm -rf ${output}
    [ -f ${output}.lock ] && rm -rf ${output}.lock
    [ -f ${output}.ok ] && rm -rf ${output}.ok
else
    [ -f ${output}.lock ] && echo -e "${RED}running${NC}: ${output}" && exit 127
    [ -f ${output}.ok ] && echo -e "${GREEN}success${NC}: ${output}" && exit 0
fi

clean_dir(){
    clean=$1
    target_dir=$2
    if [ $clean -eq 1 ];then
        rm -rf ${target_dir}
    fi
}

touch ${output}.lock

set -euo pipefail
## activate conda env
echo "source activate ${conda_env}"
set +u
source activate ${conda_env}
set -u

echo "Run ckeckv on ${input}"
checkv end_to_end -d ${db_path} -t ${threads} ${input} ${output}

echo "filter ckv results"
## 之前的过滤条件是尽可能的找到更多的病毒，现在主要是为了找到高质量病毒
# grep -v "contig >1.5x longer than expected genome length" ${output}/quality_summary.tsv | awk '($3=="No"||$1=="contig_id") && !($7/($5+1e-16)>0.5 && ($6==0 || ($6>0 && $7/($6+1e-16)>5)))' > ${output}/virus.select.ckv
head -n 1 ${output}/quality_summary.tsv > ${output}/virus.select.ckv
grep -v "contig >1.5x longer than expected genome length" ${output}/quality_summary.tsv \
    | awk -F "\t" '$3=="No" && ($9 == "High-quality" || $8=="Medium-quality")' >> ${output}/virus.select.ckv

### 如果过滤完，ckv结果为0，就退出，后续不做了
nckv=$(cat ${output}/virus.select.ckv | wc -l)
if [ $nckv -eq 1 ];then
    echo "can't find any virus from ${input}"
    rm ${output}.lock && clean_dir $clean "$output/tmp" && touch ${output}.ok
    exit 0
fi

echo "match features and gene features"
awk -F "\t" 'ARGIND==1{a[$3]=$8; next}; ARGIND==2{if(FNR==1) {printf "ctg_name\tgene_name\tstart\tend\ttype\n";next}; printf "%s\t%s\t%d\t%d\t%s\n", $1, $1"_"$2, $3, $4, a[$9]}' ${db_path}/hmm_db/checkv_hmms.tsv ${output}/tmp/gene_features.tsv > ${output}/gene.feature.tsv

echo "extract virus from ${output}/virus.select.ckv"
cut -f 1 ${output}/virus.select.ckv | seqkit grep -w 0 -f - -o ${output}/virus.select.fna ${input}

rm ${output}.lock && clean_dir $clean "$output/tmp" && touch ${output}.ok
