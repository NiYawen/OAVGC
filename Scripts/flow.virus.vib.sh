#!/usr/bin/bash

############################################
#
#           VIRBRANT
#
############################################

# 初始化默认值
threads=8
conda_env="/share/data1/software/miniconda3/envs/vibrant"
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
    -env|--conda_env    /path/to/conda_env.[default: ${GREEN}${conda_env}${NC}]
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
    rm -rf ${output}
    exit 1
}

# 捕获终止信号
trap cleanup SIGINT SIGTERM

exec > >(tee ${output}.run.log) 2> >(tee ${output}.run.err >&2)

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
echo "activate ckeckv envs ${conda_env}"
set +u && source activate ${conda_env} && set -u

echo "Run VIBRANT"
VIBRANT_run.py -i ${input} -virome -f nucl -folder ${output} -t ${threads} -no_plot

echo "extract virus  and info"
perl -ne 's/_frag.*$//;print $_;' $(find ${output} -name "*.phages_combined.txt") | cut -d " " -f 1 > ${output}/virus.name 
perl -ne 'chomp; @l=split/\t/; $l[0]=~s/(\S+).*/$1/; print "$l[0]\t$l[1]\n"' $(find . -name "VIBRANT_genome_quality*.tsv") > ${output}/virus.info

tmp_dir=$(find ${output} -maxdepth 1 -type d -name "VIBRANT_*")
rm -rf ${output}.lock ${output}.run.err && clean_dir $clean "${tmp_dir}" && touch ${output}.ok
exit 0