#!/usr/bin/bash

############################################
#
#           DeepVirFinder
#
############################################


# 初始化默认值
threads=8
db_path="/share/data1/lish/virome/2020_kideny/27.scaffolding/busco/bacteria_odb10/tot.hmm"
cut_off="/share/data1/lish/virome/2020_kideny/27.scaffolding/busco/bacteria_odb10/scores_cutoff"
bin_path=$(which hmmsearch)
force=0 # 是否强制重跑
clean=0 # 是删除临时文件

# 定义颜色
GREEN='\033[0;32m'  # 绿色
LIGHT_GRAY='\033[1;30m'   # 浅灰色
YELLOW='\033[0;33m'  # 黄色
RED='\033[0;31m'  # 红色
NC='\033[0m'        # 无颜色

help=$(cat << EOF

    Usage:  ${GREEN}$0 ${YELLOW}<input.faa> <output_dir> [Parameters] ${NC}


    Parameters:

    -p|--threads        Number of threads.[${GREEN}${threads}${NC}]
    -env|--db_path    /path/to/conda_env.[default: ${GREEN}${db_path}${NC}]
    -cutoff|--cut_off    /path/to/conda_env.[default: ${GREEN}${cut_off}${NC}]
    -bin|--bin_path    /path/to/dvf.py [default: ${GREEN}${bin_path}${NC}]
    -f|--force          overwrite [Flag]
    ${LIGHT_GRAY}-c|--clean          remove tmp directory [Flag]. (discard)${NC}
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
        -db|--db_path)
            conda_env=$2
            shift 2
            ;;
        -cutoff|--cut_off)
            conda_env=$2
            shift 2
            ;;
        -bin|--bin_path)
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
    exit 1
}

# 捕获终止信号
trap cleanup SIGINT SIGTERM


## 会将标准输出和错误输出分别定向到文件和终端
exec > >(tee ${output}.run.log) 2> >(tee ${output}.run.err >&2)

if [ $force -eq 1 ];then
    [ -d ${output} ] && rm -rf ${output}
    [ -f ${output}.lock ] && rm -rf ${output}.lock
    [ -f ${output}.ok ] && rm -rf ${output}.ok
else
    [ ! -f ${output}/busco.ok ] && { [ -d ${output} ] && rm -rf ${output}; }
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

set -euo pipefail
touch ${output}.lock
mkdir ${output}

reader="cat"
([[ $input == *.gz ]]) && reader="pigz -dc"

## 分割文件
[ ! -f ${output}/split.ok ] && { echo -e "\tsplit input" && $reader ${input} | seqkit split2 -p $threads -O ${output}/tmp && touch ${output}/split.ok || ! echo ERROR || exit 127; }


## 比对
if [ ! -f ${output}/tmp.ok ];then

    echo "Run hmmsearch"

    > ${output}/tmp.sh
    for i in `find ${output}/tmp -name "*.fasta"`
    do
        s1="${bin_path} --cpu 5 -o /dev/null --noali --tblout ${i}.busco ${db_path} ${i}"
        printf "[ ! -f $i.ok ] && { %s && touch $i.ok || exit 127; }\n" "$s1" >> ${output}/tmp.sh
    done

    parallel -j $threads < ${output}/tmp.sh;

    nsplit=`find ${output}/tmp -name "*.fasta" | wc -l`
    nsplit_ok=`find ${output}/tmp -name "*.ok" | wc -l`
    echo -e "\tmerged blast results"
    [ $nsplit -eq $nsplit_ok ] && { cat ${output}/tmp/*.busco > ${output}/busco && rm -rf ${output}/tmp ${output}/tmp.sh && touch ${output}/busco.ok || ! echo "error blastn" || exit 127; }
fi



## activate conda env
# 去除细菌基因比例大于0.05的序列

echo "filter busco by cutoff"
perl -e '%a;%b; open I,"$ARGV[0]";while(<I>){chomp;@l=split/\s+/;$a{$l[0]}=$l[1]}; close(I); open I,"$ARGV[1]"; open O, ">$ARGV[2]"; while(<I>){chomp; next unless $_!~/^#/;@l=split/\s+/; if ($l[5] >= $a{$l[2]}){ print O "$_\n"; $l[0]=~/(.*)_\d+$/; $b{$1}++; }};close(I);close(O); open I,"$ARGV[3]"; while(<I>){if(/^>(\S+)/){$n=$1; $n=~/(.*)_\d+$/; $c{$1}{"all"}++ }}; close(I); open O,">$ARGV[4]"; foreach $k(keys %c){$prok=$b{$k}//0; print O "$k\t$c{$k}{all}\t$prok\n"}; close(O)' ${cut_off}  ${output}/busco ${output}/busco.f  ${input} ${output}/xx.busco.summ

echo "Remove virus with a high proportion of bacterial genes"
awk '$3/$2>=0.05' ${output}/xx.busco.summ > ${output}/xx.busco.summ.f
cut -f 1 ${output}/xx.busco.summ.f > ${output}/need_rm
rm -rf ${output}.lock && touch ${output}.ok
exit 0
