#!/usr/bin/bash

#======================================
# https://github.com/yexianingyue
#======================================



# 初始化默认值
threads=8
db_crispr="/share/data1/Database/uhgg_20220401/csp.tot.fa"
db_genome="/share/data1/Database/uhgg_20220401/genomes_rep.fasta"
force=0
exec=0

# 定义颜色
GREEN='\033[0;32m'  # 绿色
LIGHT_GRAY='\033[1;30m'   # 浅灰色
YELLOW='\033[0;33m'  # 黄色
RED='\033[0;31m'  # 红色
NC='\033[0m'        # 无颜色

help=$(cat << EOF
    
    Usage:  ${GREEN}$0 ${YELLOW}<virus.fasta> <out_prefix> [Parameters] ${NC}


    Parameters:

    -p|--threads    Number of threads.[${GREEN}${threads}${NC}]
    -d1|--db_crispr /path/to/btn.index.csp.[default: ${GREEN}${db_crispr}${NC}]
    -d2|--db_genome /path/to/btn.index.genome.[default: ${GREEN}${db_genome}${NC}]
    -d3|--db_taxo /path/to/host.taxonomy [None]. 至少2列，第一列是sgb name， 第二列及以后是sgb信息,使用tab分割
    -e|--exec       execute. ${YELLOW} default print command.${NC}
    -f|--force      overwrite


    ${LIGHT_GRAY}recommended DB
        ## Gut
        --db_genome /share/data1/Database/uhgg_20220401/genomes_rep.fasta
        --db_crispr /share/data1/Database/uhgg_20220401/csp.tot.fa
        --db_taxo /share/data1/Database/uhgg_20220401/uhgg.taxo.group

        ## Oral
        --db_genome /share/data1/zhangy2/00.pub_data/12.oral_prok/v1.20241207/genome
        --db_crispr /share/data1/zhangy2/00.pub_data/12.oral_prok/v1.20241207/csp.tot.fa
        --db_taxo /share/data1/zhangy2/00.pub_data/12.oral_prok/v1.20241207/genomes.taxo.format
        ${NC}
        \n
EOF
)

set -e # 如果出错，就不再执行下一步
if [[ $# -lt 3 ]];then
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
        -d1|--db_crispr)
            db_crispr=$2
            shift 2
            ;;
        -d2|--db_genome)
            db_genome=$2
            shift 2
            ;;
        -d3|--db_taxo)
            db_taxo=$2
            shift 2
            ;;
        -f|--force)
            force=1
            shift 1
            ;;
        -e|--exec) ## 是执行
            exec=1
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

fa=$1
outf=$2

logger() {
    local level="$1"
    shift 1
    local message="$*"
    local timestamp
    

    timestamp="${LIGHT_GRAY}$(date '+%Y-%m-%d %H:%M:%S')\033[0m"

    local color
    case "$level" in
        INFO)
            color=$GREEN
            ;;
        WARN)
            color=$YELLOW
            ;;
        ERROR)
            color=$RED
            ;;
        *)
            color=$NC
            ;;
    esac

    message=$(echo "$message" | tr '\n' ' ' | sed 's/[ ][ ]*/ /g' | sed 's/\\t/\\\\t/g; s/\\n/\\\\n/g')
    # 输出到控制台
    echo -e "${timestamp} [${color}${level}\033[0m]\n${message}"
}

if [ $exec -eq 1 ];then
    if [ $force -eq 0 ] && [ -f ${outf}.ok ] ;then
        logger "WARN" "${outf}.ok is exists. If you wan't overwrite. please add params -f"
        exit 0
    fi
fi


run_cmd(){
    cmd=$1
    stats=$2

    ## 如果只是打印命令行，那就不判断别的
    if [ $exec -eq 0 ];then
        cmd=$(echo "$cmd" | tr '\n' ' ' | sed 's/[ ][ ]*/ /g' | sed 's/\\t/\\\\t/g; s/\\n/\\\\n/g')
        echo -e "$cmd\n"
        return 0
    fi

    logger "INFO" "$cmd"

    if [ $force -eq 0 ];then
        if [[ -f ${stats}.running ]];then
            logger "ERROR" "${stats}.running is exists. exit."
            exit 127
        elif [[ -f ${stats}.ok ]];then
            logger "WARN" "${stats}.ok; skip $stats."
            return 0
        fi
    fi

    eval $cmd 2> ${stats}.running \
        && mv ${stats}.running ${stats}.ok \
        || ! mv ${stats}.running ${stats}.err \
        || ! logger "ERROR" "details in ${stats}.err" \
        || exit 127
}

echo -e "# --- 1.By CRISPRs ---\n#===================="
# 1、比对CRISPRS
cmd="blastn -query $fa -db $db_crispr -evalue 1e-2 -out $outf.csp.bt -outfmt 6 -num_alignments 999999 -num_threads ${threads} -word_size 8"
run_cmd "$cmd" ${outf}.step1.1

cmd="filter_blast -i $outf.csp.bt -o $outf.csp.bt.f --evalue 1e-5 --score 45 --tops 20"
run_cmd "$cmd" ${outf}.step1.2

cmd=$(cat << EOF
less ${outf}.csp.bt.f | perl -e 'while(<>){ chomp;@s=split /\s+/;(\$a,\$b)=(\$s[0],\$s[1]);\$b=~s/_\d+.sp\d+\$//; push @{\$g{\$a}},\$b unless exists \$h{"\$a \$b"}; \$h{"\$a \$b"}++;} for(sort keys %g){@a=@{\$g{\$_}};print "\$_\t".(join ",",@a)."\n";}' 
> $outf.csp.list
EOF
)
run_cmd "$cmd" ${outf}.step1.3

echo -e "\n\n# --- 2.By Genome ---\n#===================="
# 2、基因组比对
cmd="blastn -query $fa -db $db_genome -evalue 1e-2 -out $outf.fa.bt -outfmt 6 -num_alignments 999999 -num_threads ${threads}"
run_cmd "$cmd" ${outf}.step2.1

cmd="connect_blast $outf.fa.bt $outf.fa.bt.conn 1"
run_cmd "$cmd" ${outf}.step2.2

cmd="get_length.pl $fa /dev/stdout |sort -rnk2 > $outf.sort.len"
run_cmd "$cmd" ${outf}.step2.3

cmd="filter_blast -i $outf.fa.bt.conn -o $outf.fa.bt.conn.f --evalue 1e-10 --qfile $outf.sort.len --qper 30 --identity 90"
run_cmd "$cmd" ${outf}.step2.4

cmd=$(cat << EOF
less $outf.fa.bt.conn.f | perl -e 'while(<>){chomp;@s=split /\s+/;(\$a,\$b)=(\$s[0],\$s[1]);\$b=~s/_\d+\$//; push @{\$g{\$a}},\$b unless exists \$h{"\$a \$b"}; \$h{"\$a \$b"}++;} for(sort keys %g){@a=@{\$g{\$_}};print "\$_\t".(join ",",@a)."\n";}' \
    > $outf.fa.list
EOF
)
run_cmd "$cmd" ${outf}.step2.5

# 3、合并两个结果
echo -e "\n\n# --- 3.Merge result ---\n#===================="
cmd=$(cat << EOF
    cat $outf.csp.list $outf.fa.list |sort | perl -ne 'chomp;@a = split /[,\t]/; foreach \$x (@a){print "\$a[0]\t\$x\n"}' |awk '\$1!=\$2'|sort -u > $outf
EOF
)
run_cmd "$cmd" ${outf}.step3.1

# 4、如果有分类文件，就顺便把分类也添加上
echo -e "\n\n$ --- 4. add taxonomy---\n#===================="
cmd=$(cat << EOF
    ([ -f "$db_taxo" ]) && { awk -F "\t" 'ARGIND==1{a[\$1]=\$0; next}; ARGIND==2{printf "%s\t%s\n", \$1, a[\$2]}' $db_taxo ${outf} > ${outf}.taxo; } || { echo "Can't find taxodb"; }
EOF
)
run_cmd "$cmd" ${outf}.step3.2

# 清除其它的标志文件
if [ $exec == 1 ];then
    touch ${outf}.ok
    logger "INFO" "rm ${outf}.step1.[123].ok ${outf}.step2.[12345].ok ${outf}.step3.[12].ok"
    rm ${outf}.step1.[123].ok ${outf}.step2.[12345].ok ${outf}.step3.[12].ok
fi



