#!/usr/bin/bash

#======================================
# https://github.com/yexianingyue
#======================================



set -e
# 定义颜色
GREEN='\033[0;32m'  # 绿色
ORANGE='\033[0;33m' # 橙色
LIGHT_GRAY='\033[1;30m'   # 浅灰色
NC='\033[0m'        # 无颜色

#################################
#        some index
db_human="/share/data1/Database/human_genome/chm13v2.0.fa"
db_other1="/share/data1/Database/phi174/phiX174.fasta"    # phi174有时候是测序添加的标准品，一般需要去除

if [ $# -lt 5 ];then
    echo -e "\n\t${LIGHT_GRAY}bbduk.sh k=23, ktrim=r, mink=11, hdist=1, tpe, tbo, ftm=5, qtrim=rl, trimq=20, and minlen=100"
    echo -e "\thttps://www.nature.com/articles/s41597-023-02622-0#Sec2"
    echo -e "\n\tbbduk.sh ktrim=r k=23 mink=11 hdist=1 tpe tbo ftm=5 qtrim=rl trimq=10 minlen=100"
    echo -e "\thttps://www.nature.com/articles/s41467-023-40726-8#Sec13 # NC"
    
    echo -e "\n\tPlease make sure to run this script using \033[31mbash instead of sh..\033[0m"
    echo -e "\n\t${LIGHT_GRAY}NOTE: add bbduk to remove low complexity sequences.\n\t=============================$NC"
    echo -e "\n\tUsage:\t${GREEN}$0${NC} ${ORANGE}<index> <out_prefix> <min_len> <nreads:0> <fq1> [<fq2>]${NC}"
    # echo -e "\n\tUsage:\t$0 <index> <out_prefix> <min_len> <nreads:0> <fq1> [<fq2>]"
    echo -e "\n\t${ORANGE}min_len$NC: remove length less than <min_len>; 150 -> 90; 100bp -> 60"
    echo -e "\n\t${ORANGE}nreads:$NC how many reads/pairs to be processed. 0 means process all reads."
    echo -e "\n\t${ORANGE}index:$NC\n\t\thuman: $db_human"
    echo ""
    exit 0
fi

#################################
#          default Configure
threads=8
mode="Single" # 默认是双端测序

index=`realpath -s $1`
out=`realpath -s $2`
len=$3
nreads=$4

fq1=`realpath -s $5`
( [ $# -eq 6 ] ) && { fq2=`realpath -s $6`; mode="Paired"; } # 如果有两个序列，改成双端


exec > $out.run.err
exec 2>&1

#################################
#        softwares
shopt -s expand_aliases
alias bowtie2="/usr/local/bin/bowtie2"
alias fastp="/share/data3/shared_bin/fastp"
alias fastp="/usr/local/bin/fastp"
alias bbduk.sh="/share/data1/software/bbmap/bbduk.sh"
alias pigz="/usr/local/bin/pigz" # /usr/bin/gzip

if [ -f ${out}.ok ];then
    echo -e "Warning:\t${out}.ok" && exit 0
fi

if [ -f $out.running ];then
    echo -e "ERROR:\tIt's running. exit.\nyou can remove '$out.running' to run." && exit 127
fi

clean_tmp(){
    echo "Accidental termination"
    rm $out.running;
    exit 127
}

trap 'clean_tmp' SIGINT SIGTERM

touch $out.running

if [ $mode == "Single" ];then
    if [ -f $out.clean_noost.1.fq.gz ];then
        echo "The file already exists, checking integrity. $out.clean_nohost.fq.gz"
        { pigz -t $out.clean_nohost.fq.gz \
            && [[ `stat --format="%s" $out.clean_nohost.fq.gz` -ne 20 ]] \
            && echo "Done" >&2; } \
            && exit 0 \
            || echo "file is incomplete. rerun"
    elif [ -f $out.clean_noost.1.fq.gz ] && [ -f $out.clean_nohost.2.fq.gz ];then
        echo "The file already exists, checking integrity. $out.clean_nohost.1.fq.gz $out.clean_nohost.2.fq.gz"
        { pigz -t $out.clean_nohost.1.fq.gz $out.clean_nohost.2.fq.gz && [[ `stat --format="%s" $out.clean_nohost.1.fq.gz` -ne 20 ]] \
          && [[ `stat --format="%s" $out.clean_nohost.2.fq.gz` -ne 20 ]] && echo "Done" >&2; } \
            || echo "file is incomplete. rerun"
    fi
fi
echo -e "\nMode:\t${mode}\n"
echo -e "nReads:\t${nreads} -> [0 means process all reads.]\n"
echo -e "input:\t${fq1}\n\t${fq2}\n"
echo -e "outputdir:\t${out%/*}\n"

echo "Remove low quality using fastp"
# if [ ! -f $out.fastp.log ] || [[ `grep ', time used: ' $out.fastp.log| wc -l` -ne 1 ]];then
if [ ! -f $out.fastp.ok ];then
    if [ $mode == "Single" ];then
        fastp --reads_to_process ${nreads} -w ${threads} -q 20 -u 30 -n 5 -y -Y 30  --trim_poly_g --trim_poly_x -j /dev/null -h /dev/null \
            -l $len \
            -o $out.fp.fq.gz \
            -i $fq1 2> $out.fastp.log \
        || ! echo "ERROR" || exit 127; 
    else
        fastp  --reads_to_process ${nreads} -w ${threads} -q 20 -u 30 -n 5 -y -Y 30  --trim_poly_g --trim_poly_x -j /dev/null -h /dev/null \
            -l $len \
            -i $fq1 -I $fq2 \
            -o $out.1.fp.fq.gz -O $out.2.fp.fq.gz 2> $out.fastp.log \
        || ! echo "ERROR" || exit 127; 
    fi
fi
touch $out.fastp.ok


echo "Remove low complexity using bbduk.sh"

# if [ ! -f $out.bbduk.log ] || [[ `grep 'Bases Processed:' $out.bbduk.log | wc -l ` -ne 1 ]];then
if [ ! -f $out.bbduk.ok ] && [ -f $out.fastp.ok ];then
    if [ $mode == "Single" ];then
        bbduk.sh entropy=0.6 entropywindow=50 entropyk=5 \
            in=$out.fp.fq.gz \
            out=$out.bbduk.fq.gz 2> $out.bbduk.log \
        && rm $out.fp.fq.gz \
        || ! echo "ERROR" || exit 127; 
    else
        bbduk.sh entropy=0.6 entropywindow=50 entropyk=5 \
            in=$out.1.fp.fq.gz in2=$out.2.fp.fq.gz \
            out=$out.1.bbduk.fq.gz out2=$out.2.bbduk.fq.gz 2> $out.bbduk.log \
        && rm $out.1.fp.fq.gz $out.2.fp.fq.gz \
        || ! echo "ERROR" || exit 127; 
    fi
fi
touch $out.bbduk.ok


echo "Remove phi174 reads using bowtie2"

if [ ! -f $out.no_phi174.ok ] && [ -f $out.bbduk.ok ];then
    if [ $mode == "Single" ];then
        bowtie2 --end-to-end --mm --fast -U $out.bbduk.fq.gz  -x $db_other1  --no-head -S /dev/null --un-gz ${out}.no_phi174.fq.gz -p ${threads} 2> $out.no_phi174.log \
        && rm $out.bbduk.fq.gz\
        || ! echo "ERROR" || exit 127

        echo -e "\nOutput:\t$out.no_phi174.fq.gz"
    else
        bowtie2 --end-to-end --mm --fast -1 $out.1.bbduk.fq.gz -2 $out.2.bbduk.fq.gz  -x $db_other1  --no-head  -p ${threads} 2> $out.no_phi174.log \
            | perl -ne 'chomp;@s=split /\s+/;if($s[1]==77){print "\@$s[0]/1\n$s[9]\n+\n$s[10]\n";}elsif($s[1]==141){print STDERR "\@$s[0]/2\n$s[9]\n+\n$s[10]\n";}' \
            > >(pigz > $out.no_phi174.1.fq.gz) 2> >(pigz > $out.no_phi174.2.fq.gz) \
        && rm $out.1.bbduk.fq.gz $out.2.bbduk.fq.gz -f || echo "ERROR" || exit 127; 
    fi
fi
touch $out.no_phi174.ok



echo "Remove host reads using bowtie2"
# if [ ! -f $out.nohost.log ] || [[ `grep '% overall alignment rate' $out.nohost.log | wc -l` -ne 1 ]];then
if [ ! -f $out.nohost.ok ] && [ -f $out.no_phi174.ok ];then
    if [ $mode == "Single" ];then
        bowtie2 --end-to-end --mm --fast -U $out.no_phi174.fq.gz  -x $index  --no-head -S /dev/null --un-gz ${out}.fq.gz -p ${threads} 2> $out.nohost.log \
        && rm $out.no_phi174.fq.gz\
        || ! echo "ERROR" || exit 127

        echo -e "\nOutput:\t$out.fq.gz"
    else
        bowtie2 --end-to-end --mm --fast -1 $out.no_phi174.1.fq.gz -2 $out.no_phi174.2.fq.gz  -x $index  --no-head  -p ${threads} 2> $out.nohost.log \
            | perl -ne 'chomp;@s=split /\s+/;if($s[1]==77){print "\@$s[0]/1\n$s[9]\n+\n$s[10]\n";}elsif($s[1]==141){print STDERR "\@$s[0]/2\n$s[9]\n+\n$s[10]\n";}' \
            > >(pigz > $out.clean_nohost.1.fq.gz) 2> >(pigz > $out.clean_nohost.2.fq.gz) \
        && rm $out.no_phi174.1.fq.gz $out.no_phi174.2.fq.gz -f || echo "ERROR" || exit 127; 

        echo -e "\nOutput:\t$out.clean_nohost.1.fq.gz\n\t$out.clean_nohost.2.fq.gz"
    fi
fi
touch $out.nohost.ok

rm $out.nohost.ok $out.fastp.ok $out.bbduk.ok $out.no_phi174.ok

mv $out.running $out.ok
mv $out.run.err $out.run.log
