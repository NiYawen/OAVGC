#!/usr/bin/bash

set -e

if [ $# -lt 2 ] || [ $# -gt 3 ];then
    echo -e "\n\n\tUsage:\t$0 <\e[32minput.fasta[.gz]\e[0m> <\e[32moutput_file\e[0m> [ncpus:\033[32m30\033[0m]\n\n"
    exit 0
fi

inf=$1
outf=$2
threads=${3:-30}

tts=6

exec > >(tee $outf.run.log) 2> >(tee $outf.run.err >&2)


([ -f $outf.ok ]) && echo "$outf.ok is already appeared." && exit 0


## 检查文件是否存在
if [ ! -f $inf ];then echo -e "\e[1;31mcan't find file: \e[32m $inf\[0m" && exit 127;fi

reader="cat"
([[ $inf == *.gz ]]) && reader="pigz -dc"

## 建库
echo -e "step1/$tts:\tcheck index"
echo -e "input file: $inf\n"
[ ! -f $outf.index.ok ] && { echo -e "\tbuild index" && less $inf | makeblastdb -dbtype nucl -out $outf.index.btn -title virus && touch $outf.index.ok || ! echo "ERROR"  || exit 127; }

## 分割文件
echo -e "\nstep2/$tts\tcheck split status."
[ ! -f $outf.split.ok ] && { echo -e "\tsplit input" && $reader $inf | seqkit split2 -p $threads -O $outf.tmp && touch $outf.split.ok || ! echo ERROR || exit 127; }


## 比对
echo -e "\nstep3/$tts\tcheck blastn results"
if [ ! -f $outf.tmp.ok ];then
    echo -e "\trun blastn"

    > $outf.tmp.ckv.sh
    for i in `find $outf.tmp -name "*.fasta"`
    do
        s1="blastn -query $i -db $outf.index.btn -out $i.btn  -outfmt '6 std qlen slen' -num_alignments 999999 -num_threads 2"
        s2="/share/data1/zhangy2/scripts/blast_cvg.v1.py $i.btn $i.cvg 70 1"
        printf "[ ! -f $i.ok ] && { %s && %s && touch $i.ok || exit 127; } || { exit 0; }\n" "$s1" "$s2" >> $outf.tmp.ckv.sh
    done

    parallel -j $threads < $outf.tmp.ckv.sh;

    nsplit=`find $outf.tmp -name "*.fasta" | wc -l`
    nsplit_ok=`find $outf.tmp -name "*.ok" | wc -l`
    echo -e "\tmerged blast results"
    [ $nsplit -eq $nsplit_ok ] && { cat $outf.tmp/*.btn > $outf.btn && cat $outf.tmp/*.cvg > $outf.cvg && rm -rf $outf.tmp && touch $outf.tmp.ok || ! echo "error blastn" || exit 127; }
fi


## 聚类
ind=95
cov=85
echo -e "\nstep4/$tts\tget length"
[ ! -f $outf.len.ok ] && { /share/data1/zhangy2/scripts/print_id_len.py $inf | sort -rnk2 > $outf.sort.len && touch $outf.len.ok || exit 127; }
echo -e "\nstep5/$tts\tcluster"
[ ! -f $outf.cvg.ok ] && {
    /share/data1/zhangy2/scripts/blast_cluster.v2.pl  $outf.sort.len $outf.cvg $outf.i${ind}_c${cov}.uniq $cov $ind \
    && less -S ${outf}.i${ind}_c${cov}.uniq.list | perl -e 'while(<>){chomp;@s=split/\s+/;push @{$h{$s[0]}},$s[1];} for(keys %h){@a=@{$h{$_}}; print "$_\t".($#a+1)."\t".(join",", @a)."\n";}' > ${outf}.i${ind}_c${cov}.uniq.list.f \
    && touch $outf.cvg.ok \
    || exit 127; }
echo -e "\nstep6/$tts\tselect represent sequences"
[ ! -f $outf.ok ] && { cut -f1 $outf.i${ind}_c${cov}.uniq.list | uniq | seqkit grep -w 0 -f - $inf > $outf.fa && touch $outf.ok || exit 127; }

## 删除临时文件
rm $outf.tmp.ckv.sh $outf.index.btn* $outf.index.ok $outf.split.ok $outf.tmp.ok $outf.len.ok $outf.cvg.ok 
pigz $outf.cvg $outf.btn
