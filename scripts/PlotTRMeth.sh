#!/usr/bin/env bash
bams=/project/mchaisso_100/mchaisso/projects/CARD/tr_methylation/HBCC/file_list.txt
inFile=$1
cwd=$PWD
while IFS= read -r line; do 
    region=`echo $line | awk '{ print $1":"$2"-"$3;}'`
    dirname=`echo $line | awk '{ print $1"_"$2"_"$3;}'`
    mkdir -p $dirname
    echo $line | awk '{ print $1"\t"$2"\t"$3;}' > $dirname/region.bed
    echo $region $dirname
    echo $region > $dirname/region.txt
    cd $cwd/$dirname && for bamfile in `tail -n +2 $bams | cut -f 3`; do bf=`echo $bamfile | tr "/" "\t" | awk '{ print $NF;}'`; echo "s=`echo $bf | tr "." "\t" | awk '{ print $1;}'`;  samtools view -b $bamfile $region -o bamfile.bam; samtools index bamfile.bam ; /project/mchaisso_100/mchaisso/projects/CARD/countme/countme/mecons bamfile.bam region.bed $s; ";  done
    
    cd $cwd/$dirname && for bamfile in `tail -n +2 $bams | cut -f 3`; do bf=`echo $bamfile | tr "/" "\t" | awk '{ print $NF;}'`; s=`echo $bf | tr "." "\t" | awk '{ print $1;}'`;  samtools view -b $bamfile $region -o bamfile.bam; samtools index bamfile.bam ; /project/mchaisso_100/mchaisso/projects/CARD/countme/countme/mecons bamfile.bam region.bed $s;  done > region.cons.txt
    cd $cwd/$dirname && Rscript ../../scripts/plot_locus_meth.R
    cd $cwd/$dirname && rm bamfile *.bai
    cd $cwd
done < "$inFile"
