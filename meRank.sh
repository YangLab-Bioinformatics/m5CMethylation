#!/bin/sh
#PBS -q tmp
#PBS -l mem=60gb,walltime=500:00:00,nodes=1:ppn=12
#HSCHED -s hschedd+maker+human

cd /pnas/yangyg_group/hanyn/2-patient-bs/Tissue_20180523_BS
core_num=10

mkdir name
cd name
mkdir 1-data
cd 1-data
gunzip -c /pnas/yangyg_group/sunbf/bladder_cancer_guangzhou/Tissue_20180523_BS/name_1.fq.gz > R1.fastq
gunzip -c /pnas/yangyg_group/sunbf/bladder_cancer_guangzhou/Tissue_20180523_BS/name_2.fq.gz > R2.fastq
cutadapt -a AGATCGGAAGAGCACACGTCTG -A AGATCGGAAGAGCGTCGTGTAG -o 1-trim_1.fastq -p 1-trim_2.fastq R1.fastq R2.fastq
java -Xmx4g -jar /software/biosoft/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 ./1-trim_1.fastq ./1-trim_2.fastq ./2-trim_qua_1.fastq ./unpaired_1.fastq ./2-trim_qua_2.fastq ./unpaired_2.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:15
mkdir ../2-fastqc
fastqc 1-trim_1.fastq -t $core_num -o ../2-fastqc
fastqc 1-trim_2.fastq -t $core_num -o ../2-fastqc
fastqc R1.fastq -t $core_num -o ../2-fastqc
fastqc R2.fastq -t $core_num -o ../2-fastqc
fastqc 2-trim_qua_1.fastq -t $core_num -o ../2-fastqc   
fastqc 2-trim_qua_2.fastq -t $core_num -o ../2-fastqc   
rm R1.fastq R2.fastq 1-trim_1.fastq 1-trim_2.fastq ../2-fastqc/*.zip unpaired_*.fastq

##DHFR#####################################################
mkdir ../5-dhfr_meRanTK
cd ../5-dhfr_meRanTK
meRanTK_index="/pnas/yangyg_group/zhangting/Reference/human/ENSEMBLE_72/meRanTK-1.2.0/Human_Dhrf/hisat2_index/"
sesFA="/pnas/yangyg_group/zhangting/Reference/human/ENSEMBLE_72/meRanTK-1.2.0/Human_Dhrf/Dhrf.fasta"
meRanGh align -o ./ -f ../1-data/2-trim_qua_2.fastq -r ../1-data/2-trim_qua_1.fastq  -t $core_num -S DHFR.sam -ds -id $meRanTK_index  -mbp -fmo
meRanCall -f $sesFA -p $core_num -bam DHFR_sorted.bam -gref -o DHFR.txt -mBQ 20 -mr 0
sed '1d'  DHFR.txt| awk '{i+=$5-$6; j+=$5}END{print i/j}' > transition.txt
rm *_unmapped.fastq 

##45S#####################################################
mkdir ../6-45S_meRanTK
cd  ../6-45S_meRanTK
meRanTK_index="/pnas/yangyg_group/zhangting/Reference/human/ENSEMBLE_72/meRanTK-1.2.0/Human_45S/hisat2_index/"
meRanGh align -o ./ -f ../1-data/2-trim_qua_2.fastq -r ../1-data/2-trim_qua_1.fastq  -t $core_num -S 45S.sam -ds -id $meRanTK_index  -mbp -fmo
rm *_unmapped.fastq

##principal#####################################################
htseq_index="/pnas/yangyg_group/zhangting/Reference/human/ENSEMBLE_72/index/genome.fa"
meRanTK_index="/pnas/yangyg_group/hanyn/data/index/human_merankYZ"
hg19_gtf="/pnas/yangyg_group/zhangting/Reference/human/ENSEMBLE_72/GTF/Homo_sapiens.GRCh37.72.chr.gtf"
mkdir ../3-meRanGhResult
cd ../3-meRanGhResult
export PATH=/pnas/yangyg_group/zhangqy/software/hisat2-2.0.4/:$PATH
meRanGh align -o ./ -f ../1-data/2-trim_qua_2.fastq -r ../1-data/2-trim_qua_1.fastq -t $core_num -S align.sam -ds -id $meRanTK_index -GTF $hg19_gtf -mbp -fmo
meRanCall -f $htseq_index -p $core_num -bam align_sorted.bam -gref -o call.txt -mBQ 20 -mr 0.01 -cr 0.99 -fdr 0.01 
rm *unmapped.fastq Homo_sapiens.*
mkdir ../4-combine_analysis
cd ../4-combine_analysis
cat ../3-meRanGhResult/call_FDR_0.01.txt|awk -v OFS="\t" '{if($12=="M")print $1,$2,$2,$3,$6,$5,$7,$12}' > methyC.bed
cat methyC.bed|awk -v OFS="\t" '$7>=0.1 && $6>=30 && $5>=5' > methyC_qua.bed
intersectBed -a methyC_qua.bed -b /pnas/yangyg_group/zhangting/Reference/human/ENSEMBLE_72/GeneList/single-transcript-chr -wa -wb|awk '$4==$12'|awk '!a[$0]++' > methyC_gene.txt
awk -v OFS='\t' 'NR==FNR{x[$1]=$2}NR>FNR{if(x[$14])print $0,x[$14]}'  /pnas/yangyg_group/zhangting/Reference/human/ENSEMBLE_72/GeneList/ensembl-genetype methyC_gene.txt> methyC_Gene.txt
cat methyC_Gene.txt|awk '$17=="protein_coding"'|sed 's/protein_coding/mRNA/'|cut -f1-14,16,17 > methy_C_mRNA.txt
cat methy_C_mRNA.txt|sed 's/_/\t/'|awk -v OFS="\t" '{print $1,$2,$3,$7,$15,$4}'|awk '!a[$1"\t"$2]++' > methylated.bed
cat methy_C_mRNA.txt|cut -f14|sort|uniq > genelist.txt
awk '{print $17}' methyC_Gene.txt> geneType.txt
mkdir ../7-mRNA_distribution
cat methylated.bed|awk -v OFS="\t" '{print $1,$2,$3,$6}' > ../7-mRNA_distribution/distribution_mRNA.bed
cd ../7-mRNA_distribution
cp /pnas/yangyg_group/zhangting/Reference/human/ENSEMBLE_72/GeneList/distribution_plot_H.pl ./
perl distribution_plot_H.pl

# mkdir ../8-motif
# cd ../8-motif
# cat ../4-combine_analysis/methylated.bed|awk -v OFS="\t" '{print $1,$2-11,$3+10,"*","*",$6}'|fastaFromBed -fi /pnas/yangyg_group/zhangting/Reference/human/ENSEMBLE_72/index/genome.fa -bed - -s -fo ./motif-20.fasta
# cat motif-20.fasta |grep -v ">"|sed 's/T/U/g' > weblogo-20.fa
# perl /pnas/yangyg_group/hanyn/data/perl/logo.pl
