
cd /pnas/yangyg_group/hanyn/14-bladder/YBX1/Rep2
mkdir 1-data
mkdir 2-fastqc
mkdir 3-bowtie
mkdir 4-paralyzer
cd 1-data
gunzip -c /pnas/yangyg_group/sunbf/bladder_cancer_guangzhou/YBX1_PARCLIP_Rep2/Y180517-1-LA-PAR_HL573CCXY_L3_1.fq.gz > R1.fastq

## data normal##
cutadapt -a AGATCGGAAGAGCACACGTCTG -o 1-trim_1.fastq R1.fastq
java -Xmx4g -jar /software/biosoft/software/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 ./1-trim_1.fastq ./2-trim_qua_1.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:15
fastqc 1-trim_1.fastq -t 10 -o ../2-fastqc/
fastqc R1.fastq -t 10 -o ../2-fastqc/
fastqc 2-trim_qua_1.fastq -t 10 -o ../2-fastqc/
## data delete A and 4bse ##
cutadapt -a AGATCGGAAGAGCACACGTCTG -o 1-trim_1.fastq R1.fastq
cutadapt -a AAAAAAAA -o 1-trim_1_1.fastq 1-trim_1.fastq
cutadapt -u 4 -o 1-trim_2_1.fastq 1-trim_1_1.fastq
java -Xmx4g -jar /software/biosoft/software/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 ./1-trim_2_1.fastq ./2-trim_qua_1.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:15
fastqc R1.fastq -t 10 -o ../2-fastqc/
fastqc 1-trim_1_1.fastq -t 10 -o ../2-fastqc/
fastqc 1-trim_2_1.fastq -t 12 -o ../2-fastqc/
fastqc 2-trim_qua_1.fastq -t 12 -o ../2-fastqc/

### mapping
cd ../3-bowtie/
genome=/pnas/yangyg_group/zhangting/Reference/human/ENSEMBLE_72/index/genome.fa
bowtie $genome -p 12 -v 2 -m 10 --best --strata ../1-data/2-trim_qua_1.fastq bowtie.sam  ### -S 输出sam格式；-v 2 -m 10（允许两个错配，最多允许10个multiple mapping）；--best，对multiple mapping，只输出得分最高的
grep -v "chrMT" bowtie.sam| awk 'length($6) < 100 {print $0}' > bowtie_rmMT.sam
#  Mycoplasma_Hyorhinis
genomeM=/pnas/yangyg_group/hanyn/data/reference/Mycoplasma_Hyorhinis/Mycoplasma_Hyorhinis.fasta
bowtie $genomeM -p 12 -v 2 -m 10 --best --strata ../1-data/2-trim_qua_1.fastq bowtie_Mycoplasma.sam  ### -S 输出sam格式；-v 2 -m 10（允许两个错配，最多允许10个multiple mapping）；--best，对multiple mapping，只输出得分最高的

### TC转换率的统计
echo 'mapped_reads:' > ct_conversion.info
wc bowtie_rmMT.sam >> ct_conversion.info
echo ''  >> ct_conversion.info
echo 'mutated reads:' >> ct_conversion.info
cat bowtie_rmMT.sam |cut -f8|awk '$1'|wc >> ct_conversion.info
echo ''  >> ct_conversion.info
echo 'mutated nucleotide:' >> ct_conversion.info
cat bowtie_rmMT.sam |cut -f8|awk '$1'|awk '{split($1,a,",");for(i=1;i<=length(a);i++){print a[i]}}'|awk -v FS=":" '{print $2}'|sort|uniq -c|sed -e 's/^ *//'|sed -e 's/ /\t/'  >> ct_conversion.info
echo ''  >> ct_conversion.info
echo 'conversion ratio:' >> ct_conversion.info
cat bowtie_rmMT.sam |cut -f8|awk '$1'|awk '{split($1,a,",");for(i=1;i<=length(a);i++){print a[i]}}'|awk -v FS=":" '{print $2}'|sort|uniq -c|sed -e 's/^ *//'|sed -e 's/ /\t/'|awk -v OFS="\t" '{num+=$1;tmp[$2]=$1}END{for(i in tmp){print i,tmp[i],num,tmp[i]/num}}'|sort  >> ct_conversion.info

##  htseq and bw
mkdir bowtie_True
cd bowtie_True
genome=/pnas/yangyg_group/zhangting/Reference/human/ENSEMBLE_72/index/genome.fa
bowtie $genome -p 12 -v 2 -m 10 --best --strata -S ../../1-data/2-trim_qua_1.fastq bowtie_true.sam  ### -S 输出sam格式；-v 2 -m 10（允许两个错配，最多允许10个multiple mapping）；--best，对multiple mapping，只输出得分最高的
hg19_gtf="/pnas/yangyg_group/zhangting/Reference/human/ENSEMBLE_72/GTF/Homo_sapiens.GRCh37.72.chr.gtf"
htseq-count -m union -s reverse bowtie_true.sam $hg19_gtf > htseq.txt
grep 'chr1'$'\t' bowtie_true.sam > bowtie_chr1.sam
samtools view -@ 12 -b -q 20 bowtie_chr1.sam > uniqmap_chr1.bam
samtools sort -@ 12 uniqmap_chr1.bam > sort_chr1.bam
samtools index sort_chr1.bam sort_chr1.bai
export LD_LIBRARY_PATH=/software/biosoft/software/openssl1.0.1/lib:$LD_LIBRARY_PATH
bamCoverage -b sort_chr1.bam -o coverage_chr1.bw
cd ..


cd ../4-paralyzer/
mkdir para1

cd para1
cp /pnas/yangyg_group/hanyn/data/parclip/PAR-CLIP_para1.ini ./
sed -i "s?MYPATH?$(pwd)?g" PAR-CLIP_para1.ini
sh /pnas/yangyg_group/hanyn/data/parclip/PARalyzer PAR-CLIP_para1.ini
sed '1d' cluster > cluster1
awk -F ',' '{OFS="\t"}{print $1,$3,$4,$2}' cluster1|sort -u > cluster.bed
intersectBed -a cluster.bed -b /pnas/yangyg_group/zhangting/Reference/human/ENSEMBLE_72/GeneList/single-transcript-chr -wa -wb|awk '$4==$8 && $11=="protein_coding"' |cut -f 10 - |sort -u > target_protein_coding.txt
rm cluster1
