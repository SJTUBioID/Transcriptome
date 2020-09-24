#!/bin/bash


mkdir Fastq
mkdir QC
mkdir Bam
mkdir Ctab
mkdir Stat
function Transcriptome()
{
infile1=$1
infile2=$2
sample=$3
ADAPTER=AGATCGGAAGAG
index=/BioID/Software/hisat2-2.0.5/index/grcm38_tran/genome_tran
gtf=/BioID/Annotation/Transcriptome/mouse/Mus_musculus.GRCm38.93.chr.gtf
rRNA_databases=/BioID/Software/SortMeRNA/rRNA_databases
rindex=/BioID/Software/SortMeRNA/index

N=8

trim1=$sample"_trim1_fq.gz"
trim2=$sample"_trim2_fq.gz"
m1=$sample"_Clean_R1.fq.gz"
u1=$sample"_unpaired_R1.fq.gz"
m2=$sample"_Clean_R2.fq.gz"
u2=$sample"_unpaired_R2.fq.gz"
c1=$sample"_rRNA-RM.R1.fq"
c2=$sample"_rRNA-RM.R2.fq"
sam=$sample".sam"
bam=$sample"_raw.bam"
Sort=$sample".bam"
stat="Stat/"$sample".stat"



# Step1 "Data Clean"
cutadapt -u 5 -a $ADAPTER -A $ADAPTER --max-n 0 --minimum-length 100 -o $trim1 -p $trim2 $infile1 $infile2
java -jar /BioID/Software/Trimmomatic-0.35/trimmomatic-0.35.jar PE -threads 6 $trim1 $trim2 $m1 $u1 $m2 $u2 SLIDINGWINDOW:3:10 LEADING:10 TRAILING:10 MINLEN:100
rm $u1
rm $u2
rm $trim1
rm $trim2



#Step 3 "m1 m2 Fastq Statistics"
echo "Sample:" $sample >> $stat
line=`zcat $m1|wc -l|awk '{print $1}' `
let reads=$line/4
echo "Total Clean Read Pairs:" $reads >> $stat
base1=`zcat $m1 |awk '{if(NR%4==2){a=length($0);s+=a}}END{print s}' `
base2=`zcat $m2 |awk '{if(NR%4==2){a=length($0);s+=a}}END{print s}' `
let totalbases=$base1+$base2
echo "R1 Bases:" $base1 >> $stat
echo "R2 Bases:" $base2 >> $stat
echo "Total Bases:" $totalbases >> $stat
echo "\n" >> $stat



# Step 4 "Remove rRNA"
pigz -p 10 -cd $m1 >> $sample"_R1.fq"
pigz -p 10 -cd $m2 >> $sample"_R2.fq"
mv $m1 Fastq
mv $m2 Fastq

bash /BioID/Software/SortMeRNA/scripts/merge-paired-reads.sh $sample"_R1.fq" $sample"_R2.fq" $sample".merged.fq"
echo "Merge"

rm -r /home/cshu/kvdb/
sortmerna --ref $rRNA_databases"/silva-bac-16s-id90.fasta",$rindex"/silva-bac-16s-db":\
$rRNA_databases"/silva-bac-23s-id98.fasta",$rindex"/silva-bac-23s-db":\
$rRNA_databases"/silva-arc-16s-id95.fasta",$rindex"/silva-arc-16s-db":\
$rRNA_databases"/silva-arc-23s-id98.fasta",$rindex"/silva-arc-23s-db":\
$rRNA_databases"/silva-euk-18s-id95.fasta",$rindex"/silva-euk-18s-db":\
$rRNA_databases"/silva-euk-28s-id98.fasta",$rindex"/silva-euk-28s-db":\
$rRNA_databases"/rfam-5s-database-id98.fasta",$rindex"/rfam-5s-db":\
$rRNA_databases"/rfam-5.8s-database-id98.fasta",$rindex"/rfam-5.8s-db" --reads $sample".merged.fq" --sam --num_alignments 1 --fastx --paired_in --aligned $sample".garbagerRNA" --other $sample".non_rRNA" -a 8 -m 32000 --log
bash /BioID/Software/SortMeRNA/scripts/unmerge-paired-reads.sh $sample".non_rRNA.fq" $c1 $c2


#Step 4 "Fastqc"
rm $sample"_R1.fq"
rm $sample"_R2.fq"
rm $sample".merged.fq"
rm $sample".non_rRNA.fq"
rm $sample.garbagerRNA*

pigz -p 10 $c1
pigz -p 10 $c2

/BioID/Software/FastQC/fastqc -t 4 $c1".gz"
/BioID/Software/FastQC/fastqc -t 4 $c2".gz"
mv *zip QC/
rm *html

#Step 5 "c1 c2 Fastq Statistics"
echo "After rRMA remove" $sample >> $stat
line=`zcat $c1".gz"|wc -l|awk '{print $1}' `
let reads=$line/4
echo "Total Clean Read Pairs:" $reads >> $stat
base1=`zcat $c1".gz" |awk '{if(NR%4==2){a=length($0);s+=a}}END{print s}' `
base2=`zcat $c2".gz" |awk '{if(NR%4==2){a=length($0);s+=a}}END{print s}' `
let totalbases=$base1+$base2
echo "R1 Bases:" $base1 >> $stat
echo "R2 Bases:" $base2 >> $stat
echo "Total Bases:" $totalbases >> $stat
echo "\n" >> $stat




#Step 4 HISAT2 mapping
/BioID/Software/hisat2-2.0.5/hisat2 --rna-strandness RF -p $N --no-softclip --dta  -x $index -1 $c1".gz" -2 $c2".gz" -S $sam 
mv $c1".gz" Fastq/
mv $c2".gz" Fastq/

#Step 5 StringTie
samtools view -bS -q 30 $sam >> $bam
samtools sort -@ 8 $bam -o $Sort
samtools index $Sort

mkdir "Ctab/"$sample
/BioID/Software/stringtie-1.3.3/stringtie  -e -b "Ctab/"$sample -p $N --rf -G $gtf -C "Ctab/"$sample"/"$sample"_Covered.gtf" -o "Ctab/"$sample"/"$sample".gtf" -A "Ctab/"$sample"/"$sample".gene.ctab" -l $sample $Sort
samtools flagstat $Sort >> $stat
rm $sam
rm $bam
mv $Sort Bam/
mv $Sort.bai Bam/
}


# Demo
# Transcriptome $R1 $R2 $sample

