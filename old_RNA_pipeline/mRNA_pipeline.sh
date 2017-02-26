#!bin/bash
project=PRJNA13758
reference=WS255
prefix=ftp://ftp.wormbase.org/pub/wormbase/releases/${reference}/species/c_elegans/${project}/

cd ~/Desktop/Pipeline/mRNA
#mkdir data

###  Download gene set
#curl ${prefix}/c_elegans.${project}.${reference}.canonical_geneset.gtf.gz > data/c_elegans.${project}.${reference}.canonical_geneset.gtf.gz

### Download reference genome
#curl ${prefix}/c_elegans.${project}.${reference}.genomic.fa.gz > data/c_elegans.${project}.${reference}.genomic.fa.gz

#gzcat data/c_elegans.${project}.${reference}.canonical_geneset.gtf.gz | python hisat2_extract_splice_sites.py - > data/${reference}.ss
#gzcat data/c_elegans.${project}.${reference}.canonical_geneset.gtf.gz | python hisat2_extract_exons.py - > data/${reference}.exon

###  Build HISAT2 index
#hisat2-build -p 1 --ss data/${reference}.ss --exon data/${reference}.exon data/c_elegans.${project}.${reference}.genomic.fa data/${reference}.hisat2_index

###  Clean up reads using Trimmomatic ()
read_dir=~/Desktop/Pipeline/test_data/
mkdir data/reads

#lane5 -> 1
trimmomatic SE -phred33 -threads 2  ${read_dir}/170112_K00242_0168_AHHF52BBXX-EA-RS8/EA-N2A_*fastq.gz data/reads/N2A_1.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
trimmomatic SE -phred33 -threads 2  ${read_dir}/170112_K00242_0168_AHHF52BBXX-EA-RS8/EA-N2B_*fastq.gz data/reads/N2B_1.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
trimmomatic SE -phred33 -threads 2  ${read_dir}/170112_K00242_0168_AHHF52BBXX-EA-RS8/EA-N2C_*fastq.gz data/reads/N2C_1.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
trimmomatic SE -phred33 -threads 2  ${read_dir}/170112_K00242_0168_AHHF52BBXX-EA-RS8/EA-N2D_*fastq.gz data/reads/N2D_1.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15

trimmomatic SE -phred33 -threads 2  ${read_dir}/170112_K00242_0168_AHHF52BBXX-EA-RS8/EA-CBA_*fastq.gz data/reads/CBA_1.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
trimmomatic SE -phred33 -threads 2  ${read_dir}/170112_K00242_0168_AHHF52BBXX-EA-RS8/EA-CBB_*fastq.gz data/reads/CBB_1.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
trimmomatic SE -phred33 -threads 2  ${read_dir}/170112_K00242_0168_AHHF52BBXX-EA-RS8/EA-CBC_*fastq.gz data/reads/CBC_1.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
trimmomatic SE -phred33 -threads 2  ${read_dir}/170112_K00242_0168_AHHF52BBXX-EA-RS8/EA-CBD_*fastq.gz data/reads/CBD_1.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15

#lane1/2 -> 2/3
trimmomatic SE -phred33 -threads 2  ${read_dir}/170118_K00242_0169_AHHCHNBBXX-EA-RS8/EA-N2A_*L001*fastq.gz data/reads/N2A_2.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
trimmomatic SE -phred33 -threads 2  ${read_dir}/170118_K00242_0169_AHHCHNBBXX-EA-RS8/EA-N2B_*L001*fastq.gz data/reads/N2B_2.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
trimmomatic SE -phred33 -threads 2  ${read_dir}/170118_K00242_0169_AHHCHNBBXX-EA-RS8/EA-N2C_*L001*fastq.gz data/reads/N2C_2.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
trimmomatic SE -phred33 -threads 2  ${read_dir}/170118_K00242_0169_AHHCHNBBXX-EA-RS8/EA-N2D_*L001*fastq.gz data/reads/N2D_2.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15

trimmomatic SE -phred33 -threads 2  ${read_dir}/170118_K00242_0169_AHHCHNBBXX-EA-RS8/EA-CBA_*L001*fastq.gz data/reads/CBA_2.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
trimmomatic SE -phred33 -threads 2  ${read_dir}/170118_K00242_0169_AHHCHNBBXX-EA-RS8/EA-CBB_*L001*fastq.gz data/reads/CBB_2.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
trimmomatic SE -phred33 -threads 2  ${read_dir}/170118_K00242_0169_AHHCHNBBXX-EA-RS8/EA-CBC_*L001*fastq.gz data/reads/CBC_2.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
trimmomatic SE -phred33 -threads 2  ${read_dir}/170118_K00242_0169_AHHCHNBBXX-EA-RS8/EA-CBD_*L001*fastq.gz data/reads/CBD_2.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15

trimmomatic SE -phred33 -threads 2  ${read_dir}/170118_K00242_0169_AHHCHNBBXX-EA-RS8/EA-N2A_*L002*fastq.gz data/reads/N2A_3.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
trimmomatic SE -phred33 -threads 2  ${read_dir}/170118_K00242_0169_AHHCHNBBXX-EA-RS8/EA-N2B_*L002*fastq.gz data/reads/N2B_3.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
trimmomatic SE -phred33 -threads 2  ${read_dir}/170118_K00242_0169_AHHCHNBBXX-EA-RS8/EA-N2C_*L002*fastq.gz data/reads/N2C_3.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
trimmomatic SE -phred33 -threads 2  ${read_dir}/170118_K00242_0169_AHHCHNBBXX-EA-RS8/EA-N2D_*L002*fastq.gz data/reads/N2D_3.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15

trimmomatic SE -phred33 -threads 2  ${read_dir}/170118_K00242_0169_AHHCHNBBXX-EA-RS8/EA-CBA_*L002*fastq.gz data/reads/CBA_3.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
trimmomatic SE -phred33 -threads 2  ${read_dir}/170118_K00242_0169_AHHCHNBBXX-EA-RS8/EA-CBB_*L002*fastq.gz data/reads/CBB_3.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
trimmomatic SE -phred33 -threads 2  ${read_dir}/170118_K00242_0169_AHHCHNBBXX-EA-RS8/EA-CBC_*L002*fastq.gz data/reads/CBC_3.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
trimmomatic SE -phred33 -threads 2  ${read_dir}/170118_K00242_0169_AHHCHNBBXX-EA-RS8/EA-CBD_*L002*fastq.gz data/reads/CBD_3.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15

#lane5/6 -> 4/5
trimmomatic SE -phred33 -threads 2  ${read_dir}/170118_K00242_0170_BHHCYYBBXX-EA-RS8/EA-N2A_*L005*fastq.gz data/reads/N2A_4.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
trimmomatic SE -phred33 -threads 2  ${read_dir}/170118_K00242_0170_BHHCYYBBXX-EA-RS8/EA-N2B_*L005*fastq.gz data/reads/N2B_4.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
trimmomatic SE -phred33 -threads 2  ${read_dir}/170118_K00242_0170_BHHCYYBBXX-EA-RS8/EA-N2C_*L005*fastq.gz data/reads/N2C_4.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
trimmomatic SE -phred33 -threads 2  ${read_dir}/170118_K00242_0170_BHHCYYBBXX-EA-RS8/EA-N2D_*L005*fastq.gz data/reads/N2D_4.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15

trimmomatic SE -phred33 -threads 2  ${read_dir}/170118_K00242_0170_BHHCYYBBXX-EA-RS8/EA-CBA_*L005*fastq.gz data/reads/CBA_4.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
trimmomatic SE -phred33 -threads 2  ${read_dir}/170118_K00242_0170_BHHCYYBBXX-EA-RS8/EA-CBB_*L005*fastq.gz data/reads/CBB_4.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
trimmomatic SE -phred33 -threads 2  ${read_dir}/170118_K00242_0170_BHHCYYBBXX-EA-RS8/EA-CBC_*L005*fastq.gz data/reads/CBC_4.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
trimmomatic SE -phred33 -threads 2  ${read_dir}/170118_K00242_0170_BHHCYYBBXX-EA-RS8/EA-CBD_*L005*fastq.gz data/reads/CBD_4.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15

trimmomatic SE -phred33 -threads 2  ${read_dir}/170118_K00242_0170_BHHCYYBBXX-EA-RS8/EA-N2A_*L006*fastq.gz data/reads/N2A_5.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
trimmomatic SE -phred33 -threads 2  ${read_dir}/170118_K00242_0170_BHHCYYBBXX-EA-RS8/EA-N2B_*L006*fastq.gz data/reads/N2B_5.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
trimmomatic SE -phred33 -threads 2  ${read_dir}/170118_K00242_0170_BHHCYYBBXX-EA-RS8/EA-N2C_*L006*fastq.gz data/reads/N2C_5.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
trimmomatic SE -phred33 -threads 2  ${read_dir}/170118_K00242_0170_BHHCYYBBXX-EA-RS8/EA-N2D_*L006*fastq.gz data/reads/N2D_5.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15

trimmomatic SE -phred33 -threads 2  ${read_dir}/170118_K00242_0170_BHHCYYBBXX-EA-RS8/EA-CBA_*L006*fastq.gz data/reads/CBA_5.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
trimmomatic SE -phred33 -threads 2  ${read_dir}/170118_K00242_0170_BHHCYYBBXX-EA-RS8/EA-CBB_*L006*fastq.gz data/reads/CBB_5.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
trimmomatic SE -phred33 -threads 2  ${read_dir}/170118_K00242_0170_BHHCYYBBXX-EA-RS8/EA-CBC_*L006*fastq.gz data/reads/CBC_5.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
trimmomatic SE -phred33 -threads 2  ${read_dir}/170118_K00242_0170_BHHCYYBBXX-EA-RS8/EA-CBD_*L006*fastq.gz data/reads/CBD_5.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15


### Align reads using HISAT2

mkdir data/align

for filename in data/reads/*.fq.gz; do
   hisat2 -p 2 -x data/WS255.hisat2_index -U "$filename" -S data/align/$(basename "$filename" .fq.gz).sam
done

for filename in data/align/*.sam; do
	samtools view -bS "$filename" > data/align/$(basename "$filename" .sam).unsorted.bam
	samtools flagstat data/align/$(basename "$filename" .sam).unsorted.bam
	samtools sort -@ 2 -o data/align/$(basename "$filename" .sam).bam data/align/$(basename "$filename" .sam).unsorted.bam
	samtools index -b data/align/$(basename "$filename" .sam).bam
done
rm data/align/*unsorted.bam
rm data/align/*sam

# Join bam files for samples split across lanes
for sample in N2A N2B N2C N2D CBA CBB CBC CBD; do
	samtools merge data/align/"$sample".unsorted.bam data/align/"$sample"_*.bam
	samtools sort -@ 2 -o data/align/"$sample".bam data/align/"$sample".unsorted.bam
	samtools index -b data/align/"$sample".bam
done
rm data/align/*unsorted.bam

### Get gene counts (stringtie) based on existing annotations

mkdir data/expression
#gzcat data/c_elegans.${project}.${reference}.canonical_geneset.gtf.gz > data/c_elegans.${project}.${reference}.canonical_geneset.gtf
Ref_GTF=data/c_elegans.${project}.${reference}.canonical_geneset.gtf

for sample in N2A N2B N2C N2D CBA CBB CBC CBD; do
	mkdir /data/expression/"$sample"
	stringtie -p 2 -G $Ref_GTF -e -B -o data/expression/"$sample"/"$sample"_expressed.gtf data/align/"$sample".bam
done

### Prepare DiffExp for R/bioconducter
mkdir data/diffexp
python prepDE.py -i data/expression/ -l 50 -g data/diffexp/gene_count_matrix.csv -t data/diffexp/transcript_count_matrix.csv

### Differential Expression with Ballgown
# creat file with sample and rep information

