# ABCC3 snippet code for TEAD4 ChIP-seq
# Data retrieved starting from ChIP-Atlas

cd /mnt/d/projects/ABCC3/
mkdir ChIP

# Download dataset from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84389
prefetch -c SRR3922075	# Input https://www.ncbi.nlm.nih.gov/sra?term=SRX1948681
prefetch -c SRR3922074 	#TEAD4 https://www.ncbi.nlm.nih.gov/sra?term=SRX1948680

# convert to fastq
fastq-dump -I --split-files SRR3922075 --outdir /mnt/d/projects/ABCC3/ChIP # split each read in separate files (Input)
fastq-dump -I --split-files SRR3922074 --outdir /mnt/d/projects/ABCC3/ChIP # (TEAD4)

cd ChIP
fastqc SRR3922075_1.fastq
fastqc SRR3922074_1.fastq

# Build hg38 Bowtie2 Index
bowtie2-build --large-index /mnt/d/genomes/Human/Homo_sapiens.GRCh38.dna.primary_assembly.fa /mnt/d/genomes/Human/bowtie/hg38

# assign index variable
index=/mnt/d/genomes/Human/bowtie/hg38
mkdir bowtie

# Align reads with bowtie2
bowtie2 -k 1 -p 10 -x $index -U SRR3922075_1.fastq | samtools view -uS - > bowtie/input.bam
bowtie2 -k 1 -p 10 -x $index -U SRR3922074_1.fastq | samtools view -uS - > bowtie/tead4.bam

# Sort and index BAMs
cd bowtie
samtools sort -m 2G -@ 10 -O BAM -o input.sorted.bam input.bam
samtools index input.sorted.bam
samtools sort -m 2G -@ 10 -O BAM -o tead4.sorted.bam tead4.bam
samtools index tead4.sorted.bam
