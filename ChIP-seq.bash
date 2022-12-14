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

# Markup and remove duplicates
 PicardCommandLine MarkDuplicates \
 REMOVE_DUPLICATES=TRUE \
 I="tead4.sorted.bam" O=tead4_NODUPS.bam \
 M=tead4_dup_metrics.txt 

 PicardCommandLine MarkDuplicates \
 REMOVE_DUPLICATES=TRUE \
 I="input.sorted.bam" O=input_NODUPS.bam \
 M=input_dup_metrics.txt 

# Blacklist regions

# ENCODE blacklists (save and unzip hg38)
# from http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/


# Blacklist with ENCODE regions
bl=/mnt/f/genomes/Human/blacklist/hg38.blacklist.bed

# Filter out black lists
mkdir filtered
bedtools intersect -v -abam input_NODUPS.bam -b $bl > filtered/input.filtered.bam
bedtools intersect -v -abam tead4_NODUPS.bam -b $bl > filtered/tead4.filtered.bam
samtools index filtered/input.filtered.bam
samtools index filtered/tead4.filtered.bam


#Call MACS peaks (one to one samples)
cd filtered
macs3 callpeak -t tead4.filtered.bam \
-c input.filtered.bam \
-f BAM -g hs -n tead4 -B -p 1e-9 --outdir macs_tead4 --verbose 2 --bdg 

# Prepare a bed file centered on transcript TSSs
gff3=/mnt/f/genomes/Human/annotations/gencode.v36.annotation.gff3
# Select only transcripts
grep -P "\ttranscript\t" $gff3 > alltranscripts_coords.bed
cut -f3 alltranscripts_coords.bed | sort | uniq

cd /mnt/d/projects/ABCC3/ChIP/bowtie/filtered/
# Remove redundant genes in custom .bed file
uniq up.bed > tmp && mv tmp genesup.bed
less -S genesup.bed | wc -l
uniq dn.bed > tmp && mv tmp genesdn.bed
uniq nonhit.bed > tmp && mv tmp genesnot.bed

# Starting the heatmap plot on TSS (center) <- center on TSS signals
# Files needed:
### BAM from ChIP-Seq (TEAD4)
### BAM from input
### BED file with TSS coordinates

cd /mnt/d/projects/ABCC3/ChIP/bowtie/filtered
# Calculate number of aligned reads
samtools flagstat "input.filtered.bam"  # 41135748
samtools flagstat "tead4.filtered.bam" # 47087786

# Convert BAM to BigWig
bamCoverage -p 6 -b "input.filtered.bam"  -o input.bw --scaleFactor 1
bamCoverage -p 6 -b "tead4.filtered.bam" -o tead4.bw --scaleFactor 0.8735 # That is, 41135748/47087786

# Subtract input from signal
bigwigCompare -b1 tead4.bw -b2 input.bw --operation subtract -p 6 -o tead4_diff.bw

# Choose reference point
bedup=/mnt/d/projects/ABCC3/ChIP/bowtie/filtered/genesup.bed
beddn=/mnt/d/projects/ABCC3/ChIP/bowtie/filtered/genesdn.bed
bednot=/mnt/d/projects/ABCC3/ChIP/bowtie/filtered/genesnot.bed
computeMatrix reference-point -p4 -S tead4_diff.bw -R $bedup -a 10000 -b 10000 --referencePoint center -o upanalysis.mat.gz
computeMatrix reference-point -p4 -S tead4_diff.bw -R $beddn -a 10000 -b 10000 --referencePoint center -o dnanalysis.mat.gz
computeMatrix reference-point -p4 -S tead4_diff.bw -R $bednot -a 10000 -b 10000 --referencePoint center -o notanalysis.mat.gz

# Plot
plotHeatmap -m upanalysis.mat.gz -out plots/001_deeptools_up.png --averageTypeSummaryPlot mean \
--heatmapHeight 13 --heatmapWidth 6 --whatToShow "heatmap and colorbar" \
--colorList 'navajowhite,orange,black' \
--refPointLabel "TEAD4 peak centers" --plotTitle "ABCC3 targets vs. TEAD4 ChIP-Seq" \
--xAxisLabel "Distance from TSS center (bp)" --yAxisLabel "peaks" \
--dpi 600

plotHeatmap -m dnanalysis.mat.gz -out plots/001_deeptools_dn.png --averageTypeSummaryPlot mean \
--heatmapHeight 13 --heatmapWidth 6 --whatToShow "heatmap and colorbar" \
--colorList 'navajowhite,orange,black' \
--refPointLabel "TEAD4 peak centers" --plotTitle "ABCC3 targets vs. TEAD4 ChIP-Seq" \
--xAxisLabel "Distance from TSS center (bp)" --yAxisLabel "peaks" \
--dpi 600

plotHeatmap -m notanalysis.mat.gz -out plots/001_deeptools_not.png --averageTypeSummaryPlot mean \
--heatmapHeight 13 --heatmapWidth 6 --whatToShow "heatmap and colorbar" \
--colorList 'navajowhite,orange,black' \
--refPointLabel "TEAD4 peak centers" --plotTitle "ABCC3 targets vs. TEAD4 ChIP-Seq" \
--xAxisLabel "Distance from TSS center (bp)" --yAxisLabel "peaks" \
--dpi 600
