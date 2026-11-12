#!/bin/bash

#### Environment setup and name maps ####

### Environment setup

eval "$(pixi shell-hook --manifest-path /mnt/share1/scratch/haodong/software/pixi_env/base)"

Software=/mnt/share1/scratch/haodong/software/
Genome=/mnt/share1/scratch/haodong/genome/
BLACKLIST=/mnt/share1/scratch/haodong/genome/hg38/ENCFF356LFX_Blacklit_hg38_sorted_chr.bed
dir=/mnt/share1/scratch/haodong/projects/EN00012472_6/
Adapters=/mnt/share1/scratch/haodong/genome/adapters/
chromSize=/mnt/share1/scratch/haodong/genome/hg38/hg38.chroms.sizes_chr
seacr=/mnt/share1/scratch/haodong/software/SEACR/SEACR_1.3.sh

export TMPDIR=/mnt/share1/scratch/haodong/software/tmp
export TEMP=/mnt/share1/scratch/haodong/software/tmp
export TMP=/mnt/share1/scratch/haodong/software/tmp
mkdir -p "${TMPDIR}" 

### Name maps 

## Biological name maps

declare -A name_map=(
    ["CUTandRUNWT1"]="WT_H3K4me3_rep1"
    ["CUTandRUNWT2"]="WT_H3K4me3_rep2"
    ["CUTandRUNWT3"]="WT_IgG_control"

    ["CUTandRUNG31"]="SETD1B_1948var_G3_H3K4me3_rep1"
    ["CUTandRUNG32"]="SETD1B_1948var_G3_H3K4me3_rep2"
    ["CUTandRUNG33"]="SETD1B_1948var_G3_IgG_control"

    ["CUTandRUN811"]="SETD1B_mutant_clone_81_H3K4me3_rep1"
    ["CUTandRUN812"]="SETD1B_mutant_clone_81_H3K4me3_rep2"
    ["CUTandRUN813"]="SETD1B_mutant_clone_81_IgG_control"

    ["CUTandRUN881"]="SETD1B_mutant_clone_88_H3K4me3_rep1"
    ["CUTandRUN882"]="SETD1B_mutant_clone_88_H3K4me3_rep2"
    ["CUTandRUN883"]="SETD1B_mutant_clone_88_IgG_control"
)

## H3K4me3 samples maps

declare -A H3K4me3_map=(
    ["WT"]="WT_H3K4me3_rep1 WT_H3K4me3_rep2"
    ["G3"]="SETD1B_1948var_G3_H3K4me3_rep1 SETD1B_1948var_G3_H3K4me3_rep2"
    ["81"]="SETD1B_mutant_clone_81_H3K4me3_rep1 SETD1B_mutant_clone_81_H3K4me3_rep2"
    ["88"]="SETD1B_mutant_clone_88_H3K4me3_rep1 SETD1B_mutant_clone_88_H3K4me3_rep2"
)

################################done##############################

#### Rename FASTQ files and first QC ####

### Rename FASTQ files to required format

for sample in "${dir}"/CRData/*_1.fastq.gz; do
    base=$(basename "$sample" "_1.fastq.gz")
    echo "${base}"
    mv "${dir}/CRData/${base}_1.fastq.gz" "${dir}/CRData/${base}_R1.fastq.gz"
    mv "${dir}/CRData/${base}_2.fastq.gz" "${dir}/CRData/${base}_R2.fastq.gz"
done 

### First QC

mkdir -p "${dir}/QC/QC_1"

fastqc -t 8 "${dir}"/CRData/*.gz -o "${dir}/QC/QC_1/"
multiqc "${dir}/QC/QC_1/." -o "${dir}/QC/QC_1/" 

echo "First_QC is finished"

################################done##############################

#### Trimming and second QC ####

### Trimming

mkdir -p "${dir}/Trim"

for sample in "${dir}"/CRData/*R1.fastq.gz; do  
    base=$(basename "$sample" "_R1.fastq.gz")
    echo "${base}"
    trimmomatic PE -Xmx16g -threads 8 -phred33 \
    "${dir}/CRData/${base}_R1.fastq.gz" "${dir}/CRData/${base}_R2.fastq.gz" \
    "${dir}/Trim/${base}_R1.trim.fastq.gz" "${dir}/Trim/${base}_R1.unpaired.trim.fastq.gz" \
    "${dir}/Trim/${base}_R2.trim.fastq.gz" "${dir}/Trim/${base}_R2.unpaired.trim.fastq.gz" \
    ILLUMINACLIP:"${Adapters}"/TruSeq3-PE-2.fa:2:30:10:8:true MINLEN:20
done
 
rm "${dir}"/Trim/*.unpaired.trim.fastq.gz

echo "Trimming is finished"

### Second QC

mkdir -p "${dir}/QC/QC_2/"

fastqc -t 8 "${dir}"/Trim/*.trim.fastq.gz -o "${dir}/QC/QC_2/" 
multiqc "${dir}/QC/QC_2/." -o "${dir}/QC/QC_2/" 

echo "Second_QC is finished"

################################done##############################

#### Alignment and deduplication ####

### Alignment

mkdir -p "${dir}/Align_hg38"

for sample in "${dir}"/Trim/*R1.trim.fastq.gz; do
    base=$(basename "$sample" "_R1.trim.fastq.gz")
    echo "${base}"
    bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 \
    -I 10 -X 700 -p 8 \
    -x "${Genome}"/Genome_Index/GRCH_38/GRCH-38_RefGenome/bowtie2/hg38 \
    -1 "${dir}/Trim/${base}_R1.trim.fastq.gz" -2 "${dir}/Trim/${base}_R2.trim.fastq.gz" \
    | samtools view -u - \
    | samtools sort -@ 8 - \
    | samtools view -@ 8 -h -b -q 30 > "${dir}/Align_hg38/${base}.bam"
done

echo "Aligning is finished"

### Remove duplication

mkdir -p "${dir}/Align_hg38/Deduplicated"

for sample in "${dir}"/Align_hg38/*.bam; do
    base=$(basename "$sample" ".bam")
    echo "${base}"
    java -Xmx16g -jar "${Software}/picard.jar" MarkDuplicates \
    I="${dir}/Align_hg38/${base}.bam" \
    O="${dir}/Align_hg38/Deduplicated/Dedup_${base}.bam" \
    M="${dir}/Align_hg38/Deduplicated/Dedup_${base}.txt" \
    REMOVE_DUPLICATES=TRUE CREATE_INDEX=FALSE VALIDATION_STRINGENCY=SILENT \
    TMP_DIR="${TMPDIR}"
done

echo "Deduplication is finished"

### Make summary of alignment and prepare input files for SEACR peak calling

mkdir -p "${dir}/UCSC_hg38"

for sample in "${dir}"/Align_hg38/Deduplicated/*.bam; do 
    base=$(basename "$sample" ".bam")
    echo "${base}"
    samtools flagstat -@ 8 "${dir}/Align_hg38/Deduplicated/${base}.bam" > "${dir}/Align_hg38/Deduplicated/${base}_flagstat.txt"
    samtools index -@ 8 "${dir}/Align_hg38/Deduplicated/${base}.bam"
    samtools sort -n -@ 8 -o "${TMPDIR}/${base}.namesorted.bam" "${dir}/Align_hg38/Deduplicated/${base}.bam"
    bedtools bamtobed -bedpe -i "${TMPDIR}/${base}.namesorted.bam" > "${dir}/Align_hg38/Deduplicated/${base}.bed"
    awk '$1==$4 && $6-$2 < 1000 {print $0}' "${dir}/Align_hg38/Deduplicated/${base}.bed" > "${dir}/Align_hg38/Deduplicated/${base}.clean.bed"
    cut -f 1,2,6 "${dir}/Align_hg38/Deduplicated/${base}.clean.bed" | sort -k1,1 -k2,2n -k3,3n > "${dir}/Align_hg38/Deduplicated/${base}.fragments.bed"
    bedtools genomecov -bg -i "${dir}/Align_hg38/Deduplicated/${base}.fragments.bed" -g "${chromSize}" > "${dir}/UCSC_hg38/${base}.fragments.raw.bedgraph"
done 

## Update bedGraph filenames by removing the 'Dedup_' prefix
for old in "${!name_map[@]}"; do 
    new=${name_map[$old]}
    old_file="${dir}/UCSC_hg38/Dedup_${old}.fragments.raw.bedgraph"
    new_file="${dir}/UCSC_hg38/${new}.fragments.raw.bedgraph"
    if [ -f "$old_file" ]; then
        mv "$old_file" "$new_file"
        echo "Renamed: $old_file → $new_file"
    else
        echo "WARNING: $old_file not found"
    fi
done

echo "preparing input bedgraph files is finished"

################################done##############################

#### SEACR peak calling ####
#### Utilize the numeric threshold 0.01 only for H3K4me3 samples

mkdir -p "${dir}/Peak_Calling"

for group in WT G3 81 88; do
    for sample in ${H3K4me3_map[$group]}; do
        echo "SEACR: ${sample}, top 1% by total signal "
        bash "$seacr" "${dir}/UCSC_hg38/${sample}.fragments.raw.bedgraph" \
             0.01 \
             non stringent "${dir}/Peak_Calling/${sample}.seacr.peaks"
    done
done

echo "SEACR peak calling is finished"

################################done##############################

#### Remove blacklist regions ####

for sample in "${dir}"/Peak_Calling/*.seacr.peaks.stringent.bed; do
    base=$(basename "$sample" ".seacr.peaks.stringent.bed")
    sort -k1,1 -k2,2n -k3,3n "${dir}/Peak_Calling/${base}.seacr.peaks.stringent.bed" | uniq | \
    bedtools intersect -v -a stdin -b "${BLACKLIST}" | \
    sort -k1,1 -k2,2n -k3,3n | uniq > "${dir}/Peak_Calling/${base}.seacr.peaks.stringent.FilteredBlacklist.bed"
done

################################done##############################

#### Consensus peaks for downstream analysis ####

 for group in WT G3 81 88; do
     rep1="${H3K4me3_map[$group]%% *}"
     rep2="${H3K4me3_map[$group]##* }"
    
     rep1_file="${dir}/Peak_Calling/${rep1}.seacr.peaks.stringent.FilteredBlacklist.bed"
     rep2_file="${dir}/Peak_Calling/${rep2}.seacr.peaks.stringent.FilteredBlacklist.bed"
    
     echo "Processing ${group}: $rep1_file x $rep2_file"
    
     bedtools intersect -a "$rep1_file" -b "$rep2_file" -wa -wb | \
         cut -f 1,2,3,4 | \
         sort -k1,1 -k2,2n -k3,3n | uniq > "${dir}/Peak_Calling/${group}_rep1.bed"
    
     bedtools intersect -a "$rep1_file" -b "$rep2_file" -wa -wb | \
         cut -f 7,8,9,10 | \
         sort -k1,1 -k2,2n -k3,3n | uniq > "${dir}/Peak_Calling/${group}_rep2.bed"
    
     sort -k1,1 -k2,2n -k3,3n "${dir}/Peak_Calling/${group}_rep1.bed" "${dir}/Peak_Calling/${group}_rep2.bed" | \
         bedtools merge -i stdin -c 4 -o mean > "${dir}/Peak_Calling/${group}_consensus.bed"
    
     rm "${dir}/Peak_Calling/${group}_rep1.bed" "${dir}/Peak_Calling/${group}_rep2.bed"
 done

### Merge consensus peaks from all four groups for downstream analysis 

 sort -k1,1 -k2,2n -k3,3n \
     "${dir}/Peak_Calling/81_consensus.bed" \
     "${dir}/Peak_Calling/88_consensus.bed" \
     "${dir}/Peak_Calling/G3_consensus.bed" \
     "${dir}/Peak_Calling/WT_consensus.bed" | \
     bedtools merge -i - > "${dir}/Peak_Calling/All_consensus.bed"

 echo "Consensus peaks and All_consensus.bed are generated"

### Count fragments within all consensus peak for each biological replicate 
### Then use DESeq2/... to identify differential H3K4me3 peaks between genotypes

################################done##############################

#### Visualization in IGV ####

for sample in "${dir}"/Align_hg38/Deduplicated/*.bam; do
    base=$(basename "$sample" ".bam")
    echo "${base}"
    
    old="${base#Dedup_}"
    new="${name_map[$old]}"
    
    bamCoverage -p 8 -bs 10 --normalizeUsing CPM --centerReads \
        --effectiveGenomeSize 2913022398 \
        -b "${dir}/Align_hg38/Deduplicated/${base}.bam" \
        -o "${dir}/UCSC_hg38/${new}.bw"
done

echo "Visualization in UCSC is finished"

################################done##############################

