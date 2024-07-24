#!/bin/bash

#### Variant Calling steps using GATK workflow
## Steps run on a single sample

# Packages/tools required pre-installed
# use BWA
# use Bowtie2
# use Java-1.8
# use Picard-Tools
# use Samtools
# use BEDTools
# use .gatk-4.1.8.1

# Define Variables
run_dir='' ## Current working directory to write all the files
gvcf_dir='' ## Directory to write GVCFs
cores='' ## Number of compute cores
hg19_ref='' ## Path to Human Genome (hg19) reference for removal of human reads
ref_path='' ## Path to reference genome (Anopheles) (.fasta)
dict='' ## Path to reference genome dictionary (.dict)
temp_dir='' ## Directory path to write temporary intermediate files
vmem='' ## Java Memory specification (could be kept as memory of the compute node)
sampleid='' ## Sample ID
bampath='' ## Path to raw BAM

# BAM to Fastq
path_to_fq="${temp_dir}/${sampleid}"
java -Xmx${vmem}G -jar $PICARD SamToFastq INPUT=${bampath} \
FASTQ=${path_to_fq}.1.fq \
SECOND_END_FASTQ=${path_to_fq}.2.fq \
VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=${temp_dir}

# Compress fastqs
gzip ${path_to_fq}.*.fq

# Choose relevant fastqs
fq1=${temp_dir}/${sampleid}.1.fq.gz
fq2=${temp_dir}/${sampleid}.2.fq.gz

# Removal of Human aligned reads
bowtie2 -x ${hg19_ref} -1 ${fq1} -2 ${fq2} | samtools view -bS -> ${temp_dir}/${sampleid}.mapAndUnmapped.hg19.bam
# filter out reads unmapped to hg19
samtools view -b -f 12 -F 256 ${temp_dir}/${sampleid}.mapAndUnmapped.hg19.bam > ${temp_dir}/${sampleid}.unmapped.bam
# sort BAM file by read name (-n) to have reads next to each other [required by bedtools]
samtools sort -n ${temp_dir}/${sampleid}.unmapped.bam -o ${temp_dir}/${sampleid}.unmapped.sorted.bam
# BAM to FASTQ files
use BEDTools
bedtools bamtofastq -i ${temp_dir}/${sampleid}.unmapped.sorted.bam -fq ${temp_dir}/${sampleid}.hostRemoved.1.fq -fq2 ${temp_dir}/${sampleid}.hostRemoved.2.fq

fq1=${temp_dir}/${sampleid}.hostRemoved.1.fq
fq2=${temp_dir}/${sampleid}.hostRemoved.2.fq


# Alignment to Reference Genome
echo "Alignment started"
bwa mem -t ${cores} \
-R "@RG\\tID:FLOWCELL_${sampleid}\\tSM:${sampleid}\\tPL:ILLUMINA\\tLB:LIB_${sampleid}" ${ref_path} \
${fq1} ${fq2} | samtools view -bS -> ${temp_dir}/${sampleid}.aligned.bam
echo "Alignment finished"

### Sort BAM
java -Xmx${vmem}G -jar $PICARD SortSam I=${temp_dir}/${sampleid}.aligned.bam \
O=${temp_dir}/${sampleid}.sorted.bam \
SO=coordinate \
TMP_DIR=${temp_dir}

### Mark Duplicates
java -Xmx${vmem}G -jar $PICARD MarkDuplicates I=${temp_dir}/${sampleid}.sorted.bam \
O=${temp_dir}/${sampleid}.marked_duplicates.bam \
M=${temp_dir}/${sampleid}.marked_duplicates.metrics \
TMP_DIR=${temp_dir}

### Re-order BAM
java -Xmx${vmem}G -jar $PICARD ReorderSam I=${temp_dir}/${sampleid}.marked_duplicates.bam \
O=${temp_dir}/${sampleid}.reordered.bam R=${ref_path} SD=${dict} \
TMP_DIR=${temp_dir}
samtools index ${temp_dir}/${sampleid}.reordered.bam

### GATK Variant calling
hapcall_input=${temp_dir}/${sampleid}.reordered.bam

# GATK HaplotypeCaller
gatk --java-options "-Xmx${vmem}G" HaplotypeCaller \
-R ${ref_path} --input ${hapcall_input} \
-ERC GVCF -ploidy 2 --interval-padding 100 --output ${gvcf_dir}/${sampleid}.g.vcf


#### Combining GVCFs into a single VCF
## Steps run on the entire batch of samples
ReferenceFile='' ## Path to reference genome (Anopheles) (.fasta)
gVCFList='' ## Path to a text file containing list of GVCFs to combine to a single VCF
IntermedFilePath='' ## Directory path to write temporary intermediate files
VCFPath='' ## Path to the final VCF output
BaseFileName='' ## Prefix filename to keep for the combined VCF output
vmem='' ## Java Memory specification (could be kept as memory of the compute node)

# Combine GVCfs (pre-requisite to joint call)
gatk --java-options "-Xmx${vmem}G" CombineGVCFs \
-R ${ReferenceFile} \
-V ${gVCFList} \
-O ${IntermedFilePath}${BaseFileName}_CombinedGVCFs.g.vcf

# Combined Genotyping
gatk --java-options "-Xmx${vmem}G" GenotypeGVCFs \
-R ${ReferenceFile} \
-V  ${IntermedFilePath}${BaseFileName}_CombinedGVCFs.g.vcf \
-O ${VCFPath}${BaseFileName}_CombinedGVCFs.vcf

# Hard filtering
# Hard filters are based on GATK best practices. See https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
# Isolate snps
gatk --java-options "-Xmx${vmem}G" SelectVariants \
-R ${ReferenceFile} \
-V ${VCFPath}${BaseFileName}_CombinedGVCFs.vcf \
-select-type SNP \
-O ${VCFPath}${BaseFileName}_CombinedGVCFs.snps.vcf

# Isolate indels
gatk --java-options "-Xmx${vmem}G" SelectVariants \
-R ${ReferenceFile} \
-V ${VCFPath}${BaseFileName}_CombinedGVCFs.vcf \
-select-type INDEL \
-O ${VCFPath}${BaseFileName}_CombinedGVCFs.indels.vcf

# Hard Filtering SNP
gatk --java-options "-Xmx${vmem}G" VariantFiltration \
-R ${ReferenceFile} \
-V ${VCFPath}${BaseFileName}_CombinedGVCFs.snps.vcf \
-filter "QD < 2.0"      --filter-name "QD2" \
-filter "QUAL < 30.0"   --filter-name "QUAL30" \
-filter "SOR > 3.0"     --filter-name "SOR3" \
-filter "FS > 60.0"     --filter-name "FS60" \
-filter "MQ < 40.0"     --filter-name "MQ40" \
-filter "MQRankSum < -12.5"     --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0"         --filter-name "ReadPosRankSum-8" \
-O ${VCFPath}${BaseFileName}.JointCall.filtered.snps.vcf

# Hard Filtering Indels
gatk --java-options "-Xmx${vmem}G" VariantFiltration \
-R ${ReferenceFile} \
-V ${VCFPath}${BaseFileName}_CombinedGVCFs.indels.vcf \
-filter "QD < 2.0"      --filter-name "QD2" \
-filter "QUAL < 30.0"   --filter-name "QUAL30" \
-filter "FS > 200.0"    --filter-name "FS200" \
-filter "ReadPosRankSum < -20.0"        --filter-name "ReadPosRankSum-20" \
-O ${VCFPath}${BaseFileName}.JointCall.filtered.indels.vcf

# Merge Hard Filtering VCFs
java -jar $PICARD MergeVcfs \
I=${VCFPath}${BaseFileName}.JointCall.filtered.snps.vcf \
I=${VCFPath}${BaseFileName}.JointCall.filtered.indels.vcf \
O=${VCFPath}${BaseFileName}.JointCall.filtered.combined.vcf
