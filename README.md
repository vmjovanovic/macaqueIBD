# macaqueIBD: a pipeline to call IBD segments in the imputed genomes of rhesus macaques

# Brief Description
This pipeline was used for calling IBD in the rhesus macaque dataset, containing the genomic sequences from 103 individuals from Cayo Santiago island population. The workflow starts with quality control of the raw data, processing and mapping the reads to the reference genome (Mmul10), proceeds with the imputation step (for the samples sequenced to lower depth) and ends with calling the IBD segments in the final recovered genotypes.

# Necessary dependencies and softwares
<li>cutadapt</li>
<li>hisat2</li>
<li>samtools</li>
<li>picard MarkDuplicates</li>
<li>bcftools</li>
<li>GLIMPSE</li>
<li>plink2</li>
<li>ancIBD</li>
