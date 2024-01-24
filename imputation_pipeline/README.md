# Script for analyzing high- and low-coverage Data

The genotyping_bcftools.sh contains script for calling genotypes from high-coverage BAM files with bcftools.

The downsample.sh contains script for downsampling high-coverage data to low coverage and then impute the low-coverage data using GLIMPSE (starting from line 71).  The same imputation pipeline is used for all the low-coverage data.

The bash script processbychr.sh contains commands used for calling IBD using ancIBD. It runs on each of tue 20 autosomes separately.

