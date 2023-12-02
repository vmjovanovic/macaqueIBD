#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N macaque_genotyping #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -q archgen.q #queue
# -m e #send an email at the end of the job
# -M annika_freudiger@eva.mpg.de #send email to this address
# -pe smp 4 #needs 8 CPU cores
#$ -l h_vmem=50G #request 4Gb of memory
#$ -V # load personal profile
#$ -o $JOB_NAME.o.$JOB_ID.$TASK_ID
#$ -t 1:5:1
#$ -l 'h=!bionode0[1-9]'


iids=(134814 134834 134854 134858 134860)
iid=${iids[$SGE_TASK_ID-1]}
echo genotyping $iid with bcftools

BAM=/mnt/archgen/users/yilei/Data/rhemac/downsample4pub/BAM/$iid.RG.bam
VCF=/mnt/archgen/users/yilei/Data/rhemac/mGAP_for_imputation/reference_panel/mGAP.allchr.sites.vcf.gz
TSV=/mnt/archgen/users/yilei/Data/rhemac/mGAP_for_imputation/reference_panel/mGAP.allchr.sites.tsv.gz
REFGEN=/mnt/archgen/users/afreudiger/rhemac/test/Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz

# test if folder $iid exists, if not, make one
if [ ! -d $iid ]; then
    mkdir $iid
fi
OUT=./$iid/$iid.vcf.gz

bcftools mpileup -f ${REFGEN} -I -E -a 'FORMAT/AD,FORMAT/DP' -Q 30 -q 30 -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20 -T ${VCF} ${BAM} -Ou | bcftools call -Aim -C alleles -T ${TSV} -Oz -o ${OUT}
bcftools index -f ${OUT}


# ############ do some basic filtering ###########
# # remove sites with missing data, with less than 15 reads, and with more than 100 reads

bcftools view -Oz -o ./$iid/$iid.filtered.vcf.gz -g ^miss -i 'FORMAT/DP>=15&&FORMAT/DP<=100&&MEDIAN(FORMAT/PL)>=30' $OUT
bcftools index ./$iid/$iid.filtered.vcf.gz