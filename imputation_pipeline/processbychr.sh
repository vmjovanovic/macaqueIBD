#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N vlad_glimpse #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -q archgen.q #queue
# -m e #send an email at the end of the job
# -M yilei_huang@eva.mpg.de #send email to this address
#$ -l h_vmem=45G #request 4Gb of memory
#$ -V # load personal profile
#$ -o $JOB_NAME.o.$JOB_ID.$TASK_ID
#$ -t 1:20:1

ch=$SGE_TASK_ID

# run the bcftools view command if the output file ./GLIMPSE/ch$ch.vcf.gz does not exist
if [ ! -f ./GLIMPSE/ch$ch.vcf.gz ]; then
    bcftools view -r $ch -m2 -M2 -v snps -Oz -o ./GLIMPSE/ch$ch.vcf.gz ./GLIMPSE/Vlad_glimpse-2023-9-15.merged.MAF5.vcf.gz
fi

ancIBD-run --vcf ./GLIMPSE/ch$ch.vcf.gz --ch $ch --af_column variants/RAF --IBD2 --ibd-in 1e-5 --ibd-out 1e-5 --post2 0.8 \
    --out test --map_path /mnt/archgen/users/yilei/Data/rhemac/mGAP_for_imputation/reference_panel/maps_new/MinMyc.snp \
    --min 4  