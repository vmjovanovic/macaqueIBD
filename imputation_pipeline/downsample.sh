#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N ds #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -q archgen.q #queue
# -m e #send an email at the end of the job
# -M yilei_huang@eva.mpg.de #send email to this address
#$ -l h_vmem=45G #request 4Gb of memory
#$ -V # load personal profile
#$ -o $JOB_NAME.o.$JOB_ID.$TASK_ID
#$ -t 1:500:1
#$ -tc 150 #max number of concurrent tasks
#$ -l 'h=!bionode0[1-9]'


# 25 replicates for each sample (5 in total) and coverage (4 types)

######### first, get sample ID and path to BAM file #########
sampleIDs=(134814 134834 134854 134858 134860)
sampleCoverage=(20.41 23.36 29.60 18.187 23.73)
id=$SGE_TASK_ID
sample_index=$((($id-1)/150))
iid=${sampleIDs[$sample_index]}
path2bam=/mnt/archgen/users/yilei/Data/rhemac/downsample4pub/BAM/$iid.RG.bam
original_coverage=${sampleCoverage[$sample_index]}
echo processing sample $iid, getting BAM file from $path2bam, original coverage $original_coverage


######### second, get target coverage of downsampling #########
covstrs=(cov1over4 cov1over2 cov1 cov2 cov4 cov6)
covfloats=(0.25 0.5 1.0 2.0 4.0 6.0)
# set id to SGE_TASK_ID's remainder when divided by 150
# this will give a number between 0 and 149
# divide by 25 to get the index of the coverage
# this will give a number between 0 and 5
# use this number to get the coverage string and float
cov_index=$((($id-1)%150/25))

covstr=${covstrs[$cov_index]}
covfloat=${covfloats[$cov_index]}
frac=$(echo "scale=6; $covfloat/$original_coverage" | bc)
echo downsampling to coverage $covstr, $covfloat, subsampling fraction $frac
# now get the batch number
# set id to SGE_TASK_ID's remainder when divided by 25
# this will give a number between 0 and 24
# add 1 to get the batch number
b=$(($id%25+1))
echo batch $b

if [ ! -d $covstr ]; then
    mkdir $covstr
fi
cd ./$covstr

if [ ! -d batch$b ]; then
    mkdir batch$b
fi
cd batch$b

if [ ! -d $iid ]; then
    mkdir $iid
fi
cd $iid

pwd

# actual work start here
samtools view --subsample-seed $b --subsample $frac -b -o $iid.batch$b.bam $path2bam
samtools index $iid.batch$b.bam
for ch in {1..20}; 
do
    echo processing chr$ch
    #prepare GL file (with bcftools) for GLIMPSE
    BAM=$iid.batch$b.bam
    VCF=/mnt/archgen/users/yilei/Data/rhemac/mGAP_for_imputation/reference_panel/mGAP.chr$ch.sites.vcf.gz
    TSV=/mnt/archgen/users/yilei/Data/rhemac/mGAP_for_imputation/reference_panel/mGAP.chr$ch.sites.tsv.gz
    REFGEN=/mnt/archgen/users/afreudiger/rhemac/test/Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz
    OUT=$iid.batch$b.chr$ch.vcf.gz

    bcftools mpileup -f ${REFGEN} --ignore-RG -I -E -a 'FORMAT/DP' -T ${VCF} -r $ch -q 30 -Q 30 ${BAM} -Ou | bcftools call -Aim -C alleles -T ${TSV} -Oz -o ${OUT}
    bcftools index -f ${OUT}

    # running GLIMPSE
    path2glimpse="/mnt/archgen/users/yilei/bin/GLIMPSE-glimpse1"
    VCF=$iid.batch$b.chr$ch.vcf.gz
    $path2glimpse/chunk/bin/GLIMPSE_chunk --input $VCF --region $ch --window-size 2000000 --buffer-size 200000 --output chunks.chr$ch.txt
    REF=/mnt/archgen/users/yilei/Data/rhemac/mGAP_for_imputation/reference_panel/mGAP.chr$ch.bcf
    MAP=/mnt/archgen/users/yilei/Data/rhemac/mGAP_for_imputation/reference_panel/maps_new/chr$ch.cleanup.gmap.gz

    if [ ! -d GLIMPSE_imputed ]; then
        mkdir GLIMPSE_imputed
    fi

    while IFS="" read -r LINE || [ -n "$LINE" ];
    do
        printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
        IRG=$(echo $LINE | cut -d" " -f3)
        ORG=$(echo $LINE | cut -d" " -f4)
        OUT=GLIMPSE_imputed/$iid.batch$b.chr$ch.${ID}.bcf
        $path2glimpse/phase/bin/GLIMPSE_phase --input ${VCF} --reference ${REF} --map ${MAP} --input-region ${IRG} --output-region ${ORG} --output ${OUT}
        bcftools index -f ${OUT}
    done < chunks.chr$ch.txt

    if [ ! -d GLIMPSE_ligated ]; then
        mkdir GLIMPSE_ligated
    fi
    LST=GLIMPSE_ligated/list.chr$ch.txt
    ls GLIMPSE_imputed/$iid.batch$b.chr$ch.*.bcf > ${LST}
    OUT=GLIMPSE_ligated/$iid.batch$b.chr$ch.merged.bcf
    $path2glimpse/ligate/bin/GLIMPSE_ligate --input ${LST} --output $OUT
    bcftools index -f ${OUT}
    rm -r ./GLIMPSE_imputed

    if [ ! -d GLIMPSE_phased ]; then
        mkdir GLIMPSE_phased
    fi
    VCF=GLIMPSE_ligated/$iid.batch$b.chr$ch.merged.bcf
    OUT=GLIMPSE_phased/$iid.batch$b.chr$ch.phased.bcf
    $path2glimpse/sample/bin/GLIMPSE_sample --input ${VCF} --solve --output ${OUT}
    bcftools index -f ${OUT}

    # merge the phasing and imputation file together
    bcftools annotate -a ./GLIMPSE_ligated/$iid.batch$b.chr$ch.merged.bcf -c FORMAT/GP,FORMAT/DS ./GLIMPSE_phased/$iid.batch$b.chr$ch.phased.bcf -Oz -o $iid.chr$ch.merged.vcf.gz
    bcftools index $iid.chr$ch.merged.vcf.gz
done
rm chunks.chr*.txt
bcftools concat -Oz -o $iid.batch$b.phased.vcf.gz $iid.chr{1..20}.merged.vcf.gz
bcftools index $iid.batch$b.phased.vcf.gz
rm $iid.batch$b.chr*.vcf.gz
rm $iid.batch$b.chr*.vcf.gz.csi
rm $iid.chr{1..20}.merged.vcf.gz
rm $iid.chr{1..20}.merged.vcf.gz.csi
rm -r GLIMPSE_ligated
rm -r GLIMPSE_phased
# delete BAM file because they are too large, takes too much space
# I used seed when downsampling, so it's easy to reproduce
rm $iid.batch$b.bam
rm $iid.batch$b.bam.bai
