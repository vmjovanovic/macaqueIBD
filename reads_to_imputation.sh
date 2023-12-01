#!/bin/bash
#The raw reads should be in the subdirectory RAW in the project home

#configuration
nthreads=25
projecthomepath=/PATH/TO/PROJECT/HOME
picardpath=/PATH/TO/PICARD/JAR
genome=/PATH/TO/GENOME/FASTA
refpanelpath=/PATH/TO/REFERENCE/PANEL
mappath=/PATH/TO/REFERENCE/RECOMBINATION/MAP

cd $projecthomepath

#setting up the directories
mkdir LOGS
mkdir TRIMMED
mkdir CLEAN_READS
mkdir GENOME
mkdir MAPPING
mkdir DEDUPBAMS


#adapter and quality trimming
cd RAW
for FILE in *R1_001.fastq.gz; do 
	NAME=$(basename $FILE _R1_001.fastq.gz)
	cutadapt -a file:adapters1.fa -A file:adapters2.fa -j $nthreads -m20 -q20,20 -o ../TRIMMED/$NAME\_R1_001.trimmed.fastq.gz -p ../TRIMMED/$NAME\_R2_001.trimmed.fastq.gz $NAME\_R1_001.fastq.gz $NAME\_R2_001.fastq.gz
done &>> ../LOGS/trimming.log


# combining the runs from the same sample

cd ../TRIMMED/
ls *.fastq.gz | cut -d"_" -f2 | sort -u > list_of_samples.txt
while IFS= read -r LINE; do 
	echo $LINE
	cat *$LINE*R1*.gz > ../CLEAN_READS/$LINE.R1.fastq.gz; cat *$LINE*R2*.gz > ../CLEAN_READS/$LINE.R2.fastq.gz
done < list_of_samples.txt

# indexing the genome 
cd ../GENOME/
hisat2-build -p $nthreads $genome Macgenome

# mapping
cd ../CLEAN_READS/
for FILE in *R1.fastq.gz*; do 
	NAME=$(basename $FILE .R1.fastq.gz)
	echo $NAME &>>../LOGS/mapping.log
	hisat2 -p $nthreads -x ../GENOME/Macgenome -1 $NAME.R1.fastq.gz -2 $NAME.R2.fastq.gz -S ../MAPPING/$NAME.sam &>>../LOGS/mapping.log
done

# deduplication

cd ../MAPPING/
for FILE in *.sam; do 
	NAME=$(basename $FILE .sam)
	samtools sort -o $NAME.bam -@ $nthreads $FILE
	samtools index $NAME.bam
	java -Xmx50G -jar $picardpath/picard.jar MarkDuplicates I=$NAME.bam M=../LOGS/$NAME.deduplication_metrics.txt REMOVE_DUPLICATES=true O=../DEDUPBAMS/$NAME.RG.bam
done

#indexing the final bam files, doing imputation

cd .../DEDUPBAMS/
for FILE in *.RG.bam; do 
	id=$(basename $BAM .RG.bam)
	echo doing$id
	samtools index $FILE
	
	for ch in {1..20}; do
		echo processing chr$ch
		#prepare GL file (with bcftools) for GLIMPSE
		VCF=$refpanelpath/mGAP.chr$ch.sites.vcf.gz
		TSV=$refpanelpath/mGAP.chr$ch.sites.tsv.gz
		VCF_OUT=$id.chr$ch.vcf.gz
		bcftools mpileup -f ${genome} --ignore-RG -I -E -a 'FORMAT/DP,FORMAT/AD' -T ${VCF} -r $ch -q 30 -Q 30 ${BAM} -Ou | bcftools call -Aim -C alleles -T ${TSV} -Oz -o ${VCF_OUT}
		bcftools index -f ${VCF_OUT}

		# running GLIMPSE

		GLIMPSE_chunk --input ${VCF_OUT} --region $ch --window-size 2000000 --buffer-size 200000 --output chunks.chr$ch.txt
		REF=$refpanelpath/mGAP.chr$ch.bcf
		MAP=$mappath/Mmul10_GeneticMap_Chr$ch.gmap.gz

		if [ ! -d GLIMPSE_imputed ]; then
			mkdir GLIMPSE_imputed
		fi

		while IFS="" read -r LINE || [ -n "$LINE" ]; do
			printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
			IRG=$(echo $LINE | cut -d" " -f3)
			ORG=$(echo $LINE | cut -d" " -f4)
			OUT=GLIMPSE_imputed/$id.chr$ch.${ID}.bcf
			GLIMPSE_phase --input ${VCF_OUT} --reference ${REF} --map ${MAP} --input-region ${IRG} --output-region ${ORG} --output ${OUT}
			bcftools index -f ${OUT}
		done < chunks.chr$ch.txt

		if [ ! -d GLIMPSE_ligated ]; then
			mkdir GLIMPSE_ligated
		fi
		
		LST=GLIMPSE_ligated/list.chr$ch.txt
		ls GLIMPSE_imputed/$id.chr$ch.*.bcf > ${LST}
		OUT=GLIMPSE_ligated/$id.chr$ch.merged.bcf
		GLIMPSE_ligate --input ${LST} --output $OUT
		bcftools index -f ${OUT}
		rm -r ./GLIMPSE_imputed

		if [ ! -d GLIMPSE_phased ]; then
			mkdir GLIMPSE_phased
		fi

		VCF=GLIMPSE_ligated/$id.chr$ch.merged.bcf
		OUT=GLIMPSE_phased/$id.chr$ch.phased.bcf
		GLIMPSE_sample --input ${VCF} --solve --output ${OUT}
		bcftools index -f ${OUT}

		# merge the phasing and imputation file together
		bcftools annotate -a ./GLIMPSE_ligated/$id.chr$ch.merged.bcf -c FORMAT/GP ./GLIMPSE_phased/$id.chr$ch.phased.bcf -Oz -o $id.chr$ch.merged.vcf.gz
		bcftools index $id.chr$ch.merged.vcf.gz
	done

	rm chunks.chr*.txt
	bcftools concat -Oz -o $id.phased.vcf.gz $id.chr{1..20}.merged.vcf.gz
	bcftools index $id.phased.vcf.gz
	rm $id.chr*.vcf.gz
	rm $id.chr*.vcf.gz.csi
	rm -r GLIMPSE_ligated
	rm -r GLIMPSE_phased	
done

cd $projecthomepath
rm -r TRIMMED
rm -r CLEAN_READS
rm -r MAPPING

