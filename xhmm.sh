#!/bin/bash

#source activate namhee21

bam_DIR=/data/analysis/project/231211_WXS_CNV/xhmm_GIAB_1000G_test/bam
ref_DIR=/data/analysis/project/231211_WXS_CNV/xhmm_GIAB_1000G_test/reference
out_DIR=/data/analysis/project/231211_WXS_CNV/xhmm_GIAB_1000G_test/output

interval_list_to_pseq_reg=/data/home/namhee21/tools_nh/xhmm/sources/scripts/interval_list_to_pseq_reg


########################################
#### step 1. GATK depth of coverage ####
########################################

## 1) reference 

function picard_interval(){

docker run --rm -v /data:/data \
	broadinstitute/picard:2.27.5 \
	java -jar /usr/picard/picard.jar BedToIntervalList \
	I=${ref_DIR}/Agilent_V5.bed \
	O=${ref_DIR}/Agilent_V5.interval_list \
	SD=/data/analysis/project/231211_WXS_CNV/public_data/reference/hs37d5.dict
}

picard_interval


## 2) gatk Depth of coverate

function gatk_run(){

## gatk 3.8 ver
java -jar /data/home/namhee21/tools_nh/GenomeAnalysisTK-3.8.1/GenomeAnalysisTK.jar \
	-T DepthOfCoverage \
	-R /data/analysis/project/231211_WXS_CNV/public_data/reference/hs37d5.fa \
	-L ${ref_DIR}/Agilent_V5.bed \
	-I /data/analysis/project/231211_WXS_CNV/public_data/BAM/bam.list \
	-dt BY_SAMPLE -dcov 5000 -l INFO \
	--omitLocusTable \
	--minBaseQuality 0 --minMappingQuality 20 \
	--start 1 --stop 5000 --nBins 200 \
	--includeRefNSites \
	--countType COUNT_FRAGMENTS \
	-o ${out_DIR}/GIAB_1000g_DepthOfCoverage

#-dt BY_SAMPLE -dcov 5000 -l INFO --omitDepthOutputAtEachBase --omitLocusTable \

}

gatk_run

#########################################################
#### step 2. xhmm -> combines GATK depth of coverage ####
#########################################################

function xhmm_RD(){

grep -v "NaN" ${out_DIR}/GIAB_1000g_DepthOfCoverage.sample_interval_summary \
> ${out_DIR}/GIAB_1000g_DepthOfCoverage.ext_NaN.sample_interval_summary

xhmm --mergeGATKdepths \
-o ${out_DIR}/DATA.RD.txt \
--GATKdepths ${out_DIR}/GIAB_1000g_DepthOfCoverage.ext_NaN.sample_interval_summary

}

xhmm_RD


##################################################################
#### step 3. calculate GC content -> GATK GCContentByInterval ####
##################################################################

function gatk_calc_GC(){

## 1) GC contents

java -jar /data/home/namhee21/tools_nh/GenomeAnalysisTK-3.8.1/GenomeAnalysisTK.jar \
	-T GCContentByInterval \
	-L ${ref_DIR}/Agilent_V5.bed \
	-R /data/analysis/project/231211_WXS_CNV/public_data/reference/hs37d5.fa \
	-o ${out_DIR}/total_sample.locus_GC.txt


awk '{if($2 < 0.1 || $2 > 0.9) print $1}' ${out_DIR}/total_sample.locus_GC.txt \
> ${out_DIR}/total_sample.locus_GC.filter.txt

}

gatk_calc_GC


#####################################################################
#### step 4. SEQDB - reference database  (plinkseq installation) ####
#####################################################################


function gen_seqdb(){

## 1) calculate the fraction of repeat-masked bases in each target

${interval_list_to_pseq_reg} \
${ref_DIR}/Agilent_V5.interval_list \
> ${out_DIR}/EXOME.targets.reg

## 2) locdb loc-load 
pseq . loc-load --locdb ${out_DIR}/EXOME.targets.LOCDB \
--file ${out_DIR}/EXOME.targets.reg \
--group targets \
--out ${out_DIR}/EXOME.targets.LOCDB.loc-load \
--noweb to skip

## 3) SeqDB loc-stats
pseq . loc-stats --locdb ${out_DIR}/EXOME.targets.LOCDB --group targets --seqdb ${out_DIR}/seqdb |\
awk '{if (NR > 1) print $_}' | sort -k1 -g | awk '{print $10}' | \
paste ${ref_DIR}/Agilent_V5.bed - | awk '{print $1"\t"$2}' \
> ${out_DIR}/DATA.locus_complexity.txt

cat ${out_DIR}/DATA.locus_complexity.txt | awk '{if ($2 > 0.25) print $1}' \
> ${out_DIR}/low_complexity_targets.txt

}

gen_seqdb

###############################################################################
#### step 5. Filters samples and targets and then mean-centers the targets ####
###############################################################################

function xhmm_matrix(){

xhmm --matrix \
	-r ${out_DIR}/DATA.RD.txt \
	--centerData \
	--centerType target \
	-o ${out_DIR}/DATA.filtered_centered.RD.txt \
	--outputExcludedTargets ${out_DIR}/DATA.filtered_centered.RD.txt.filtered_targets.txt \
	--outputExcludedSamples ${out_DIR}/DATA.filtered_centered.RD.txt.filtered_samples.txt \
	--excludeTargets ${out_DIR}/total_sample.locus_GC.filter.txt \
	--minTargetSize 10 --maxTargetSize 10000 --minMeanTargetRD 10 --maxMeanTargetRD 500 \
	--minMeanSampleRD 25 --maxMeanSampleRD 200 --maxSdSampleRD 150
}

xhmm_matrix

################################################
#### step 6. Runs PCA on mean-centered data ####
################################################

function xhmm_pca(){

xhmm --PCA -r ${out_DIR}/DATA.filtered_centered.RD.txt --PCAfiles ${out_DIR}/DATA.RD_PCA

}

xhmm_pca

####################################################################
#### step 7. Normalize mean-centered data using PCA information ####
####################################################################

function xhmm_norm(){

xhmm --normalize -r ${out_DIR}/DATA.filtered_centered.RD.txt \
	--PCAfiles ${out_DIR}/DATA.RD_PCA \
	--normalizeOutput ${out_DIR}/DATA.PCA_normalized.txt \
	--PCnormalizeMethod PVE_mean \
	--PVE_mean_factor 0.7
}

xhmm_norm

#################################################################################
#### step 8. Filters and z-score centers (by sample) the PCA-normalized data ####
#################################################################################

function xhmm_z_score(){

xhmm --matrix -r ${out_DIR}/DATA.PCA_normalized.txt \
	--centerData --centerType sample --zScoreData \
	-o ${out_DIR}/DATA.PCA_normalized.filtered.sample_zscores.RD.txt \
	--outputExcludedTargets ${out_DIR}/DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
	--outputExcludedSamples ${out_DIR}/DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
	--maxSdTargetRD 30
}


xhmm_z_score

##############################################################################################
#### step 9. Filters original read-depth data to be the same as filtered, normalized data ####
##############################################################################################

function xhmm_filter_data(){

xhmm --matrix -r ${out_DIR}/DATA.RD.txt \
	--excludeTargets ${out_DIR}/DATA.filtered_centered.RD.txt.filtered_targets.txt \
	--excludeTargets ${out_DIR}/DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
	--excludeSamples ${out_DIR}/DATA.filtered_centered.RD.txt.filtered_samples.txt \
	--excludeSamples ${out_DIR}/DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
	-o ${out_DIR}/DATA.sample_filtered.RD.txt

}

xhmm_filter_data

#################################################################
#### step 10. Genotypes discovered CNVs after normalize data ####
#################################################################


function discover_CNVs(){

### 1) Discovers CNVs in normalized data

xhmm --discover -p /data/home/namhee21/tools_nh/xhmm/params.txt \
	-r ${out_DIR}/DATA.PCA_normalized.filtered.sample_zscores.RD.txt \
	-R ${out_DIR}/DATA.sample_filtered.RD.txt \
	-c ${out_DIR}/DATA.xcnv \
	-a ${out_DIR}/DATA.aux_xcnv \
	-s ${out_DIR}/DATA


### 2) Genotypes discovered CNVs in all samples

xhmm --genotype -p /data/home/namhee21/tools_nh/xhmm/params.txt \
	-r ${out_DIR}/DATA.PCA_normalized.filtered.sample_zscores.RD.txt \
	-R ${out_DIR}/DATA.sample_filtered.RD.txt \
	-g ${out_DIR}/DATA.xcnv \
	-F /data/analysis/project/231211_WXS_CNV/public_data/reference/hs37d5.fa \
	-v ${out_DIR}/DATA.vcf

}

discover_CNVs

