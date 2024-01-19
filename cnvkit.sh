#!/bin/bash

DIR=/data/analysis/project/231211_WXS_CNV/GIAB_1000G_test/bam
out_DIR=/data/analysis/project/231211_WXS_CNV/GIAB_1000G_test/cnvkit_output
ref_DIR=/data/analysis/project/231211_WXS_CNV/GIAB_1000G_test/reference
#ref_DIR=/data/analysis/project/231211_WXS_CNV/public_data/reference

function cnvkit_step1()
{
cnvkit.py batch \
    --normal ${DIR}/HG002.bam \
             ${DIR}/HG003.bam \
             ${DIR}/HG004.bam \
    --targets ${ref_DIR}/Agilent_V5.bed \
    --antitargets ${ref_DIR}/Agilent_V5.antitarget.bed \
    --fasta /data/analysis/project/231211_WXS_CNV/public_data/reference/hs37d5.fa \
    --access ${out_DIR}/access.hs37d5.bed \
    --output-reference ${out_DIR}/reference.cnn\
    --output-dir ${out_DIR}/ \
    --diagram --scatter

#--annotate /home/gordeeva/tools/cnvkit/refFlat.txt \

}



function cnvkit_step2()
{

reference_cnn=${out_DIR}/reference.cnn

for antitarget_cnn in `ls -1 ${out_DIR}/* | grep "antitargetcoverage.cnn"`
do

target_cnn=`echo ${antitarget_cnn} | sed 's/antitargetcoverage.cnn/targetcoverage.cnn/'`
sample=`echo ${antitarget_cnn} | awk -F"." '{print $1}'`

echo ${target_cnn}
echo ${antitarget_cnn}
echo ${sample}

## normalization
cnvkit.py fix ${target_cnn} ${antitarget_cnn} ${reference_cnn} -o ${sample}.cnr

## segment
cnvkit.py segment ${sample}.cnr -m cbs -o ${sample}.cns

## cnv_call
cnvkit.py call ${sample}.cns -y -m threshold -t=-1.1,-0.4,0.3,0.7 -o ${sample}.call.cns

## scatter plot
cnvkit.py scatter ${sample}.cnr --y-min -2.0 --y-max 2.0 -s ${sample}.call.cns -o ${sample}.scatter.png

done
}

#cnvkit_step1
cnvkit_step2

