#!/bin/bash

cd ./TTS

echo 'start converting bam to bed pseudomonas'
bedtools bamtobed -i trimmed_0_pseudo.sorted.bam > pseudo_0.bed
bedtools bamtobed -i trimmed_5_pseudo.sorted.bam > pseudo_5.bed
bedtools bamtobed -i trimmed_10_pseudo.sorted.bam > pseudo_10.bed
bedtools bamtobed -i trimmed_15_pseudo.sorted.bam > pseudo_15.bed

echo 'start converting bam to bed LUZ7'
bedtools bamtobed -i trimmed_0.sorted.bam > vir_0.bed
bedtools bamtobed -i trimmed_5.sorted.bam > vir_5.bed
bedtools bamtobed -i trimmed_10.sorted.bam > vir_10.bed
betoold bamtobed -i trimmed_15.sorted.bam > vir_15.bed

echo 'start python count pseudomonas'
python2 TSS_analysis.py count --input pseudo_0.bed --output pseudo_0.output
python2 TSS_analysis.py count --input pseudo_5.bed --output pseudo_5.output
python2 TSS_analysis.py count --input pseudo_10.bed --output pseudo_10.output
python2 TSS_analysis.py count --input pseudo_15.bed --output pseudo_15.output

echo 'start python count LUZ7'
python2 TSS_analysis.py count --input vir_0.bed --output vir_0.output
python2 TSS_analysis.py count --input vir_5.bed --output vir_5.output
python2 TSS_analysis.py count --input vir_10.bed --output vir_10.output
python2 TSS_analysis.py count --input vir_15.bed --output vir_15.output

echo 'GIVE CLUSTERING VARIABLE (1-10)'
echo 'This determines which ends are clustered together, default is 5'
read variableClustering

python2 TSS_analysis.py cluster --input pseudo_0.output --output clustered_pseudo_0.output --control cluster_info_pseudo_0.txt --cutoff $variableClustering
python2 TSS_analysis.py cluster --input pseudo_5.output --output clustered_pseudo_5.output --control cluster_info_pseudo_5.txt --cutoff $variableClustering
python2 TSS_analysis.py cluster --input pseudo_10.output --output clustered_pseudo_10.output --control cluster_info_pseudo_10.txt --cutoff $variableClustering
python2 TSS_analysis.py cluster --input pseudo_15.output --output clustered_pseudo_15.output --control cluster_info_pseudo_15.txt --cutoff $variableClustering

python2 TSS_analysis.py cluster --input vir_0.output --output clustered_vir_0.output --control cluster_info_vir_0.txt --cutoff $variableClustering
python2 TSS_analysis.py cluster --input vir_5.output --output clustered_vir_5.output --control cluster_info_vir_5.txt --cutoff $variableClustering
python2 TSS_analysis.py cluster --input vir_10.output --output clustered_vir_10.output --control cluster_info_vir_10.txt --cutoff $variableClustering
python2 TSS_analysis.py cluster --input vir_15.output --output clustered_vir_15.output --control cluster_info_vir_15.txt --cutoff $variableClustering

echo 'STATISTICS R pseudomonas'
R binomialtest.R --no-save clustered_pseudo_0.output stats_pseudo_0.bed 0.2 < binomialtest.R
R binomialtest.R --no-save clustered_pseudo_5.output stats_pseudo_5.bed 0.2 < binomialtest.R
R binomialtest.R --no-save clustered_pseudo_10.output stats_pseudo_10.bed 0.2 < binomialtest.R
R binomialtest.R --no-save clustered_pseudo_15.output stats_pseudo_15.bed 0.2 < binomialtest.R

echo 'STATISTICS R LUZ7'
R binomialtest.R --no-save clustered_vir_0.output stats_vir_0.bed 0.2 < binomialtest.R
R binomialtest.R --no-save clustered_vir_5.output stats_vir_5.bed 0.2 < binomialtest.R
R binomialtest.R --no-save clustered_vir_10.output stats_vir_10.bed 0.2 < binomialtest.R
R binomialtest.R --no-save clustered_vir_15.output stats_vir_15.bed 0.2 < binomialtest.R
