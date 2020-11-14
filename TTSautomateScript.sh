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

echo 'STATISTICS R pseudomonas'
R binomialtest.R --no-save pseudo_0.output stats_pseudo_0.bed 0.2 < binomialtest.R
R binomialtest.R --no-save pseudo_5.output stats_pseudo_5.bed 0.2 < binomialtest.R
R binomialtest.R --no-save pseudo_10.output stats_pseudo_10.bed 0.2 < binomialtest.R
R binomialtest.R --no-save pseudo_15.output stats_pseudo_15.bed 0.2 < binomialtest.R

echo 'STATISTICS R LUZ7'
R binomialtest.R --no-save vir_0.output stats_vir_0.bed 0.2 < binomialtest.R
R binomialtest.R --no-save vir_5.output stats_vir_5.bed 0.2 < binomialtest.R
R binomialtest.R --no-save vir_10.output stats_vir_10.bed 0.2 < binomialtest.R
R binomialtest.R --no-save vir_15.output stats_vir_15.bed 0.2 < binomialtest.R
