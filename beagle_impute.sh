#!/bin/bash

#################################################
# File: beagle_impute.sh
#################################################
# Author: Rajiv McCoy
# This file provides the commands used for
# genotype imputation for subset of patients
# of European ancestry using the 1000 Genomes
# reference panel.
#################################################

for i in {1..4}
do
    echo "Running imputation on chromosome "
    echo $i
    java -Xmx8g -jar b4.r1230.jar gt=conform_chr$i.vcf.gz ref=chr$i.1kg.ref.phase1_release_v3.20101123.vcf.gz excludesamples=non_eur.excl chrom=$i out=imputed_chr$i nthreads=5 &
done
wait

for i in {5..8}
do
    echo "Running imputation on chromosome "
    echo $i
    java -Xmx8g -jar b4.r1230.jar gt=conform_chr$i.vcf.gz ref=chr$i.1kg.ref.phase1_release_v3.20101123.vcf.gz excludesamples=non_eur.excl chrom=$i out=imputed_chr$i nthreads=5 &
done
wait

for i in {9..12}
do
    echo "Running imputation on chromosome "
    echo $i
    java -Xmx8g -jar b4.r1230.jar gt=conform_chr$i.vcf.gz ref=chr$i.1kg.ref.phase1_release_v3.20101123.vcf.gz excludesamples=non_eur.excl chrom=$i out=imputed_chr$i nthreads=5 &
done
wait

for i in {13..16}
do
    echo "Running imputation on chromosome "
    echo $i
    java -Xmx8g -jar b4.r1230.jar gt=conform_chr$i.vcf.gz ref=chr$i.1kg.ref.phase1_release_v3.20101123.vcf.gz excludesamples=non_eur.excl chrom=$i out=imputed_chr$i nthreads=5 &
done
wait

for i in {17..20}
do
    echo "Running imputation on chromosome "
    echo $i
    java -Xmx8g -jar b4.r1230.jar gt=conform_chr$i.vcf.gz ref=chr$i.1kg.ref.phase1_release_v3.20101123.vcf.gz excludesamples=non_eur.excl chrom=$i out=imputed_chr$i nthreads=5 &
done
wait

for i in {21..23}
do
    echo "Running imputation on chromosome "
    echo $i
    java -Xmx8g -jar b4.r1230.jar gt=conform_chr$i.vcf.gz ref=chr$i.1kg.ref.phase1_release_v3.20101123.vcf.gz excludesamples=non_eur.excl chrom=$i out=imputed_chr$i nthreads=5 &
done
wait
