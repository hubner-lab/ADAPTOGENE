#!/bin/bash
#########################
VCF=$1 ; VCF_name=$(basename $(basename $VCF .vcf.gz) .vcf) # in VCF or VCFGZ format
SAMPLES=$2
maf=$3
miss=$4
winsize=$5 #100
stepsize=$6 #20
r2threshold=$7 #0.1
#########################
VCF_f=${VCF_name}_MAF${maf}_MISS${miss}.vcf
VCF_f_name=$(basename $VCF_f .vcf)

# Filter MAF MISS
plink --vcf $VCF \
--const-fid \
--allow-extra-chr \
--set-missing-var-ids @:# \
--keep $SAMPLES \
--maf $maf  \
--geno $miss \
--recode vcf \
--out ${VCF_f_name}

# Calculate LD
plink --vcf $VCF_f \
--const-fid \
--allow-extra-chr \
--indep-pairwise $winsize $stepsize $r2threshold \
--out ${VCF_f_name}.LD${r2threshold} \
--set-missing-var-ids @:#

# Make PCA
plink --vcf $VCF_f \
--const-fid \
--allow-extra-chr \
--set-missing-var-ids @:# \
--extract ${VCF_f_name}.LD${r2threshold}.prune.in \
--make-bed \
--pca \
--recode vcf \
--out ${VCF_f_name}.LD${r2threshold}

VCF_f_LD="${VCF_f_name}.LD${r2threshold}.vcf"

echo $VCF_f
echo $VCF_f_LD
# Edit OUT vcf - remove FAMID
sed -i '/^#CHROM/s/\t0_/\t/g' ${VCF_f}
sed -i '/^#CHROM/s/\t0_0_/\t/g' ${VCF_f_LD}

echo 'INFO: Process vcf completed successfully'
