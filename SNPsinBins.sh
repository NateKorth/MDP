#!/bin/sh
#SBATCH --partition=benson
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=120:00:00          # Run time in hh:mm:ss
#SBATCH --mem=10000       # --mem-per-cpu Maximum memory required per CPU (in megabytes)
#SBATCH --job-name=genKinship
#SBATCH --error=./bins.err
#SBATCH --output=./bins.out

ml anaconda 
conda activate mypython
#python -m schnablelab.SNPcalling.base SubsamplingSNPs MADP_geno_imputed.vcf SNPsSetV2/MADP_BIN1_SNPs.csv 
#mv MADP_geno_imputed_subSNPs.vcf MADP_BIN1_SNPs.vcf

#python -m schnablelab.SNPcalling.base SubsamplingSNPs MADP_geno_imputed.vcf SNPsSetV2/MADP_BIN2_SNPs.csv 
#mv MADP_geno_imputed_subSNPs.vcf MADP_BIN2_SNPs.vcf

#python -m schnablelab.SNPcalling.base SubsamplingSNPs MADP_geno_imputed.vcf SNPsSetV2/MADP_BIN3_SNPs.csv 
#mv MADP_geno_imputed_subSNPs.vcf MADP_BIN3_SNPs.vcf

#python -m schnablelab.SNPcalling.base SubsamplingSNPs MADP_geno_imputed.vcf SNPsSetV2/MADP_BIN4.2_SNPs.csv 
#mv MADP_geno_imputed_subSNPs.vcf MADP_BIN4_SNPs.vcf

#python -m schnablelab.SNPcalling.base SubsamplingSNPs MADP_geno_imputed.vcf SNPsSetV2/MADP_BIN5_SNPs.csv 
#mv MADP_geno_imputed_subSNPs.vcf MADP_BIN5_SNPs.vcf

python -m schnablelab.SNPcalling.base SubsamplingSNPs MADP_geno_imputed.vcf SNPsSetV2/MADP_BIN6.2_SNPs.csv 
mv MADP_geno_imputed_subSNPs.vcf MADP_BIN6_SNPs.vcf

#python -m schnablelab.SNPcalling.base SubsamplingSNPs MADP_geno_imputed.vcf SNPsSetV2/MADP_BIN7_SNPs.csv 
#mv MADP_geno_imputed_subSNPs.vcf MADP_BIN7_SNPs.vcf
