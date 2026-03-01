#!/bin/bash -l
#SBATCH -A sens2023005         
#SBATCH -M bianca              
#SBATCH -p core                
#SBATCH -n 2                   
#SBATCH -t 04:00:00           
#SBATCH -J split_vcf           
#SBATCH --mail-type=ALL   
#SBATCH --mail-user=

# Load required modules
module load bioinfo-tools
module load bcftools

# Directory for outputs
VCF=""
OUTDIR=""
mkdir -p "$OUTDIR"

# Get sample names
bcftools query -l "$VCF" > samples.txt

# Split each sample
while read -r SAMPLE; do
    echo "Processing $SAMPLE"
    bcftools view -s "$SAMPLE" -O z -o "${OUTDIR}/${SAMPLE}.vcf.gz" "$VCF"
    bcftools index "${OUTDIR}/${SAMPLE}.vcf.gz"
done < samples.txt

