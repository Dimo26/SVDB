#!/bin/bash -l
#SBATCH -A sens2023005
#SBATCH -M bianca
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 8:00:00
#SBATCH -J build_new_plat_db
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH -o plat_db_%j.out
#SBATCH -e plat_db_%j.err

# Load modules
module load python/3.9.5
module load bioinfo-tools

# Install required packages
python3 -m pip install --user --no-index \
    --find-links=/ \
    numpy

VCF_DIR=""
DB_DIR=""

echo "Building G2 parent"
python3 -m svdb --build \
    --files \
        ${VCF_DIR}/NA12877.vcf.gz \
        ${VCF_DIR}/NA12878.vcf.gz \
    --prefix ${DB_DIR}/plat_G2_parents.db


echo "Building G3 siblings"
python3 -m svdb --build \
    --files \
        ${VCF_DIR}/NA12879.vcf.gz \
        ${VCF_DIR}/NA12881.vcf.gz \
        ${VCF_DIR}/NA12882.vcf.gz \
        ${VCF_DIR}/NA12885.vcf.gz \
        ${VCF_DIR}/NA12886.vcf.gz \
    --prefix ${DB_DIR}/plat_G3_siblings.db


echo "Building combined G2 + G3"
python3 -m svdb --build \
    --files \
        ${VCF_DIR}/NA12877.vcf.gz \
        ${VCF_DIR}/NA12878.vcf.gz \
        ${VCF_DIR}/NA12879.vcf.gz \
        ${VCF_DIR}/NA12881.vcf.gz \
        ${VCF_DIR}/NA12882.vcf.gz \
        ${VCF_DIR}/NA12885.vcf.gz \
        ${VCF_DIR}/NA12886.vcf.gz \
    --prefix ${DB_DIR}/plat_G2_G3_parents_kids.db

python3 -m svdb --build --files ${VCF_DIR}/NA12877.vcf.gz --prefix ${DB_DIR}/NA12877.db
python3 -m svdb --build --files ${VCF_DIR}/NA12878.vcf.gz --prefix ${DB_DIR}/NA12878.db
python3 -m svdb --build --files ${VCF_DIR}/NA12879.vcf.gz --prefix ${DB_DIR}/NA12879.db
python3 -m svdb --build --files ${VCF_DIR}/NA12881.vcf.gz --prefix ${DB_DIR}/NA12881.db
python3 -m svdb --build --files ${VCF_DIR}/NA12882.vcf.gz --prefix ${DB_DIR}/NA12882.db
python3 -m svdb --build --files ${VCF_DIR}/NA12885.vcf.gz --prefix ${DB_DIR}/NA12885.db
python3 -m svdb --build --files ${VCF_DIR}/NA12886.vcf.gz --prefix ${DB_DIR}/NA12886.db

