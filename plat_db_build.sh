#!/bin/bash -l
#SBATCH -A sens2023005
#SBATCH -M bianca
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 8:00:00
#SBATCH -J build_new_plat_db
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dima.mohsin.1749@student.uu.se
#SBATCH -o plat_db_%j.out
#SBATCH -e plat_db_%j.err

# Load modules
module load python/3.9.5
module load bioinfo-tools

# Install required packages
echo "Installing required Python packages..."
python3 -m pip install --user --no-index \
    --find-links=/proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/python_packages \
    numpy

VCF_DIR="/proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/Degree_project/SVDB/platinum_samples/platinum_samples_reloaded"
DB_DIR="/proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/Degree_project/SVDB/Platinum_experiments"

echo ""
echo "--- Building G2 parents SVDB ---"
python3 -m svdb --build \
    --files \
        ${VCF_DIR}/NA12877.vcf.gz \
        ${VCF_DIR}/NA12878.vcf.gz \
    --prefix ${DB_DIR}/plat_G2_parents.db

echo ""
echo "--- Building G3 siblings SVDB ---"
python3 -m svdb --build \
    --files \
        ${VCF_DIR}/NA12879.vcf.gz \
        ${VCF_DIR}/NA12881.vcf.gz \
        ${VCF_DIR}/NA12882.vcf.gz \
        ${VCF_DIR}/NA12885.vcf.gz \
        ${VCF_DIR}/NA12886.vcf.gz \
    --prefix ${DB_DIR}/plat_G3_siblings.db

echo ""
echo "--- Building combined G2 + G3 SVDB ---"
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

