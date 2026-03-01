#!/bin/bash -l
#SBATCH -A sens2023005
#SBATCH -M bianca
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 17:00:00
#SBATCH -J build_pairwise_dbs
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH -o plat_pairwise_db_%j.out
#SBATCH -e plat_pairwise_db_%j.err

# Load modules
module load python/3.9.5
module load bioinfo-tools

# Install required packages
echo "Installing required Python packages..."
python3 -m pip install --user --no-index \
    --find-links=
    numpy

VCF_DIR=""
DB_DIR=""

# Navigate to SVDB directory
cd 

# Define all NA samples
SAMPLES=(NA12877 NA12878 NA12879 NA12881 NA12882 NA12885 NA12886)


echo "Total samples: ${#SAMPLES[@]}"
echo "Pairwise combinations: $((${#SAMPLES[@]} * (${#SAMPLES[@]} - 1) / 2))"
echo ""

# Counter for progress
counter=0
total=$((${#SAMPLES[@]} * (${#SAMPLES[@]} - 1) / 2))

# Build all pairwise combinations
for ((i=0; i<${#SAMPLES[@]}; i++)); do
    for ((j=i+1; j<${#SAMPLES[@]}; j++)); do
        sample1=${SAMPLES[$i]}
        sample2=${SAMPLES[$j]}
        
        counter=$((counter + 1))
        
        echo ""
        echo "[$counter/$total] Building: ${sample1}_${sample2}"

        
        python3 -m svdb --build \
            --files \
                ${VCF_DIR}/${sample1}.vcf.gz \
                ${VCF_DIR}/${sample2}.vcf.gz \
            --prefix ${DB_DIR}/${sample1}_${sample2}
        
        if [ $? -eq 0 ]; then
            echo "✓ Successfully built ${sample1}_${sample2}.db"
        else
            echo "✗ Failed to build ${sample1}_${sample2}.db"
        fi
    done
done

echo "Databases created in: $DB_DIR"
echo "Total databases: $(ls -1 ${DB_DIR}/*_NA12*.db 2>/dev/null | wc -l)"
echo ""
echo "Database list:"
ls -1 ${DB_DIR}/*_NA12*.db 2>/dev/null | while read db; do
    size=$(du -h "$db" | cut -f1)
    echo "  $(basename $db) ($size)"
done


# echo ""
# echo "Building G2 parents"
# python3 -m svdb --build \
#     --files \
#         ${VCF_DIR}/NA12877.vcf.gz \
#         ${VCF_DIR}/NA12878.vcf.gz \
#     --prefix ${DB_DIR}/plat_G2_parents

# echo ""
# echo "Building G3 siblings"
# python3 -m svdb --build \
#     --files \
#         ${VCF_DIR}/NA12879.vcf.gz \
#         ${VCF_DIR}/NA12881.vcf.gz \
#         ${VCF_DIR}/NA12882.vcf.gz \
#         ${VCF_DIR}/NA12885.vcf.gz \
#         ${VCF_DIR}/NA12886.vcf.gz \
#     --prefix ${DB_DIR}/plat_G3_siblings

# echo ""
# echo "Building combined G2 + G3"
# python3 -m svdb --build \
#     --files \
#         ${VCF_DIR}/NA12877.vcf.gz \
#         ${VCF_DIR}/NA12878.vcf.gz \
#         ${VCF_DIR}/NA12879.vcf.gz \
#         ${VCF_DIR}/NA12881.vcf.gz \
#         ${VCF_DIR}/NA12882.vcf.gz \
#         ${VCF_DIR}/NA12885.vcf.gz \
#         ${VCF_DIR}/NA12886.vcf.gz \
#     --prefix ${DB_DIR}/plat_G2_G3_parents_kids


# python3 -m svdb --build --files ${VCF_DIR}/NA12877.vcf.gz --prefix ${DB_DIR}/NA12877
# python3 -m svdb --build --files ${VCF_DIR}/NA12878.vcf.gz --prefix ${DB_DIR}/NA12878
# python3 -m svdb --build --files ${VCF_DIR}/NA12879.vcf.gz --prefix ${DB_DIR}/NA12879
# python3 -m svdb --build --files ${VCF_DIR}/NA12881.vcf.gz --prefix ${DB_DIR}/NA12881
# python3 -m svdb --build --files ${VCF_DIR}/NA12882.vcf.gz --prefix ${DB_DIR}/NA12882
# python3 -m svdb --build --files ${VCF_DIR}/NA12885.vcf.gz --prefix ${DB_DIR}/NA12885
# python3 -m svdb --build --files ${VCF_DIR}/NA12886.vcf.gz --prefix ${DB_DIR}/NA12886

