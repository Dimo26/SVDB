#!/bin/bash
#SBATCH -A sens2023005
#SBATCH -M bianca
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 4:00:00
#SBATCH -J svdb_1000genomes_database
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH -o build_db_%j.out
#SBATCH -e build_db_%j.err

module load bioinfo-tools
module load python/3.9.5
module load bcftools

# Directories
Out_dir=""
VCF_DIR=""

# Create output directory if it doesn't exist
mkdir -p "$Out_dir"
cd "$Out_dir"

# Find VCF files and create a shuffled list
all_vcfs_file="${Out_dir}/all_vcfs_list.txt"
find "$VCF_DIR" -name "*.sites.vcf.gz" > "$all_vcfs_file"

# Shuffle for random sampling
shuf "$all_vcfs_file" > "${all_vcfs_file}.shuffled"

total_vcfs=$(wc -l < "${all_vcfs_file}.shuffled")
echo "Total VCFs found: $total_vcfs"

for size in 10 50 100 500 1000; do
    if [[ $size -gt $total_vcfs ]]; then
        echo "Skipping size $size (only $total_vcfs samples available)"
        continue
    fi
    
    selected_vcfs=$(head -n $size "${all_vcfs_file}.shuffled" | tr '\n' ' ')
    
    python -m svdb --build \
        --files $selected_vcfs \
        --prefix "${Out_dir}/svdb_${size}samples"
    
    if [[ $? -eq 0 ]]; then
        echo "Completed: svdb_${size}samples.db"
        # Show database info
        if [[ -f "${Out_dir}/svdb_${size}samples.db" ]]; then
            db_size=$(du -h "${Out_dir}/svdb_${size}samples.db" | cut -f1)
            echo "  Database size: $db_size"
        fi
    else
        echo "Failed to build database with $size samples"
    fi
    echo ""
done

all_vcfs=$(cat "${all_vcfs_file}.shuffled" | tr '\n' ' ')

python -m svdb --build \
    --files $all_vcfs \
    --prefix "${Out_dir}/svdb_all_samples"

if [[ $? -eq 0 ]]; then
    echo "Completed: svdb_all_samples.db"
    if [[ -f "${Out_dir}/svdb_all_samples.db" ]]; then
        db_size=$(du -h "${Out_dir}/svdb_all_samples.db" | cut -f1)
        echo "  Database size: $db_size"
    fi
else
    echo "Failed to build database with all samples"
fi

# Cleanup temporary files
rm -f "$all_vcfs_file" "${all_vcfs_file}.shuffled"

echo "Database files created:"
ls -lh "${Out_dir}"/svdb_*.db 2>/dev/null | awk '{print $9, $5}'


