#!/bin/bash -l
#SBATCH -A sens2023005
#SBATCH -M bianca
#SBATCH -p core
#SBATCH -n 4                    
#SBATCH -t 26:00:00              
#SBATCH -J spatial_lev_comparison
#SBATCH --mail-type=ALL
<<<<<<< Updated upstream
#SBATCH --mail-user=
#SBATCH -o spatial_lev_comparison_%j.out
#SBATCH -e spatial_lev_comparison_%j.err
=======
#SBATCH --mail-user=dima.mohsin.1749@student.uu.se
#SBATCH -o plat_stat_lev_ham_comparison_%j.out
#SBATCH -e plat_stat_lev_ham_comparison_%j.err
>>>>>>> Stashed changes

# Load required modules
module load python/3.9.5

# Change to SVDB directory
cd 

# Install required packages
echo "Installing required Python packages..."
python3 -m pip install --user --no-index \
    --find-links=your_dir \
    numpy matplotlib psutil

# Set PYTHONPATH
export PYTHONPATH="${PWD}:${PYTHONPATH}"

# Verify SVDB modules

python3 -c "from svdb import database; print('✓ SVDB modules loaded')" || {
    echo "✗ Error: Cannot import SVDB modules"
    exit 1
}



DB_DIR="./Platinum_experiments"
REQUIRED_DBS=(
    "NA12877_NA12878.db"
    "NA12877_NA12879.db"
    "NA12879_NA12881.db"
)

MISSING=0
for db in "${REQUIRED_DBS[@]}"; do
    if [ -f "${DB_DIR}/${db}" ]; then
        echo "   ${db}"
    else
        echo "   ${db} - NOT FOUND"
        MISSING=$((MISSING + 1))
    fi
done

if [ $MISSING -gt 0 ]; then
    echo ""
    echo "✗ Error: $MISSING required database file(s) not found in ${DB_DIR}"
    exit 1
fi



# Run comparison
python3 cluster_ratio_hamming_vs_lv.py

# Check exit status
if [ $? -eq 0 ]; then

    echo "Generated files:"
    ls -lh clustering_ratio_*.png processing_time_*.png 2>/dev/null || echo "  No plot files found"
else

    echo "✗ Error: Comparison failed"
    exit 1
fi

<<<<<<< Updated upstream
echo "Job completed: $(date)"
=======
echo ""
echo "Job completed: $(date)"
>>>>>>> Stashed changes
