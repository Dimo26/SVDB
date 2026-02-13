#!/bin/bash -l
#SBATCH -A sens2023005
#SBATCH -M bianca
#SBATCH -p core
#SBATCH -n 4                    
#SBATCH -t 26:00:00              
#SBATCH -J hamming_lev_comparison
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dima.mohsin.1749@student.uu.se
#SBATCH -o hamming_lev_comparison_%j.out
#SBATCH -e hamming_lev_comparison_%j.err

# Load required modules
module load python/3.9.5

# Change to SVDB directory
cd /proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/Degree_project/SVDB

# Install required packages
echo "Installing required Python packages..."
python3 -m pip install --user --no-index \
    --find-links=/proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/python_packages \
    numpy matplotlib

# Set PYTHONPATH
export PYTHONPATH="${PWD}:${PYTHONPATH}"

# Verify SVDB modules
echo ""
echo "Verifying SVDB installation..."
python3 -c "from svdb import database; print('✓ SVDB modules loaded')" || {
    echo "✗ Error: Cannot import SVDB modules"
    exit 1
}

# Check for databases
echo ""
echo "Checking for pairwise databases..."
DB_COUNT=$(find ./Platinum_experiments -type f -name "NA12*_NA12*.db" 2>/dev/null | wc -l)

if [ $DB_COUNT -eq 0 ]; then
    echo "✗ Error: No pairwise database files found in ./Platinum_experiments"
    echo ""
    echo "Please run pairwise database creation first"
    exit 1
fi

echo "✓ Found $DB_COUNT pairwise databases"

# Run comparison
echo ""
echo "Starting Hamming vs Levenshtein comparison..."
echo "========================================================================"
python3 cluster_ratio_hamming_vs_lv.py

# Check exit status
if [ $? -eq 0 ]; then
    echo ""
    echo "=================================================="
    echo "✓ Comparison completed successfully"
    echo ""
    echo "Generated files:"
    ls -lh comparison_*.png 2>/dev/null || echo "  No plot files found"
else
    echo ""
    echo "✗ Error: Comparison failed"
    exit 1
fi

echo ""
echo "Job completed: $(date)"
