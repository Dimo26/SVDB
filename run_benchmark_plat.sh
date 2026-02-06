#!/bin/bash -l
#SBATCH -A sens2023005
#SBATCH -M bianca
#SBATCH -p core
#SBATCH -n 4                    
#SBATCH -t 12:00:00              
#SBATCH -J svdb_benchmark_final
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dima.mohsin.1749@student.uu.se
#SBATCH -o stat_plat_benchmark_%j.out
#SBATCH -e err_benchmark_%j.err

# Load required modules
module load python/3.9.5
module load scikit-learn/1.4.2

# Change to SVDB directory
cd /proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/Degree_project/SVDB

# Install required packages
echo "Installing required Python packages..."
python3 -m pip install --user --no-index \
    --find-links=/proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/python_packages \
    numpy psutil matplotlib

# Set PYTHONPATH
export PYTHONPATH="${PWD}:${PYTHONPATH}"

# Verify SVDB modules
echo ""
echo "Verifying SVDB installation..."
python3 -c "from svdb import database, export_module, optics_clustering, interval_tree_overlap; print('✓ All SVDB modules loaded')" || {
    echo "✗ Error: Cannot import SVDB modules"
    exit 1
}

# Find all databases
echo ""
echo "Searching for databases..."
DB_FILES=$(find . -maxdepth 1 -type f -name "*.db" | \
  grep -E './(NA12877|NA12878|NA12879|NA12881|NA12882|NA12885|NA12886|plat_G2_G3_parents_kids|plat_G2_parents|plat_G3_siblings)\.db$' | \
  sort)
DB_COUNT=$(echo "$DB_FILES" | grep -v '^$' | wc -l)

if [ $DB_COUNT -eq 0 ]; then
    echo "✗ Error: No database files found"
    exit 1
fi

echo "Found $DB_COUNT database file(s):"
echo "$DB_FILES" | sed 's/^/  /'
echo ""

# Run benchmark
echo "========================================"
echo "RUNNING SVDB BENCHMARK - FINAL VERSION"
echo "========================================"
echo "Features:"
echo "  • All SVs plotted together (one plot per algorithm+state)"
echo "  • Before/After Hamming comparison for insertions"
echo "  • Stacked bar charts for time/memory"
echo "  • Legend at bottom of plots"
echo "  • Detailed statistics in console"
echo "  • Noise analysis"
echo ""
echo "Output per database:"
echo "  • 1 performance comparison plot"
echo "  • 6 algorithm plots (3 algos × 2 states)"
echo "  • 3 insertion comparison plots"
echo "  Total: 10 plots per database"
echo "========================================"
echo ""

# Run with default parameters (distance=500bp, min_samples=2)
python3 algorithm_benchmark_plat.py $DB_FILES

# Check exit status
if [ $? -eq 0 ]; then
    echo ""
    echo "========================================"
    echo "✓ BENCHMARK COMPLETED SUCCESSFULLY"
    echo "========================================"
    echo ""
    echo "Generated plots:"
    ls -1 *chr1*.png *performance*.png 2>/dev/null | wc -l | awk '{print "  Total: " $1 " plot files"}'
    echo ""
    echo "Breakdown by database:"
    for db_file in $DB_FILES; do
        db_name=$(basename "$db_file" .db)
        plot_count=$(ls -1 ${db_name}*.png 2>/dev/null | wc -l)
        echo "  $db_name: $plot_count plots"
    done
    echo ""
else
    echo ""
    echo "✗ BENCHMARK FAILED - Check error log"
    exit 1
fi

# Summary
echo ""
echo "========================================"
echo "JOB SUMMARY"
echo "========================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Databases processed: $DB_COUNT"
echo "Output files:"
echo "  - Plots: *.png"
echo "  - Stdout: benchmark_${SLURM_JOB_ID}.out"
echo "  - Stderr: benchmark_${SLURM_JOB_ID}.err"
echo ""
echo "Check benchmark_${SLURM_JOB_ID}.out for:"
echo "  • Detailed SV statistics"
echo "  • Clustering results per algorithm"
echo "  • Noise analysis"
echo "  • Performance metrics"
echo "========================================"
