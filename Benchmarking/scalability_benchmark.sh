#!/bin/bash -l
#SBATCH -A sens2023005
#SBATCH -M bianca
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 36:00:00
#SBATCH -J svdb_scalability
#SBATCH --mail-type=ALL
#SBATCH --mail-user= 
#SBATCH -o scalability_%j.out
#SBATCH -e scalability_%j.err

module load python/3.9.5

#find SVDB dir
cd 

python3 -m pip install --user --no-index \
    --find-links=your_dir \
    numpy psutil matplotlib 2>&1 | grep -v "Requirement already satisfied" > /dev/null

export PYTHONPATH="${PWD}:${PYTHONPATH}"

mkdir -p scalability_results

python3 -c "from svdb.database import DB; from svdb.export_module import DBSCAN; from svdb.optics_clustering import optics_cluster; from svdb.interval_tree_overlap import interval_tree_cluster" 2>&1

if [ $? -ne 0 ]; then
    echo "Module import failed"
    exit 1
fi

DB_COUNT=$(ls -1 svdb_*samples.db 2>/dev/null | grep -v "all_samples" | wc -l)

if [ $DB_COUNT -eq 0 ]; then
    echo "No sample databases found"
    exit 1
fi

echo "Found $DB_COUNT databases"
ls -1 svdb_*samples.db 2>/dev/null | grep -v "all_samples"

python3 scalability_benchmark.py

if [ $? -eq 0 ]; then
    mv scalability_*.png scalability_results/ 2>/dev/null
    echo "Complete"
    ls scalability_results/ 2>/dev/null
else
    echo "Failed"
    exit 1
fi
