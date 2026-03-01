#!/bin/bash -l
#SBATCH -A sens2023005
#SBATCH -M bianca
#SBATCH -p core
#SBATCH -n 10                    
#SBATCH -t 6:00:00              
#SBATCH -J filter_samples_optimized
#SBATCH --mail-type=ALL
#SBATCH --mail-user=

# Setup paths
PROJECT_BASE=""
VCF_DIR="${PROJECT_BASE}/1000genomes/20190825_Yale_cnvnator"
METADATA_FILE="${PROJECT_BASE}/1000genomes/igsr_samples.tsv"
ANALYSIS_DIR="${PROJECT_BASE}/"
OUTPUT_LISTS="${ANALYSIS_DIR}"
OUTPUT_VCFS="${ANALYSIS_DIR}/"

# Target populations and sample limits
declare -A POPULATIONS=(
    ["GBR"]="British"
    ["FIN"]="Finnish"
    ["CHS"]="Southern Han Chinese"
    ["CLM"]="Colombian"
    ["ACB"]="African Caribbean"
    ["IBS"]="Iberian"
    ["PEL"]="Peruvian"
    ["KHV"]="Kinh Vietnamese"
    ["PUR"]="Puerto Rican"
    ["CDX"]="Dai Chinese"
)
SAMPLES_PER_POP=200

# Load modules
module load bioinfo-tools
module load bcftools


echo "Target: ${SAMPLES_PER_POP} samples × ${#POPULATIONS[@]} populations = $((SAMPLES_PER_POP * ${#POPULATIONS[@]})) total samples"
# Create directories
mkdir -p "$OUTPUT_LISTS" "$OUTPUT_VCFS"

# Extract and limit samples per population
echo "Step 1: Extracting sample lists..."
for pop_code in "${!POPULATIONS[@]}"; do
    echo "  Processing ${pop_code} (${POPULATIONS[$pop_code]})..."
    
    tail -n +2 "$METADATA_FILE" | \
        awk -v pop="$pop_code" -F'\t' '$4 == pop {print $1}' | \
        head -n "$SAMPLES_PER_POP" > "${OUTPUT_LISTS}/${pop_code}_samples.txt"
    
    sample_count=$(wc -l < "${OUTPUT_LISTS}/${pop_code}_samples.txt")
    
    if [ "$sample_count" -lt "$SAMPLES_PER_POP" ]; then
        echo "    WARNING: Only ${sample_count} samples available (wanted ${SAMPLES_PER_POP})"
    fi
done

# Create combined sample list for the populations
echo "Creating combined sample list..."
cat "${OUTPUT_LISTS}"/*_samples.txt > "${OUTPUT_LISTS}/all_target_samples.txt"
total_samples=$(wc -l < "${OUTPUT_LISTS}/all_target_samples.txt")
echo "Total samples to extract: ${total_samples}"

# Step 2: Filter VCFs in parallel

filter_vcf() {
    local vcf="$1"
    local sample_list="$2"
    local output_dir="$3"
    
    local vcf_name=$(basename "$vcf")
    local output_vcf="${output_dir}/${vcf_name}"
    
    if [ -f "$output_vcf" ]; then
        echo "  [SKIP] ${vcf_name} already exists"
        return 0
    fi
    
    echo "  [PROCESSING] ${vcf_name}..."
    
    bcftools view \
        -S "$sample_list" \
        --force-samples \
        -O z \
        -o "$output_vcf" \
        "$vcf" 2>/dev/null
    
    # Index if successful
    if [ $? -eq 0 ]; then
        bcftools index "$output_vcf" 2>/dev/null
        echo "  [DONE] ${vcf_name}"
    else
        echo "  [ERROR] Failed to filter ${vcf_name}"
        rm -f "$output_vcf"
        return 1
    fi
}

export -f filter_vcf


mkdir -p "${OUTPUT_VCFS}/filtered"

vcf_count=$(ls "$VCF_DIR"/*.vcf.gz 2>/dev/null | wc -l)
echo "Found ${vcf_count} VCF files to process"
echo ""

# Process VCFs 
if command -v parallel &> /dev/null; then
    # Use GNU parallel if available (much better)
    echo "Using GNU parallel..."
    ls "$VCF_DIR"/*.vcf.gz | \
        parallel -j 8 filter_vcf {} "${OUTPUT_LISTS}/all_target_samples.txt" "${OUTPUT_VCFS}/filtered"
else
    # Fallback to xargs
    echo "Using xargs for parallel processing..."
    ls "$VCF_DIR"/*.vcf.gz | \
        xargs -P 8 -I {} bash -c "filter_vcf '{}' '${OUTPUT_LISTS}/all_target_samples.txt' '${OUTPUT_VCFS}/filtered'"
fi


echo "Sample counts by population:"
for pop_code in "${!POPULATIONS[@]}"; do
    count=$(wc -l < "${OUTPUT_LISTS}/${pop_code}_samples.txt")
    echo "  ${pop_code}: ${count} samples"
done

