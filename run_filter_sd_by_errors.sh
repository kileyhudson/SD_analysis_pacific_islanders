#!/bin/bash

TEST_MODE=false
TEST_SAMPLE="UKS17D00107"  # Change this to test different samples

# Path to python script
PYTHON_SCRIPT="/projects/standard/hsiehph/shared/DIR_homes/hudso501/analysisPI/scripts/filter_sd_by_errors.py"

# Path to sample list  
SAMPLE_LIST="/projects/standard/hsiehph/shared/DIR_homes/hudso501/analysisPI/data/metadata/PI_sample_names.txt"

# Base paths
WGAC_BASE="/scratch.global/hudso501/projects/wgac/pacificIslander"
ERROR_BASE="/projects/standard/hsiehph/shared/globus-incoming/assembly_qc_files"

# Read sample names and loop through them
while read -r sample; do
    # Skip other samples during testing
    if [[ "$TEST_MODE" == "true" && "$sample" != "$TEST_SAMPLE" ]]; then
        continue
    fi
    echo "=== Processing sample: ${sample} ==="

    # Find the error directory for this sample
    ERROR_DIR=$(find ${ERROR_BASE} -maxdepth 1 -type d -name "*_${sample}")
        
    if [[ -z "$ERROR_DIR" ]]; then
        echo "ERROR: Could not find error directory for ${sample}"
        continue
    fi
    
    echo "Found error directory: ${ERROR_DIR}"
    
    for hap in hap1 hap2; do
        echo "  Processing ${sample} ${hap}..."
        
        # Build file paths
        WGAC_FILE="${WGAC_BASE}/${sample}/${hap}/data/GenomicSuperDup.tab"
        SMALL_ERROR="${ERROR_DIR}/${hap}/small_scale_error.bed"
        STRUCT_ERROR="${ERROR_DIR}/${hap}/structural_error.bed"
        
        # Check if files exist
        if [[ ! -f "$WGAC_FILE" ]]; then
            echo "    WARNING: WGAC file not found: $WGAC_FILE"
            continue
        fi
        
        if [[ ! -f "$SMALL_ERROR" ]]; then
            echo "    WARNING: Small error file not found: $SMALL_ERROR"
            continue
        fi
        
        if [[ ! -f "$STRUCT_ERROR" ]]; then
            echo "    WARNING: Structural error file not found: $STRUCT_ERROR"
            continue
        fi
        
        # Run the Python script
        python ${PYTHON_SCRIPT} \
            --szGenomicSuperDup "${WGAC_FILE}" \
            --szSmallScaleErrors "${SMALL_ERROR}" \
            --szStructuralErrors "${STRUCT_ERROR}" \
            --szSampleName "${sample}" \
            --szHaplotype "${hap}"
            
    done
done < "$SAMPLE_LIST"

# Define the summary file
SUMMARY="/projects/standard/hsiehph/shared/DIR_homes/hudso501/analysisPI/data/filter_by_asm_errors/filtering_errors_summary.tsv"

# Create the header
echo -e "Sample\tHaplotype\tSD_pairs\tError_Overlap_pairs\tFiltered_pairs\tPercent_Removed" > "$SUMMARY"

# Find all individual summary files and combine them
find /projects/standard/hsiehph/shared/DIR_homes/hudso501/analysisPI/data/filter_by_asm_errors/ -name "*.filtering_summary.txt" | while read summary_file; do
    # Skip the header line and append the data line to master file
    tail -n +2 "$summary_file" >> "$SUMMARY"
done

# Delete individual summary files since we have the master file
find /projects/standard/hsiehph/shared/DIR_homes/hudso501/analysisPI/data/filter_by_asm_errors/ -name "*.filtering_summary.txt" -delete

