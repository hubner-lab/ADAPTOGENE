#!/bin/bash

# VCF.GZ Validation Report Generator
# Usage: ./vcf_validate.sh <vcf.gz> [output_table.tsv]
# Returns: TSV table with validation statistics

VCF_GZ=$1
OUTPUT=${2:-/dev/stdout}

# Function to safely count lines in compressed file
count_variants() {
    bcftools view -H "$1" 2>/dev/null | wc -l | awk '{print $1}'
}

count_samples() {
    bcftools query -l "$1" 2>/dev/null | wc -l | awk '{print $1}'
}

get_chromosomes() {
    bcftools view -H "$1" 2>/dev/null | cut -f1 | sort -u | wc -l | awk '{print $1}'
}

get_first_chr() {
    bcftools view -H "$1" 2>/dev/null | head -1 | cut -f1
}

get_last_chr() {
    bcftools view -H "$1" 2>/dev/null | tail -1 | cut -f1
}

check_index() {
    if [ -f "${1}.csi" ]; then
        echo "CSI_PRESENT"
    elif [ -f "${1}.tbi" ]; then
        echo "TBI_PRESENT" 
    else
        echo "NO_INDEX"
    fi
}

# Validation function
validate_vcf() {
    local file=$1
    local filename=$(basename "$file")
    local file_exists="NO"
    local file_readable="NO"
    local file_size="0"
    local is_bgzip="NO"
    local index_status="NO_INDEX"
    local header_valid="NO"
    local vcf_version="UNKNOWN"
    local sample_count="0"
    local variant_count="0"
    local chromosome_count="0"
    local first_chr="UNKNOWN"
    local last_chr="UNKNOWN"
    local header_lines="0"
    local has_format="NO"
    local has_info="NO"
    local bcftools_readable="NO"
    local validation_status="FAIL"
    local error_message="NONE"
    
    # Basic file checks
    if [ -f "$file" ]; then
        file_exists="YES"
        file_size=$(stat -f%z "$file" 2>/dev/null || stat -c%s "$file" 2>/dev/null || echo "0")
        
        if [ -r "$file" ]; then
            file_readable="YES"
        fi
        
        # Check if bgzipped
        if file "$file" | grep -q "gzip"; then
            is_bgzip="YES"
        fi
        
        # Check index
        index_status=$(check_index "$file")
        
        # Try to read with bcftools
        if bcftools view -h "$file" >/dev/null 2>&1; then
            bcftools_readable="YES"
            
            # Header validation
            if bcftools view -h "$file" 2>/dev/null | head -1 | grep -q "^##fileformat=VCF"; then
                header_valid="YES"
                vcf_version=$(bcftools view -h "$file" 2>/dev/null | head -1 | cut -d= -f2)
            fi
            
            # Count header lines
            header_lines=$(bcftools view -h "$file" 2>/dev/null | wc -l | awk '{print $1}')
            
            # Check for FORMAT and INFO in header
            if bcftools view -h "$file" 2>/dev/null | grep -q "^##FORMAT="; then
                has_format="YES"
            fi
            
            if bcftools view -h "$file" 2>/dev/null | grep -q "^##INFO="; then
                has_info="YES"
            fi
            
            # Count samples and variants (only if readable)
            sample_count=$(count_samples "$file")
            variant_count=$(count_variants "$file")
            chromosome_count=$(get_chromosomes "$file")
            
            if [ "$variant_count" -gt 0 ]; then
                first_chr=$(get_first_chr "$file")
                last_chr=$(get_last_chr "$file")
            fi
            
            # Overall validation
            if [ "$header_valid" = "YES" ] && [ "$variant_count" -gt 0 ] && [ "$sample_count" -gt 0 ]; then
                validation_status="PASS"
            elif [ "$variant_count" -eq 0 ]; then
                validation_status="FAIL"
                error_message="NO_VARIANTS"
            elif [ "$sample_count" -eq 0 ]; then
                validation_status="FAIL"
                error_message="NO_SAMPLES"
            else
                validation_status="FAIL"
                error_message="HEADER_INVALID"
            fi
        else
            error_message="BCFTOOLS_CANNOT_READ"
        fi
    else
        error_message="FILE_NOT_FOUND"
    fi
    
    # Output results
    echo -e "${filename}\t${file_exists}\t${file_readable}\t${file_size}\t${is_bgzip}\t${index_status}\t${bcftools_readable}\t${header_valid}\t${vcf_version}\t${header_lines}\t${has_format}\t${has_info}\t${sample_count}\t${variant_count}\t${chromosome_count}\t${first_chr}\t${last_chr}\t${validation_status}\t${error_message}"
}

# Header for output table
{
    echo -e "FILENAME\tFILE_EXISTS\tFILE_READABLE\tFILE_SIZE_BYTES\tIS_BGZIP\tINDEX_STATUS\tBCFTOOLS_READABLE\tHEADER_VALID\tVCF_VERSION\tHEADER_LINES\tHAS_FORMAT\tHAS_INFO\tSAMPLE_COUNT\tVARIANT_COUNT\tCHROMOSOME_COUNT\tFIRST_CHR\tLAST_CHR\tVALIDATION_STATUS\tERROR_MESSAGE"
    validate_vcf "$VCF_GZ"
} > "$OUTPUT"
