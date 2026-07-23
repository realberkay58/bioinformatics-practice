import os
import pandas as pd
import vcfpy

# ==========================================
# CONFIGURATION & SETUP
# ==========================================
RESULTS_DIR = "results"
VCF_FILE_PATH = os.path.join(RESULTS_DIR, "sample_variants.vcf")
QUALITY_THRESHOLD = 30.0

os.makedirs(RESULTS_DIR, exist_ok=True)

# Generate a mock VCF file if it doesn't exist (ensures the script runs flawlessly)
def generate_mock_vcf():
    if not os.path.exists(VCF_FILE_PATH):
        print(f"[*] Generating sample VCF file at {VCF_FILE_PATH}...")
        vcf_content = """##fileformat=VCFv4.2
##source=SimulatedData
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chr3,length=198295559>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=GENE,Number=1,Type=String,Description="Overlapping Gene">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t10177\trs367896724\tA\tAC\t50.5\tPASS\tDP=100;GENE=DDX11L1
chr1\t10352\trs555500075\tT\tTA\t20.0\tFAIL\tDP=80;GENE=DDX11L1
chr1\t69091\trs532656360\tA\tC\t99.9\tPASS\tDP=120;GENE=OR4F5
chr2\t15000\t.\tG\tA\t15.5\tPASS\tDP=40;GENE=WASH7P
chr2\t23040\t.\tC\tT\t45.0\tPASS\tDP=200;GENE=WASH7P
chr3\t55000\t.\tT\tG\t85.2\tPASS\tDP=150;GENE=TRIM71
chr3\t56000\t.\tC\tA\t10.0\tFAIL\tDP=30;GENE=TRIM71
chr3\t57000\t.\tG\tC\t65.4\tPASS\tDP=90;GENE=TRIM71
"""
        with open(VCF_FILE_PATH, "w") as f:
            f.write(vcf_content)

# ==========================================
# MAIN ANALYSIS PIPELINE
# ==========================================
def run_snp_analysis():
    print("="*50)
    print(" 🧬 SNP Detection & Genetic Variation Pipeline")
    print("="*50)
    
    generate_mock_vcf()
    
    print(f"\n[*] Reading VCF file: {VCF_FILE_PATH}")
    print(f"[*] Applying Quality Threshold: QUAL > {QUALITY_THRESHOLD}\n")
    
    # Initialize counters and storage
    snp_records = []
    chromosome_counts = {}
    
    try:
        reader = vcfpy.Reader.from_path(VCF_FILE_PATH)
        
        for record in reader:
            # Check if it's an SNP (single base substitution)
            is_snp = all(len(alt.value) == 1 for alt in record.ALT) and len(record.REF) == 1
            
            if is_snp:
                chrom = record.CHROM
                qual = record.QUAL if record.QUAL is not None else 0.0
                gene = record.INFO.get('GENE', 'Unknown')
                
                # Count total SNPs per chromosome
                chromosome_counts[chrom] = chromosome_counts.get(chrom, 0) + 1
                
                # Filter by Quality Score
                if qual > QUALITY_THRESHOLD:
                    snp_records.append({
                        "Chromosome": chrom,
                        "Position": record.POS,
                        "Reference": record.REF,
                        "Alternate": record.ALT[0].value,
                        "Quality": qual,
                        "Gene_Region": gene
                    })
                    
    except Exception as e:
        print(f"[!] Error parsing VCF: {e}")
        return

    # Process and Export Results
    df_snps = pd.DataFrame(snp_records)
    
    print("📊 TOTAL SNPs PER CHROMOSOME (Unfiltered):")
    for chrom, count in sorted(chromosome_counts.items()):
        print(f"   - {chrom}: {count} SNPs")
        
    print(f"\n✅ HIGH-QUALITY SNPs FOUND (QUAL > {QUALITY_THRESHOLD}): {len(df_snps)}")
    
    if not df_snps.empty:
        print("\n🧬 Sample of High-Quality SNPs in Gene Regions:")
        print(df_snps[['Chromosome', 'Position', 'Quality', 'Gene_Region']].head().to_string(index=False))
        
        # Save to CSV
        output_csv = os.path.join(RESULTS_DIR, "filtered_snps_report.csv")
        df_snps.to_csv(output_csv, index=False)
        print(f"\n📁 Full report saved successfully to: {output_csv}")

if __name__ == "__main__":
    run_snp_analysis()