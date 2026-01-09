#!/usr/bin/env python3
"""
APOE Haplotype Caller from VCF Files

Identifies APOE haplotypes (ε2, ε3, ε4) based on two SNPs:
- rs429358 (chr19:45411941) - T>C
- rs7412 (chr19:45412079) - C>T
"""

import argparse
import sys
from pathlib import Path
import pandas as pd
import pysam


# APOE haplotype lookup table
# Key: (rs429358_genotype, rs7412_genotype)
# Value: APOE haplotype
# Based on: ε2=T,T; ε3=T,C; ε4=C,C (for rs429358, rs7412 respectively)
APOE_LOOKUP = {
    ('T/T', 'T/T'): 'ε2/ε2',  # rs429358=TT, rs7412=TT
    ('T/T', 'C/T'): 'ε2/ε3',  # rs429358=TT, rs7412=CT
    ('C/T', 'T/T'): 'ε2/ε4',  # rs429358=CT, rs7412=TT (cis phasing)
    ('C/T', 'C/T'): 'ε2/ε4',  # rs429358=CT, rs7412=CT (trans phasing)
    ('T/T', 'C/C'): 'ε3/ε3',  # rs429358=TT, rs7412=CC
    ('C/T', 'C/C'): 'ε3/ε4',  # rs429358=CT, rs7412=CC
    ('C/C', 'C/C'): 'ε4/ε4',  # rs429358=CC, rs7412=CC
}


def parse_genotype(variant, sample_name):
    """
    Parse genotype from a pysam VariantRecord for a specific sample.
    
    Args:
        variant: pysam VariantRecord
        sample_name: Sample name to extract genotype for
    
    Returns:
        Tuple of (allele1, allele2) or (None, None) if missing
    """
    if sample_name not in variant.samples:
        return None, None
    
    sample = variant.samples[sample_name]
    gt = sample.get('GT', None)
    
    if gt is None or None in gt:
        return None, None
    
    # Get alleles based on genotype indices
    alleles = [variant.ref] + list(variant.alts) if variant.alts else [variant.ref]
    
    try:
        allele1 = alleles[gt[0]]
        allele2 = alleles[gt[1]]
        return allele1, allele2
    except (IndexError, TypeError):
        return None, None


def genotype_to_string(allele1, allele2):
    """
    Convert allele pair to genotype string.
    
    Args:
        allele1: First allele
        allele2: Second allele
    
    Returns:
        Genotype string (e.g., 'T/T', 'C/T') or None if missing
    """
    if allele1 is None or allele2 is None:
        return None
    return f"{allele1}/{allele2}"


def determine_apoe_haplotype(rs429358_geno, rs7412_geno):
    """
    Determine APOE haplotype from two SNP genotypes.
    
    Args:
        rs429358_geno: Genotype string for rs429358
        rs7412_geno: Genotype string for rs7412
    
    Returns:
        APOE haplotype string or None if missing/ambiguous
    """
    if rs429358_geno is None or rs7412_geno is None:
        return None
    
    # Normalize genotype order (e.g., T/C -> C/T for lookup)
    def normalize_geno(geno):
        alleles = geno.split('/')
        return '/'.join(sorted(alleles))
    
    rs429358_norm = normalize_geno(rs429358_geno)
    rs7412_norm = normalize_geno(rs7412_geno)
    
    return APOE_LOOKUP.get((rs429358_norm, rs7412_norm), 'Unknown')


def call_apoe_haplotypes(vcf_file, output_file=None, chrom='19', 
                         rs429358_pos=45411941, rs7412_pos=45412079):
    """
    Call APOE haplotypes from tabix-indexed VCF file.
    
    Args:
        vcf_file: Path to tabix-indexed VCF file (.vcf.gz or .bcf)
        output_file: Optional output file path for results
        chrom: Chromosome name (default: '19')
        rs429358_pos: Position of rs429358 (default: 45411941)
        rs7412_pos: Position of rs7412 (default: 45412079)
    """
    vcf_path = Path(vcf_file)
    
    # Check if VCF file exists
    if not vcf_path.exists():
        print(f"Error: {vcf_file} not found", file=sys.stderr)
        sys.exit(1)
    
    # Check for index
    index_extensions = ['.tbi', '.csi']
    has_index = any((vcf_path.parent / f"{vcf_path.name}{ext}").exists() 
                    for ext in index_extensions)
    
    if not has_index:
        print(f"Error: No index found for {vcf_file}", file=sys.stderr)
        print("Please create an index with: tabix -p vcf {vcf_file}", file=sys.stderr)
        sys.exit(1)
    
    try:
        # Open VCF file
        vcf = pysam.VariantFile(str(vcf_path))
        
        # Get sample names
        samples = list(vcf.header.samples)
        print(f"Found {len(samples)} samples in VCF")
        
        # Fetch rs429358 variant
        rs429358_variant = None
        try:
            for variant in vcf.fetch(chrom, rs429358_pos - 1, rs429358_pos):
                if variant.pos == rs429358_pos:
                    rs429358_variant = variant
                    break
        except Exception as e:
            print(f"Error fetching rs429358 at {chrom}:{rs429358_pos}: {e}", file=sys.stderr)
            sys.exit(1)
        
        if rs429358_variant is None:
            print(f"Error: No variant found at {chrom}:{rs429358_pos} (rs429358)", file=sys.stderr)
            sys.exit(1)
        
        # Fetch rs7412 variant
        rs7412_variant = None
        try:
            for variant in vcf.fetch(chrom, rs7412_pos - 1, rs7412_pos):
                if variant.pos == rs7412_pos:
                    rs7412_variant = variant
                    break
        except Exception as e:
            print(f"Error fetching rs7412 at {chrom}:{rs7412_pos}: {e}", file=sys.stderr)
            sys.exit(1)
        
        if rs7412_variant is None:
            print(f"Error: No variant found at {chrom}:{rs7412_pos} (rs7412)", file=sys.stderr)
            sys.exit(1)
        
        # Print variant information
        rs429358_alts = ','.join(rs429358_variant.alts) if rs429358_variant.alts else '.'
        rs7412_alts = ','.join(rs7412_variant.alts) if rs7412_variant.alts else '.'
        
        print(f"\nFound rs429358 at {chrom}:{rs429358_variant.pos}")
        print(f"  REF={rs429358_variant.ref}, ALT={rs429358_alts}")
        print(f"  ID={rs429358_variant.id if rs429358_variant.id else '.'}")
        
        print(f"\nFound rs7412 at {chrom}:{rs7412_variant.pos}")
        print(f"  REF={rs7412_variant.ref}, ALT={rs7412_alts}")
        print(f"  ID={rs7412_variant.id if rs7412_variant.id else '.'}")
        print()
        
        # Extract genotypes for all samples
        results = []
        for sample in samples:
            # Get genotypes
            rs429358_alleles = parse_genotype(rs429358_variant, sample)
            rs7412_alleles = parse_genotype(rs7412_variant, sample)
            
            rs429358_geno = genotype_to_string(rs429358_alleles[0], rs429358_alleles[1])
            rs7412_geno = genotype_to_string(rs7412_alleles[0], rs7412_alleles[1])
            
            # Determine APOE haplotype
            apoe_haplotype = determine_apoe_haplotype(rs429358_geno, rs7412_geno)
            
            results.append({
                'Sample': sample,
                'rs429358': rs429358_geno if rs429358_geno else 'missing',
                'rs7412': rs7412_geno if rs7412_geno else 'missing',
                'APOE_haplotype': apoe_haplotype if apoe_haplotype else 'missing'
            })
        
        vcf.close()
        
        # Create results dataframe
        results_df = pd.DataFrame(results)
        
        # Print summary
        print("APOE Haplotype Summary:")
        print("=" * 50)
        print(results_df['APOE_haplotype'].value_counts().to_string())
        print()
        
        # Save or display results
        if output_file:
            results_df.to_csv(output_file, sep='\t', index=False)
            print(f"Results written to: {output_file}")
        else:
            print("\nFirst 10 samples:")
            print(results_df.head(10).to_string(index=False))
            if len(results_df) > 10:
                print(f"\n... ({len(results_df)} total samples)")
        
        return results_df
    
    except Exception as e:
        print(f"Error processing VCF: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description='Call APOE haplotypes from tabix-indexed VCF file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --vcf mydata.vcf.gz
  %(prog)s --vcf mydata.vcf.gz --output apoe_results.txt
  %(prog)s --vcf mydata.vcf.gz --chrom chr19
  %(prog)s --vcf mydata.vcf.gz --rs429358-pos 45411941 --rs7412-pos 45412079

Note: VCF file must be bgzip-compressed and tabix-indexed.
  To index: bgzip mydata.vcf && tabix -p vcf mydata.vcf.gz
        """
    )
    
    parser.add_argument(
        '--vcf',
        required=True,
        help='Path to tabix-indexed VCF file (.vcf.gz or .bcf)'
    )
    
    parser.add_argument(
        '--output', '-o',
        help='Output file for results (default: print to stdout)'
    )
    
    parser.add_argument(
        '--chrom',
        default='19',
        help='Chromosome name (default: 19)'
    )
    
    parser.add_argument(
        '--rs429358-pos',
        type=int,
        default=45411941,
        help='Position of rs429358 (default: 45411941)'
    )
    
    parser.add_argument(
        '--rs7412-pos',
        type=int,
        default=45412079,
        help='Position of rs7412 (default: 45412079)'
    )
    
    args = parser.parse_args()
    
    try:
        call_apoe_haplotypes(
            args.vcf,
            args.output,
            args.chrom,
            args.rs429358_pos,
            args.rs7412_pos
        )
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
