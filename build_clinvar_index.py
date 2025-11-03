#!/usr/bin/env python3
"""
Convert ClinVar VCF to optimized JSON index for browser-based processing.
Only includes pathogenic/likely pathogenic variants with rsIDs.
"""

import json
import re
import sys
from pathlib import Path


def parse_info_field(info_str):
    """Parse VCF INFO field into a dictionary"""
    info_dict = {}
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
        else:
            info_dict[item] = True
    return info_dict


def parse_clnsig(clnsig_text):
    """Parse CLNSIG field and map to numeric class (4 or 5)"""
    if not clnsig_text:
        return None
    
    txt = clnsig_text.lower().replace('_', ' ')
    
    exclude_terms = [
        'conflicting', 'uncertain', 'benign', 'not provided',
        'no classification', 'association', 'protective',
        'affects', 'confers sensitivity'
    ]
    
    for term in exclude_terms:
        if term in txt:
            return None
    
    # Pathogenic (class 5)
    if txt.startswith('pathogenic'):
        if 'likely' not in txt:
            return '5'
        elif txt.startswith('pathogenic/likely pathogenic') or txt.startswith('pathogenic|likely pathogenic'):
            return '5'
    
    # Likely pathogenic (class 4)
    if txt.startswith('likely pathogenic'):
        return '4'
    
    # Handle pipe-separated values
    if '|' in txt:
        first_part = txt.split('|')[0].strip()
        if first_part == 'pathogenic':
            return '5'
        elif first_part == 'likely pathogenic':
            return '4'
    
    return None


def normalize_ref_alt(ref, alt):
    """Normalize REF/ALT values to single-nucleotide alleles"""
    def tokens(s):
        return [t for t in re.split(r"[,/;|\s]+", str(s).strip().upper()) if t]
    
    valid = {"A", "C", "G", "T"}
    ref_parts = tokens(ref) if ref is not None else []
    alt_parts = tokens(alt) if alt is not None else []
    
    if ref_parts and not alt_parts and len(ref_parts) >= 2:
        ref_norm = ref_parts[0]
        alt_parts = ref_parts[1:]
    else:
        ref_norm = ref_parts[0] if ref_parts else None
    
    if not ref_norm or ref_norm not in valid:
        return None, None
    
    alt_parts = [a for a in alt_parts if a in valid]
    if not alt_parts:
        return ref_norm, None
    
    alt_norm = ",".join(sorted(set(alt_parts)))
    return ref_norm, alt_norm


def build_clinvar_index(vcf_path, output_json):
    """Build optimized JSON index from ClinVar VCF"""
    print(f"Reading ClinVar VCF: {vcf_path}")
    index = {}
    total_lines = 0
    processed = 0
    
    # Count total lines for progress
    with open(vcf_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if not line.startswith('#'):
                total_lines += 1
    
    print(f"Processing {total_lines:,} variants...")
    
    with open(vcf_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
            
            chrom, pos, var_id, ref, alt, qual, filt, info = parts[:8]
            
            info_dict = parse_info_field(info)
            
            # Extract RS ID
            rs_raw = info_dict.get('RS')
            if not rs_raw or rs_raw.strip() in ('', '.', '-1'):
                continue
            
            rs_ids = [f"rs{r.strip()}".lower() for r in rs_raw.split(',') if r.strip()]
            if not rs_ids:
                continue
            
            # Get clinical significance
            clnsig_raw = info_dict.get('CLNSIG', '')
            cls = parse_clnsig(clnsig_raw)
            if cls not in ('4', '5'):
                continue
            
            # Extract gene
            geneinfo = info_dict.get('GENEINFO', '')
            gene = 'Unknown'
            if geneinfo:
                gene_parts = geneinfo.split(':')
                if gene_parts:
                    gene = gene_parts[0].split('|')[0]
            
            # Extract disease
            disease = info_dict.get('CLNDN', 'Unknown').replace('_', ' ')
            if disease:
                disease = disease.split('|')[0]
            
            # Normalize REF/ALT
            ref_norm, alt_norm = normalize_ref_alt(ref, alt)
            if alt_norm is None:
                continue
            
            # Extract additional info
            allele_id = info_dict.get('ALLELEID', '')
            hgvs = info_dict.get('CLNHGVS', '')
            review_status = info_dict.get('CLNREVSTAT', '').replace('_', ' ')
            
            # Molecular consequence
            mc_raw = info_dict.get('MC', '')
            molecular_consequence = ''
            if mc_raw:
                mc_parts = mc_raw.split('|')
                if len(mc_parts) > 1:
                    molecular_consequence = mc_parts[1].replace('_', ' ')
            
            # Origin
            origin_code = info_dict.get('ORIGIN', '')
            origin_map = {
                '0': 'unknown', '1': 'germline', '2': 'somatic',
                '4': 'inherited', '8': 'paternal', '16': 'maternal',
                '32': 'de-novo', '64': 'biparental', '128': 'uniparental'
            }
            origin = origin_map.get(origin_code, 'unknown')
            
            # Create compact record
            record = {
                'g': gene,  # gene
                'd': disease,  # disease
                'c': chrom,  # chromosome
                'p': pos,  # position
                'r': ref_norm,  # ref
                'a': alt_norm,  # alt
                's': cls,  # significance class
                'v': var_id,  # variation_id
                'al': allele_id,  # allele_id
                'h': hgvs,  # hgvs
                'rv': review_status,  # review_status
                'mc': molecular_consequence,  # molecular_consequence
                'o': origin  # origin
            }
            
            # Add to index for each RS ID
            for rsid in rs_ids:
                index[rsid] = record
            
            processed += 1
            if processed % 10000 == 0:
                print(f"  Processed {processed:,} / {total_lines:,} ({processed*100//total_lines}%)")
    
    print(f"\n✓ Indexed {len(index):,} pathogenic/likely pathogenic variants with rsIDs")
    
    # Write JSON
    print(f"Writing JSON index: {output_json}")
    with open(output_json, 'w', encoding='utf-8') as f:
        json.dump(index, f, separators=(',', ':'))  # Compact format
    
    # Report file size
    size_mb = Path(output_json).stat().st_size / (1024 * 1024)
    print(f"✓ JSON index size: {size_mb:.2f} MB")
    
    return len(index)


def main():
    if len(sys.argv) < 2:
        print("Usage: python build_clinvar_index.py <clinvar.vcf> [output.json]")
        print("\nExample:")
        print("  python build_clinvar_index.py clinvar.vcf clinvar_index.json")
        sys.exit(1)
    
    vcf_path = sys.argv[1]
    output_json = sys.argv[2] if len(sys.argv) > 2 else 'clinvar_index.json'
    
    if not Path(vcf_path).exists():
        print(f"Error: VCF file not found: {vcf_path}")
        sys.exit(1)
    
    build_clinvar_index(vcf_path, output_json)
    print("\n✓ Done! You can now use this JSON file with the web application.")


if __name__ == '__main__':
    main()
