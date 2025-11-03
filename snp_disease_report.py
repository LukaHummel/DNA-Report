#!/usr/bin/env python3
"""
ClinVar VCF Report Generator
Parses ClinVar VCF data and 23andMe results
"""

import csv
import re
import sys
from pathlib import Path
from datetime import datetime
import time


class SNPReportGenerator:
    def __init__(self, clinvar_vcf_path, raw_data_path):
        self.clinvar_vcf = clinvar_vcf_path
        self.raw_data_path = raw_data_path
        self.user_genotypes = {}
        self.findings = []
        self.cv_index = None  # for ClinVar VCF index
        
    class ProgressBar:
        """Minimal stdout progress bar without external deps"""
        def __init__(self, total, prefix=""):
            self.total = max(1, int(total))
            self.prefix = prefix
            self.last_pct = -1
            self.start_time = time.time()

        def update(self, current):
            current = min(current, self.total)
            pct = int((current * 100) / self.total)
            if pct != self.last_pct or current == self.total:
                bar_len = 30
                filled = int(bar_len * current / self.total)
                bar = '█' * filled + '-' * (bar_len - filled)
                elapsed = time.time() - self.start_time
                sys.stdout.write(f"\r{self.prefix} [{bar}] {pct}% ({current}/{self.total})  Elapsed: {int(elapsed)}s")
                sys.stdout.flush()
                self.last_pct = pct
            if current >= self.total:
                sys.stdout.write("\n")

    @staticmethod
    def _genotype_matches_ref_alt(user_geno: str, ref: str, alt: str) -> bool:
        """Return True if the user's genotype is compatible with REF/ALT alleles.
        - Accepts match when both alleles are in {REF} ∪ {ALTs}
        - If REF/ALT missing or ambiguous, returns False (strict filtering)
        - Handles multi-ALT values separated by comma/semicolon/slash/space
        - Normalizes case and allele order (AG == GA)
        """
        if not user_geno:
            return False
        if not ref or not alt:
            # Strict: require both REF and ALT to be present
            return False

        ref = str(ref).strip().upper()
        # ALT may contain multiple values (e.g., "A,G")
        raw_alt = str(alt).strip().upper()
        alt_list = [a for a in re.split(r"[,/;|\s]+", raw_alt) if a]

        # Only allow single-nucleotide alleles (A,C,G,T) for validation
        valid_tokens = {"A","C","G","T"}
        if ref not in valid_tokens:
            return False
        alt_list = [a for a in alt_list if a in valid_tokens]
        if not alt_list:
            return False

        # Build valid allele and alt sets
        alt_set = set(alt_list)
        valid = set([ref]) | alt_set

        # Extract user alleles
        # Typical 23andMe genotypes are two characters like 'AG'.
        # Keep only letters to be safe.
        user_alleles = re.findall(r"[A-Z]", user_geno.upper())
        if len(user_alleles) != 2:
            # Strict: require a diploid SNP call like 'AG'
            return False

        a1, a2 = user_alleles
        both_valid = (a1 in valid) and (a2 in valid)
        any_alt = (a1 in alt_set) or (a2 in alt_set)
        # Require at least one ALT allele present, and both alleles must be valid
        return both_valid and any_alt

    @staticmethod
    def _get_value_case_insensitive(d: dict, candidates):
        if not isinstance(d, dict):
            return None
        keys_map = {str(k).lower(): k for k in d.keys()}
        for name in candidates:
            k = keys_map.get(str(name).lower())
            if k is not None:
                return d[k]
        return None

    @staticmethod
    def _normalize_ref_alt(ref, alt):
        """Normalize REF/ALT values from various formats.
        - Splits tokens on comma/semicolon/slash/pipe/whitespace
        - If REF contains multiple tokens and ALT is empty, assumes first is REF and rest are ALT
        - Filters to single-nucleotide alleles only (A,C,G,T)
        Returns (ref_str, alt_str or None). If normalization fails, returns (None, None).
        """
        def tokens(s):
            return [t for t in re.split(r"[,/;|\s]+", str(s).strip().upper()) if t]

        valid = {"A", "C", "G", "T"}
        ref_parts = tokens(ref) if ref is not None else []
        alt_parts = tokens(alt) if alt is not None else []

        if ref_parts and not alt_parts and len(ref_parts) >= 2:
            # Case like "C/A" in REF and ALT empty
            ref_norm = ref_parts[0]
            alt_parts = ref_parts[1:]
        else:
            ref_norm = ref_parts[0] if ref_parts else None

        if not ref_norm or ref_norm not in valid:
            return None, None

        alt_parts = [a for a in alt_parts if a in valid]
        if not alt_parts:
            return ref_norm, None

        # Deduplicate and stable order
        alt_norm = ",".join(sorted(set(alt_parts)))
        return ref_norm, alt_norm

    @staticmethod
    def _match_status(user_geno: str, ref: str, alt_csv: str) -> str:
        """Return human-friendly match type: Ref/Ref, Ref/Alt, Alt/Alt, or Unknown."""
        if not user_geno or not ref or not alt_csv:
            return "Unknown"
        alleles = re.findall(r"[A-Z]", user_geno.upper())
        if len(alleles) != 2:
            return "Unknown"
        a1, a2 = alleles
        alt_set = set(alt_csv.split(','))
        if a1 in alt_set and a2 in alt_set:
            return "Alt/Alt"
        if (a1 in alt_set and a2 == ref) or (a2 in alt_set and a1 == ref):
            return "Ref/Alt"
        if a1 == ref and a2 == ref:
            return "Ref/Ref"
        return "Other"
        
    def connect_database(self):
        """Connect to ClinVar VCF database"""
        print("✓ Using ClinVar VCF source")
        return True

    def _load_clinvar_vcf(self):
        """Load ClinVar VCF file into an index by rsid"""
        path = self.clinvar_vcf
        print("Indexing ClinVar VCF (this may take a moment)...")
        self.cv_index = {}
        
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
            """Parse CLNSIG field and map to numeric class (4 or 5)
            
            Handles ClinVar CLNSIG values:
            - Pathogenic (pure or with modifiers like low_penetrance)
            - Likely_pathogenic (pure or with modifiers)
            - Pathogenic/Likely_pathogenic (combined)
            
            Excludes:
            - Conflicting_classifications_of_pathogenicity
            - Uncertain_significance
            - Benign/Likely_benign
            - Risk factors, drug response, associations (not disease-causing)
            """
            if not clnsig_text:
                return None
            
            # Normalize: lowercase and remove underscores
            txt = clnsig_text.lower().replace('_', ' ')
            
            # Explicitly exclude non-pathogenic categories
            exclude_terms = [
                'conflicting',
                'uncertain',
                'benign',
                'not provided',
                'no classification',
                'association',  # genetic association, not disease-causing
                'protective',
                'affects',
                'confers sensitivity'
            ]
            
            for term in exclude_terms:
                if term in txt:
                    return None
            
            # Match pathogenic classifications
            # Use exact matching on normalized text to avoid false positives
            
            # Pathogenic (class 5) - including variants like "pathogenic, low penetrance"
            if txt.startswith('pathogenic'):
                # Check it's not "pathogenic/likely pathogenic" (handled separately)
                if 'likely' not in txt:
                    return '5'
                # Handle "pathogenic/likely pathogenic" - treat as pathogenic
                elif txt.startswith('pathogenic/likely pathogenic') or txt.startswith('pathogenic|likely pathogenic'):
                    return '5'
            
            # Likely pathogenic (class 4)
            if txt.startswith('likely pathogenic'):
                return '4'
            
            # Handle pipe-separated values (e.g., "pathogenic|drug_response")
            # Take the first classification if it's pathogenic
            if '|' in txt:
                first_part = txt.split('|')[0].strip()
                if first_part == 'pathogenic':
                    return '5'
                elif first_part == 'likely pathogenic':
                    return '4'
            
            return None
        
        with open(path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                # Skip header lines
                if line.startswith('#'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 8:
                    continue
                
                chrom, pos, var_id, ref, alt, qual, filt, info = parts[:8]
                
                # Parse INFO field
                info_dict = parse_info_field(info)
                
                # Extract RS ID
                rs_raw = info_dict.get('RS')
                if not rs_raw or rs_raw.strip() in ('', '.', '-1'):
                    continue
                
                # Handle multiple RS IDs (comma-separated)
                rs_ids = [f"rs{r.strip()}".lower() for r in rs_raw.split(',') if r.strip()]
                if not rs_ids:
                    continue
                
                # Get clinical significance
                clnsig_raw = info_dict.get('CLNSIG', '')
                cls = parse_clnsig(clnsig_raw)
                if cls not in ('4', '5'):
                    continue
                
                # Extract gene name from GENEINFO (format: "GENE:ID")
                geneinfo = info_dict.get('GENEINFO', '')
                gene = 'Unknown'
                if geneinfo:
                    gene_parts = geneinfo.split(':')
                    if gene_parts:
                        gene = gene_parts[0].split('|')[0]  # Handle multiple genes
                
                # Extract disease name from CLNDN
                disease = info_dict.get('CLNDN', 'Unknown').replace('_', ' ')
                if disease:
                    disease = disease.split('|')[0]  # Take first disease if multiple
                
                # Normalize REF/ALT
                ref_norm, alt_norm = self._normalize_ref_alt(ref, alt)
                if alt_norm is None:
                    continue
                
                # Extract additional clinical information
                allele_id = info_dict.get('ALLELEID', '')
                hgvs = info_dict.get('CLNHGVS', '')
                review_status = info_dict.get('CLNREVSTAT', '').replace('_', ' ')
                
                # Extract molecular consequence (format: "SO:ID|consequence_term")
                mc_raw = info_dict.get('MC', '')
                molecular_consequence = ''
                if mc_raw:
                    mc_parts = mc_raw.split('|')
                    if len(mc_parts) > 1:
                        molecular_consequence = mc_parts[1].replace('_', ' ')
                
                # Extract origin (germline=1, somatic=2, etc.)
                origin_code = info_dict.get('ORIGIN', '')
                origin_map = {'0': 'unknown', '1': 'germline', '2': 'somatic', 
                              '4': 'inherited', '8': 'paternal', '16': 'maternal',
                              '32': 'de-novo', '64': 'biparental', '128': 'uniparental'}
                origin = origin_map.get(origin_code, 'unknown')
                
                # Create record for each RS ID
                for rsid in rs_ids:
                    self.cv_index[rsid] = {
                        'rsid': rsid,
                        'gene': gene,
                        'disease': disease,
                        'chromosome': chrom,
                        'position': pos,
                        'ref': ref_norm,
                        'alt': alt_norm,
                        'clnsig': cls,
                        'significance': 'Pathogenic' if cls == '5' else 'Likely pathogenic',
                        'variation_id': var_id,
                        'allele_id': allele_id,
                        'hgvs': hgvs,
                        'review_status': review_status,
                        'molecular_consequence': molecular_consequence,
                        'origin': origin
                    }
        
        print(f"✓ Indexed {len(self.cv_index)} ClinVar variants with rsIDs from VCF")
        return True
    
    def load_23andme_data(self):
        """Load and parse 23andMe raw data file"""
        print(f"Loading 23andMe data...")
        
        with open(self.raw_data_path, 'r') as f:
            lines = [line for line in f if not line.startswith('#')]
        
        # Detect delimiter
        delimiter = '\t' if '\t' in lines[0] else ','
        
        reader = csv.DictReader(lines, delimiter=delimiter)
        # Estimate total record count (exclude header, blank lines)
        total_rows = sum(1 for l in lines[1:] if l.strip()) if lines else 0
        pb = self.ProgressBar(total_rows or 1, prefix="Parsing 23andMe data")
        processed = 0
        
        for row in reader:
            rsid = row.get('rsid') or row.get('# rsid')
            genotype = row.get('genotype')
            
            if rsid and genotype:
                rsid_lower = rsid.lower()
                # Accept RSIDs starting with 'rs' or 'i' (insertion variants)
                if rsid_lower.startswith('rs') or rsid_lower.startswith('i'):
                    # Normalize genotype (remove hyphens for no-calls)
                    if genotype not in ['--', 'II', 'DD', 'NN']:
                        # Normalize RSID to lowercase for case-insensitive matching
                        self.user_genotypes[rsid_lower] = {
                            'genotype': genotype,
                            'chromosome': row.get('chromosome'),
                            'position': row.get('position')
                        }
            processed += 1
            pb.update(processed)
        
        print(f"✓ Loaded {len(self.user_genotypes)} SNPs from your data")
        return len(self.user_genotypes) > 0
    
    def query_clinvar(self):
        """Query ClinVar VCF for matching variants"""
        return self._query_clinvar_vcf()

    def _query_clinvar_vcf(self):
        """Query ClinVar VCF index for matches"""
        if self.cv_index is None:
            self._load_clinvar_vcf()
        print("Scanning ClinVar VCF for matches...")
        total = len(self.user_genotypes)
        pb = self.ProgressBar(total or 1, prefix="Scanning VCF          ")
        processed = 0
        for rsid, user_data in self.user_genotypes.items():
            rec = self.cv_index.get(rsid)
            if rec:
                # genotype match
                if self._genotype_matches_ref_alt(user_data.get('genotype'), rec['ref'], rec['alt']):
                    finding = {
                        'rsid': rec['rsid'],
                        'user_genotype': user_data.get('genotype'),
                        'chromosome': user_data.get('chromosome') or rec.get('chromosome'),
                        'position': user_data.get('position') or rec.get('position'),
                        'gene': rec.get('gene', 'Unknown'),
                        'disease': rec.get('disease', 'Unknown'),
                        'clnsig': rec['clnsig'],
                        'significance': rec['significance'],
                        'ref': rec['ref'],
                        'alt': rec['alt'],
                        'match': self._match_status(user_data.get('genotype'), rec['ref'], rec['alt']),
                        'variation_id': rec.get('variation_id'),
                        'allele_id': rec.get('allele_id'),
                        'hgvs': rec.get('hgvs'),
                        'review_status': rec.get('review_status'),
                        'molecular_consequence': rec.get('molecular_consequence'),
                        'origin': rec.get('origin')
                    }
                    self.findings.append(finding)
            processed += 1
            pb.update(processed)
        print(f"✓ Found {len(self.findings)} disease-associated variants")
        return len(self.findings)
    
    def generate_html_report(self, output_path):
        """Generate comprehensive HTML report"""
        print(f"Generating report...")
        
        # Sort by clinical significance
        pathogenic = [f for f in self.findings if f['clnsig'] == '5']
        likely_pathogenic = [f for f in self.findings if f['clnsig'] == '4']
        
        html_header = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Genomic Disease Report</title>
    <style>
        :root {{
            --color-pathogenic: #c0152f;
            --color-likely: #e67e22;
            --color-border: #ecf0f1;
            --color-text: #2c3e50;
            --color-bg: #ecf0f1;
            --color-surface: #ffffff;
        }}
        
        * {{ box-sizing: border-box; }}
        
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            line-height: 1.6;
            color: var(--color-text);
            background: var(--color-bg);
            margin: 0;
            padding: 20px;
        }}
        
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background: var(--color-surface);
            padding: 40px;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        
        h1 {{
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 15px;
            margin-bottom: 30px;
        }}
        
        h2 {{
            color: #34495e;
            margin-top: 40px;
            margin-bottom: 20px;
            font-size: 20px;
        }}
        
        .summary {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-bottom: 40px;
        }}
        
        .summary-card {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            text-align: center;
        }}
        
        .summary-card.pathogenic {{
            background: linear-gradient(135deg, #c0152f 0%, #8b0000 100%);
        }}
        
        .summary-card.likely {{
            background: linear-gradient(135deg, #e67e22 0%, #d35400 100%);
        }}
        
        .summary-number {{
            font-size: 32px;
            font-weight: bold;
            margin-bottom: 5px;
        }}
        
        .summary-label {{
            font-size: 14px;
            opacity: 0.9;
        }}
        
        .finding-card {{
            background: white;
            border: 1px solid var(--color-border);
            border-left: 5px solid #3498db;
            padding: 20px;
            margin-bottom: 20px;
            border-radius: 6px;
        }}
        
        .finding-card.pathogenic {{
            border-left-color: var(--color-pathogenic);
            background: rgba(192, 21, 47, 0.02);
        }}
        
        .finding-card.likely {{
            border-left-color: var(--color-likely);
            background: rgba(230, 126, 34, 0.02);
        }}
        
        .finding-header {{
            display: flex;
            justify-content: space-between;
            align-items: start;
            margin-bottom: 15px;
        }}
        
        .finding-title {{
            font-size: 18px;
            font-weight: 600;
        }}
        
        .badge {{
            padding: 6px 12px;
            border-radius: 4px;
            font-size: 12px;
            font-weight: 600;
            text-transform: uppercase;
        }}
        
        .badge-pathogenic {{
            background: var(--color-pathogenic);
            color: white;
        }}
        
        .badge-likely {{
            background: var(--color-likely);
            color: white;
        }}
        
        .finding-details {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
            gap: 15px;
            margin: 15px 0;
        }}
        
        .detail-item {{
            font-size: 13px;
        }}
        
        .detail-label {{
            color: #7f8c8d;
            font-weight: 600;
            text-transform: uppercase;
            font-size: 11px;
            margin-bottom: 3px;
        }}
        
        .detail-value {{
            color: var(--color-text);
            font-family: 'Courier New', monospace;
        }}
        
        .disease-info {{
            background: #f8f9fa;
            padding: 15px;
            border-radius: 4px;
            margin: 15px 0;
        }}
        
        .disease-name {{
            font-weight: 600;
            color: var(--color-pathogenic);
            margin-bottom: 8px;
        }}
        
        .disease-description {{
            font-size: 13px;
            line-height: 1.5;
        }}
        
        .disclaimer {{
            background: #fff3cd;
            border: 1px solid #ffc107;
            padding: 20px;
            border-radius: 6px;
            margin-top: 40px;
        }}
        
        .disclaimer h3 {{
            color: #856404;
            margin-top: 0;
        }}
        
        .disclaimer p {{
            color: #856404;
            margin: 10px 0;
        }}
        
        .external-link {{
            color: #3498db;
            text-decoration: none;
            font-size: 12px;
        }}
        
        .external-link:hover {{
            text-decoration: underline;
        }}
        
        .timestamp {{
            text-align: right;
            font-size: 12px;
            color: #7f8c8d;
            margin-top: 30px;
            border-top: 1px solid var(--color-border);
            padding-top: 15px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 Genomic Disease Risk Report</h1>
        
        <div class="disclaimer">
            <h3>⚠️ Important Medical Disclaimer</h3>
            <p><strong>This report is for educational and informational purposes only.</strong></p>
            <p>This is NOT medical advice and should NOT be used for self-diagnosis or treatment decisions. 
            Many genetic variants have incomplete penetrance or require specific environmental factors to manifest. 
            Please consult with a qualified healthcare provider, genetic counselor, or medical geneticist to:</p>
            <ul>
                <li>Properly interpret these findings in your personal context</li>
                <li>Determine clinical relevance and actionability</li>
                <li>Discuss screening, prevention, or management options</li>
                <li>Understand inheritance patterns and family implications</li>
            </ul>
        </div>
        
        <div class="summary">
            <div class="summary-card pathogenic">
                <div class="summary-number">{len(pathogenic)}</div>
                <div class="summary-label">Pathogenic Variants</div>
            </div>
            <div class="summary-card likely">
                <div class="summary-number">{len(likely_pathogenic)}</div>
                <div class="summary-label">Likely Pathogenic Variants</div>
            </div>
            <div class="summary-card">
                <div class="summary-number">{len(self.findings)}</div>
                <div class="summary-label">Total Disease Associations</div>
            </div>
        </div>
"""
        
        # Add pathogenic findings
        html_body = ""
        if pathogenic:
            html_body += '<h2>🔴 Pathogenic Variants</h2>\n'
            for finding in sorted(pathogenic, key=lambda x: x['disease']):
                html_body += self._format_finding_card(finding, 'pathogenic')
        
        # Add likely pathogenic findings
        if likely_pathogenic:
            html_body += '<h2>🟠 Likely Pathogenic Variants</h2>\n'
            for finding in sorted(likely_pathogenic, key=lambda x: x['disease']):
                html_body += self._format_finding_card(finding, 'likely')
        
        html_footer = f"""
        <div class="timestamp">
            Report generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} CET
            <br>
            Based on ClinVar VCF database
        </div>
    </div>
</body>
</html>
"""
        
        html = html_header + html_body + html_footer
        
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html)
        
        print(f"✓ Report saved to: {output_path}")
        return True
    
    def _format_finding_card(self, finding, priority):
        """Format a disease finding card"""
        rsid = finding['rsid']
        gene = finding['gene']
        disease = finding['disease']
        genotype = finding['user_genotype']
        ref = finding.get('ref', '?')
        alt = finding.get('alt', '?')
        match = finding.get('match', 'Unknown')
        variation_id = finding.get('variation_id', '')
        molecular_consequence = finding.get('molecular_consequence', '')
        origin = finding.get('origin', '')
        review_status = finding.get('review_status', '')
        
        snpedia_link = f'https://www.snpedia.com/index.php/{rsid}'
        clinvar_link = f'https://www.ncbi.nlm.nih.gov/clinvar/variation/{variation_id}/' if variation_id else ''
        
        # Build optional fields
        mol_effect_html = f'''<div class="detail-item">
                    <div class="detail-label">Molecular Effect</div>
                    <div class="detail-value">{molecular_consequence}</div>
                </div>''' if molecular_consequence else ''
        
        origin_html = f'''<div class="detail-item">
                    <div class="detail-label">Origin</div>
                    <div class="detail-value">{origin.capitalize()}</div>
                </div>''' if origin else ''
        
        review_html = f'''<div class="detail-item">
                    <div class="detail-label">Review Status</div>
                    <div class="detail-value" style="font-size: 10px;">{review_status}</div>
                </div>''' if review_status else ''
        
        clinvar_link_html = f' | <a href="{clinvar_link}" target="_blank" class="external-link">View on ClinVar →</a>' if clinvar_link else ''
        clinvar_id_html = f' | ClinVar ID: {variation_id}' if variation_id and not clinvar_link else ''
        
        html = f"""    <div class="finding-card {priority}">
            <div class="finding-header">
                <div>
                    <div class="finding-title">{gene}: {disease}</div>
                </div>
                <span class="badge badge-{priority}">{finding['significance']}</span>
            </div>
            
            <div class="finding-details">
                <div class="detail-item">
                    <div class="detail-label">SNP ID</div>
                    <div class="detail-value">{rsid}</div>
                </div>
                <div class="detail-item">
                    <div class="detail-label">Your Genotype</div>
                    <div class="detail-value">{genotype}</div>
                </div>
                <div class="detail-item">
                    <div class="detail-label">Reference/Alt</div>
                    <div class="detail-value">{ref} / {alt.replace(',', ', ') if isinstance(alt, str) else alt}</div>
                </div>
                <div class="detail-item">
                    <div class="detail-label">Position</div>
                    <div class="detail-value">Chr {finding['chromosome']}:{finding['position']}</div>
                </div>
                <div class="detail-item">
                    <div class="detail-label">Genotype Match</div>
                    <div class="detail-value">{match}</div>
                </div>
                {mol_effect_html}
                {origin_html}
                {review_html}
            </div>
            
            <div class="disease-info">
                <div class="disease-name">Associated Disease/Condition:</div>
                <div class="disease-description">{disease}</div>
            </div>
            
            <div>
                <a href="{snpedia_link}" target="_blank" class="external-link">View on SNPedia →</a>{clinvar_link_html}{clinvar_id_html}
            </div>
        </div>
"""
        return html
    
    
    
    def run(self, output_path='genomics_disease_report.html'):
        """Run complete analysis"""
        print("=" * 70)
        print("ClinVar VCF Disease Report Generator")
        print("=" * 70)
        
        if not self.connect_database():
            return False
        
        if not self.load_23andme_data():
            return False
        
        self.query_clinvar()
        
        if self.findings:
            self.generate_html_report(output_path)
        else:
            print("\n✓ No disease-associated variants found in your data")
        
        print("\n" + "=" * 70)
        print("✓ Analysis complete!")
        print("=" * 70)
        return True


def main():
    if len(sys.argv) < 3:
        print("Usage: python snp_disease_report.py <clinvar.vcf> <23andme_raw_data.txt> [output.html]")
        print("\nExample:")
        print("  python snp_disease_report.py clinvar.vcf genome_results.txt report.html")
        sys.exit(1)
    
    clinvar_vcf = sys.argv[1]
    raw_data = sys.argv[2]
    output = sys.argv[3] if len(sys.argv) > 3 else 'genomics_disease_report.html'
    
    if not Path(clinvar_vcf).exists():
        print(f"Error: VCF file not found: {clinvar_vcf}")
        sys.exit(1)
    
    if not Path(raw_data).exists():
        print(f"Error: Data file not found: {raw_data}")
        sys.exit(1)
    
    generator = SNPReportGenerator(clinvar_vcf, raw_data)
    generator.run(output)


if __name__ == '__main__':
    main()
