# SNPedia Disease Report Generator - User Guide

## Quick Start

### Installation Requirements
```bash
pip install -r requirements.txt
```

Or just Python 3.7+ with standard library (no external dependencies needed).

### Running the Script

```bash
python snp_disease_report.py <snpedia.db> <23andme_raw_data.txt> [output.html]
```

**Example:**
```bash
python snp_disease_report.py snpedia.db genome_results.txt my_report.html
```

### Input Files

#### 1. SNPedia Database (.db file)
- SQLite database with table structure: `rsid | content | scraped_at | attribution`
- Content field contains wiki-markup with:
  - `{{Rsnum}}` section: Basic SNP information
  - `{{ClinVar}}` section: Clinical significance data

#### 2. 23andMe Raw Data (.txt file)
- Tab or comma-separated format
- Required columns: `rsid`, `genotype`, `chromosome`, `position`
- Example:
  ```
  # rsid  chromosome  position  genotype
  rs4988235  2  136608646  AA
  rs797044442  16  47502989  AG
  ```

### Output

Generates an HTML report (`genomics_disease_report.html` by default) with:

- **Summary Statistics**: Count of pathogenic and likely pathogenic variants
- **Pathogenic Variants Section**: Disease-causing variants (CLNSIG = 5)
- **Likely Pathogenic Section**: Probably disease-causing variants (CLNSIG = 4)
- **For Each Finding**:
  - SNP ID (rsid)
  - Your genotype
  - Reference/Alternative alleles
  - Chromosomal position
  - Associated disease/condition
  - Links to SNPedia and ClinVar for detailed information

## Clinical Significance Levels (CLNSIG)

The script filters for pathogenic findings using ClinVar classifications:

- **5 = Pathogenic**: Strong evidence for disease causation
- **4 = Likely Pathogenic**: Probable disease causation
- **3 = Likely Benign**: Probably not disease-causing (filtered out)
- **2 = Benign**: Not associated with disease (filtered out)
- **0 = Uncertain Significance**: Unknown effect (filtered out)
- **1 = Not Provided**: Insufficient data (filtered out)

## Understanding the Report

### Pathogenic Variants (Red 🔴)
These variants have strong evidence of causing disease. If you carry these:
- Consider contacting a genetic counselor
- May warrant medical screening
- Family members may be affected

### Likely Pathogenic Variants (Orange 🟠)
These variants probably cause disease, but evidence is less definitive. Same recommendations as above apply.

## Important Notes

### Medical Disclaimer
**This report is for educational purposes only and is NOT medical advice.** 

- Not all carriers of pathogenic variants will develop disease (incomplete penetrance)
- Environmental and other genetic factors may be required
- Many variants have age-dependent or variable expression
- Consult a healthcare provider before making any decisions

### Limitations

1. **Database Coverage**: Only shows variants in your SNPedia database
2. **Population Bias**: Most genomic databases overrepresent European ancestry
3. **Incomplete Data**: Not all disease variants are yet characterized
4. **Strand Orientation**: May miss variants reported on opposite DNA strand
5. **Genotype Matching**: Requires exact match between your data and database

### What's Not Included

- Carrier status for recessive conditions
- Pharmacogenomics (drug response)
- Ancestry and population markers
- Complex disease susceptibility (requires multiple variants)
- Risk scores
- Polygenic predictions

## Troubleshooting

### No matches found
1. Verify SNPedia database contains your rsids
2. Check 23andMe file is properly formatted
3. Ensure rsid format: `rs` followed by numbers (e.g., `rs4988235`)

### Few matches compared to expected
1. Your SNPedia database might be incomplete
2. Many variants may have "uncertain significance" (filtered out)
3. Most variants are not disease-causing

### Database Connection Error
1. Verify `.db` file path is correct
2. Confirm file is a valid SQLite database: `sqlite3 snpedia.db ".tables"`

## Advanced Usage

### Checking Database Structure
```bash
sqlite3 snpedia.db ".schema"
```

### Exporting specific SNPs
Modify the script to include specific gene filtering:
```python
# In query_snpedia method, add gene filter:
if finding['gene'] in ['BRCA1', 'BRCA2', 'TP53']:
    self.findings.append(finding)
```

### Filtering by disease category
Add custom disease filtering to focus on specific conditions:
```python
disease_keywords = ['cancer', 'cardiovascular', 'metabolic']
for keyword in disease_keywords:
    if keyword.lower() in disease.lower():
        self.findings.append(finding)
```

## Data Files Provided

- `snp_disease_report.py`: Main script to generate reports
- This guide file

## Support

For issues with:
- **SNPedia database**: https://www.snpedia.com
- **ClinVar data**: https://www.ncbi.nlm.nih.gov/clinvar/
- **23andMe format**: Contact 23andMe support or check their documentation
