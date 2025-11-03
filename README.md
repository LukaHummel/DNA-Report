# DNA Report - Privacy-First Genomic Analysis

A browser-based tool for analyzing 23andMe genetic data against the ClinVar database to identify disease-associated variants. **All processing happens locally in your browser** - your genetic data never leaves your device.

## 🔒 Privacy First

- **100% Client-Side Processing**: Your genetic data is processed entirely in your browser
- **No Data Upload**: Files are read locally and never sent to any server
- **No Tracking**: No analytics, cookies, or third-party scripts (If using the GitHub page; GitHub saves your IP)
- **Open Source**: Full transparency - review the code yourself

## ✨ Features

- 🧬 Analyze 23andMe raw data files
- 🔍 Match against ClinVar pathogenic/likely pathogenic variants
- 📊 Beautiful, comprehensive HTML reports
- 🖨️ Print and download your results
- 📱 Responsive design for mobile and desktop

## 🚀 Quick Start

### Option 1: Use GitHub Pages (Recommended)

Visit: **[https://lukahummel.github.io/DNA-Report](https://lukahummel.github.io/DNA-Report)**

1. Click "Select your 23andMe file"
2. Upload your raw genetic data file
3. Wait for analysis to complete
4. Review your personalized report

### Option 2: Run Locally

1. Clone this repository:
   ```bash
   git clone https://github.com/LukaHummel/DNA-Report.git
   cd DNA-Report
   ```

2. Generate the ClinVar index (one-time setup):
   ```bash
   # Download ClinVar VCF (if you don't have it)
   # Then run:
   python build_clinvar_index.py clinvar.vcf clinvar_index.json
   ```

3. Serve locally:
   ```bash
   # Using Python
   python -m http.server 8000
   
   # Or using Node.js
   npx http-server
   ```

4. Open `http://localhost:8000` in your browser

## 📋 Requirements

### For Users (Browser-Based):
- Modern web browser (Chrome, Firefox, Safari, Edge)
- Your 23andMe raw data file

### For Developers (Building ClinVar Index):
- Python 3.7+
- ClinVar VCF file (download from [NCBI](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/))

## 🏗️ How It Works

1. **Pre-processing** (one-time): The large ClinVar VCF file is converted to an optimized JSON index containing only pathogenic/likely pathogenic variants with rsIDs
2. **Upload**: User uploads their 23andMe file (never sent anywhere)
3. **Parse**: JavaScript parses the genetic data in-browser
4. **Match**: Variants are matched against the ClinVar index
5. **Report**: Results are displayed in a comprehensive, printable report

## 📁 File Structure

```
DNA-Report/
├── index.html              # Main application interface
├── app.js                  # Core JavaScript logic
├── styles.css              # Styling and responsive design
├── build_clinvar_index.py  # Python script to build ClinVar index
├── clinvar_index.json      # Pre-built ClinVar database (generated)
└── README.md               # This file
```

## ⚠️ Medical Disclaimer

**This tool is for educational and informational purposes only.**

This is NOT medical advice and should NOT be used for:
- Self-diagnosis
- Treatment decisions
- Clinical decision-making

Many genetic variants have:
- Incomplete penetrance (may never cause disease)
- Variable expressivity (different severity in different people)
- Environmental dependencies (require other factors to manifest)

**Always consult with a qualified healthcare provider, genetic counselor, or medical geneticist** to properly interpret any genetic findings.

## 🔧 Building the ClinVar Index

The ClinVar VCF file is too large for browsers (~1-2GB). We convert it to an optimized JSON index:

```bash
# Download latest ClinVar VCF
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
gunzip clinvar.vcf.gz

# Build index (takes a few minutes)
python build_clinvar_index.py clinvar.vcf clinvar_index.json
```

The resulting JSON file is much smaller (~10-50MB) and contains only:
- Pathogenic and likely pathogenic variants
- Variants with rsIDs (matching 23andMe format)
- Essential clinical information

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

### Development

```bash
# Clone and setup
git clone https://github.com/LukaHummel/DNA-Report.git
cd DNA-Report

# Make changes to HTML/CSS/JS
# Test locally with a web server

# Submit PR
```

## 📊 Data Sources

- **ClinVar**: NCBI's public archive of genetic variants and their clinical significance
- **SNPedia**: Community-curated wiki for SNP information (linked in reports)

## 📝 License

MIT License - see LICENSE file for details

## 🔗 Links

- [Live Demo](https://lukahummel.github.io/DNA-Report)
- [Report Issues](https://github.com/LukaHummel/DNA-Report/issues)
- [ClinVar Database](https://www.ncbi.nlm.nih.gov/clinvar/)
- [23andMe Raw Data](https://customercare.23andme.com/hc/en-us/articles/360004944654-What-s-In-Your-Account-Settings)
