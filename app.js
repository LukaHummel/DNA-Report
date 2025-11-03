/**
 * DNA Report - Browser-based Genomic Analysis
 * All processing happens client-side for complete privacy
 */

class DNAReportApp {
    constructor() {
        this.clinvarIndex = null;
        this.userGenotypes = new Map();
        this.findings = [];
        
        this.initEventListeners();
    }
    
    initEventListeners() {
        const fileInput = document.getElementById('fileInput');
        const uploadArea = document.getElementById('uploadArea');
        
        fileInput.addEventListener('change', (e) => this.handleFileSelect(e));
        
        // Drag and drop support
        uploadArea.addEventListener('dragover', (e) => {
            e.preventDefault();
            uploadArea.classList.add('drag-over');
        });
        
        uploadArea.addEventListener('dragleave', () => {
            uploadArea.classList.remove('drag-over');
        });
        
        uploadArea.addEventListener('drop', (e) => {
            e.preventDefault();
            uploadArea.classList.remove('drag-over');
            const file = e.dataTransfer.files[0];
            if (file) {
                this.processFile(file);
            }
        });
    }
    
    handleFileSelect(event) {
        const file = event.target.files[0];
        if (file) {
            this.processFile(file);
        }
    }
    
    async processFile(file) {
        // Show file info
        const fileInfo = document.getElementById('fileInfo');
        fileInfo.innerHTML = `
            <strong>File loaded:</strong> ${file.name} (${(file.size / 1024 / 1024).toFixed(2)} MB)
        `;
        fileInfo.style.display = 'block';
        
        // Show processing section
        document.getElementById('uploadSection').style.display = 'none';
        document.getElementById('processingSection').style.display = 'block';
        
        try {
            // Load ClinVar index
            await this.loadClinvarIndex();
            
            // Parse 23andMe file
            await this.parse23andMeFile(file);
            
            // Match variants
            await this.matchVariants();
            
            // Generate report
            this.generateReport();
            
        } catch (error) {
            this.showError(error.message);
        }
    }
    
    async loadClinvarIndex() {
        this.updateProgress(10, 'Loading ClinVar database...');
        
        try {
            const response = await fetch('clinvar_index.json');
            if (!response.ok) {
                throw new Error('Failed to load ClinVar database. Please ensure clinvar_index.json is in the same directory.');
            }
            
            this.clinvarIndex = await response.json();
            this.updateProgress(30, `✓ Loaded ${Object.keys(this.clinvarIndex).length.toLocaleString()} ClinVar variants`);
            
        } catch (error) {
            throw new Error('Error loading ClinVar database: ' + error.message);
        }
    }
    
    async parse23andMeFile(file) {
        this.updateProgress(40, 'Parsing your genetic data...');
        
        return new Promise((resolve, reject) => {
            const reader = new FileReader();
            
            reader.onload = (e) => {
                try {
                    const text = e.target.result;
                    const lines = text.split('\n');
                    
                    // Find header line: handle both "rsid..." and "# rsid..."
                    let headerIndex = -1;
                    let headerLine = '';
                    for (let i = 0; i < lines.length; i++) {
                        const raw = lines[i];
                        if (!raw) continue;
                        const trimmed = raw.replace(/\r$/, '').trim();
                        if (!trimmed) continue;
                        const noHash = trimmed.replace(/^#+\s*/, '');
                        const lower = noHash.toLowerCase();
                        if (lower.startsWith('rsid') && (lower.includes('genotype') || lower.includes('chromosome'))) {
                            headerIndex = i;
                            headerLine = noHash; // header without leading '#'
                            break;
                        }
                    }
                    
                    if (headerIndex === -1) {
                        throw new Error('No valid header found. Expecting a line starting with "rsid" (possibly prefixed by #).');
                    }
                    
                    // Detect delimiter: prefer tab, else comma, else whitespace
                    let delimiter = '\t';
                    if (headerLine.includes('\t')) {
                        delimiter = '\t';
                    } else if (headerLine.includes(',')) {
                        delimiter = ',';
                    } else {
                        delimiter = null; // use regex split on whitespace
                    }
                    
                    // Normalize headers: strip quotes, '#', trim and lowercase
                    const headers = (delimiter ? headerLine.split(delimiter) : headerLine.split(/\s+/))
                        .map(h => h.replace(/^#+\s*/, '').replace(/^["']|["']$/g, '').trim().toLowerCase());
                    
                    // Find column indices (support '# rsid' already normalized)
                    const rsidIndex = headers.findIndex(h => h === 'rsid' || h.endsWith('rsid'));
                    const genotypeIndex = headers.findIndex(h => h.includes('genotype'));
                    const chromosomeIndex = headers.findIndex(h => h.includes('chromosome'));
                    const positionIndex = headers.findIndex(h => h.includes('position'));
                    
                    if (rsidIndex === -1 || genotypeIndex === -1) {
                        throw new Error('Invalid file format: missing required columns (rsid, genotype). Make sure you uploaded the unzipped 23andMe raw data file.');
                    }
                    
                    // Parse data rows
                    let validCount = 0;
                    for (let i = headerIndex + 1; i < lines.length; i++) {
                        let line = lines[i];
                        if (!line) continue;
                        line = line.replace(/\r$/, '').trim();
                        if (!line) continue;
                        if (line.startsWith('#')) continue; // skip comment lines after header if any
                        
                        const cols = (delimiter ? line.split(delimiter) : line.split(/\s+/));
                        const rsidRaw = cols[rsidIndex] ?? '';
                        const genotypeRaw = cols[genotypeIndex] ?? '';
                        const rsid = rsidRaw.trim().toLowerCase();
                        const genotype = genotypeRaw.trim();
                        
                        if (rsid && genotype && (rsid.startsWith('rs') || rsid.startsWith('i'))) {
                            // Skip no-calls
                            if (!['--', 'II', 'DD', 'NN'].includes(genotype)) {
                                const chromosome = (chromosomeIndex !== -1 ? (cols[chromosomeIndex] || '').trim() : '');
                                const position = (positionIndex !== -1 ? (cols[positionIndex] || '').trim() : '');
                                this.userGenotypes.set(rsid, { genotype, chromosome, position });
                                validCount++;
                            }
                        }
                        
                        // Update progress periodically
                        if ((i - headerIndex) % 10000 === 0) {
                            const progress = 40 + Math.floor(((i - headerIndex) / Math.max(1, (lines.length - headerIndex))) * 20);
                            this.updateProgress(progress, `Parsing... ${validCount.toLocaleString()} SNPs loaded`);
                        }
                    }
                    
                    this.updateProgress(60, `✓ Loaded ${validCount.toLocaleString()} SNPs from your data`);
                    resolve();
                    
                } catch (error) {
                    reject(error);
                }
            };
            
            reader.onerror = () => reject(new Error('Failed to read file'));
            reader.readAsText(file);
        });
    }
    
    async matchVariants() {
        this.updateProgress(70, 'Matching against ClinVar database...');
        
        let matched = 0;
        let processed = 0;
        const total = this.userGenotypes.size;
        
        for (const [rsid, userData] of this.userGenotypes) {
            const clinvarRecord = this.clinvarIndex[rsid];
            
            if (clinvarRecord) {
                // Check if genotype matches ref/alt
                if (this.genotypeMatchesRefAlt(userData.genotype, clinvarRecord.r, clinvarRecord.a)) {
                    this.findings.push({
                        rsid: rsid,
                        userGenotype: userData.genotype,
                        chromosome: userData.chromosome || clinvarRecord.c,
                        position: userData.position || clinvarRecord.p,
                        gene: clinvarRecord.g,
                        disease: clinvarRecord.d,
                        clnsig: clinvarRecord.s,
                        significance: clinvarRecord.s === '5' ? 'Pathogenic' : 'Likely pathogenic',
                        ref: clinvarRecord.r,
                        alt: clinvarRecord.a,
                        match: this.getMatchStatus(userData.genotype, clinvarRecord.r, clinvarRecord.a),
                        variationId: clinvarRecord.v,
                        alleleId: clinvarRecord.al,
                        hgvs: clinvarRecord.h,
                        reviewStatus: clinvarRecord.rv,
                        molecularConsequence: clinvarRecord.mc,
                        origin: clinvarRecord.o
                    });
                    matched++;
                }
            }
            
            processed++;
            if (processed % 1000 === 0) {
                const progress = 70 + Math.floor((processed / total) * 20);
                this.updateProgress(progress, `Scanning... ${matched} matches found`);
            }
        }
        
        this.updateProgress(95, `✓ Found ${matched} disease-associated variants`);
    }
    
    genotypeMatchesRefAlt(userGeno, ref, alt) {
        if (!userGeno || !ref || !alt) return false;
        
        const validTokens = new Set(['A', 'C', 'G', 'T']);
        
        // Validate ref
        if (!validTokens.has(ref.toUpperCase())) return false;
        
        // Parse alt (may be comma-separated)
        const altList = alt.split(',').filter(a => validTokens.has(a.toUpperCase()));
        if (altList.length === 0) return false;
        
        const altSet = new Set(altList.map(a => a.toUpperCase()));
        const validAlleles = new Set([ref.toUpperCase(), ...altSet]);
        
        // Extract user alleles
        const userAlleles = userGeno.toUpperCase().match(/[ACGT]/g);
        if (!userAlleles || userAlleles.length !== 2) return false;
        
        const [a1, a2] = userAlleles;
        const bothValid = validAlleles.has(a1) && validAlleles.has(a2);
        const anyAlt = altSet.has(a1) || altSet.has(a2);
        
        return bothValid && anyAlt;
    }
    
    getMatchStatus(userGeno, ref, alt) {
        if (!userGeno || !ref || !alt) return 'Unknown';
        
        const alleles = userGeno.toUpperCase().match(/[ACGT]/g);
        if (!alleles || alleles.length !== 2) return 'Unknown';
        
        const [a1, a2] = alleles;
        const altSet = new Set(alt.split(',').map(a => a.toUpperCase()));
        const refUpper = ref.toUpperCase();
        
        if (altSet.has(a1) && altSet.has(a2)) return 'Alt/Alt';
        if ((altSet.has(a1) && a2 === refUpper) || (altSet.has(a2) && a1 === refUpper)) return 'Ref/Alt';
        if (a1 === refUpper && a2 === refUpper) return 'Ref/Ref';
        return 'Other';
    }
    
    generateReport() {
        this.updateProgress(100, '✓ Analysis complete!');
        
        // Hide processing, show report
        setTimeout(() => {
            document.getElementById('processingSection').style.display = 'none';
            document.getElementById('reportSection').style.display = 'block';
            
            const reportContent = document.getElementById('reportContent');
            
            if (this.findings.length === 0) {
                reportContent.innerHTML = `
                    <div class="no-findings">
                        <h2>✓ No Disease-Associated Variants Found</h2>
                        <p>Your genetic data shows no pathogenic or likely pathogenic variants in the ClinVar database.</p>
                        <p class="note">This doesn't mean you have no genetic risks - it only reflects variants documented in ClinVar with rsIDs present in your 23andMe data.</p>
                    </div>
                `;
                return;
            }
            
            // Separate by severity
            const pathogenic = this.findings.filter(f => f.clnsig === '5');
            const likelyPathogenic = this.findings.filter(f => f.clnsig === '4');
            
            reportContent.innerHTML = `
                <h1>🧬 Genomic Disease Risk Report</h1>
                
                <div class="summary">
                    <div class="summary-card pathogenic">
                        <div class="summary-number">${pathogenic.length}</div>
                        <div class="summary-label">Pathogenic Variants</div>
                    </div>
                    <div class="summary-card likely">
                        <div class="summary-number">${likelyPathogenic.length}</div>
                        <div class="summary-label">Likely Pathogenic Variants</div>
                    </div>
                    <div class="summary-card">
                        <div class="summary-number">${this.findings.length}</div>
                        <div class="summary-label">Total Disease Associations</div>
                    </div>
                </div>
                
                ${pathogenic.length > 0 ? `
                    <h2>🔴 Pathogenic Variants</h2>
                    ${pathogenic.sort((a, b) => a.disease.localeCompare(b.disease))
                        .map(f => this.formatFindingCard(f, 'pathogenic')).join('')}
                ` : ''}
                
                ${likelyPathogenic.length > 0 ? `
                    <h2>🟠 Likely Pathogenic Variants</h2>
                    ${likelyPathogenic.sort((a, b) => a.disease.localeCompare(b.disease))
                        .map(f => this.formatFindingCard(f, 'likely')).join('')}
                ` : ''}
                
                <div class="timestamp">
                    Report generated: ${new Date().toLocaleString('en-US', { 
                        year: 'numeric', 
                        month: 'long', 
                        day: 'numeric', 
                        hour: '2-digit', 
                        minute: '2-digit' 
                    })}
                    <br>
                    Based on ClinVar database
                    <br><br>
                    <button onclick="window.print()" class="print-button">🖨️ Print Report</button>
                    <button onclick="app.downloadReport()" class="download-button">💾 Download HTML</button>
                </div>
            `;
        }, 500);
    }
    
    formatFindingCard(finding, priority) {
        const snpediaLink = `https://www.snpedia.com/index.php/${finding.rsid}`;
        const clinvarLink = finding.variationId 
            ? `https://www.ncbi.nlm.nih.gov/clinvar/variation/${finding.variationId}/` 
            : '';
        
        const molEffectHtml = finding.molecularConsequence 
            ? `<div class="detail-item">
                <div class="detail-label">Molecular Effect</div>
                <div class="detail-value">${finding.molecularConsequence}</div>
              </div>` 
            : '';
        
        const originHtml = finding.origin && finding.origin !== 'unknown'
            ? `<div class="detail-item">
                <div class="detail-label">Origin</div>
                <div class="detail-value">${finding.origin.charAt(0).toUpperCase() + finding.origin.slice(1)}</div>
              </div>`
            : '';
        
        const reviewHtml = finding.reviewStatus
            ? `<div class="detail-item">
                <div class="detail-label">Review Status</div>
                <div class="detail-value" style="font-size: 10px;">${finding.reviewStatus}</div>
              </div>`
            : '';
        
        const clinvarLinkHtml = clinvarLink 
            ? ` | <a href="${clinvarLink}" target="_blank" class="external-link">View on ClinVar →</a>`
            : '';
        
        return `
            <div class="finding-card ${priority}">
                <div class="finding-header">
                    <div>
                        <div class="finding-title">${finding.gene}: ${finding.disease}</div>
                    </div>
                    <span class="badge badge-${priority}">${finding.significance}</span>
                </div>
                
                <div class="finding-details">
                    <div class="detail-item">
                        <div class="detail-label">SNP ID</div>
                        <div class="detail-value">${finding.rsid}</div>
                    </div>
                    <div class="detail-item">
                        <div class="detail-label">Your Genotype</div>
                        <div class="detail-value">${finding.userGenotype}</div>
                    </div>
                    <div class="detail-item">
                        <div class="detail-label">Reference/Alt</div>
                        <div class="detail-value">${finding.ref} / ${finding.alt}</div>
                    </div>
                    <div class="detail-item">
                        <div class="detail-label">Position</div>
                        <div class="detail-value">Chr ${finding.chromosome}:${finding.position}</div>
                    </div>
                    <div class="detail-item">
                        <div class="detail-label">Genotype Match</div>
                        <div class="detail-value">${finding.match}</div>
                    </div>
                    ${molEffectHtml}
                    ${originHtml}
                    ${reviewHtml}
                </div>
                
                <div class="disease-info">
                    <div class="disease-name">Associated Disease/Condition:</div>
                    <div class="disease-description">${finding.disease}</div>
                </div>
                
                <div>
                    <a href="${snpediaLink}" target="_blank" class="external-link">View on SNPedia →</a>${clinvarLinkHtml}
                </div>
            </div>
        `;
    }
    
    downloadReport() {
        const reportContent = document.getElementById('reportContent').innerHTML;
        const fullHtml = `<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>DNA Report - ${new Date().toISOString().split('T')[0]}</title>
    <link rel="stylesheet" href="styles.css">
</head>
<body>
    <div class="app-container">
        <div class="disclaimer">
            <h3>⚠️ Important Medical Disclaimer</h3>
            <p><strong>This report is for educational and informational purposes only.</strong></p>
            <p>This is NOT medical advice. Please consult with a qualified healthcare provider.</p>
        </div>
        ${reportContent}
    </div>
</body>
</html>`;
        
        const blob = new Blob([fullHtml], { type: 'text/html' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = `dna-report-${new Date().toISOString().split('T')[0]}.html`;
        a.click();
        URL.revokeObjectURL(url);
    }
    
    updateProgress(percent, message) {
        const progressFill = document.getElementById('progressFill');
        const progressText = document.getElementById('progressText');
        
        progressFill.style.width = percent + '%';
        progressText.textContent = message;
    }
    
    showError(message) {
        const processingSection = document.getElementById('processingSection');
        processingSection.innerHTML = `
            <h2>❌ Error</h2>
            <div class="error-message">
                <p>${message}</p>
                <button onclick="location.reload()" class="retry-button">Try Again</button>
            </div>
        `;
    }
}

// Initialize app
const app = new DNAReportApp();
