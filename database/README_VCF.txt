Last updated January 30, 2024

README file for: 

ftp.ncbi.nih.gov/pub/clinvar/vcf_GRCh37/

ftp.ncbi.nih.gov/pub/clinvar/vcf_GRCh38/

This document contains the following sections:

1. GENERAL INFO
2. SCOPE
3. FILE FORMAT
4. ARCHIVE
5. VCF FILES FROM OTHER NCBI DATABASES

We welcome your feedback. Please email us at clinvar@ncbi.nlm.nih.gov.


==================================================================================
1. GENERAL INFO
==================================================================================

These directories contain VCF files for the ClinVar dataset.
* Each row in the file represents a single variant that has been reported to ClinVar.
* Files are provided for the GRCh37 (hg19) and GRCh38 (hg38) human genome assemblies. 
* The 4.1 version of the VCF specification is used (http://www.internationalgenome.org/wiki/Analysis/variant-call-format).
* Files are generated weekly, in the weekly sub-directory.
* Only the monthly release (first Thursday of the month, in each VCF directory) is archived as part of ClinVar's regular release.

As an example, VCF files for the GRCh37 assembly are:
* clinvar_vcf.GRCh37.vcf.gz
* clinvar_vcf.GRCh37_papu.vcf.gz - variants with precise locations in the PAR region on the Y chromosome and on alternate loci, patches, 
  or unplaced/unlocalized contigs on NT_ and NW_ accessions (papu: patch, alternate, PAR, unplaced). This file pairs with clinvar_vcf.GRCh37.vcf.gz.

Ambiguous bases: ClinVar accepts all IUPAC ambiguity codes for nucleotides. However, the 4.1 version of the VCF specification only allows 
ambiguity code N. Thus, ClinVar XML files retain the actual ambiguous bases, but all ambiguous values are converted to N in the VCF files.


==================================================================================
2. SCOPE
==================================================================================

The ClinVar VCF files include variants with the following characteristics:

* simple allele (not a haplotype or genotype)
* precise endpoints
* mapped to the GRCh37 and/or GRCh38 genome assembly
* <10 kb in length

All other variants are not in scope, including cytogenetic variants, copy number variants with inner and/or outer start and stop coordinates, 
and variants >10 kb. 

Inconsistencies between the VCF and VCV XML files: ClinVar's VCF files are generated directly from the ClinVar VCV XML (ClinVarVariationRelease) file, 
which is a comprehensive extraction of the ClinVar database. Because some variants are not in scope for the ClinVar VCF file, it is expected that some 
variants are in the VCV XML file but not in the VCF file.

Included variants: Classifications in ClinVar may be made for a single variant or a set of variants, such as a haplotype. Variants classified only as 
part of a set of variants (i.e. no direct classification for the variant itself) are considered "included" variants. The VCF files include both variants 
with a direct classification and included variants. Included variants do not have an associated disease (CLNDN, CLNDISDB) or a classification (CLNSIG). 
Instead, there are three INFO tags specific to included variants - CLNDNINCL, CLNDISDBINCL, and CLNSIGINCL (see below).


==================================================================================
3. FILE FORMAT
==================================================================================

The eight fixed fields of the ClinVar VCF files include:

CHROM		chromosome number, from either the GRCh37 or GRCh38 human genome reference assembly
POS		position
ID		the ClinVar Variation ID (https://www.ncbi.nlm.nih.gov/clinvar/docs/identifiers/#variation)
REF		reference base(s)
ALT		alternate base(s)
QUAL		quality. This field is always '.' (no value)
FILTER		filter status. This field is always '.' (no value)
INFO		additional information. INFO tags include:

AF_ESP		allele frequencies from GO-ESP
AF_EXAC		allele frequencies from ExAC 
AF_TGP		allele frequencies from TGP (1000 Genomes Project)
ALLELEID	the ClinVar Allele ID (https://www.ncbi.nlm.nih.gov/clinvar/docs/identifiers/#allele)
CLNDN		ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB
CLNDNINCL	for included variants only. ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDBINCL
CLNDISDB	tag-value pairs of disease database name and identifier, e.g. OMIM:NNNNNN
CLNDISDBINCL	for included variants only. Tag-value pairs of disease database name and identifier, e.g. OMIM:NNNNNN
CLNHGVS		top-level (primary assembly, alt, or patch) HGVS expression
CLNREVSTAT	ClinVar's review status of the germline classfication for the Variation ID
CLNSIG		aggregate germline classification for this single variant; multiple values are separated by a vertical bar
CLNSIGCONF	conflicting germline classification for this single variant; multiple values are separated by a vertical bar
CLNSIGINCL	germline classification for a haplotype or genotype that includes this variant. Reported as 
		pairs of Variation ID:classification; multiple values are separated by a vertical bar
CLNVC		variant type
CLNVCSO		Sequence Ontology ID for the variant type (www.sequenceontology.org)
CLNVI		identifiers for the variant from other databases, e.g. OMIM Allelic variant IDs
DBVARID		nsv accessions from dbVar for the variant
GENEINFO	gene(s) for the variant reported as gene symbol:NCBI GeneID. The gene symbol and ID are delimited by a colon 
 		and each pair is delimited by a vertical bar.
MC		comma separated list of molecular consequence in the form of Sequence Ontology ID|molecular_consequence
ONCDN		ClinVar's preferred disease name for the concept specified by disease identifiers in ONCDISDB
ONCDNINCL	For included variant: ClinVar's preferred disease name for the concept specified by disease identifiers in ONCDISDBINCL
ONCDISDB	Tag-value pairs of disease database name and identifier submitted for oncogenicity classifications, e.g. MedGen:NNNNNN
ONCDISDBINCL	For included variant: Tag-value pairs of disease database name and identifier for oncogenicity classifications, e.g. OMIM:NNNNNN
ONC		Aggregate oncogenicity classification for this single variant; multiple values are separated by a vertical bar
ONCINCL		Oncogenicity classification for a haplotype or genotype that includes this variant. Reported as pairs of VariationID:classification; 
		multiple values are separated by a vertical bar
ONCREVSTAT	ClinVar review status of oncogenicity classification for the Variation ID
ONCCONF		Conflicting oncogenicity classifications for this single variant; multiple values are separated by a vertical bar
ORIGIN		 allele origin reported to ClinVar. One or more of the following values: 0 - unknown; 1 - germline; 2 - somatic; 4 - inherited;
 		 8 - paternal; 16 - maternal; 32 - de-novo; 64 - biparental; 128 - uniparental; 256 - not-tested; 512 - tested-inconclusive;
		 1073741824 - other
RS		 dbSNP ID (i.e. rs number) from dbSNP build 155 
SCIDN		ClinVar's preferred disease name for the concept specified by disease identifiers in SCIDISDB
SCIDNINCL	For included variant: ClinVar's preferred disease name for the concept specified by disease identifiers in SCIDISDBINCL
SCIDISDB	Tag-value pairs of disease database name and identifier submitted for somatic clinial impact classifications, e.g. MedGen:NNNNNN
SCIDISDBINCL	For included variant: Tag-value pairs of disease database name and identifier for somatic clinical impact classifications, e.g. OMIM:NNNNNN
SCIREVSTAT	ClinVar review status of somatic clinical impact for the Variation ID
SCI		Aggregate somatic clinical impact for this single variant; multiple values are separated by a vertical bar
SCIINCL		Somatic clinical impact classification for a haplotype or genotype that includes this variant. Reported as pairs of VariationID:classification; 
		multiple values are separated by a vertical bar


==================================================================================
4. ARCHIVE
==================================================================================

VCF files are archived in the following directories, with subdirectories for each year:

ftp.ncbi.nih.gov/pub/clinvar/vcf_GRCh37/archive_2.0/
ftp.ncbi.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/

VCF files in the old format for the Sept 2017 release and all prior releases are archived also:

ftp.ncbi.nih.gov/pub/clinvar/vcf_GRCh37/archive_1.0/
ftp.ncbi.nih.gov/pub/clinvar/vcf_GRCh38/archive_1.0/


==================================================================================
5. VCF FILES FROM OTHER NCBI DATABASES
==================================================================================

Didn't find what you were looking for? Other NCBI products generate VCF files independently of ClinVar, including:

dbSNP (ftp.ncbi.nih.gov/snp/latest_release/), which includes all small variants regardless of clinical relevance
dbVAR (www.ncbi.nlm.nih.gov/dbvar/content/ftp_manifest/), which includes all large variants such as copy number variants and variants >10 kb in length.