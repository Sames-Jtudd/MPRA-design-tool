# MPRA design tool
An R script for designing probes for massively parallel reporter assay (MPRA). MPRA investigates the influence of sequence variation on enhancer/promoter activity, simultaneously assessing hundreds or thousands of DNA sequences. MPRA can be used to assess the influence small genetic polymorphisms (SNPs / indels) on enhancer activity for genetic variants linked to trait, for example from a GWAS study. The script designs probes to be used in an MRPA experiment such as (https://www.nature.com/articles/s41596-020-0333-5). 
The user submits a list of genetic variants which are used for probe design, creating REF, ALT alleles and Fwd, Rev orientations, outputing a text file of sequences that can be sent directly for oligo synthesis. 
Users can set their own probe length, adapter sequence and other features such as removing/including indels. 
A feature of the script is that it will refine user provided risk loci by reference to linkage from LDlink (https://ldlink.nci.nih.gov). Users must register here for an API token before use, see instructions for further information. Users provide a summary SNP file containing variants tagging regions to be investigated. This is used to extract variants within a user defined genomic space. Within each region loci are recursively split based on an optional R2 threshold into independent loci. Independent loci are then filtered using a P value threshold based on the lead variant in that loci. 
The script is designed to work with the output from META (https://mathgen.stats.ox.ac.uk/genetics_software/meta/meta.html) a frequently used tool for GWAS meta-analysis. 

Functions performed by the script ;

Part 1 : extract variants within a defined genetic interval around sentinel/lead SNPs from which  probe sequences will be created. Sentinel/lead SNPs are defined by **--summary_SNP_file** and full SNP list by **--all_SNP_file**.

Look up SNPs against a reference database file to ensure mappings and alleles are correct.

Determine reference and alternate alleles.

Optionally, remove SNPs that are known sequencing artifacts.

Part 2 : split SNPs into linkage based loci by extracting R2 information from NIH LDlink https://ldlink.nci.nih.gov/?tab=apiaccess.
Apply a P-value filter based on the lead SNP in each loci.
Optionally add proxies that are not in the original SNP list (which may have been filtered eg due to info score).

Part 3 : Extract the sequences for the REF and ALT alleles.
Add adapters and then optionally alter or filter sequences that have homopolymers.
Generate control sequences.
Perform some final checks (remap sequences to verify consistency). 

Part 4 : Optionally generate figures of the SNP filtering

# Instructions

**Before running** 
The script will connect to the NIH LDlink website to extract linkage information. You will need to go to this link (https://ldlink.nci.nih.gov/?tab=apiaccess) and register for an API token and parse this, using **-T, --API_token**.

**Quick start** 
When running the script ensure you have at least 40Gb of RAM as the script reads both genome and reference SNP files into memory. 
To execute, in its simplest form
> Rscript /path_to_script/design_library_v0.22.R \ \
-i / path_to/summary_SNPs.txt \ \
-a /path_to/all_SNPs.txt \ \
-d /path_to/SNPdb.vcf.gz \ \
-g /path_to/reference_genome.fa.gz \ \
-o /path_to/outdir \ \
-T {API_token}

**Libraries**

The script has the following dependencies, they should install automatically on first run, if they don’t please install manually.
From CRAN - stringr, optparse, LDlinkR, data.table, dplyr, remotes, BiocManager, taRifx, RColorBrewer.
From BiocManager – rtracklayer, VariantAnnotation, Biostrings.
From Github – sarlacc, SparseSummarizedExperiment.

**Run environment information**

R version 4.0.3 (2020-10-10)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: CentOS Linux 8 (Core)

Matrix products: default
BLAS/LAPACK: /opt/software/applications/anaconda/3/envs/r-essentials4.0/lib/libopenblasp-r0.3.12.so

attached base packages:
parallel, stats4,stats, graphics, grDevices, utils, datasets,methods, base

other attached packages:
SparseSummarizedExperiment_0.0.0.9021, sarlacc_1.0.0, VariantAnnotation_1.36.0, Rsamtools_2.6.0, Biostrings_2.58.0, XVector_0.30.0, SummarizedExperiment_1.20.0           Biobase_2.50.0, MatrixGenerics_1.2.1, matrixStats_0.58.0, rtracklayer_1.50.0                    GenomicRanges_1.42.0, GenomeInfoDb_1.26.7,  IRanges_2.24.1, S4Vectors_0.28.1                      BiocGenerics_0.36.1, RColorBrewer_1.1-2, taRifx_1.0.6.2, BiocManager_1.30.15                   remotes_2.3.0, dplyr_1.0.5, data.table_1.14.0, LDlinkR_1.1.2, optparse_1.7.1, stringr_1.4.0

loaded via a namespace (and not attached):
httr_1.4.2, bit64_4.0.5, assertthat_0.2.1, askpass_1.1, BiocFileCache_1.14.0, latticeExtra_0.6-29, blob_1.2.1, BSgenome_1.58.0, GenomeInfoDbData_1.2.4, progress_1.2.2, pillar_1.6.0, RSQLite_2.2.3, lattice_0.20-41, glue_1.4.2, Matrix_1.3-2, plyr_1.8.6, XML_3.99-0.6, pkgconfig_2.0.3, ShortRead_1.48.0, biomaRt_2.46.3, zlibbioc_1.36.0, purrr_0.3.4, jpeg_0.1-8.1, getopt_1.20.3, BiocParallel_1.24.1, tibble_3.1.1, openssl_1.4.3, generics_0.1.0, ellipsis_0.3.1, cachem_1.0.4, GenomicFeatures_1.42.3, magrittr_2.0.1, crayon_1.4.1, memoise_2.0.0, fansi_0.4.2, hwriter_1.3.2, xml2_1.3.2, tools_4.0.3, prettyunits_1.1.1, hms_0.5.3, lifecycle_1.0.0, DelayedArray_0.16.3, AnnotationDbi_1.52.0, compiler_4.0.3, rlang_0.4.10, grid_4.0.3, RCurl_1.98-1.3, rstudioapi_0.13, rappdirs_0.3.3, bitops_1.0-7, curl_4.3, DBI_1.1.0, reshape2_1.4.4, R6_2.5.0, GenomicAlignments_1.26.0 fastmap_1.0.1, bit_4.0.4, utf8_1.1.4, stringi_1.5.3, Rcpp_1.0.7, png_0.1-7, vctrs_0.3.6, dbplyr_2.0.0, tidyselect_1.1.0

**INPUT FILES** 

For both the summary_SNP_file and main_SNP_file;
chromosome notation eg (chr1 or 1) should be consistent in both files.
P_value notation should be exponential format eg 2.3755e-64

*Required inputs*
1)	**-i, --summary_SNP_file** : "/path_to/summary_SNPs.txt"
A headered text file, containing a list of lead SNPs used to define loci for probe design. Tab delimited (can be gzipped). 

Header must contain the columns “CHR”, “POS” and “RSID”, can contain others.

eg
CHR	POS	RSID	P \
1	22503282	rs2807367	3.27E-08


2)	**-a, --all_SNP_file** : "/path_to/meta_results_header.txt"
Headered text file containing the full list of genetic variants. May be the output from META. Tab delimited (can be gzipped).

Header must contain the columns “chr”, “pos”, “rsid”, "P_value",”allele_A”,” allele_B”;

optional headers;

“RAF” (relative allele freq) if present variants can be filtered to remove those with a RAF below a value specified by **--min_RAF**.
“P_heterogeneity” if present and **--filter_phet**=T variants with significant heterogeneity in a region where the lead SNP doesn’t exhibit heterogeneity are filtered. (Intended for studies from multiple populations) 
Can contain others columns.

3)	**-d, --snp_db_file** : ”/path_to /00-common_all.vcf.gz" 
Reference data base SNP file for the species and genomic build used for generating SNP lists/GWAS.
Must be vcf file format.
Human dbSNP reference file (Hg37) can be download from https://ftp.ncbi.nih.gov/snp/organisms/ 

4)	**-g, --genome_file** : ”/path_to /hg19.fa.gz"
Fasta file of the relevant species and genome build, must be the same as the SNPdb file.
Must be fasta file format.


*Optional inputs*

5)	**-n, --negative_control_file** : "/path_to /H3K27me3.broadPeak.gz" 
A bed file use to design control probes. Enriched regions are defined by **--ChIP_enrichment**. Regions with an enrichment score below this value are discarded. Retained regions are overlapped with variants in **--snp_db_file** and randomly selected for probe design. 
Must be headerless bed file with the columns arranged as per;
"**chr**","**start**","**end**","rank","score","strand","**enrichment**","p-val","q-val"
columns in bold are critical other can contain dummy/sham data. However the order of the columns and total number must be preserved.  
Eg https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/broadPeak/


6)	**-u, --unfiltered_file_path** : "/path_to/unfiltered_rsid_{CHROM}.txt.gz"
Directory containing the per chr output from meta, without any P value filtering. If **--add_proxies** is true, these files are used to prevent importing variants which were not in **--all_SNP_file** due to P-value filtering but are linked to risk variants. Chromosome names are inferred from the SNPs in the **--summary_SNP_file**. These are then substituted for “CHROM” in the path to specify the relevant file path. Eg for variants on chromosome 1 the script will attempt to import unfiltered_rsid_1.txt.gz. The “CHROM” in the path must be present to allow substitution with chr names. The remainder of the path must match the relevant file name. File format same as **--all_SNP_file**.

7)	**-b, --black_list_snp_file** : "/path_to/b151_rs_without_GRCh38_mapping.bcp"
Text file.
Must be headerless file with the columns arranged as per;
"chr",”pos”,"rsid”

Other columns can after these can exist. Any consistent separator can be used. 
Some SNPs in build hg37 are artefacts and have been removed in hg38 dbSNP builds they are listed here   # https://gatk.broadinstitute.org/hc/en-us/articles/360035890951-Human-genome-reference-builds-GRCh38-or-hg38-b37-hg19

**OPTIONS**
1)	**-o, --outdir** : ”/path_to/output_dir”
Path to output dir where output files are saved.
2)	**-s, --suffix** : string.
Suffix added to output file names.
3)	**-T, --API_token** : string
Token to access LDlink database see 'https://ldlink.nci.nih.gov/?tab=apiaccess' to obtain.
4)	**-F, --Forward_adapter** : character string,  default="AGGACCGGATCAACT".
Sequence added to the 5' end of each oligo.
5)	**-R, --Reverse_adapter** : character string, default="CATTGCGTGAACCGA". 
Sequence added to the 3' end of each oligo.
6)	**-w,--loci_window** : positive integer, default=250000. 
Variants less than this number of base pairs from a lead SNP will be used for designing probes.
7)	**-z, --insert_size** :  positive integer, default=199 eg 99bp/SNP/99bp.
Desired length of probes/cis-regulatory sequences containing a SNP; if indels are retained probe length will differ from selected value.
8)	**-c, --n_controls** : positive integer, default=100. 
The number of control sequences to design. One control = 4 oligos; Fwd/Rev : Allele1/Allele 2.
9)	**-e, --ChIP_enrichment** : positive integer, default=4.
If providing a control region file only design control probes/CRSs in regions with an enrichment > than this value.
10)	**--prop_control** : numeric, default=0.5 , range 0-1.
If providing a control region file what proportion of controls sequences should be taken from this. If less than 1 remainder are randomly generated.
11)	**-m, --minor_allele_frequency** : numeric, default=0.01, range 0-1.
Variants with a MAF/RAF lower than this value are removed.
12)	**-P, --P_value_cuttoff** : numeric, default=5e-8, format must be scientific notation Xe-y.
Hard P-value cut-off. Variants with a P-value above this value are removed.
13)	**-f, --P_cuttoff_factor** : numeric, default=0.7, range 0-1, higher the number the more variants filtered.
Used to generate locus specific P_value thresholds. Calculated as -log(lead_Pval,10)*P_cuttoff_factor. Variants within a locus with a P value above this are removed. lead_Pval = P value of the lead SNP within a defined loci.
14)	**-S, --SNVs_only** : logical, default=TRUE, values TRUE|FALSE.
If true retain only single nucleotide variants.
15)	**-L, --max_indel_length** : positive integer, default=10. 
Indels with an alt allele length > this value are removed; irrelevent if -S is TRUE.
16)	**-U,--keep_unmapped_variants** : logical, default=TRUE, values TRUE|FALSE. 
If TRUE user provided variants not present in reference SNP file (snp_db_file) are retained.
17)	**--update_alleles** : logical, default=FALSE, values TRUE|FALSE.
If TRUE user provided variants not matching the reference SNP file (snp_db_file) alleles are reverted to the ref db alleles.
18)	**--output_figures** : logical, default=TRUE, values TRUE|FALSE.
Generate figures showing filtered and retained variants.
19)	**--split_by_r2** : logical, default=TRUE, values TRUE|FALSE.
If TRUE variants in a region surrounding a sentinel SNP will be split into separate loci based on linkage. The SNP with the lowest p-value in a new loci will be the new lead and used for p-value filtering.
20)	**-r, --r2_cuttoff** : numeric, default=0.3, range 0-1. 
When splitting loci variants with an r2, to the lead SNP, below this value are split into separate loci.
21)	**-M, --mutate_homopolyers** : logical, default=TRUE, values TRUE|FALSE.
Probes with homopolymers (specified by --max_hp_len) will be mutated to remove homopolymers.
22)	**--filter_homopolymers** : , logical, default=TRUE, values TRUE|FALSE.
Probes with homopolymers (specified by --max_hp_len) will be filtered. If --mutate_homopolyers is TRUE successfully mutated homopolymers are retained.
23)	**-l, --min_dist_mut2snp** : positive integer, default=5.
Bases within a homoploymer < this number of base pairs from the reference or alternative allele are not mutated.
24)	**-p, --populations** : character, default="EUR", values 'EUR|SAS|EAS|AMR|AFR'.
Populations used for r2 extraction. Multiple populations separate by a comma.
25)	**-x, --add_proxies** : logical, default=FALSE, values TRUE|FALSE. 
Import variants linked to a lead SNP but absent from user provided variants. If TRUE necessary to specify --unfiltered_file_path.
26)	**--proxy_maf** : numeric, default=0.01, value 0-1  
If --add_proxies is TRUE, proxies with a MAF less than this value are filtered.
27)	**--proxy_r2** : numeric, default=0.7, range 0-1 
Filter proxies with an r2, to loci lead SNP, less than this value. 
28)	**--filter_phet** : logical, default=FALSE, values TRUE|FALSE. 
Where lead SNPs has phet > value specified in --P_het_cuttoff option (ie not heterogenous) filter variants within that loci with phet < P_het_cuttoff. If FALSE no filtering is performed.
29)	**--P_het_cuttoff** : numeric, default=0.01, value 0-1.
If --filter_phet = TRUE. This value defines the cut off. 
30)	**-C, --cycle_limit** : Positive integer, default=5.
Risk loci are split into loci containing only linked variants, this value sets the upper limit on the number of cycles used to split loci. 
31)	**--strict** : logical, default=FALSE, values TRUE|FALSE.
If TRUE script will exit where it encounters issues that while may indicate some errors, but do not prevent execution.
 
**OUTPUT FILES**

**MPRA_library_sequences.txt**  - contains the final probe set names and sequences that can be submitted directly to an oligo synthesis company for manufacture.

**MPRA_summary.txt** - details the number of oligos, SNPs, and loci in the final probe set.
MPRA_filtered_SNPs.txt - contains all the filtered SNPs and the filter which was applied to them.

**MPRA_summary_filtered.txt** – contains summary information about the numbers of SNPs filtered due to various constraints.

**MPRA_summary_stats.txt** – describes each loci after filtering, provides information on the number of SNPs pre and post filtering and the P-value threshold with was applied to each. 

**MPRA_library_all_variants.txt** – contains additional information on all test SNPs in the final probe set. Including, lead SNP, original probe sequence and loci P-value cut-off.

**MPRA_library_included_variants_and_cotrols.txt** – contains additional information on all, test and control SNPs, in the final probe set. Including, lead SNP, original probe sequence and loci P-value cut-off.

**MPRA_prefilter_sequences.txt** – contains all variants extracted for probe design before filtering.
 in all output files unless otherwise specified RSID names are those updated based on the supplied SNP reference file.   
