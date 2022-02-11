# MPRA-design-tool
An R script for designing probes for massively parallel reporter assays, investigating the influence of SNP genotype on enhancer/promoter activity. 
The scripts accepts lists of genetic variants, designs probes ( Allele A / Allele B / Fwd / Rev) and outputs a fasta file of sequences that can be sent directly for oligo synthesis. 

Massively parallel reporter assays (MPRA) simultaneously assess the enhancer / promoter activity of hundreds/thousands of DNA sequences. MPRA can be used to assess the influence small genetic polymorphisms on enhancer activity of all genetic variants linked to trait, for example from a GWAS study. This script accepts a list of SNP/Indels and designs probes to be used in an MRPA experiment such as (https://www.nature.com/articles/s41596-020-0333-5). Users can set their own probes length, adapter sequence and other features such as removing/including indels. 

A feature of the script is that it will refine user provided risk loci by reference to linkage scores from LDlink (https://ldlink.nci.nih.gov). Users must register here for an API token before use, see instructions for further information. Users provide a summary SNP file containing variants tagging regions to be investigated. This is used to extract variants within a user defined genomic space. Within each region loci are recursively split based on an optional R2 threshold into independent loci. Independent loci are then filtered using a P values threshold based on the lead variant in that loci. 
The script is designed to work with the output from META (https://mathgen.stats.ox.ac.uk/genetics_software/meta/meta.html) a frequently used tool for GWAS meta-analysis. 
