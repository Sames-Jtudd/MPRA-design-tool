#!/usr/bin/env Rscript
##########################################################
# design MPRA oligos 
# the scripts has three parts 
# part 1 - this will look up all the snps against a ref db file 
# 			this will ensure all the mappings and alleles are correct 
# 			it will also remove snps that are known sequencing artifacts
# part 2 - this will extract r2 inforation for all the snps and define loci based on linkage
# 			it will then apply a p-value filter based on the lead snp in that loci and an optional factor
#			finally it will add proxies that are not in the oringal snp list back in 
#			as these may have been filtered for info score or other reasons
# part 3 - extract the sequences for the ref and alt alleles ; 
#			add the adapters and then alter sequences that have homopolmers
#			generate control alleles
#			perform some final checks 
# part 4 - generate some figures of the snp filtering
#
##########################################################


#set up the libraries
sink("/dev/null", type = c("output","message"))

required_libs_CRAN=c("stringr","optparse")
for (lib in required_libs_CRAN) {
	if ( ! require(lib,character.only=T)) {
		install.packages(lib) 
	} else {
		 suppressPackageStartupMessages(library(lib,character.only=T))
		 #suppressMessages(library(lib,character.only=T),"all")
	}
}
sink(NULL)
sink(stdout(), type = c("message"))


############# parse options 

option_list = list(
	make_option(c("-i", "--summary_SNP_file"), type="character", default=NULL, 
		help="path to file containing lead variants, required input [default= %default]"),
	make_option(c("-a", "--all_SNP_file"), type="character", default=NULL, 
		help="path to file containing all variants, required input [default= %default]"),
	make_option(c("-d", "--snp_db_file"), type="character", default=NULL, 
		help="path to SNP db VCF file, required input [default= %default]"),
	make_option(c("-g", "--genome_file"), type="character", default=NULL, 
		help="path to genome fasta file, required input [default= %default]"),
	make_option(c("-n", "--negative_control_file"), type="character", default=NULL, 
		help="path to bed file contaning regions for designing control regions, optional input [default= %default]"),
	make_option(c("-u", "--unfiltered_file_path"), type="character", default=NULL, 
		help="path to files contaning unfilteted variants, see intructions for more information relating to format, to add proxies this must be provided, optional input [default= %default]"),
	make_option(c("-b", "--black_list_snp_file"), type="character", default=NULL, 
		help="path to file containing known techincal artefacts or other SNPs to remove, if provided user provided variants are filtered against this list, optional input [default= %default]"),

	make_option(c("-o", "--outdir"), type="character", default=".", 
		help="directory to save output [default= %default]") ,
	make_option(c("-s", "--suffix"), type="character", default="", 
		help="suffix added to ouput file names [default= %default]", metavar="character"),
	make_option(c("-T", "--API_token"), type="character", default=NULL, 
		help="token to access LDlink database see 'https://ldlink.nci.nih.gov/?tab=apiaccess' to obtain", metavar="character"),
	make_option(c("-F", "--Forward_adapter"), type="character", default="AGGACCGGATCAACT", 
		help="sequence added to the 5\\' end of each oligo [default= %default]", metavar="character") ,
	make_option(c("-R", "--Reverse_adapter"), type="character", default="CATTGCGTGAACCGA", 
		help="sequence added to the 3\\' end of each oligo [default= %default]", metavar="character") ,
	make_option(c("-w", "--loci_window"), type="integer", default=250000, 
		help="all variants less than this number of basepairs from a lead SNP will be considered for designing probes. Value - postive integer [default= %default]", metavar="integer"),
	make_option(c("-z", "--insert_size"), type="integer", default=199, 
		help="length of a probe/ cis-regulatory sequences that contains a SNP ; indels if retained will produces sequences of different length. Value - postive integer [default= %default]", metavar="integer"),
	make_option(c("-c", "--n_controls"), type="integer", default=100, 
		help="The number of control sequences to design. One control = 4 oligos; Fwd/Rev : Allele1/Allele 2. Value postive integer [default= %default]", metavar="integer") ,
	make_option(c("-e", "--ChIP_enrichment"), type="numeric", default=4, 
 		help="If providing a control region file only design control probes/CRSs in regions with an enrichemnt > than this value. [default= %default]", metavar="numeric"),	
	make_option(c("--prop_control"), type="numeric", default=0.5, 
 		help="If providing a control region file what proportion of controls sequences should be from this. Value - 0-1 [default= %default]", metavar="numeric"),		
	make_option(c("-m", "--minor_allele_frequency"), type="numeric", default=0.01, 
		help="variants with a MAF/RAF lower than this value are removed. Value 0-1 [default= %default]", metavar="numeric"),
	make_option(c("-P", "--P_value_cuttoff"), type="numeric", default=5e-8, 
		help="any varient with a P-value above this vaule is removed [default= %default]", metavar="numeric"),
	make_option(c("-f", "--P_cuttoff_factor"), type="numeric", default=0.7, 
		help="P_value thresholds for each risk loci are caluated as -log(lead_Pval,10)*P_cuttoff_factor. Value 0-1  [default= %default]", metavar="numeric") ,
	make_option(c("-S", "--SNVs_only"), type="logical", default=T, 
		help="retain only single nucleotide varients [default= %default]", metavar="logical"),
	make_option(c("-L", "--max_indel_length"), type="integer", default=10, 
		help="indels with an alt allele > this value are removed; irrelevent if -S is T. Value - postive integer [default= %default]", metavar="integer"),
	make_option(c("-U", "--keep_unmapped_variants"), type="logical", default=T, 
		help="retain user provided variants not present in reference SNP file [default= %default]", metavar="logical"),
	make_option(c("--output_figures"), type="logical", default=T, 
		help="generate figures showing filtered and retained variants [default= %default]", metavar="logical"),
	make_option(c("--split_by_r2"), type="logical", default=T, 
		help="if T varinats in a region surrounding a sentinal SNP will be split into separate loci based on linkage. The SNP with the lowest p-value in a new loci will be the new lead and used for p-value filtering. [default= %default]", metavar="logical"),
	make_option(c("-r", "--r2_cuttoff"), type="numeric", default=0.3, 
		help="when splittting loci variants with an r2 , to the lead SNP, below this value are split into separate loci. Value - 0-1 [default= %default]", metavar="numeric"),
	make_option(c("--max_hp_len"), type="integer", default=8, 
		help="mono nucleotide runs > than is value are considered as homopolymers. Value - positive integer [default= %default]", metavar="integer"),
	make_option(c("-M", "--mutate_homopolyers"), type="logical", default=T, 
		help="probes with homopolymers (specified by --max_hp_len) will be mutated to remove homopolymers. If F any homopolymers will be filtered [default= %default]", metavar="logical") ,
	make_option(c("-l", "--min_dist_mut2snp"), type="integer", default=5, 
		help="bases within a homoploymer < this values are not mutated. Instead the probe/CRS is removed. Value - positive integer [default= %default]", metavar="integer"),
	make_option(c("-p", "--populations"), type="character", default="EUR", 
		help="populations used for r2 extraction one or more of 'EUR,SAS,EAS,AMR,AFR' If multiple populations separate by a comma [default= %default]", metavar="character"),
	make_option(c("-x", "--add_proxies"), type="logical", default=F, 
		help="import variants linked to a lead SNP but absent from user provided variants [default= %default]", metavar="logical"),
	make_option(c("--proxy_maf"), type="numeric", default=0.01, 
		help="filter proxies with a maf less than this value. Value 0-1  [default= %default]", metavar="numeric"),
	make_option(c("--proxy_r2"), type="numeric", default=0.7, 
		help="filter proxies with a r2, to loci lead, less than this value. Value 0-1 [default= %default]", metavar="numeric"),
	make_option(c("--filter_phet"), type="logical", default=F, 
		help="where lead SNPs has phet > value specified in P_het_cuttoff option (ie not heterogenous) filter variants that have phet < P_het_cuttoff. If F no filtering is performed [default= %default]", metavar="logical"),
	make_option(c("--P_het_cuttoff"), type="numeric", default=0.01, 
		help="if --filter_phet = T. The phet cutoff for filtering set by this value. Value 0-1 [default= %default]", metavar="numeric"),
	make_option(c("-C", "--cycle_limit"), type="integer", default=5, 
		help="risk loci are split into loci containing only linked variants, this value sets the upper limit on the number of cycles used to split loci. Value - positve integer [default= %default]", metavar="integer"),
	make_option(c("--strict"), type="logical", default=F, 
		help="If TRUE script will exit where it encounters issues that while may indicate some errors, but do not prevent execution [default= %default]", metavar="logical")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#input files
#required
summary_SNP_file=opt$summary_SNP_file
main_SNP_file=opt$all_SNP_file
snp_db_file=opt$snp_db_file
genome_file=opt$genome_file

#optional
negative_control_region_file=opt$negative_control_file
unfiltered_file_path=opt$unfiltered_file_path
black_list_snp_file=opt$black_list_snp_file


#set variables
outdir=opt$outdir
suffix=opt$suffix
CRS_adpater_F<-opt$Forward_adapter
CRS_adpater_R<-opt$Reverse_adapter
loci_window<-opt$loci_window *2
insert_size=opt$insert_size -1
min_RAF=opt$minor_allele_frequency
hard_P_val_cutoff=opt$P_value_cuttoff
remove_indels=opt$SNVs_only
max_indel_len=opt$max_indel_length
keep_unmapped_variants=opt$keep_unmapped_variants
API_token=opt$API_token

if ( is.null(API_token) == T) stop("please supply an API token; to obtain a token please visit 'https://ldlink.nci.nih.gov/?tab=apiaccess'")  

cycle_limit=opt$cycle_limit
lower_r2_cuttoff=opt$r2_cuttoff
max_hp_len=opt$max_hp_len
mutate_homopolyers= opt$mutate_homopolyers
min_dist_mut2snp=opt$min_dist_mut2snp
split_by_r2=opt$split_by_r2
p_val_cuttoff_factor<-opt$P_cuttoff_factor
populations<-unlist(str_split(opt$populations,","))

ld_pops=c("EUR","SAS","EAS","AMR","AFR")
if ( ! all(populations %in% ld_pops)) stop ( "one or more populations do not appear in ldlink populations")

add_proxies<-opt$add_proxies
if ( add_proxies == T & is.null(unfiltered_file_path)) stop ( "you have choosen to add_proxies but have not provided an unfiltered_file_path")

MAF_threshold=opt$proxy_maf
r2_cuttoff=opt$proxy_r2
n_controls=opt$n_controls
ChIP_enrichment=opt$ChIP_enrichment
prop_control=opt$prop_control
filter_phet=opt$filter_phet
P_het_cuttoff=opt$P_het_cuttoff
output_figures=opt$output_figures
strict=opt$strict


# check values of variables are sensible
pos_ints<-c(cycle_limit,loci_window,insert_size,max_hp_len,max_indel_len,min_dist_mut2snp,n_controls)
names(pos_ints)<-c("cycle_limit","loci_window","insert_size","max_hp_len","max_indel_len","min_dist_mut2snp","n_controls")

for (i in 1:length(pos_ints)) {
if ( ! ( pos_ints[i] %% 1 == 0) &  pos_ints[i] > 0) stop ( names(pos_ints[i])," must be a postive integer")
}

factors <-c(min_RAF,lower_r2_cuttoff,p_val_cuttoff_factor,MAF_threshold,prop_control)
names(factors)<-c("min_RAF","lower_r2_cuttoff","p_val_cuttoff_factor","MAF_threshold","prop_control")

for (i in 1:length(factors)) {
if ( ! (( factors[i] >= 0) &  factors[i] <= 1)) stop ( names(factors[i])," must be a postive integer")
}

#internal variables
swap_base_list=list("A"="T","T"="A","C"="G","G"="C") # when mutating hps bases will be changed according to this pattern
updating_proxies=F # internal dont change unless clear on logic

#make ouput dir 
if ( ! dir.exists(outdir)) dir.create(outdir)
setwd(outdir)

# check input files
if ( is.null ( "summary_SNP_file")) stop ( "please specify the summary_SNP_file")
if ( is.null ( "main_SNP_file")) stop ( "please specify the main_SNP_file")
if ( is.null ( "snp_db_file")) stop ( "please specify the snp_db_file")
if ( is.null ( "genome_file")) stop ( "please specify the genome_file")

required_input_files=list(summary_SNP_file,main_SNP_file,snp_db_file,genome_file)

for ( file in required_input_files ) {
	if ( ! file.exists(file)) stop ( "cant located file")
}

optional_input_files=c()

if ( ! is.null( "negative_control_region_file")) {
	optional_input_files=c(optional_input_files,negative_control_region_file)
}
if ( ! is.null( "black_list_snp_file")) {
	optional_input_files=c(optional_input_files,black_list_snp_file)
}

for ( file in optional_input_files) {
	if ( ! file.exists(file)) stop ( "cant located file")
} 

sink(NULL,type="message")
sink("/dev/null", type = c("output","message"))
# load the remaining libs
required_libs_CRAN=c("LDlinkR","data.table","dplyr","remotes","BiocManager","taRifx","RColorBrewer")
required_libs_BiocManager=c("rtracklayer","VariantAnnotation","Biostrings")
required_libs_other=c("sarlacc","SparseSummarizedExperiment")

for (lib in required_libs_CRAN) {
	if ( ! require(lib,character.only=T)) {
		install.packages(lib) 
	} else {
		 suppressPackageStartupMessages(library(lib,character.only=T))
	}
}

for (lib in required_libs_BiocManager) {
	if ( ! require(lib,character.only=T)) {
		BiocManager::install(lib)
	} else {
		 suppressPackageStartupMessages(library(lib,character.only=T))
		 #suppressMessages(library(lib,character.only=T),"all")
	}
}

for (lib in required_libs_other) {
	if ( ! require(lib,character.only=T)) {
		remotes::install_github(lib)
	} else {
		 suppressPackageStartupMessages(library(lib,character.only=T))
		 #suppressMessages(library(lib,character.only=T),"all")
	}
}

sink(NULL, type = c("output","message"))
sink(stdout(), type = c("message"))

#calc the insert size
if ( (insert_size)/2 == round((insert_size)/2))
{
	upstream=(insert_size)/2
	downstream <- upstream
} else {
	upstream<- (insert_size/2)-0.5
	downstream <- upstream +1 
}

# read in the input files
summary_SNPs<-fread(summary_SNP_file,header=T)
nec_cols=c("CHR","POS","RSID")
if ( ! all(nec_cols %in% colnames(summary_SNPs))) stop ( "summary snp file doesnt contain the necessary columns")

main_SNPs<-fread(main_SNP_file,header=T)
nec_cols=c("chr", "pos", "rsid", "P_value")
if ( ! all(nec_cols %in% colnames(main_SNPs))) stop ( "main snp file doesnt contain the necessary columns")

if ( ! is.null( black_list_snp_file )) { 
	non_mapping_snps<-fread(black_list_snp_file)
	colnames(non_mapping_snps)[1:3]<-c("chr","pos","rsid") #,"unknown","unknown2")
	rsid<-! is.na(non_mapping_snps$rsid)
	is_rsid<- ! str_sub(non_mapping_snps$rsid[rsid],1,2) == "rs"
	add_rs<-which(rsid)[is_rsid]
	non_mapping_snps$rsid[add_rs]<-paste0("rs",non_mapping_snps$rsid[add_rs])
}


###################
#
# FUNCTIONS
#
###################

read_head <- function(...) read.delim(..., header=T , stringsAsFactors=F,  sep="\t")
write_simp<- function (...) write.table (...,quote=F,row.names=F,sep="\t")

####################################################
# take a nested list of data frames and return a list of dataframe with no nesting

flatten_list=function(list_obj) {	
	y <- list()
	for (i in 1:length(list_obj)) y <- if (class(list_obj[[i]]) == "data.frame") c(y, list_obj[i]) else c(y, list_obj[[i]])
	return(y)
	}
#####################

#this will change a vector without chrs to be all chrs even if mixed 
add_chr<-function(chr_vector) {
	chr_ok<-substr(chr_vector,1,3)
	not_chr<-which(! chr_ok == "chr")
	if(length(not_chr) >0 ) {
		chrom2<-sapply(not_chr, function(row) {
			gsub( "^","chr",chr_vector[row])
		})
		chr_vector[not_chr]<-chrom2
	}
	return(chr_vector)
}

#####################################################################################
# return the indexes for snps where the snp sequence differes out side the snp after changing due to homoployers

ref_alt_diff<-function(indexes) {
	diff_ref_alt<-sapply(as.numeric(indexes), function(ind) {
		ref_len<- nchar(all_variants$REF[ind])	
		alt_len<- nchar(all_variants$ALT[ind])
		ref_no_snp<-paste0(str_sub(all_variants$REF_mod[ind],1,downstream),str_sub(all_variants$REF_mod[ind],downstream+1+ref_len))
		alt_no_snp<-paste0(str_sub(all_variants$ALT_mod[ind],1,downstream),str_sub(all_variants$ALT_mod[ind],downstream+1+alt_len))
		if ( ref_no_snp !=alt_no_snp ) return(ind)
	})

is_null<-sapply(1:length(diff_ref_alt) , function(i) is.null(diff_ref_alt[[i]]))
return(unlist(diff_ref_alt[! is_null]))
}

####################################################################################
# mutate bases with a homopolymer 

mutate_homopolyer<-function(homopolymer_index,has_homo ,allele) {
	if ( ! allele %in% colnames(all_variants)) return ( "allele must be a column in the all_variants table")

	bases_to_change=c()
	mutated_seqs<-lapply(homopolymer_index, function(ind) {
		#message(ind)
		seq_col=paste0(allele,"_seq")
		seq<-all_variants[ind,..seq_col] # $REF_seq
		homos<-has_homo[[ind]]
		homos_revised=c()
		for ( row in 1:length(homos)) {
			polybase=mcols(homos)$base[row] # the base in the homopolymer
			hom_width<-width(homos[row]) # homopolymer width
			n_breaks<-floor(hom_width/(max_hp_len +1))
			spacing<-ceiling(hom_width/(n_breaks+1)) # evenly spaced every n bases
			ref_len=nchar(all_variants[ind,..allele])
			rsID=all_variants[ind,rsid]
			snp_end=downstream+ref_len
			homo_SNP_ol<-end(homos)[row] >= downstream +1-min_dist_mut2snp & start(homos)[row] <= downstream+ ref_len + min_dist_mut2snp
			if (homo_SNP_ol) {
				# senario 1 - must alter bases upsteam of snp - Full overlap and distance from the start homopolymer to the snp is less than the "min_dist_mut2snp" (min distance of an introduced muation from a SNP)
				if (  downstream +1 >= start(homos)[row]  & downstream +1 - start(homos)[row] < min_dist_mut2snp ) { 
				#message(ind)
					first_break=downstream+ref_len + min_dist_mut2snp
					if (first_break-start(homos)[row] <= max_hp_len) { 
						#bases_to_change=start(homos)[row]+max_hp_len
						bases_to_change<-first_break
						hom_width<-end(homos)[row]-bases_to_change
						n_breaks<-floor(hom_width/(max_hp_len +1))
						ad_bases_to_change=c()
						if ( n_breaks > 0) {
							spacing<-ceiling(hom_width/(n_breaks+1)) # evenly spaced every n bases
							ad_bases_to_change<-sapply(1:n_breaks, function(jumps) {
								bases_to_change+(spacing*jumps)
							})
						}
					} else { 
						message("SNP ",rsID, " unable to modify oligos within constratints")
						return()
					}
					bases_to_change<-c(bases_to_change,ad_bases_to_change)
				###########################################################	
				# senario 2 - must alter bases downstream of snp - Full overlap and distance from the SNP to the end of the homopolymer is less than the "min_dist_mut2snp" (min distance of an introduced muation from a SNP)
				} else if ( end(homos)[row] >= snp_end & end(homos)[row] - snp_end < min_dist_mut2snp ) {  
					#message(ind)
					first_break=downstream + 1 - min_dist_mut2snp
					if (end(homos)[row] - first_break <= max_hp_len) { 
						#bases_to_change=end(homos)[row]-max_hp_len
						bases_to_change<-first_break
						hom_width<-bases_to_change-start(homos)[row]
						n_breaks<-floor(hom_width/(max_hp_len +1))
						ad_bases_to_change=c()
						if ( n_breaks > 0) {
							spacing<-ceiling(hom_width/(n_breaks+1)) # evenly spaced every n bases
							ad_bases_to_change<-sapply(1:n_breaks, function(jumps) {
								bases_to_change-(spacing*jumps)
							})
						}
					} else { 
						message( "SNP ",rsID, " unable to modify oligos within constratints")
						return()
					}
					bases_to_change<-c(bases_to_change,ad_bases_to_change)
				##########################################################	
				#} else if ( start(homos)[row] < 101 & end(homos)[row] < 101 & 101 - end(homos)[row] < min_dist_mut2snp) { # senario 4 (like 2 above) - no actual overlap but homoployer is 
				} else if ( start(homos)[row] < downstream + 1 & end(homos)[row] < downstream + 1 & downstream + 1 - end(homos)[row] <= min_dist_mut2snp) { # senario 4 (like 2 above) - no actual overlap but homoployer is 
				#message(ind)
					bases_to_change=min(c(downstream + 1-min_dist_mut2snp,downstream + 1-spacing))
					hom_width<-bases_to_change-start(homos)[row]
					n_breaks<-floor(hom_width/(max_hp_len +1))
					ad_bases_to_change=c()
						if ( n_breaks > 0) {
							spacing<-ceiling(hom_width/(n_breaks+1)) # evenly spaced every n bases
							ad_bases_to_change<-sapply(1:n_breaks, function(jumps) {
								bases_to_change-(spacing*jumps)
							})
						}
					bases_to_change<-c(bases_to_change,ad_bases_to_change)
				##################################################################
				} else if ( start(homos)[row] > snp_end & end(homos)[row] > snp_end & start(homos)[row] - snp_end <= min_dist_mut2snp ) { # senario 3 (like 1 above) - no actual overlap but homoployer is 
					bases_to_change=max(c(downstream+ ref_len+min_dist_mut2snp,downstream + ref_len+spacing))
					hom_width<- end(homos)[row]-bases_to_change
					n_breaks<-floor(hom_width/(max_hp_len +1))
					ad_bases_to_change=c()
						if ( n_breaks > 0) {
							spacing<-ceiling(hom_width/(n_breaks+1)) # evenly spaced every n bases
							ad_bases_to_change<-sapply(1:n_breaks, function(jumps) {
								bases_to_change-(spacing*jumps)
							})
						}
					bases_to_change<-c(bases_to_change,ad_bases_to_change)
				}
				#if ( end(homos)[row] <= 101 ) 
			} else {		
				bases_to_change<-sapply(1:(n_breaks), function(jumps) {
					start(homos)[row]-1+(spacing*jumps)
				}) # starting from the first base this means changing these bases
			}
			if (any( bases_to_change > (downstream + 1 - min_dist_mut2snp) & bases_to_change <  snp_end + min_dist_mut2snp) | length(bases_to_change) == 0) {
				message( "SNP ",rsID, " error modifing oligos within constratints, returning NA")
				seq<-(NA)
			} else { 
				homos_revised=c(homos_revised,row)
				for ( base in bases_to_change) {
					str_sub(seq,base,base)<-swap_base_list[[polybase]]
				}
			} # alter the declared bases
		}
		return(list(seq,homos_revised))
		})
	return(mutated_seqs)
}
###########################################
#concatenate a granges list 

c_granges <- function(grl, 
                      use.names = TRUE, 
                      sep = ".", 
                      save.names = FALSE) {
  grl <- as(grl, "CompressedGRangesList")
  gr <- unlist(grl, use.names = FALSE)
  list.names <- rep(names(grl), elementNROWS(grl)) 
  if (use.names) {
    names(gr) <- paste(list.names, names(gr), sep = sep) 
  }
  if (save.names != FALSE) {
    attr.names <- ifelse(save.names == TRUE, 
                         "list.names", 
                         save.names) 
    mcols(gr)[[attr.names]] <- list.names
  }
  
  return(gr)
}
#########################################

################################# FUNCTION #############################################
########################################################################################
# merge the proxy snps from different populations 

merge_proxies<- function( proxyList, both_missing2) {
#merge_proxies<- function( proxyList, missing_snps, both_missing2) {
	if (any(is.na(names(proxyList)))) {
		message("proxyList must be a named list")
		return()
	}
	populations=names(proxyList)
	proxyList<-lapply(populations, function(pop) {
		rowname_TF<-sapply(proxyList[[pop]], function(slot) {
			if( is.null(slot)) return (NA)
			ncol(slot) == 11
		})
		rowname_slots<-which(rowname_TF)
		for ( ind in rowname_slots ) {
			#print(ind)
			#proxyList[[pop]][[ind]]<-proxyList[[pop]][[ind]][,2:ncol(proxyList[[pop]][[ind]])]
			proxyList[[pop]][[ind]]<-proxyList[[pop]][[ind]]
		}
		return(proxyList[[pop]])
	})
	names(proxyList)<-populations
###########################################################################
#  merge the proxies lists for each population

	proxyList_merge<-lapply(populations, function(pop) {
		orig_lead_snps<-names(proxyList[[pop]])
		proxyList[[pop]]<-lapply(orig_lead_snps, function(snp) { 
			proxyList[[pop]][[snp]]$original_lead_snp<- snp
			return(proxyList[[pop]][[snp]])
		})
		names(proxyList[[pop]])<-orig_lead_snps
		proxies<-proxyList[[pop]][ ! names(proxyList[[pop]]) %in% rownames(both_missing2)[ both_missing2[,pop]]]
		proxies<-rbindlist(proxies)
		proxies$population = pop
		proxies<-filter(proxies, MAF > {{MAF_threshold}})
		return(proxies)
	})
	names(proxyList_merge)<-populations

	all_proxy_snps<-rbindlist(proxyList_merge)$RS_Number
	all_proxy_original_lead_snp<-rbindlist(proxyList_merge)$original_lead_snp
	all_proxy_lead_snp<-rbindlist(proxyList_merge)$lead_snp
	
	merged_proxies<-data.frame(RS_Number=all_proxy_snps,Alleles=NA,Coord=NA,lead_snp=all_proxy_lead_snp,original_lead_snp=all_proxy_original_lead_snp,stringsAsFactors=F)
	merged_proxies<-unique(merged_proxies)
	
	all_proxy_snps<-merged_proxies$RS_Number
	all_proxy_original_lead_snp<-merged_proxies$original_lead_snp
	all_proxy_lead_snp<-merged_proxies$lead_snp
	
	col_names<-colnames(proxyList_merge[[1]])

	match_merged<-paste(all_proxy_snps,all_proxy_original_lead_snp,all_proxy_lead_snp,sep="_")

	match_unmerged<-list()
	for ( pop in populations) {
	match_unmerged[[pop]]<-paste(proxyList_merge[[pop]]$RS_Number,proxyList_merge[[pop]]$original_lead_snp,proxyList_merge[[pop]]$lead_snp,sep="_")
	}
	# this generates a table of T/F showing which populations a snp exists in at a maf > 0.01
	pop_snps<-sapply(populations, function(pop) {
		match_merged %in% match_unmerged[[pop]] 
	})

	both_tf<-sapply(1:nrow(pop_snps), function(i) all(pop_snps[i,]))
	rownames(pop_snps) <- match_merged

	#pop_snps<-t(pop_snps)
	var_col_names<-c("Alleles","MAF","Dprime","R2","Correlated_Alleles")
	var_col_names=var_col_names[var_col_names %in% col_names]

	#merged_proxies<-data.frame(RS_Number=all_proxy_snps,Alleles=NA,Coord=NA,lead_snp=NA,original_lead_snp=NA,stringsAsFactors=F)

	start_cols<-c("RS_Number","Alleles","Coord","lead_snp","original_lead_snp")
	cols_edit=c("RS_Number","MAF","Dprime","R2")

	for ( pop in populations ) {
		add_cols<-paste0(pop,"_",cols_edit)
		merged_proxies[,add_cols] <-NA
	}

	#highly optimsied method for merge the proxies for any number of populations
	for ( pop in  populations) {
		match_ind<-match(match_merged[pop_snps[,pop]],match_unmerged[[pop]])
		pop_prox<-proxyList_merge[[pop]][match_ind]
		#merged_proxies[pop_snps[,pop],start_cols]<-pop_prox[,..start_cols]
		merged_proxies[pop_snps[,pop],start_cols]<-pop_prox[,..start_cols]
		add_cols<-paste0(pop,"_",cols_edit)
		colnames(pop_prox) [colnames(pop_prox) %in% cols_edit] <- add_cols
		merged_proxies[pop_snps[,pop],add_cols]<-pop_prox[,..add_cols]
	}
	####
	rs_cols<-sapply(populations, function(pop) {
		col<-paste0(pop,"_RS_Number")
		merged_proxies[,col]
	})

	#check these have been appropaitely combined
	rs_cols<-as.data.frame(rs_cols, stringsAsFactors=F)
	rs_u_cols<-sapply(1:nrow(rs_cols),function(i) unique(rs_cols[i,]))

	if (! all(match_new_old<-sapply(1:length(all_proxy_snps), function(i) all(all_proxy_snps[i]== unlist(rs_u_cols[i]), na.rm=T)))) {
		stop ( "merging the proxy snps was not successful")
	}

	#remve the redundant colums 
	for ( pop in populations ) {
		col_rm<-paste0(pop,"_RS_Number")
		merged_proxies[,col_rm]<-NULL
	}
return(merged_proxies)
}
################################### END FUNCTION ###################################

##########################################################################################
################################### FUNCTION ############################################
#find clumps/clusters with no proxies and if requested look up additional proxies

find_replace_empty_proxies<- function( proxyList) {

	if (any(is.na(names(proxyList)))) {
		message("proxyList must be a named list")
		return()
	}
	
	populations<-names(proxyList)
	no_ld_snps<-lapply ( populations, function(pop) {
		nvars<-sapply(1:length(proxyList[[pop]]), function(i) nrow(proxyList[[pop]][[i]]))
		one_var=which(nvars==1)
		error_vars<-grepl("error",proxyList[[pop]][one_var])
		error_vars<-one_var[error_vars]
		missing_vars=names(proxyList[[pop]])[error_vars]
		#missing_vars=p_val_lead_snp[error_vars]
		return(missing_vars)
	})
	names(no_ld_snps)<-c(populations)

	if ( length(unlist(no_ld_snps)) == 0 ) {
		both_missing2<-matrix(ncol=length(populations))
		colnames(both_missing2)<-populations
		return(both_missing2)
	}

	unique_missing_ld_snps<-unique(unlist(no_ld_snps))

	both_missing<-sapply(populations, function(pop) {
		unique_missing_ld_snps %in% no_ld_snps[[pop]]
	})

	both_missing2<-both_missing

	both_missing<-unique_missing_ld_snps[sapply(1:nrow(both_missing), function(i) all(both_missing[i,]))]
	# when the lead snp has no proxies try the next highest pval snp
	loci_SNPs2$proxy_lead<-NA

	both_missing_round<-list(rep(1,times=length(both_missing)))
	names(both_missing_round)<- both_missing

	n_pval_list<-lapply(both_missing, function(snp) {
		loc_ind<-which(loci_SNPs2$lead == snp)
		loc_snps<-loci_SNPs2[loc_ind,]
		if ( nrow(loc_snps) == 1 ) {
			warning (snp, "; is an orphan SNP ; having no other significiant SNPs in the region")
			#remove_ind=c(remove_ind,error_vars)
			return(NA)
		}
		lead_p<-unlist(filter(loc_snps,rsid == {{snp}}) %>% dplyr::select(P_value))
		pval_cuttoff<- -log(lead_p,10)/p_val_cuttoff_factor
		min_pval_cuttoff<- 10^-pval_cuttoff
		pval_cuttoff<- -log(lead_p,10)*p_val_cuttoff_factor
		max_pval_cuttoff<- 10^-pval_cuttoff
		p_order<-unlist(order(abs(loc_snps$P_value -  lead_p),decreasing=F))
		return(list(p_order,min_pval_cuttoff,max_pval_cuttoff))
		
	})
	names(n_pval_list)<-both_missing

	min_p<-sapply(n_pval_list, function(slot) slot[2])
	max_p<-sapply(n_pval_list, function(slot) slot[3])
	n_pval_list<-sapply(n_pval_list, function(slot) slot[1])

	#remove from both missing those with no other snps in the regions
	both_missing<-both_missing[!is.na(n_pval_list)]

	remove_ind<-c()
	next_pval_snps_rounds=list()
	while (length(both_missing) > 0) {
		# for any snp where we were unable to extract linked snps look up the nearest p-val snp in the same locus
		next_pval_snps<- sapply(both_missing, function(snp) { 
			lead_snp=unique(loci_SNPs2$lead[loci_SNPs2$rsid== snp ])
			loc_snps<-filter(loci_SNPs2 , lead == {{ lead_snp }})
			lead_p<-loci_SNPs2$P_value[loci_SNPs2$rsid == lead_snp ]
			next_snp<-loc_snps[n_pval_list[[snp]][both_missing_round[[snp]]+1],]
			next_pval<-next_snp$P_value

			if ( next_pval < min_p[[lead_snp]] | next_pval > max_p[[lead_snp]] ) { 
				warning(paste("lead snp ", snp, " has no other snps in the loci within the defined threshold"))
				#and what
			}
			next_pval_snp<-next_snp$rsid
			
			if ( str_sub (next_pval_snp,1,2) != "rs" ) {
				snp_pos<-unlist(str_split(next_pval_snp,":"))
				next_pval_snp<-paste0("chr",snp_pos[1],":",snp_pos[2])
			}
			return(next_pval_snp)
		})

		for ( i in 1:length(next_pval_snps)) {
			next_snp=next_pval_snps[[i]]
			old_snp<-both_missing[[i]]
			
			if ( is.na(next_snp) ) { 
				message(i," might be in an infinite loop")
				next
			}
			#update the object that tracks the snp name and order 
			both_missing_round[[old_snp]]<-both_missing_round[[old_snp]]+1
			names(both_missing_round) [ names(both_missing_round) == old_snp ] <-next_snp
			names(n_pval_list)[ names(n_pval_list) == old_snp ] <-next_snp
		}
		#create a list to keep track of the sequence of events
		next_pval_snps_rounds=append(next_pval_snps_rounds,next_pval_snps)
			
		#rextract the r2 info
		new_proxies<-lapply( populations , function(pop) {
			suppressMessages(LDproxy_batch(next_pval_snps, 
				pop = pop, 
				r2d = "r2", 
				token = API_token
			))
			new_proxies<-lapply(next_pval_snps, function(snp) {
				#new_proxy<-fread(paste0(snp,".txt"))
				new_proxy<-read.delim(paste0(snp,".txt"),header=T, stringsAsFactors=F)
				new_proxy$lead_snp<-snp
				return(new_proxy)
			})
			names(new_proxies)<-next_pval_snps
		})
		names(new_proxies)<-populations
		
		#update the proxy list
		#for ( x in 1:length(next_pval_snps) {
		proxyList<-lapply(populations, function(pop) { 
			proxyList[[pop]][names(next_pval_snps)]<-new_proxies[[pop]][next_pval_snps]
			return(proxyList[[pop]])
		})
		names(proxyList)<-populations
		
		no_ld_snps<-lapply ( populations, function(pop) {
			nvars<-sapply(1:length(proxyList[[pop]]), function(i) nrow(proxyList[[pop]][[i]]))
			one_var=which(nvars==1)
			error_vars<-grepl("error",proxyList[[pop]][one_var])
			error_vars<-one_var[error_vars]
			missing_vars=names(proxyList[[pop]])[error_vars]
			return(missing_vars)
			})

		names(no_ld_snps)<-populations
		unique_missing_ld_snps<-unique(unlist(no_ld_snps))
		both_missing<-sapply(populations, function(pop) {
			unique_missing_ld_snps %in% no_ld_snps[[pop]]
		})
		both_missing2<-both_missing
		both_missing<-unique_missing_ld_snps[sapply(1:nrow(both_missing), function(i) all(both_missing[i,]))]
	}
	rownames(both_missing2)<-unique_missing_ld_snps

return(list(both_missing2,new_proxies))

}
#################################### END FUNCTION ##############################################
###############################################################################################

##################################### FUNCTION ####################################################

##################################################################################
# adjust the ref alt format of the proxies and gwas snps before mergeing
# first I cant be bothered to make sence of highly complex loci, ie any multiallelic snps
# lets remove these
# this will also filter out any variants below the maf threshold

adjust_alleles<-function (mergedProxies,ref_genome) {

	split_alleles<-sapply(mergedProxies$Alleles, function(allele) str_split(allele,"/"))
	split_alleles<-lapply(split_alleles, function(allele) gsub("(\\(|\\))","",allele))
	not_2_alleles<-which( ! sapply(1:length(split_alleles), function(i) length(split_alleles[[i]])) == 2 )

	if ( length(not_2_alleles) > 0 ) { 
		cat ("the following variants have a complex structure and will be filtered\n")
		cat( "RS_Number","Coord","\n",sep="\t")
		for ( i in not_2_alleles) { 
			cat( mergedProxies$RS_Number[i],mergedProxies$Coord[i],"\n",sep="\t")
		}
		mergedProxies_old<-mergedProxies
		split_alleles_old<-split_alleles
		mergedProxies<-mergedProxies[-not_2_alleles,]
		split_alleles<-split_alleles[-not_2_alleles]
	}

	chr_pos<-str_split(mergedProxies$Coord, ":")
	mergedProxies$chr=sapply(chr_pos, function(loc) loc[1])
	mergedProxies$pos=sapply(chr_pos, function(loc) loc[2])
	mergedProxies$chr<-gsub("chr","",mergedProxies$chr)

	mergedProxies$REF<-unname(sapply(split_alleles,function(allele) allele[1]))
	mergedProxies$ALT<-unname(sapply(split_alleles,function(allele) allele[2]))

	# label as snp / indel 
	indel_TF=mergedProxies$REF == "-" | mergedProxies$ALT == "-"

	mergedProxies$snp_type <- NA
	mergedProxies$snp_type[indel_TF] <- "InDel"
	mergedProxies$snp_type[ ! indel_TF] <- "SNP"

	ref_snp<-mergedProxies$REF[ ! indel_TF]
	alt_snp<-mergedProxies$ALT[ ! indel_TF]
	pos=as.numeric(mergedProxies$pos[! indel_TF])
	chrom=mergedProxies$chr[! indel_TF]

	chr_ok<-substr(chrom,1,3)
	not_chr<-which(! chr_ok == "chr")
	if(length(not_chr) >0 ) {
		chrom2<-sapply(not_chr, function(row) {
			gsub( "^","chr",chrom[row])
			})
		chrom[not_chr]<-chrom2
	}

	snps<-subseq(ref_genome[chrom],start=pos,end=pos)
	snps_ref<-sapply(as.character(snps), function(snp) substr(snp,1,1))

	if (all(snps_ref == ref_snp)) message("first allele is always ref allele")

	# get the addtional base for indels from the ref genome
	# the mergedProxies have indels recored as eg (-/AA) for ins AA however gwas snps are recorded as eg (C/CAA)
	# for consistency update the former to the later 

	ref_snp<-mergedProxies$REF[ indel_TF]
	alt_snp<-mergedProxies$ALT[ indel_TF]
	pos=as.numeric(mergedProxies$pos[ indel_TF])
	chrom=mergedProxies$chr[ indel_TF]

	chr_ok<-substr(chrom,1,3)
	not_chr<-which(! chr_ok == "chr")
	if(length(not_chr) >0 ) {
		chrom2<-sapply(not_chr, function(row) {
			gsub( "^","chr",chrom[row])
			})
		chrom[not_chr]<-chrom2
	}

	snps<-subseq(ref_genome[chrom],start=pos,end=pos)
	snps_ref<-sapply(as.character(snps), function(snp) substr(snp,1,1))

	#changing the format of the ref alt allelels in the proxy file
	new_ref_snp<-sapply(1:length(snps_ref), function(i) paste0(snps_ref[i],gsub("-","",ref_snp[i])))
	new_alt_snp<-sapply(1:length(snps_ref), function(i) paste0(snps_ref[i],gsub("-","",alt_snp[i])))

	startpos<- pos
	ref_length<-nchar(new_ref_snp)-1
	endpos=pos+ref_length

	check_snps_seq<-subseq(ref_genome[chrom],start=pos,end=endpos)
	check_snps_seq<-as.character(check_snps_seq)

	not_ref_seq_ind<-which( ! check_snps_seq == new_ref_snp)
	if (length(not_ref_seq_ind) >0 ) {
		cat("oops unable to extract the reference base for converting annotation format for the following SNPs\n")
		cat(mergedProxies$RS_Number[ indel_TF][not_ref_seq_ind],sep="\n")
		stop ()
	}

	mergedProxies$REF[indel_TF]<-new_ref_snp
	mergedProxies$ALT[indel_TF]<-new_alt_snp

	chr_pos<-str_split(mergedProxies$Coord, ":")
	mergedProxies$chr=sapply(chr_pos, function(loc) loc[1])
	mergedProxies$pos=sapply(chr_pos, function(loc) loc[2])
	mergedProxies$chr<-gsub("chr","",mergedProxies$chr)
	return(mergedProxies)
}
############################################################################################
##################################### END FUNCTION #########################################

######################################## FUNCTION ########################################

add_r2<- function( snp_list, variants, proxies, populations) {

	if ( ! "data.frame" %in% class(variants)) stop ("variants must be a single data frame") 
	if ( ! "data.frame" %in% class(proxies)) stop ("proxies must be a single data frame")

	variants_r2<-lapply( snp_list, function(snp) {
	#message(snp)
	loc_snps<-filter(variants , lead == {{snp}})
	lead_p<-loc_snps$P_value[loc_snps$rsid == snp]
	match_snps<-paste(loc_snps$chr,loc_snps$pos,loc_snps$REF,loc_snps$ALT,sep="_")
	prox_snps<-filter(proxies, original_lead_snp == {{snp}})
	if ( nrow(prox_snps ) == 0) { 
		proxy_lead <- NA 
		pop_r2<-sapply(populations, function (pop) {
			return(rep(NA,times=nrow(loc_snps)))
		})
	} else {
		proxy_lead <-unique(prox_snps$lead_snp)
		pop_r2<-sapply(populations, function (pop) {
			r2_col<-paste0(pop,"_R2")
			snp_proxies_match_ind<-match(loc_snps$rsid,prox_snps$RS_Number)
			match_proxies<-paste(prox_snps$chr,prox_snps$pos,prox_snps$REF,prox_snps$ALT,sep="_")
			unmatched_ind<-which( is.na(snp_proxies_match_ind))
			snp_proxies_rematch_ind<-match(match_snps,match_proxies)
			check_TF<- (!is.na(snp_proxies_match_ind)) & (! is.na(snp_proxies_rematch_ind))
			if ( ! all(snp_proxies_match_ind[check_TF] == snp_proxies_rematch_ind[check_TF])) stop ( "mismatching between SNPs and proxies")
			return(prox_snps[snp_proxies_match_ind,r2_col])	
		})
	}
	pop_cols<-paste0(populations,"_R2")
	if ( nrow( loc_snps) == 1 )	pop_r2<-t(as.data.frame(pop_r2, stringsAsFactors =F))
	colnames(pop_r2) <-pop_cols
	
	pop_r22<-lapply(pop_cols, function(col) pop_r2[,col])
	names(pop_r22)<-pop_cols
	loc_snps[, (pop_cols) := pop_r22]
	
	loc_snps$proxy_lead<-proxy_lead
	if (proxy_lead != snp | is.na(proxy_lead))  warning ( paste("no proxy lead",snp , proxy_lead ))
	return(loc_snps)
	})
return(variants_r2)
}
	
############################################################################

######################################## FUNCTIONS #################################

find_unlinked_snps<- function(snp_list,lead_snps,populations,lower_r2_cuttoff) {
	both_NA_tf<-lapply( lead_snps, function(snp) {
		loc_snps<-filter(loci_SNPs2 , lead == {{snp}})
		pop_r2<-sapply(populations, function (pop) {
		r2_col<-paste0(pop,"_R2")
		(is.na(loc_snps[,r2_col]) | loc_snps[,r2_col] < lower_r2_cuttoff) 
		})
		if ( nrow( loc_snps) == 1 )	pop_r2<-t(as.data.frame(pop_r2, stringsAsFactors =F))
		# return a T / F for each other snps in the loci depending on whether there is any linkage to the lead snp
		sapply(1:nrow(pop_r2),function(i) all(pop_r2[i,]))
		# if all snps are linked 
	})
	names(both_NA_tf)<-lead_snps
	return(both_NA_tf)
}


#########################################################################

append_to_column<- function(input_vector, string,index_to_append_to=NULL) {
	if (is.null(index_to_append_to)) index_to_append_to=c(1:length(input_vector)) 
	is_na_TF<-is.na(input_vector[index_to_append_to])
	cols_update<-sapply(1:length(is_na_TF), function(i) {
		TF=is_na_TF[i]
		main_ind<-index_to_append_to[i]
		ifelse(TF,string,paste(input_vector[main_ind],string,sep=";"))
		})
		input_vector[index_to_append_to]<- cols_update
		return(input_vector)
}

#############################################
#
# END FUNCTIONS
#
#############################################


###################################################################################################################
################################################## part 1 #########################################################

#delete duplicated columns
rsid_cols<-which(colnames(main_SNPs) == "rsid")
if ( length(rsid_cols) > 1) {
main_SNPs[,rsid_cols[2]] <- NULL
}

message("standardising SNP format")
message("	checking variants against dbSNP")

#extract all snps around lead snp based on the r2 values provided 
lead_snps<-unlist(filter(summary_SNPs, P < hard_P_val_cutoff) %>%  dplyr::select(RSID))

missing_leads<-lead_snps[! lead_snps %in% main_SNPs$rsid]

if (length(missing_leads) > 0 ) {
	message(cat(missing_leads,sep="\n"),"are missing from the all snps file")
	stop("script exiting")
}

#filter snps on maf 
if ( "RAF" %in% colnames(main_SNPs)) {
	main_SNPs<- main_SNPs[(main_SNPs$RAF > min_RAF & main_SNPs$RAF < 1-min_RAF) | is.na(main_SNPs$RAF),]
}

missing_leads<-lead_snps[! lead_snps %in% main_SNPs$rsid]

if (length(missing_leads) > 0 ) {
	warning("the above SNPs have been filtered due to '-m, --minor_allele_frequency' threshold")
	cat(missing_leads,sep="\n")
}

lead_snps<-lead_snps [ ! lead_snps %in% missing_leads]
summary_SNPs<-summary_SNPs[! summary_SNPs$RSID %in% missing_leads,]

main_SNPs.gr<-makeGRangesFromDataFrame(main_SNPs, start.field= "pos" , end.field="pos" ,keep.extra.columns=T)

# make a granges of the lead snps
lead_snps.gr<-GRanges(summary_SNPs$CHR,IRanges(start=summary_SNPs$POS,end=summary_SNPs$POS))
lead_snps.gr$rsid=summary_SNPs$RSID

#message(loci_window)

loci_snps_position<-lapply(lead_snps , function(snp) {
	snp_row<-main_SNPs[main_SNPs$rsid == snp ]
	loci.gr<-GRanges(snp_row$chr, IRanges(start=max(c(snp_row$pos - loci_window/2,0)),end=snp_row$pos +loci_window/2))
	loci_snp.gr<-subsetByOverlaps(main_SNPs.gr,loci.gr)
	loci_snp.gr$lead=snp
	return(loci_snp.gr)
})
names(loci_snps_position)<-lead_snps

loci_snps_position.gr<-c_granges(loci_snps_position)
names(loci_snps_position.gr)<- NULL
loci_SNPs2<-as.data.frame(loci_snps_position.gr,stringsAsFactors=F)

colnames(loci_SNPs2)[colnames(loci_SNPs2) == "seqnames"]<-"chr"
colnames(loci_SNPs2)[colnames(loci_SNPs2) == "start"]<-"pos"

loci_SNPs2<-remove.factors(loci_SNPs2)

loci_SNPs2$end= NULL
loci_SNPs2$width= NULL

# first of all we will check the snps have the correct alleles and occur in snpdb
#read in the hg37 vcf 
message("	reading in SNPdb file")
dbSNPs<-readVcf(snp_db_file)
dbSNPs.gr<-rowRanges(dbSNPs)
rm(dbSNPs)

loci_SNPs2$REF<-NA
loci_SNPs2$ALT<-NA


matched_rsid<-match(loci_SNPs2$rsid,names(dbSNPs.gr))
names(matched_rsid)<-1:length(matched_rsid)
unmatched_ids<-names(matched_rsid)[is.na(matched_rsid)]
if (any (is.na(matched_rsid))) {
	missing_id_SNPs<-loci_SNPs2[as.numeric(unmatched_ids),]
	missing_id_SNPs.gr<-GRanges(missing_id_SNPs$chr,IRanges(start=missing_id_SNPs$pos,end=missing_id_SNPs$pos))
	dbsnps.ol<-findOverlaps(missing_id_SNPs.gr,dbSNPs.gr)
	number_mappings<-table(queryHits(dbsnps.ol))
	if ( any ( number_mappings > 1)) {
		warning( for (var in as.numeric(names(number_mappings )[number_mappings > 1])) {
			cat (missing_id_SNPs$rsid[as.numeric(var)],"\n") 
		},"the above snps have more than one variant at the same pos")
	cat("\n")
	}

	ind<-as.numeric(unmatched_ids)[queryHits(dbsnps.ol)]

	missing_id_allele_A<-loci_SNPs2$EA[ind]
	missing_id_allele_B<-loci_SNPs2$NEA[ind]

	ref<-as.character(dbSNPs.gr$REF[subjectHits(dbsnps.ol)])
	alt_list<-as.list(dbSNPs.gr$ALT[subjectHits(dbsnps.ol)])
	alt_list<-lapply(alt_list, function(alleles) as.character(alleles))

	same_alleles<-((missing_id_allele_A == ref | sapply(1:length(ind), function(i) any(missing_id_allele_A[i] == alt_list[[i]]))) &
	(missing_id_allele_B == ref | sapply(1:length(ind), function(i) any(missing_id_allele_B[i] == alt_list[[i]]))))

	new_names<-names(dbSNPs.gr)[subjectHits(dbsnps.ol)][same_alleles]
	new_REF=ref[same_alleles]
	
	#new_ind<-as.numeric(unmatched_ids)[as.numeric(names(number_mappings))]
	new_ind<-ind[same_alleles]

	new_ALT<-sapply( 1:length(new_ind), function(i) { 
		ind=new_ind[i] 
		ref_alt<-c(loci_SNPs2$EA[ind],loci_SNPs2$NEA[ind])
		ref_alt[! ref_alt %in% new_REF[i]]
	})

	ref_alt_ind<-new_ind[new_REF != new_ALT]

	if( any (new_REF == new_ALT) ) {
		no_ref_alt_ind<-new_ind[new_REF == new_ALT]
		cat ( "unable to extract ref alt for ;", "\n")
		cat (sep="\n", loci_SNPs2$rsid[no_ref_alt_ind])
		cat("\n")
	}	

	loci_SNPs2$old_rsid<-NA
	loci_SNPs2$old_rsid[new_ind]<-loci_SNPs2$rsid[new_ind]
	loci_SNPs2$rsid[new_ind]<-new_names
	loci_SNPs2$REF[ref_alt_ind]<-new_REF[new_REF != new_ALT]
	loci_SNPs2$ALT[ref_alt_ind]<-new_ALT[new_REF != new_ALT]

	all_missing<-1:length(unmatched_ids) 
	not_remapped<-all_missing [ ! all_missing %in% queryHits(dbsnps.ol)]
	not_remapped_ind<-as.numeric(unmatched_ids[not_remapped])
	# now remove no match snps if desired
	if ( keep_unmapped_variants == T ) 
	{
		message( "the following snps are not in the reference file but will be retained") 
		cat("\n")
		cat(sep="\n",loci_SNPs2$rsid[not_remapped_ind])
		remove_ind=c()
		cat("\n")
	} else if (keep_unmapped_variants == F ){
		message( "the following snps are not in the reference file and will be filtered") 
		cat("\n")
		cat(sep="\n",loci_SNPs2$rsid[not_remapped_ind])
		remove_ind<-not_remapped_ind
		cat("\n")
	}
	
}

if ( length(new_names) > 0 ) {
	message ("	updating old rsIDs")
	for ( i in ind ) {
		message("swapped ", loci_SNPs2$old_rsid[i], " for ", loci_SNPs2$rsid[i])	
	}	
	cat("\n")
}

# this next piece of text will extract the ref and alt columns from the snpdb file
names(matched_rsid)= 1:length(matched_rsid)
matched_rsid<-matched_rsid[! is.na(matched_rsid)]

# get the allele from the input snp file
allele_A<-loci_SNPs2$EA[as.numeric(names(matched_rsid))] 
allele_B<-loci_SNPs2$NEA[as.numeric(names(matched_rsid))] 

# get the ref and alt from the refdb file
ref<-as.character(dbSNPs.gr$REF[matched_rsid[names(matched_rsid)]])
alt_list<-as.list(dbSNPs.gr$ALT[matched_rsid[names(matched_rsid)]])
alt_list<-lapply(alt_list, function(alleles) as.character(alleles))

# check if the aleles match
same_alleles<-((allele_A == ref | sapply(1:length(allele_A), function(i) allele_A[i] %in% alt_list[[i]])) &
	(allele_B == ref | sapply(1:length(allele_A), function(i) allele_B[i] %in% alt_list[[i]])) & 
	(! sapply(1:length(allele_A), function(i) allele_A[i] %in% alt_list[[i]] & allele_B[i] %in% alt_list[[i]]))
	)

#add ref alt columns allele columns 
new_REF=ref[same_alleles]
new_ind<-as.numeric(names(matched_rsid))[same_alleles]

new_ALT<-sapply( 1:length(new_ind), function(i) { 
	ind=new_ind[i] 
	ref_alt<-c(loci_SNPs2$EA[ind],loci_SNPs2$NEA[ind])
	ref_alt[! ref_alt %in% new_REF[i]]
})

ref_alt_ind<-new_ind[new_REF != new_ALT]

if( any (new_REF == new_ALT) ) {
	no_ref_alt_ind<-new_ind[new_REF == new_ALT]
	cat ( "unable to extract ref alt for ;", "\n")
	cat (sep="\n", loci_SNPs2$rsid[no_ref_alt_ind])
}	

loci_SNPs2$REF[ref_alt_ind]<-new_REF
loci_SNPs2$ALT[ref_alt_ind]<-new_ALT

#check all the alleles match
all_mismatches<- ! sapply(ref_alt_ind , function(ind) {
	alleleA=loci_SNPs2$EA[ind]
	alleleB=loci_SNPs2$NEA[ind]
	old_alleles=c(alleleA,alleleB)
	ref=loci_SNPs2$REF[ind]
	alt=loci_SNPs2$ALT[ind]
	ref_alt_alleles=c(ref,alt)
	all(ref_alt_alleles %in% old_alleles)
})

if ( any (all_mismatches) ) { 
	#loci_SNPs2[ref_alt_ind[all_mismatches],]
	warning( for (i in which(all_mismatches)) {
		ind=ref_alt_ind[i]
		cat (loci_SNPs2$rsid[ind],"\n") 
		},"The provided REF and ALT alleles for the above SNP dont match the dbSNP file. The SNP will be filtered"
	)
	remove_ind<-c(remove_ind,ref_alt_ind[which(all_mismatches)])
	cat("\n")
}
##########################

if ( any (! same_alleles)) { 
	message( "the following snps have miss-matched alleles") 
	remove_ind<- c(remove_ind,as.numeric(names(matched_rsid))[! same_alleles])
	cat("\n")
	cat(sep="\n",loci_SNPs2$rsid[remove_ind])
	cat("\n")
}

##############################################
# some snps legit dont have the correct alleles 
# so take the ref db alleles 
#some alleles are not recorded consistently eg AT : ATT -> A : AT 

can_fix_ind<-remove_ind[ loci_SNPs2$rsid[remove_ind] %in% names(dbSNPs.gr)]

gwas_rsid=loci_SNPs2$rsid[can_fix_ind]
dbsnp_rsid<-names(dbSNPs.gr)[matched_rsid[as.character(can_fix_ind)]]

if ( ! all( gwas_rsid == dbsnp_rsid )) stop( "the object matching gwas and db snps is not opperating correctly") 

ref_sub<-as.character(dbSNPs.gr$REF[matched_rsid[as.character(can_fix_ind)]])
alt_sub<-as.list(dbSNPs.gr$ALT[matched_rsid[as.character(can_fix_ind)]])
alt_sub<-lapply(alt_sub, function(alleles) as.character(alleles))

bi_allelic_TF<-sapply(alt_sub, function(alt) length(alt) == 1)

can_fix_ind<-can_fix_ind[bi_allelic_TF]
remove_ind<-c(remove_ind[ ! remove_ind %in% can_fix_ind])

loci_SNPs2$REF[can_fix_ind]<-ref_sub[bi_allelic_TF]
loci_SNPs2$ALT[can_fix_ind]<-unlist(alt_sub[bi_allelic_TF])

#########################################################################
# for those snps without not in snpdb we can ref alt we can deduce this from the sequence 
# extract sequence for both alleles for both EA nad NEA and then compare to the ref sequenece

#loci_SNPs2<-remove.factors(loci_SNPs2)
no_ref_alt_ind<-which(is.na(loci_SNPs2$REF))
no_ref_alt_ind<-no_ref_alt_ind[! no_ref_alt_ind %in% remove_ind]

if ( length( no_ref_alt_ind) > 0) {
	#read in the refgenome
	message( "	reading in the ref genome")
	ref_genome<-readDNAStringSet(genome_file, format = "fasta")
	#lapply(no_ref_alt_ind, function(ind) {
	pos=loci_SNPs2$pos[no_ref_alt_ind]
	pos=as.numeric(pos)
	chrom=loci_SNPs2$chr[no_ref_alt_ind]
	chr_ok<-substr(chrom,1,3)
	not_chr<-which(! chr_ok == "chr")
	if(length(not_chr) >0 ) {
		chrom2<-sapply(not_chr, function(row) {
			gsub( "^","chr",chrom[row])
			})
		chrom[not_chr]<-chrom2
	}

	snps<-subseq(ref_genome[chrom],start=pos,end=pos)
	snps_ref<-sapply(as.character(snps), function(snp) substr(snp,1,1))

	EA_is_ref<-snps_ref == loci_SNPs2$EA[no_ref_alt_ind]
	NEA_is_ref<-snps_ref == loci_SNPs2$NEA[no_ref_alt_ind]

	ok_ref<-which(EA_is_ref | NEA_is_ref)
	ok_ref_ind<-no_ref_alt_ind[ok_ref]

	missing_ref<-no_ref_alt_ind[ ! no_ref_alt_ind %in% ok_ref_ind ]

	if ( length(missing_ref) > 0 ) {
		cat("unable to find matching allele info in ref genome for\n")
		cat (loci_SNPs2$rsid[missing_ref],sep="\n") 
		remove_ind<-c(remove_ind,missing_ref)
	}

	EA_is_ref_ind<-no_ref_alt_ind[EA_is_ref]
	NEA_is_ref_ind<-no_ref_alt_ind[NEA_is_ref]

	loci_SNPs2$REF[EA_is_ref_ind]<-loci_SNPs2$EA[EA_is_ref_ind]
	loci_SNPs2$ALT[EA_is_ref_ind]<-loci_SNPs2$NEA[EA_is_ref_ind]

	loci_SNPs2$REF[NEA_is_ref_ind]<-loci_SNPs2$NEA[NEA_is_ref_ind]
	loci_SNPs2$ALT[NEA_is_ref_ind]<-loci_SNPs2$EA[NEA_is_ref_ind]
}

wtf<-which(nchar(loci_SNPs2$REF) > 1 & nchar(loci_SNPs2$ALT) > 1)
if ( length(wtf) > 0) {
	cat (" both alleles of the following variants are greater then 1 base in lenght\n")
	cat (" this may lead to unusual behavior and they will be removed\n")
	cat (loci_SNPs2$rsid[wtf],sep="\n") 
	remove_ind<-unique(c(remove_ind,wtf))
	cat("\n")
}

#remove any unmatched snps
if ( length(remove_ind) > 0) {
	cat ("the following variants have been removed\n")
	cat (loci_SNPs2$rsid[remove_ind],sep="\n") 
	loci_SNPs2<-loci_SNPs2[-remove_ind,]
	cant_map_filt<-loci_SNPs2[remove_ind,]
	cant_map_filt$comments<-"cant_map_SNP"
	cat("\n")
}

## any remaining snps with no allele information
# just aribarilly assign the ref alt
no_ref_alt_ind<-which(is.na(loci_SNPs2$REF))

if ( length(no_ref_alt_ind) >  0) {
	cat (" unable to extract the ref alt alleles for the following SNPs\n")
	cat (loci_SNPs2$rsid[no_ref_alt_ind],sep="\n")
	cat (" this may lead to unusual behavior\n")
	loci_SNPs2$REF[no_ref_alt_ind]<- loci_SNPs2$EA[no_ref_alt_ind]
	loci_SNPs2$ALT[no_ref_alt_ind]<- loci_SNPs2$NEA[no_ref_alt_ind]
	cat("\n")
}

# add the snp type field 
indels_TF<-nchar(loci_SNPs2$REF) > 1 | nchar(loci_SNPs2$ALT) > 1

# add the snp_type
loci_SNPs2$snp_type<-"SNP"
loci_SNPs2$snp_type[indels_TF]<- "InDel"

message ("finished mapping snps")
cat("\n")

############################################# finished part1 mapping snps #################################################

##########################################################################################################################
################################################## part 2 ################################################################

message ("Part 2") 
message ("	splitting loci based on linkage and p-value filtering")
cat("\n")

# download r2 info for the lead snps
# split regions into linked loci 
# filter by pvalue 
# now re import linked snps that may have been filtered in gwas processing
# get all the proxies 
# repeat for each population used

# try to define the locus snps based on ld 
#fist extract ld scores
wk_dir=getwd()

if (split_by_r2 == T) {
	if(! exists("proxy_file_list")) { 
		proxylist<-lapply(populations, function(pop) {
			prox_dir=file.path(wk_dir,pop)
			if (! dir.exists(prox_dir))  { 
				dir.create(prox_dir)
			}
			setwd(prox_dir)
			suppressMessages(LDproxy_batch(lead_snps,pop=pop,token=API_token))
			pop_proxies<-lapply(lead_snps,function(snp) {
				#print(snp)
				if (file.exists(paste0(snp,".txt"))) {
					snp_proxies<-read.delim(paste0(snp,".txt"),header=T,stringsAsFactors=F)
					snp_proxies$lead_snp<-snp
					return(list(snp_proxies,snp))
				}
			})
			pop_proxies<- lapply(1:length(pop_proxies) , function(i) pop_proxies[[i]][[1]])
			names(pop_proxies)<-lead_snps
			return(pop_proxies)
		})
	} else { 
		proxylist<- lapply(populations, function(pop) {
			load(proxy_file_list[[pop]])
			return(proxylist)
		})
	}
	setwd(wk_dir)
	names(proxylist)<-populations

	#######################################################
	#find clumps/clusters with no proxies

	no_ld_snps<-lapply ( populations, function(pop) {
		nvars<-sapply(1:length(proxylist[[pop]]), function(i) nrow(proxylist[[pop]][[i]]))
		missing_vars=names(proxylist[[pop]])[sapply(nvars, is.null)]
		return(missing_vars)
	})
	names(no_ld_snps)<-c(populations)

	unique_missing_ld_snps<-unique(unlist(no_ld_snps))

	both_missing<-sapply(populations, function(pop) {
		unique_missing_ld_snps %in% no_ld_snps[[pop]]
	})

	both_missing2<-both_missing

	both_missing<-unique_missing_ld_snps[sapply(1:nrow(both_missing), function(i) all(both_missing[i,]))]
	# when the lead snp has no proxies try the next highest pval snp
	loci_SNPs2$proxy_lead<-NA

	both_missing_round<-list()
	for ( i in 1:length(both_missing)) { 
	both_missing_round[i] <- 1
	}
	names(both_missing_round)<- both_missing

	n_pval_list<-lapply(both_missing, function(snp) {
		loc_ind<-which(loci_SNPs2$lead == snp)
		loc_snps<-loci_SNPs2[loc_ind,]
		if ( nrow(loc_snps) == 1 ) {
			warning (snp, "; is an orphan SNP ; having no other significiant SNPs in the region")
			#remove_ind=c(remove_ind,error_vars)
			return(NA)
		}
		lead_p<-unlist(filter(loc_snps,rsid == {{snp}}) %>% dplyr::select(P_value))
		pval_cuttoff<- -log(lead_p,10)/p_val_cuttoff_factor
		min_pval_cuttoff<- 10^-pval_cuttoff
		pval_cuttoff<- -log(lead_p,10)*p_val_cuttoff_factor
		max_pval_cuttoff<- 10^-pval_cuttoff
		p_order<-unlist(order(abs(loc_snps$P_value -  lead_p),decreasing=F))
		return(list(p_order,min_pval_cuttoff,max_pval_cuttoff))
	})
	names(n_pval_list)<-both_missing

	min_p<-sapply(n_pval_list, function(slot) slot[2])
	max_p<-sapply(n_pval_list, function(slot) slot[3])
	n_pval_list<-sapply(n_pval_list, function(slot) slot[1])

	#remove from both missing those with no other snps in the regions
	both_missing<-both_missing[!is.na(n_pval_list)]

	remove_ind<-c()
	next_pval_snps_rounds=list()
	while (length(both_missing) > 0) {
		# for any snp where we were unable to extract linked snps look up the nearest p-val snp in the same locus
		next_pval_snps<- sapply(both_missing, function(snp) { 
			lead_snp=unique(loci_SNPs2$lead[loci_SNPs2$rsid== snp ])
			loc_snps<-filter(loci_SNPs2 , lead == {{ lead_snp }})
			lead_p<-loci_SNPs2$P_value[loci_SNPs2$rsid == lead_snp ]
			#loc_snps<-filter(loci_SNPs2 , lead == {{ snp }})
			nround=both_missing_round[[snp]]+1
			if ( nround > nrow(loc_snps)) return (NA)
			next_snp<-loc_snps[n_pval_list[[snp]][nround],]
			next_pval<-next_snp$P_value
			next_pval_snp<-next_snp$rsid
			
			if ( next_pval < min_p[[lead_snp]] | next_pval > max_p[[lead_snp]] ) { 
				warning(paste("lead snp ", snp, " has no other snps in the loci within the defined threshold"))
			}
			if ( str_sub (next_pval_snp,1,2) != "rs" ) {
				snp_pos<-unlist(str_split(next_pval_snp,":"))
				next_pval_snp<-paste0("chr",snp_pos[1],":",snp_pos[2])
			}
			return(next_pval_snp)
		})

		na_TF<-is.na(next_pval_snps)
		if(all(na_TF)) {
			for ( i in which(na_TF)) {
				warning(paste("None of the SNPs in the loci defined by the lead snp", names(next_pval_snps)[na_TF], "has any ld proxy data, are you sure you have choosen the correct population(s)"))
			}
			if ( strict == T ) {
				stop("as strict is TRUE, the script will now exit")
			} else {
				break
			}
		}

		for ( i in 1:length(next_pval_snps)) {
			if ( is.na(next_pval_snps[i])) next
			next_snp=next_pval_snps[[i]]
			old_snp<-both_missing[[i]]

			#update the object that tracks the snp name and order 
			both_missing_round[[old_snp]]<-both_missing_round[[old_snp]]+1
			next_pval_snps_rounds[[old_snp]]=append(next_pval_snps_rounds[[old_snp]],next_snp)
		}
		#create a list to keep track of the sequence of events

		#rextract the r2 info
		new_proxies<-lapply( populations , function(pop) {
			prox_dir=file.path(wk_dir,pop)
			setwd(prox_dir)
			got_data_tf<-sapply(next_pval_snps, function(snp) file.exists(paste0(snp,".txt")))
			suppressMessages(LDproxy_batch(next_pval_snps[ ! got_data_tf], 
				pop = pop, 
				r2d = "r2", 
				token = API_token
			))
			new_proxies<-lapply(next_pval_snps, function(snp) {
				if ( file.exists(paste0(snp,".txt"))) {
					#new_proxy<-suppressWarnings(fread(paste0(snp,".txt")))
					new_proxy<-read.delim(paste0(snp,".txt"),header=T,stringsAsFactors=F)
					new_proxy$lead_snp<-snp
					return(new_proxy)
				}
			})
			names(new_proxies)<-next_pval_snps
			return(new_proxies)
		})
		setwd(wk_dir)
		names(new_proxies)<-populations

		#update the proxy list
		proxylist<-lapply(populations, function(pop) { 
			not_null<-sapply(1:length(new_proxies[[pop]]), function(i) ! is.null(unlist(new_proxies[[pop]][i])))
			if ( any(not_null)) {  
				snps_add<-new_proxies[[pop]][ sapply(1:length(new_proxies[[pop]]), function(i) ! is.null(unlist(new_proxies[[pop]][i])))]
				#snps_add<-lapply(1:length(snps_add), function(i) snps_add[[i]][,2:12])
				names_add<-next_pval_snps[sapply(1:length(new_proxies[[pop]]), function(i) ! is.null(unlist(new_proxies[[pop]][i])))]
				proxylist[[pop]][names(names_add)]<-snps_add
			}
			return(proxylist[[pop]])
		})
		names(proxylist)<-populations

		no_ld_snps<-lapply ( populations, function(pop) {
			nvars<-sapply(1:length(proxylist[[pop]]), function(i) nrow(proxylist[[pop]][[i]]))	
			missing_vars=names(proxylist[[pop]])[sapply(nvars, is.null)]
			return(missing_vars)
			})

		names(no_ld_snps)<-populations
		unique_missing_ld_snps<-unique(unlist(no_ld_snps))
		both_missing<-sapply(populations, function(pop) {
			unique_missing_ld_snps %in% no_ld_snps[[pop]]
		})
		both_missing2<-both_missing
		if ( length(populations) > 1) {
			both_missing<-unique_missing_ld_snps[sapply(1:nrow(both_missing), function(i) all(both_missing[i,]))]
		} else {
			both_missing<-unique_missing_ld_snps[sapply(1:length(both_missing), function(i) all(both_missing[i]))]
		}
	}

	if ( typeof(both_missing2) == "logical") {
		both_missing3<-as.data.frame(both_missing2)
		colnames(both_missing3)<-names(both_missing2)
		both_missing2<-both_missing3
		rm(both_missing3)
	}
	rownames(both_missing2)<-unique_missing_ld_snps
	##########################################################################

	#merge the proxies
	merged_proxies<-merge_proxies(proxylist,both_missing2)
	#adjust the ref_alt_allele 
	merged_proxies<-adjust_alleles(merged_proxies,ref_genome)

	#create a ref for the snp that were subbed as proxies
	proxy_swaps<-sapply(lead_snps,function(snp) {
		unique(filter(merged_proxies, original_lead_snp == {{snp}}) %>% dplyr::select(lead_snp))
	})
	names(proxy_swaps)<- lead_snps

	# make any leads without proxies NA
	proxy_swaps[sapply(proxy_swaps, length) == 0]<-NA

	##########################################################
	# extract r2 values from proxy lists 
	# this next bit will attempt to define loci that have unlinked snps in them
	# where these can be attributed to other loci they are dropped else a new loci is produced 
	# this happens recussivly untill they are all separated out or the max number of itterations is reached
	regions_to_split=T
	lead_snps_list=list()
	loci_snps_list=list()
	cycle=0

	##########
	## start the while loop
	## no proxy lead
	#######
	loci_snps_orig <-loci_SNPs2 # store the original snp for later if ness
	merged_proxies_orig<-merged_proxies
	lead_snps_orig<-lead_snps

	loci_snps_r2<-lapply( lead_snps, function(snp) {
		#message(snp)
		loc_snps<-filter(loci_SNPs2 , lead == {{snp}})
		loc_snps$prev_lead=NA
		if ( snp %in% merged_proxies$original_lead_snp) {
			lead_p<-loc_snps$P_value[loc_snps$rsid == snp]
			match_snps<-paste(loc_snps$chr,loc_snps$pos,loc_snps$REF,loc_snps$ALT,sep="_")
			prox_snps<-filter(merged_proxies, original_lead_snp == {{snp}})
			proxy_lead<-unique(prox_snps$lead_snp)
			pop_r2<-sapply(populations, function (pop) {
				r2_col<-paste0(pop,"_R2")
				loc_snps[,r2_col] <- NA
				snp_proxies_match_ind<-match(loc_snps$rsid,prox_snps$RS_Number)
				match_proxies<-paste(prox_snps$chr,prox_snps$pos,prox_snps$REF,prox_snps$ALT,sep="_")
				unmatched_ind<-which( is.na(snp_proxies_match_ind))
				snp_proxies_rematch_ind<-match(match_snps,match_proxies)
				check_TF<- (!is.na(snp_proxies_match_ind)) & (! is.na(snp_proxies_rematch_ind))
				if ( ! all(snp_proxies_match_ind[check_TF] == snp_proxies_rematch_ind[check_TF])) stop ( "mismatching between SNPs and proxies")
				return(prox_snps[snp_proxies_match_ind,r2_col])	
			})
			
			if ( nrow( loc_snps) == 1 )	pop_r2<-t(as.data.frame(pop_r2, stringsAsFactors =F))
		} else {
			pop_r2<-matrix(ncol=length(populations),nrow=nrow(loc_snps))
			proxy_lead=NA
		}
		pop_cols<-paste0(populations,"_R2")
		colnames(pop_r2) <-pop_cols
		loc_snps[,pop_cols]<-pop_r2

		loc_snps$proxy_lead<-proxy_lead
		#if (proxy_lead != snp)  print ( paste("no proxy lead",snp , proxy_lead ))
		return(loc_snps)
	})
	names(loci_snps_r2) <-lead_snps

		# any snps without any proxy data or with an ld score lower then the 'lower_r2_cuttoff' are separated in to a separte locus
		# the lead snp of this seconary locus used to extract more proxy data
	cat("\n")
	while ( regions_to_split==T ) {
	# this next piece of code is very inefficient. Especially i should specifiy which loci have been optimsied (as in split as much as poss)
	# these could then be skipped and not reassessed
		cycle=cycle+1
		message("on cycle ", cycle , " of ", cycle_limit)
		loci_snps_r2<-lapply( lead_snps, function(snp) {
			#print(snp)
			loc_snps<-as.data.frame(loci_snps_r2[[snp]])
			if ( cycle == 1 ) { 
				loc_snps$orig_lead<-snp
				loc_snps$assessed=F
			}
			# if this loci cant be split just return
			if ( all(loc_snps$assessed)) { 
				#print ("wtf")
				return(list(list(loc_snps),c(snp),F))
			}

			pop_r2<-sapply(populations, function (pop) {
				r2_col<-paste0(pop,"_R2")
				(is.na(loc_snps[,r2_col]) | loc_snps[,r2_col] < lower_r2_cuttoff) 
			})
			if ( nrow( loc_snps) == 1 )	{
				loc_snps$assessed=T
				return(list(list(loc_snps),c(snp),F))
			}
			# return a T / F for each other snps in the loci depending on whether there is any linkage to the lead snp in any population
			both_na_tf<-sapply(1:nrow(pop_r2),function(i) all(pop_r2[i,]))
			# if all snps are linked | if none are
			if ( (! any(both_na_tf)) | all(both_na_tf) ) {
				loc_snps$assessed=T
				return(list(list(loc_snps),c(snp),F))
			}
			isSNP<-loc_snps$rsid ==  snp
			linked_snps<-loc_snps[ (! both_na_tf) | isSNP, ]
			unlinked_snps<- loc_snps[ both_na_tf & ( ! isSNP),]

			if ( nrow(unlinked_snps) == 0 ) {
				loc_snps$assessed=T
				return(list(list(loc_snps),c(snp),F))	
			}
			#####################
			#unlinked_snps<-find_unlinked_snps(loci_SNPs2,lead_snps,populations,lower_r2_cuttoff)
			################
			filter_unlinked_snps=T
			if (nrow(unlinked_snps) == 1 ) {
				#loc_snps$assessed=T
				linked_snps$assessed=T
				unlinked_snps$assessed=T
				return(list(list(linked_snps,unlinked_snps),c(snp,unlinked_snps$rsid),F))
			}
			#if threre are unlinked snps if they are part of another locus then 
			if ( ! any( unlinked_snps$rsid %in% lead_snps)) {
				# as we are going to look up snps in LDlinkR we need to make sure they exist and are biallelic 
				# this will test both until we get a hit
				snp_ok=FALSE
				round=1
				while ( snp_ok ==F) {
					p_order<-order(unlinked_snps$P_value,decreasing=F)
					if ( round > nrow(unlinked_snps)) {
						loc_snps$assessed=T
						unlinked_snps$assessed=T
						return(list(list(loc_snps,unlinked_snps),c(snp,unlinked_snps$rsid[p_order[1]]),F))
					}
					if ( round > 10 ) message( snp , " none of the top 10 snps are bi allelic or are in the snpdb file")
					min_p_ind<-p_order[round]
					min_p_snp<-unlinked_snps[min_p_ind,]
					snp_ind= match(min_p_snp$rsid,names(dbSNPs.gr))
					if ( is.na(snp_ind)) {
						round=round+1
						next()
					}
					alt=as.character(unlist(mcols(dbSNPs.gr[snp_ind])$ALT))
					if 	( length(alt) != 1) {
						round=round+1
						next()
					}
					snp_ok =T
				}
				min_p_snp_wide.gr<-GRanges(min_p_snp$chr, IRanges(start=max(min_p_snp$pos-loci_window/2,0),end=min_p_snp$pos+loci_window/2))
				unlinked_lead_ol<-findOverlaps(lead_snps.gr[lead_snps.gr$rsid != snp],min_p_snp_wide.gr)
				if (length(unlinked_lead_ol) == 0) { 
					filter_unlinked_snps=F
				} else {
					pos_linked_leads<-lead_snps.gr$rsid[queryHits(unlinked_lead_ol)]
					r2_res="NULL"
					while ( r2_res == "NULL" & (!is.na(r2_res))) {
						r2_res=NA
						r2_out<-tryCatch({ LDpair(min_p_snp$rsid,"rs145407148",pop=pop,token=API_token)},
						error = function(e) { min_p_ind<-p_order[round+1] 
									round=round+1
									return(list("NULL",round,min_p_ind))
								})
						if ( class(r2_out) == "list") {
							r2_res=r2_out[[1]]
							round=r2_out[[2]]
							min_p_ind<-r2_out[[3]]
							min_p_snp<-unlinked_snps[min_p_ind,]
							if ( round >= nrow(unlinked_snps)) {
								loc_snps$assessed=T
								unlinked_snps$assessed=T
								return(list(list(loc_snps,unlinked_snps),c(snp,unlinked_snps$rsid[p_order[1]]),F))
							}
						} else {
							r2_res=NA
						}
					}
					unlinked_other_lead_r2<-sapply(populations, function(pop) {
						r2<-sapply(pos_linked_leads, function(snpy) {
							r2<-tryCatch({ LDpair(min_p_snp$rsid,snpy,pop=pop,token=API_token)},
								error = function(e) {
									return(NA)
									})
							r2<-ifelse(class(r2) == "data.frame",r2$r2,NA)
							return(r2)
						})
						return(r2)
					})
					if ( is.null( nrow(unlinked_other_lead_r2))) unlinked_other_lead_r2<-t(as.data.frame(unlinked_other_lead_r2, stringsAsFactors =F))
					any_pop_link_tf<-sapply(1:nrow(unlinked_other_lead_r2),function(i) any(unlinked_other_lead_r2[i,] > lower_r2_cuttoff, na.rm=T))
					if ( ! any(any_pop_link_tf)) filter_unlinked_snps=F
				}
			}
			if ( filter_unlinked_snps ) { 
				linked_snps$assessed=T
				return(list(list(linked_snps),c(snp),F))
			} else {
				if ( nrow(unlinked_snps) == 1 ) return(list(list(linked_snps,unlinked_snps),c(snp,min_p_snp$rsid),F))
				return(list(list(linked_snps,unlinked_snps),c(snp,min_p_snp$rsid),T))
			} #data
		### next extract ld based on the next snps 
		})
		names(loci_snps_r2)<-lead_snps

		# the object loci_snps_r2 has 3 slots
		#1) the snp table for the old loci and any new loci containing unlinked snps
		#2) the lead snp in each of the above
		#3) T/F whether this loci has been split
		
		#########################
		## now have to go through the loci_snps_r2 object and split into componant parts
		## for those with unlinked variants we need to re extact the proxies and again redefine the loci based on ld 
		## this will have to be changed to a while loop. untill all loci have been separated. 

		variants <-lapply(lead_snps, function(snp) loci_snps_r2[[snp]][[1]])
		names(variants)<-lead_snps
		leads <-sapply(lead_snps, function(snp) loci_snps_r2[[snp]][[2]])
		names(leads)<-lead_snps
		new_loci <- sapply(lead_snps, function(snp) loci_snps_r2[[snp]][[3]])
		names(new_loci)<-lead_snps
		
		if ( all(unlist(new_loci)== F | cycle == cycle_limit)) regions_to_split=F  # exit the loop next time

		#save.image(file=file.path(wk_dir,".RData"))

		if ( any(unlist(new_loci)== T)) {
		# for those loci that have been split now extract the 
			split<-sapply(leads ,function(lead) length(lead) > 1)

			new_leads_lookup<-sapply(leads[new_loci], function(i) i[2]) # lead snps in new split loci that need to have thier proxies looked up 
			new_leads<-sapply(leads[split], function(i) i[2])

			#extract the proxies for the new leads
			new_loci=list()
			old_loci=list()
			for ( i in 1:length(variants)) {
				if ( length(variants[[i]]) == 2) {
					old_loci<- c(old_loci, variants[[i]][1])
					new_loci<- c(new_loci, variants[[i]][2]) 
				} else if (length(variants[[i]]) == 1) {
					old_loci<- c(old_loci, variants[[i]][1])
				} else stop ("split loci has more than two regions")
			}
			names(new_loci) <-new_leads
			names(old_loci) <-lead_snps

			#add a old lead column 
			new_loci<-lapply(1:length(new_loci), function(i) { 
				snp<-names(new_loci)[i]
				new_loci[[snp]]$prev_lead<-new_loci[[snp]]$lead
				new_loci[[snp]]$lead<-snp
				if ( snp %in% new_leads_lookup ) {
					new_loci[[snp]]$proxy_lead<-snp
				}

			return(new_loci[[snp]])
			})
			new_loci2<-rbindlist(new_loci)

			## get the proxies for the split loci 
			new_proxies<-lapply( populations , function(pop) {
				prox_dir=file.path(wk_dir,pop)
				setwd(prox_dir)
				got_data_tf<-sapply(new_leads_lookup, function(snp) file.exists(paste0(snp,".txt")))

				suppressMessages(LDproxy_batch(new_leads_lookup[ ! got_data_tf], 
					pop = pop, 
					r2d = "r2", 
					token = API_token
				))
				got_data_tf<-sapply(new_leads_lookup, function(snp) file.exists(paste0(snp,".txt")))
				new_proxies<-lapply(new_leads_lookup[got_data_tf], function(snp) {
					new_proxy<-read.delim(paste0(snp,".txt"), header=T, sep="\t", stringsAsFactors=F)
					#new_proxy<-data.table(new_proxy)
					new_proxy$lead_snp<-snp
					return(new_proxy)
				})
				names(new_proxies)<-new_leads_lookup[got_data_tf]
				return(new_proxies)
				setwd(wk_dir)
			})
			names(new_proxies)<-populations

			# now find duplications and
			sub_proxies<-find_replace_empty_proxies(new_proxies)
			snps_no_proxies= if( class(sub_proxies[1]) == "list") sub_proxies[[1]] else sub_proxies
			new_proxies= if( class(sub_proxies[1]) == "list") sub_proxies[[2]] else new_proxies

			merged_new_proxies<-merge_proxies(new_proxies,snps_no_proxies)
			merged_new_proxies<-adjust_alleles(merged_new_proxies,ref_genome)
			new_loci2<-add_r2(new_leads_lookup,new_loci2 ,merged_new_proxies,populations)

			#merge the split loci with the original
			old_loci<-lapply(old_loci, function(df) {
				return(df)
			}) 

			new_old_loci<- list(old_loci,new_loci2)
			new_old_loci<-flatten_list(new_old_loci)

			loci_snps_r2<-new_old_loci
			lead_snps<-sapply(1:length(loci_snps_r2), function(i) unique(loci_snps_r2[[i]]$lead))
			names(loci_snps_r2)<-lead_snps
			# merge the new proxies with the original 
			merged_proxies<-rbind(merged_proxies,merged_new_proxies)

		} else {

			split<-sapply(leads ,function(lead) length(lead) > 1) # T/F those loci that have been split
			#new_leads_lookup<-sapply(leads[new_loci], function(i) i[2]) # lead snps in new split loci that need to have thier proxies looked up 
			new_leads<-sapply(leads[split], function(i) i[2]) # lead snps in the new split loci

			#merge the new and old split loci
			new_loci=list()
			old_loci=list()
			for ( i in 1:length(variants)) {
				if ( length(variants[[i]]) == 2) {
					old_loci<- c(old_loci, variants[[i]][1])
					new_loci<- c(new_loci, variants[[i]][2]) 
				} else if (length(variants[[i]]) == 1) {
					old_loci<- c(old_loci, variants[[i]][1])
				} else stop ("split loci has more than two regions")
			}
			names(new_loci) <-new_leads
			names(old_loci) <-lead_snps

			#add a old lead column 
			new_loci<-lapply(1:length(new_loci), function(i) { 
				snp<-names(new_loci)[i]
				new_loci[[snp]]$prev_lead<-new_loci[[snp]]$lead
				new_loci[[snp]]$lead<-snp
				if ( snp %in% new_leads_lookup ) {
					new_loci[[snp]]$proxy_lead<-snp
				}
			return(new_loci[[snp]])
			})

			new_old_loci<- list(old_loci,new_loci)
			new_old_loci<-flatten_list(new_old_loci)
			loci_snps_r2<-new_old_loci
			lead_snps<-sapply(1:length(loci_snps_r2), function(i) unique(loci_snps_r2[[i]]$lead))
			names(loci_snps_r2)<-lead_snps		
		}

		loci_snps_list[[cycle]]<-loci_snps_r2
		lead_snps_list[[cycle]]<-lead_snps
		cat("\n")
	} # END while loop

	loci_SNPs3<-rbindlist(loci_snps_r2)
} else {
	loci_SNPs3<-loci_SNPs2
}
cat("\n")

if ( ! "comments" %in% colnames( loci_SNPs3) ) loci_SNPs3$comments <- NA

# redefine lead_snps and filter on P_value
lead_snps2<-unique(loci_SNPs3$lead)

#################################################################
# filter based on loci p-value
################################################################
leading_snps<-lapply(lead_snps2, function(lead) {
	#print (lead)
	loc_snps<-filter(loci_SNPs3, lead == {{ lead }})
	loc_snps$pre_filt_lead= loc_snps$lead
	loc_snps$lead <- NULL
	if ( nrow(loc_snps) == 1  | all( loc_snps$P_value == 0)) { 
		loc_snps$p_cuttoff <- NA
		loc_snps[,P_lead:=loc_snps$rsid[1]]
		loc_snps$dist_to_P_lead<-NA
		return(list(loc_snps,loc_snps$rsid[1],NULL))
	}
	P_order<-order(loc_snps$P_value,decreasing=F)
	p_val_lead_snp<-loc_snps[P_order[which(loc_snps$P_value[P_order] != 0)[1]],]

	pval_cuttoff<- -log(p_val_lead_snp$P_value,10)*p_val_cuttoff_factor
	pval_cuttoff<- 10^-pval_cuttoff
	loc_snps$p_cuttoff <- pval_cuttoff
	p_filt_ind<-which(! (loc_snps$P_value < pval_cuttoff & loc_snps$P_value < hard_P_val_cutoff))
	loc_snps$comments<-append_to_column(loc_snps$comments,"exceeds p_val filter",p_filt_ind)

	#if the phet of the lead snp is not significant remove snps with a high P-het
	if ("P_heterogeneity" %in% colnames(loc_snps) & filter_phet) {
		if ( p_val_lead_snp$P_heterogeneity > P_het_cuttoff ) {
			p_het_ind<-which(loc_snps$P_heterogeneity < P_het_cuttoff)
			loc_snps$comments<-append_to_column(loc_snps$comments,"P_het_filter",p_het_ind)
		}
	}

	suppressWarnings(loc_snps[,P_lead := p_val_lead_snp$rsid])
	pos=loc_snps$pos[loc_snps$rsid == p_val_lead_snp$rsid]
	loc_snps$dist_to_P_lead<-loc_snps$pos-pos
	#return(list(loc_snps,p_val_lead_snp$rsid,p_val_filt_snps))
	return(list(loc_snps,p_val_lead_snp$rsid))
})

######################################################################
update_snps<-lapply(1:length(leading_snps), function(i) leading_snps[[i]][[1]])
p_val_lead_snp<-sapply(1:length(leading_snps), function(i) leading_snps[[i]][[2]])

retained_loci<-sapply (update_snps, function(slot) nrow(slot) > 0)

p_val_lead_snp<-p_val_lead_snp[retained_loci]
lead_snps2<-lead_snps2[retained_loci]
names(p_val_lead_snp)=lead_snps2
names(lead_snps2)<-p_val_lead_snp

diff_snps<-p_val_lead_snp[ ! p_val_lead_snp %in% lead_snps2]

if ( length(diff_snps) > 0 ) {
	message("the lead snps in these clusters from the summary file are not the lowest p_val variants in the cluster:")
	cat(diff_snps, sep="\n")
	message("reverting to lowest p variant")
	cat("\n")
}
#save.image(file=file.path(wk_dir,".RData"))

loci_SNPs3<-rbindlist(update_snps)
loci_SNPs3$chr<-add_chr(loci_SNPs3$chr)

# columns to extract from the gwas snps
if( add_proxies == T  & (! is.null(unfiltered_file_path)) & split_by_r2 == T) { 
	cols_keep_gwas<-c("rsid","chr","pos","REF","ALT","snp_type","P_value","proxy","P_lead","proxy_lead","orig_lead","pre_filt_lead","p_cuttoff","comments")
	loci_SNPs3$proxy= "False"
	all_proxies<-merged_proxies
	#filter proxies for r2 and MAF cuttoff
	r2_cols<-paste0(populations,"_R2")
	r2<-all_proxies[,r2_cols] # add a column for each population
	r2_tf <- r2 > r2_cuttoff # keep only the variants exceeding the r2 threshold
	keep <- sapply(1:nrow(r2_tf), function(i) any(r2_tf[i,] > r2_cuttoff,na.rm=T))

	maf_cols<-paste0(populations,"_MAF")
	maf<-all_proxies[,maf_cols]
	maf_tf<- maf > MAF_threshold
	keep_maf<-sapply(1:nrow(maf_tf), function(i) any(maf_tf[i,] > MAF_threshold,na.rm=T))

	keep_all<-sapply(1:length(keep_maf), function(i) all( keep[i], keep_maf[i]))
	all_proxies_r2<-all_proxies[keep_all,]

	# filter out the rsid in the merged_proxies that are in the gwas snps
	all_proxies_r2<-all_proxies_r2[ ! all_proxies_r2$RS_Number %in% loci_SNPs3$rsid, ]

	#add some necessary cols
	all_proxies_r2$P_value=NA
	all_proxies_r2$proxy="True"

	#on the fly edit not required if running from scratch 
	#snps_without_proxy_slot<-p_val_lead_snp[ ! p_val_lead_snp %in% names(proxylist[[pop]])]
	#snps_w_proxy_slot<-p_val_lead_snp[ p_val_lead_snp %in% names(proxylist[[pop]])]
	#snps_without_proxy_slot_clumps<-sapply(snps_without_proxy_slot, function(snp) filter(loci_SNPs2, rsid == {{snp}}) %>% dplyr::select(clump))

	if ( updating_proxies == T ) {
		if ( all( snps_without_proxy_slot_clumps %in% changed_lead_snps_clumps )) {
			message ( "updating some of the proxy lists")
			new_proxies<-lapply( populations, function(pop) {
				#proxylist[[pop]]<-proxylist[[pop]][ snps_w_proxy_slot ] 
				my_proxies<-lapply(snps_without_proxy_slot, function(snp) {
					prox<-suppressMessages(LDproxy(snp = snp, 
					pop = pop, 
					r2d = "r2", 
					token = API_token
					))
				prox$lead_snp<-snp
				return(prox)
				})
				names(my_proxies)<-snps_without_proxy_slot
				return(my_proxies)
			})
			names(new_proxies)<-populations

			proxylist<-lapply(populations, function(pop) { 
				pop_proxylist<-append(proxylist[[pop]][ snps_w_proxy_slot ],new_proxies[[pop]])
				#pop_proxylist<-pop_proxylist[ names(pop_proxylist) %in% p_val_lead_snp]
				return(pop_proxylist)
			})
			names(proxylist)<-populations
		}
	}

	# for those loci with with an updated lead snp based on pvalue need to update the correspnding column in the proxies
	#update_diff_snps<-diff_snps[lead_snps2[diff_snps] %in% unique(all_proxies_r2$lead_snp)]

	####################
	if ( ! all(all(p_val_lead_snp %in% unique(loci_SNPs3$P_lead)), all(unique(loci_SNPs3$P_lead[!is.na(loci_SNPs3$P_lead)]) %in% p_val_lead_snp))) stop ( "lead snps dont match")

	#### change the column names and extract the information we need for library design
	cols_keep <-c("RS_Number","chr","pos","REF","ALT","snp_type","P_value","proxy","lead_snp")

	# filter for low R2
	r2_TF<-sapply( 1:nrow(all_proxies_r2) , function(i) {
	min(c(all_proxies_r2$EUR_R2[i],all_proxies_r2$EAS_R2[i]),na.rm=T) < r2_cuttoff
	})

	proxies_merge<-all_proxies_r2[r2_TF,cols_keep]
	colnames(proxies_merge)<-c("rsid","chr","pos","REF","ALT","snp_type","P_value","proxy","lead_snp")

	gwas_snps_merge<-loci_SNPs3[,..cols_keep_gwas]
	
	#add any columns in the gwas snps not present in the proxies
	add_cols<-colnames(gwas_snps_merge)[ ! colnames(gwas_snps_merge) %in% colnames(proxies_merge) ]
	for ( col in add_cols ) proxies_merge[,col]<-NA
	
	if (str_sub(gwas_snps_merge$chr[1],0,3) == "chr")
	gwas_snps_merge$chr<-gsub("chr","",gwas_snps_merge$chr)

	gwas_match<-paste(gwas_snps_merge$chr,gwas_snps_merge$pos,gwas_snps_merge$REF,gwas_snps_merge$ALT,sep="_")
	proxy_match<-paste(proxies_merge$chr,proxies_merge$pos,proxies_merge$REF,proxies_merge$ALT,sep="_")

	#############################################
	# filter the proxies for low p-val snps previously removed
	chr_proxies_filt<-lapply(unique(proxies_merge$chr), function(chr) {
		chr=as.numeric(chr)
		chr_leads<-unlist(unique(filter(gwas_snps_merge, chr == {{chr}}) %>% dplyr::select(proxy_lead)))
		chr<-sub("chr","",chr)
		unfilt_file=sub("CHROM",chr,unfiltered_file_path)
		unfilt_snps<-fread(unfilt_file)
		chr_proxies_filt<-sapply(chr_leads, function(lead) {
			proxies_SNP<-dplyr::filter(proxies_merge, lead_snp == {{lead}})
			if ( nrow(proxies_SNP) == 0) return()
			prox_min<-as.numeric(min(proxies_SNP$pos))
			prox_max<-as.numeric(max(proxies_SNP$pos))
			proxies_loci_match<-paste(proxies_SNP$chr,proxies_SNP$pos,proxies_SNP$REF,proxies_SNP$ALT,sep="_")
			unfilt_loci_snps<-filter(unfilt_snps, pos >= {{prox_min}} & pos <= {{prox_max}})
			unfilt_loci_match1<-paste(unfilt_loci_snps$chr,unfilt_loci_snps$pos,unfilt_loci_snps$EA,unfilt_loci_snps$NEA,sep="_")
			unfilt_loci_match2<-paste(unfilt_loci_snps$chr,unfilt_loci_snps$pos,unfilt_loci_snps$NEA,unfilt_loci_snps$EA,sep="_")
			proxy_gwas_hits<-proxies_loci_match %in% unfilt_loci_match1 | proxies_loci_match %in% unfilt_loci_match2 
			match(proxies_loci_match[! proxy_gwas_hits],proxy_match)			
		})
		return(chr_proxies_filt)
	})

	#combine the set of all proxies not filtered due to p-val
	proxies_remove_ind<-sort(unname(unlist(chr_proxies_filt)))
	proxies_merge<-proxies_merge[proxies_remove_ind,]
	proxy_match<-paste(proxies_merge$chr,proxies_merge$pos,proxies_merge$REF,proxies_merge$ALT,sep="_")

	# check for duplicates
	dups_snps<-match(gwas_match,proxy_match)
	if (any(! is.na(dups_snps))) stop ( "some variants appear to be in both gwas and proxy sets")

	#generate the P_lead col for the proxies 
	proxies_merge$P_lead<-gwas_snps_merge$P_lead[match(proxies_merge$lead_snp,gwas_snps_merge$proxy_lead)]

	colnames(proxies_merge)[colnames(proxies_merge) == "lead_snp"] <- "proxy_lead"

	missing_cols<- colnames(gwas_snps_merge)[ ! colnames(gwas_snps_merge) %in% colnames(proxies_merge)]
	for ( col in missing_cols) proxies_merge[,col]<- NA

	if ( ! all(colnames(proxies_merge) %in% colnames(gwas_snps_merge))) stop("cols dont match")

	#all(colnames(proxies_merge) %in% colnames(gwas_snps_merge))
	all_variants<-rbind(proxies_merge[,cols_keep_gwas],gwas_snps_merge)
} else { 
	cols_keep_gwas<-c("rsid","chr","pos","REF","ALT","snp_type","P_value","P_lead","proxy_lead","orig_lead","pre_filt_lead","p_cuttoff","comments")
	gwas_snps_merge<-loci_SNPs3[,..cols_keep_gwas]
	if (str_sub(gwas_snps_merge$chr[1],0,3) == "chr")
	gwas_snps_merge$chr<-gsub("chr","",gwas_snps_merge$chr)
	all_variants<-gwas_snps_merge
}
all_variants<-all_variants[order(sub("chr","",all_variants$chr),all_variants$pos),]

#remove duplicates save prefiltered file
remove_dups<- ! duplicated(all_variants$rsid)
all_variants<-all_variants[remove_dups,]

message("finished defining loci")
cat("\n")

############################################## finished defining loci ####################################################

##########################################################################################################################
################################################## part 3 ################################################################

# remove any snps that are known sequencing artifacts 
remove_snps=c()
#artefact_filt_snps=c()
message("Part 3")
message("	filtering artefacts")
if ( ! is.null( black_list_snp_file )) {

	remove_snps<-which(all_variants$rsid %in% non_mapping_snps$rsid)
	match_lociSNPs<-paste(all_variants$chr,all_variants$pos,sep="_")
	match_nonmap_snps<-paste(non_mapping_snps$chr,non_mapping_snps$pos,sep="_")
	remove_snps_pos<-which(! is.na(match(match_lociSNPs,match_nonmap_snps)))
	remove_snps<-c(remove_snps,remove_snps_pos)

	if ( length(remove_snps) > 0) { 
		message( "the following snps have have been removed as they are known sequencing artifacts") 
		for ( i in remove_snps)  {
		cat(sep="\n",all_variants$rsid[i])
		}
		all_variants$comments<-append_to_column(all_variants$comments,"SNP is a sequencing artifact",remove_snps)
	}
}

#remove indels if desired
# else remove indels above the threshold
if ( remove_indels==T) { 
	message ("	removing indels")
	indels_ind<-which(all_variants$snp_type != "SNP")
	all_variants$comments <- append_to_column(all_variants$comments,"indel",indels_ind)

} else { 
	long_indels_ind<-nchar(all_variants$REF) > max_indel_len | nchar(all_variants$ALT) > max_indel_len
	if ( any ( long_indels_ind)) {
		cat("flagging the following Indels as they exceed the max_indel_len threshold\n")
		# these are subsequently removed
		cat (all_variants$rsid[long_indels_ind],sep="\n")
		all_variants$comments <- append_to_column(all_variants$comments,"exceeds max indel length",long_indels_ind)
	}
}

########### removed SNP with comments and combine unmappable SNPs into filtered_snps object ############

#harmonise columns
cant_map_filt$p_cuttoff<-NA
cant_map_filt$P_lead<-NA
cant_map_filt$orig_lead<-NA
merge_cols<- colnames(all_variants) [colnames(all_variants) %in% colnames(cant_map_filt)]

# to do just make all filter comments only and filter at the appropiate time
# filtered_snps<-rbindlist(list(p_val_filt_snps[,..merge_cols],cant_map_filt[,merge_cols]))
# now add  the variants from the commented rows 

filtered_snps<-all_variants[! is.na(all_variants$comments),]
filtered_snps<-rbindlist(list(filtered_snps[,..merge_cols],cant_map_filt[,merge_cols]))
all_variants<-all_variants[is.na(all_variants$comments),]

############################### finished removing snps ############################################

###################################### add controls ################################
# add controls +ve not necessary  
# -ve take from random snps? non gwas or repressed region??
# add a few scrambled 

####################################  extract the sequence for each snp and generate the alt allele ############
snp_ind<-which(all_variants$snp_type == "SNP")
#all_snps<-all_variants[snp_ind,]

message("	extracting sequence around variants")
snp_seqs<-lapply(snp_ind , function(i) {
	pos=as.numeric(all_variants$pos[i])
	chrom=paste0("chr",all_variants$chr[i])
	ref_seq<-subseq(ref_genome[chrom],start=pos-downstream,end=pos+upstream)
	ref_seq=unname(as.character(ref_seq))
	alt_seq=ref_seq
	str_sub(alt_seq,downstream+1,downstream+1)<-all_variants$ALT[i]
	return(c("ref"=ref_seq,"alt"=alt_seq))
})
names(snp_seqs)<-all_variants$rsid[snp_ind]

## now deal with the indels

indel_ind<-which(all_variants$snp_type == "InDel")
indel_seq<-lapply( indel_ind, function(i) {
	ref_len<-nchar(all_variants$REF[i])
	pos=as.numeric(all_variants$pos[i])
	start_pos=pos - downstream  
	end_pos= (pos-1 + ref_len) + upstream
	chrom=paste0("chr",all_variants$chr[i])
	ref_seq<-subseq(ref_genome[chrom],start=start_pos,end=end_pos)
	ref_seq=unname(as.character(ref_seq))
	alt_seq=ref_seq
	ref_end=downstream+1+ref_len-1
	if ( ! all_variants$REF[i] == str_sub(alt_seq,downstream+1,ref_end)) stop ( "error unable to recover the reference allele") 
	str_sub(alt_seq,downstream+1,ref_end)<-all_variants$ALT[i]
	return(c("ref"=ref_seq,"alt"=alt_seq))
})
names(indel_seq)<-all_variants$rsid[indel_ind]

all_variants$REF_seq<-NA
all_variants$ALT_seq<-NA

#merge the snps into the all_variants table
all_variants$REF_seq[snp_ind]<- sapply(1:length(snp_seqs), function(i) snp_seqs[[i]]["ref"])
all_variants$ALT_seq[snp_ind]<- sapply(1:length(snp_seqs), function(i) snp_seqs[[i]]["alt"])

if ( length(indel_ind) > 0) {
all_variants$REF_seq[indel_ind]<- sapply(1:length(indel_ind), function(i) indel_seq[[i]]["ref"])
all_variants$ALT_seq[indel_ind]<- sapply(1:length(indel_ind), function(i) indel_seq[[i]]["alt"])
}

# other libarary considerations 
###############################
# remove homopolymers (>8 bases) due to higher synthesis error rate.

message("	identifying homopolymers")
homo_REF<- homopolymerFinder(DNAStringSet(all_variants$REF_seq))
homo_ALT<- homopolymerFinder(DNAStringSet(all_variants$ALT_seq))

has_homo_ref<- lapply(1:length(homo_REF), function(i) {
	return(homo_REF[[i]][width(homo_REF[[i]]) > max_hp_len])
	})
has_homo_ref_ind<-which(sapply(has_homo_ref,function(row) length(row)) > 0)


has_homo_alt<- lapply(1:length(homo_ALT), function(i) {
	return(homo_ALT[[i]][width(homo_ALT[[i]]) > max_hp_len])
	})
has_homo_alt_ind<-which(sapply(has_homo_alt,function(row) length(row)) > 0)

 
if ( mutate_homopolyers == T ) {
	message("	mutating homopolymers") 
	ref_alter<-mutate_homopolyer(has_homo_ref_ind,has_homo_ref,"REF")
	alt_alter<-mutate_homopolyer(has_homo_alt_ind,has_homo_alt,"ALT")
	cat("\n")

	#######################################################################
	# add the new sequences to the master table

	all_variants$REF_mod<-NA
	all_variants$ALT_mod<-NA

	ref_changed_homos<-sapply(ref_alter, function(row) row[[2]])
	alt_changed_homos<-sapply(alt_alter, function(row) row[[2]])

	names(ref_changed_homos)<-has_homo_ref_ind
	names(alt_changed_homos)<-has_homo_alt_ind

	ref_alter<-sapply(ref_alter, function(row) row[[1]])
	alt_alter<-sapply(alt_alter, function(row) row[[1]])

	zero_ind<-which(sapply(ref_alter, function(i) length(i) == 0))
	ref_alter[zero_ind]<-NA
	zero_ind<-which(sapply(alt_alter, function(i) length(i) == 0))
	alt_alter[zero_ind]<-NA

	all_variants$REF_mod[has_homo_ref_ind]<-unlist(ref_alter)
	all_variants$ALT_mod[has_homo_alt_ind]<-unlist(alt_alter)

	# for those with both homopolymer in the ref and alt  
	###########################################################
	## now reconsile ref/alt 

	comb_ind<-c(has_homo_ref_ind,has_homo_alt_ind)
	count_ind<-table(comb_ind) 

	# go though the homoployers and separate out those in only one allele
	remove_revised_seqs<-c()
	ref_only=c()
	alt_only=c()
	for ( i in names(count_ind)) {
		only_ref<-c()
		only_alt<-c()
		ref_hom_change<-ref_changed_homos[[i]]
		alt_hom_change<-alt_changed_homos[[i]]
		if ( all( ! ref_hom_change %in% alt_hom_change )) {
			#message(i)
			only_ref<-ref_hom_change[ ! ref_hom_change %in% alt_hom_change ]
			only_alt<-alt_hom_change[ ! alt_hom_change %in% ref_hom_change ]
			if ( length(only_ref) > 0 & length(only_alt) > 0 ) {
				remove_revised_seqs<-c(remove_revised_seqs,i)			
			} else {
				if ( length(only_ref) > 0 ) {
				ref_only<-c(ref_only,as.numeric(i))
				} else if ( length(only_alt) > 0 ) {
				alt_only<-c(alt_only,as.numeric(i))
				}
			}	
		}
		if ( length(remove_revised_seqs) > 0) {
			remove_inds<-as_numeric(remove_revised_seqs)
			cat("unable to reconsile differences between ref and alt allele revised sequences for\n")
			cat(all_variants$rsid[remove_inds],sep="\n")
			message("removing the revised sequences, these snps will be dropped")
			all_variants$REF_mod[remove_inds]<-NA
			all_variants$ALT_mod[remove_inds]<-NA
		}
	}

	# this will for those snps where the homopolyer was only on one allele make the same change on the other allele and check it doesnt 
	# violate any constratins

	alt_mod<-sapply(ref_only, function(ind) {
		if ( is.na(all_variants$REF_mod[ind])) { 
			message( "SNP " , ind , " failed to modify ref")
			return(NA)
		}
		rsID=all_variants$rsid[ind]
		ref_len<- nchar(all_variants$REF[ind])	
		alt_len<- nchar(all_variants$ALT[ind])
		orig<-unlist(str_split(all_variants$REF_seq[ind],""))
		new<-unlist(str_split(all_variants$REF_mod[ind],""))
		seq<-all_variants$ALT_seq[ind]
		diff_base<-which(orig != new)
		snp_dist=diff_base - (ref_len + downstream)
		if ( all(snp_dist >= min_dist_mut2snp)) {
			#able to map diff base to other allele
			base_on_other=downstream+alt_len+snp_dist
			for (base in base_on_other) {
				str_sub(seq,base,base)<-swap_base_list[[str_sub(seq,base,base)]]
			}
		return(seq)
		} else {
		message("SNP ", rsID ," unable to generate same homoploymer mutations on the boths alleles, within constratins, SNP will be filtered")
		return(NA)
		}
	})	

	ref_mod<-sapply(alt_only, function(ind) {
		if ( is.na(all_variants$ALT_mod[ind])) { 
			message( "SNP" , ind , "failed to modify alt") 
			return(NA)
		}
		rsID=all_variants$rsid[ind]
		ref_len<- nchar(all_variants$REF[ind])	
		alt_len<- nchar(all_variants$ALT[ind])
		orig<-unlist(str_split(all_variants$ALT_seq[ind],""))
		new<-unlist(str_split(all_variants$ALT_mod[ind],""))
		seq<-all_variants$REF_seq[ind]
		diff_base<-which(orig != new)
		snp_dist=diff_base - (alt_len + downstream)
		if ( all(snp_dist >= min_dist_mut2snp)) {
			#able to map diff base to other allele
			base_on_other=downstream +ref_len+snp_dist
			for (base in base_on_other) {
				str_sub(seq,base,base)<-swap_base_list[[str_sub(seq,base,base)]]
			}
		return(seq)
		} else {
		message("SNP ", rsID ," unable to generate same homoploymer mutations on the boths alleles, within constratins, SNP will be filtered")
		return(NA)
		}
	})	
	
	cat("\n")
	all_variants$ALT_mod[ref_only]<-alt_mod
	all_variants$REF_mod[alt_only]<-ref_mod

	## get the index of all that were changed 
	any_homo_ind<-unique(c(has_homo_ref_ind,has_homo_alt_ind))
	was_changed<- which((! is.na(all_variants$REF_mod)) | ( ! is.na(all_variants$ALT_mod)))
	any_changed_ind<-any_homo_ind[any_homo_ind %in% was_changed]

	diff_change<-sapply(any_changed_ind, function(i) {
		ref_mod<-all_variants$REF_mod[i]
		alt_mod<-all_variants$ALT_mod[i]
		ref<-all_variants$REF[i]
		alt<-all_variants$ALT[i]
		ref_no_SNP<-paste0(str_sub(ref_mod,1,downstream),str_sub(ref_mod,downstream +1+nchar(ref),nchar(ref_mod)))
		alt_no_SNP<-paste0(str_sub(alt_mod,1,downstream),str_sub(alt_mod,downstream +1+nchar(alt),nchar(alt_mod)))
		ref_no_SNP==alt_no_SNP
	})

	diff_change_ind<-any_changed_ind[! diff_change]

	#remove those where we failed to map the ref to alt or vice versa
	remove_inds<-is.na(all_variants$REF_mod) & (! is.na(all_variants$ALT_mod)) | ( ! is.na(all_variants$REF_mod)) & is.na(all_variants$ALT_mod)
	remove_inds<-which(remove_inds)
	diff_change_ind<-diff_change_ind[!diff_change_ind %in% remove_inds] 

	revise_change<-lapply(diff_change_ind, function(ind) {
		#print(ind)
		rsID=all_variants$rsid[ind]
		ref_mod<-all_variants$REF_mod[ind]
		alt_mod<-all_variants$ALT_mod[ind]
		ref<-all_variants$REF[ind]
		alt<-all_variants$ALT[ind]
		ref_len<-nchar(ref)
		alt_len<-nchar(alt)

		#move ref_2_alt 
		alt_revise<- unlist(ref_mod)
		str_sub(alt_revise,downstream +1,downstream +nchar(ref)) <-alt

		new<-unlist(str_split(alt_revise,""))
		orig<-unlist(str_split(all_variants$ALT_seq[ind],""))
		if ( !length( new) == length(orig) ) stop ( "the length of the doner and recpiant sequence must be the same")

		diff_base<-which(orig != new)
		hp_TF<-any(width(unlist(homopolymerFinder(DNAStringSet(alt_revise)))) > 8)
		snp_dist=sapply(diff_base, function(pos) {
			min(c(abs(pos - (alt_len + downstream)),abs(downstream +1 -pos)))
		})
		if ( all(snp_dist >= min_dist_mut2snp) & ( ! hp_TF)) {
			return(c("ALT",alt_revise))
		}

		#move alt_2_ref
		ref_revise<- unlist(alt_mod)
		str_sub(ref_revise,downstream +1,downstream +nchar(alt)) <-ref

		new<-unlist(str_split(ref_revise,""))
		orig<-unlist(str_split(all_variants$REF_seq[ind],""))
		if ( ! length( new) == length(orig) ) stop ( "the length of the doner and recpiant sequence must be the same")

		diff_base<-which(orig != new)
		hp_TF<-any(width(unlist(homopolymerFinder(DNAStringSet(ref_revise)))) > 8)
		snp_dist=sapply(diff_base, function(pos) {
			min(c(abs(pos - (alt_len + downstream)),abs(downstream +1-pos)))
		})
		if ( all(snp_dist >= min_dist_mut2snp) & ( ! hp_TF)) {
			return(c("REF",ref_revise))
		} else {
			return(c(NA,NA))
			message("SNP ", rsID ,"unable to generate same homoploymer mutations on the boths alleles, within constratins, SNP will be filtered")
		}
	})	

	## add these into all_variants 
	col<-sapply(revise_change,function(slot) slot[1])
	seq<-sapply(revise_change,function(slot) slot[2])

	keep=! is.na(col)
	keep_change_ind<-diff_change_ind[keep]

	for ( i in 1:length(keep_change_ind)) {
		ind=keep_change_ind[i]
		mod_col=paste0(revise_change[[i]][1],"_mod")
		new_seq<-revise_change[[i]][2]
		all_variants[ind,mod_col] <- new_seq
	}

	remove_change_ind<-unique(c(diff_change_ind[! keep],remove_inds))

} else { # just remove any homopolymers 

	homopolymer_filter_ind<-unique(c(has_homo_ref_ind,has_homo_alt_ind))

	if ( length(homopolymer_filter_ind) > 0) { 
		message( "the following snps have have been removed as they contain homopolymers") 
		cat(sep="\n",all_variants$rsid[homopolymer_filter_ind])
		cat("\n")
		all_variants$comments<-append_to_column(all_variants$comments,"contains homoployer",homopolymer_filter_ind)
		homopolymer_filter<-all_variants[homopolymer_filter_ind,]
	}
	all_variants$REF_mod<-NA
	all_variants$ALT_mod<-NA
}

############################################## done ###############################################

#############################################combine the sequences passing filtering #####################################################
#find any remaining homoployers and filter 
mod_ref=! is.na(all_variants$REF_mod) #get those that have been modified
mod_alt=! is.na(all_variants$ALT_mod) 

if (any(mod_ref)) {
	homo_REF<- homopolymerFinder(DNAStringSet(unlist(all_variants$REF_mod[mod_ref]))) #re extract homos
	has_homo_ref<- lapply(1:length(homo_REF), function(i) {
		return(homo_REF[[i]][width(homo_REF[[i]]) > max_hp_len])
	})
	has_homo_ref_ind<-which(sapply(has_homo_ref,function(row) length(row)) > 0)
}

if (any(mod_ref)) {
	homo_ALT<- homopolymerFinder(DNAStringSet(unlist(all_variants$ALT_mod[mod_alt])))
	has_homo_alt<- lapply(1:length(homo_ALT), function(i) {
		return(homo_ALT[[i]][width(homo_ALT[[i]]) > max_hp_len])
	})
	has_homo_alt_ind<-which(sapply(has_homo_alt,function(row) length(row)) > 0)
}

still_homo<-unique(c(which(mod_ref)[has_homo_ref_ind],which(mod_alt)[has_homo_alt_ind]))

################
# merge the modified and unmodified oligos

if ( mutate_homopolyers == T ) {
	unchanged_ind <- which(is.na(all_variants$REF_mod) & is.na(all_variants$ALT_mod)) # get the unchanged oligos
	cant_use_homo <-c (unchanged_ind[unchanged_ind %in% names(count_ind)],still_homo,remove_change_ind)
	cant_use_homo <- unique(cant_use_homo)
	cant_use_homo <- cant_use_homo[! is.na(cant_use_homo)]

	all_variants$comments<-append_to_column(all_variants$comments,"contains homoployer",cant_use_homo)
	homopolymer_filter<-all_variants[cant_use_homo,]

	use_orig<- unchanged_ind [ ! unchanged_ind %in% names(count_ind)] # filter for those with a homo in either ref or alt
	use_orig <- use_orig[ ! use_orig %in% cant_use_homo]

	use_mod<- which(mod_ref & mod_alt )
	use_mod<-use_mod[ ! use_mod %in% cant_use_homo] # filter those oligos that were not sucessfully filtered 

	all_vars<-c(use_orig,use_mod,cant_use_homo) # index of all_variants

	all_variants$use_modified<-NA
	all_variants$use_modified[use_mod]<- "yes"

	if(any(duplicated(all_vars))) stop ( "unable to determine whether snp modification was successful for all alleles")
	if ( ! all(1:nrow(all_variants) %in% all_vars)) stop ( "some varinats are unaccounted for")

	if ( any(duplicated(all_variants$rsid)) ) stop ( "there are duplicated rsid in the variants file")
} else {
	all_variants$use_modified<-NA
	use_orig<-which(is.na(all_variants$comments))
}

cols<-colnames(filtered_snps) # because data.table is shit 
filtered_snps<-rbindlist(list(filtered_snps,homopolymer_filter[,..cols]))

#######################
# find and mutate any occurances of the restriction site SceI

#make a dna set and find any occurnaces for the restriction site

ref_seqs<-c(unlist(all_variants$REF_seq[use_orig]),unlist(all_variants$REF_mod[use_mod])) # get the sequences either modified or not
alt_seqs<-c(unlist(all_variants$ALT_seq[use_orig]),unlist(all_variants$ALT_mod[use_mod]))
names(ref_seqs)<-c(use_orig,use_mod) # name the seqs with the original index
names(alt_seqs)<-c(use_orig,use_mod)

#all(rle(as.numeric(names(alt_seqs)[order(as.numeric(names(alt_seqs)),decreasing =F)]))$length == 1)
ref_seqs<-ref_seqs[order(as.numeric(names(ref_seqs)),decreasing =F)]
alt_seqs<-alt_seqs[order(as.numeric(names(alt_seqs)),decreasing =F)]

ref_seqs_DNA<- DNAStringSet(ref_seqs)
alt_seqs_DNA<- DNAStringSet(alt_seqs)

SceI="TAGGGATAACAGGGTAAT"
SceI<- DNAString(SceI)

SceI_ref_matches<-vmatchPattern(SceI,ref_seqs_DNA)
SceI_alt_matches<-vmatchPattern(SceI,alt_seqs_DNA)

ref_sceI_TF<-sapply(SceI_ref_matches, function(row) length(row)) > 0
if ( any(ref_sceI_TF) ) {
	stop ( "some ref oloigos have an SceI site in them, mine didnt so I didnt code for this")
}

alt_sceI_TF<-sapply(SceI_alt_matches, function(row) length(row)) > 0
if ( any(alt_sceI_TF) ) {
	stop ( "some alt oloigos have an SceI site in them, mine didnt so I didnt code for this")
}

#############################
#
# make CONTROLS
#
#######################
# now take the top repressed regions based on H3K27me3 marks in colon extract SNPs and use these as negative controls

message("	generating controls")

if ( ! (is.null(negative_control_region_file) & prop_control > 0)) {

	message( "using user provided regions to design controls")
	neg_cont_region<-fread(negative_control_region_file)
	colnames(neg_cont_region) <- c("chr","start","stop","rank","score","strand","enrichment","p-val","q-val")

	#get the enriched regions from the repressed bed
	repressed<-filter(neg_cont_region,enrichment > ChIP_enrichment)
	repressed$chr<-gsub("chr","",repressed$chr)
	repressed.gr<-makeGRangesFromDataFrame(repressed, keep.extra.columns=T)
	repressed_SNPs.gr<-subsetByOverlaps(dbSNPs.gr,repressed.gr)

	#get SNPs in respresed regions 
	#provide column names
	colnam<-c("rsid","chr","pos","REF","ALT","snp_type","P_value","proxy","lead","REF_seq","ALT_seq","REF_mod","ALT_mod")
	#create data frame with 0 rows and 3 columns
	repressed_snp_df<- data.frame(matrix(ncol = length(colnam), nrow = 0))
	colnames(repressed_snp_df) <- colnam
	
	if ( length(repressed.gr) < n_controls * prop_control) { 
		warning ("Using the provided ChIP_enrichment yields ", length(repressed.gr) ," of " , n_controls * prop_control , " desired control regions\nThe remaining controls will be generated from scrambled sequence, plesae use a different control region file or lower the ChIP_enrichment threshold.")
	}

	done=c()
	while ( nrow(repressed_snp_df) < n_controls * prop_control ) {
		pick<-sample(c(1:length(repressed.gr))[! 1:length(repressed.gr) %in% done],1)
		regions_search.gr<-repressed.gr[pick]
		region_SNPs.gr<-subsetByOverlaps(repressed_SNPs.gr,regions_search.gr)
		if (length(region_SNPs.gr) == 0 ) next 
		seq_OK=F
		list_picks<-sample(1:length(region_SNPs.gr),length(region_SNPs.gr))

		j=0
		while( seq_OK == F & j <= length(list_picks)) {
			j=j+1
			pick<-list_picks[j]
			if (is.na(pick)) stop ("pick is NA")
			snp<-region_SNPs.gr[pick]
			nalts<-length(unlist(mcols(snp)$ALT))
			if ( nalts > 1 ) next 
			chrom<-paste0("chr",seqnames(snp))
			start_pos<-start(snp) - downstream
			ref_len<-nchar(mcols(snp)$REF)
			alt_len<-nchar(as.character(unlist(mcols(snp)$ALT)))
			end_pos<-start(snp)-1 + ref_len + upstream
			ref_seq<-subseq(ref_genome[chrom],start=start_pos,end=end_pos)
			end_pos<-start(snp)-1 + alt_len + upstream
			alt_seq<-subseq(ref_genome[chrom],start=start_pos,end=end_pos)
			str_sub(alt_seq,downstream+1,downstream+ref_len)<-unlist(mcols(snp)$ALT)
			snp_alleles<-DNAStringSet(c(ref_seq,alt_seq))
			homos<-homopolymerFinder(snp_alleles)
			snps_homo_TF<-sapply(homos, function(homo) {
				any(width(homo) > max_hp_len)
			})
			if (any (snps_homo_TF) ) next 
			SceI_match<-vmatchPattern(SceI,snp_alleles)
			any_SceI_TF<-sapply(SceI_match, function(sceI_hit) {
				length(sceI_hit) > 0
			})
			if (any (any_SceI_TF )) next 
			seq_OK=T
		}
		if ( seq_OK==T ) {
			snp_type<-ifelse ( ref_len > 1 | alt_len > 1, "InDel", "SNP" )
			if ( remove_indels == T & snp_type == "InDel" ) next
			repressed_snp_df<-rbind(repressed_snp_df,data.frame(rsid=names(snp),chr=chrom,pos=start(snp),REF=mcols(snp)$REF,ALT=unlist(mcols(snp)$ALT),snp_type=snp_type,P_value=NA,test_control="repressed_control",lead=NA,REF_seq=ref_seq,ALT_seq=alt_seq,REF_mod=NA,ALT_mod=NA))
			done=c(done,pick)
		}
	}
	if ( n_controls - nrow(repressed_snp_df) >  n_controls * prop_control ) message( "unable to extract ", n_controls/2, "control sequences from the supplied region file. Consider lowering the enrichment threshold currently hard coded in the line repressed<-filter(neg_cont_region,enrichment > 4)")
	add_controls = n_controls - nrow(repressed_snp_df) 
} else {
	add_controls=n_controls
}
# create n scrabled sequences
scrambled_seqs=c()
scrambled_alt_seqs=c()
scrambled_ref=c()
scrambled_alt=c()
bases=c("A","T","C","G")

while ( length(scrambled_seqs) < add_controls ) {
	scram_seq<-sample(bases,size=insert_size+1,replace = T)
	scram_seq<-paste0(scram_seq,collapse="")
	scram_seq_DNA<-DNAStringSet(scram_seq)
	homo<-homopolymerFinder(scram_seq_DNA)
	if (any(width(homo) > max_hp_len)) next 
	SceI_match<-vmatchPattern(SceI,scram_seq_DNA)
	if (length(SceI_match[[1]]) > 0 ) next 
	if ( scram_seq %in% scrambled_seqs ) next 
	ref=str_sub(scram_seq,downstream + 1,downstream + 1)
	swap<-bases [ ! bases %in% str_sub(scram_seq,downstream + 1,downstream + 1)]
	alt_seq<-scram_seq
	str_sub(alt_seq,downstream + 1,downstream + 1)<- sample(swap,1)
	alt=str_sub(alt_seq,downstream + 1,downstream + 1)
	alt_seq_DNA<-DNAStringSet(alt_seq)
	homo<-homopolymerFinder(alt_seq_DNA)
	if (any(width(homo) > max_hp_len)) next 
	SceI_match<-vmatchPattern(SceI,alt_seq_DNA)
	if (length(SceI_match[[1]]) > 0 ) next 
	if ( alt_seq %in% scrambled_seqs ) next 
	scrambled_seqs <- c(scrambled_seqs,scram_seq)
	scrambled_alt_seqs <- c(scrambled_alt_seqs,alt_seq)
	scrambled_ref<-c(scrambled_ref,ref)
	scrambled_alt=c(scrambled_alt,alt)
}

if ( add_controls > 0) {
	scrambled_controls<-data.frame(rsid=paste0("scram",1:(add_controls)),chr=NA,pos=NA,REF=scrambled_ref,ALT=scrambled_alt,snp_type="SNP",P_value=NA,test_control="scrambled_control",lead=NA,REF_seq=scrambled_seqs,ALT_seq=scrambled_alt_seqs,REF_mod=NA,ALT_mod=NA)
	if (exists("repressed_snp_df")) {
		controls<-rbind(repressed_snp_df,scrambled_controls)
	} else {
		controls<-scrambled_controls
	}
} else {
	controls<-repressed_snp_df
}
# add the differnt SNPs together (test + control)
# test + control table 

keep_rows<-sort(c(use_orig,use_mod))
keep_cols<-c("rsid","snp_type","REF","ALT","use_modified","test_control")

all_variants$test_control<- "test"

all_variants<-all_variants[keep_rows,] # filter out any SNPs with homopolymers

test_snps<-all_variants[,..keep_cols]
#test_snps$test_control="test"
test_snps$REF_seq=ref_seqs
test_snps$ALT_seq=alt_seqs

control_snps<-data.frame(rsid=controls$rsid, snp_type=controls$snp_type, REF=controls$REF, ALT=controls$ALT ,use_modified= NA, test_control="control", REF_seq= controls$REF_seq,ALT_seq=controls$ALT_seq)
test_control_SNPs<-rbind(test_snps,control_snps)


message( "	finished generating controls")
#test_snps
#####################################
# identify those snps in coding or promoter regions ?  Not implimented

# add the adapters in both directions 
CRS_adpater_F_DNA<-DNAString(CRS_adpater_F)
CRS_adpater_R_DNA<-DNAString(CRS_adpater_R)
CRS_adpater_F_rev_comp<-as.character(reverseComplement(CRS_adpater_F_DNA))
CRS_adpater_R_rev_comp<-as.character(reverseComplement(CRS_adpater_R_DNA))

library_adapters<-lapply(1:nrow(test_control_SNPs) , function(i) {
	ref_F<-paste0(CRS_adpater_F,test_control_SNPs$REF_seq[i],CRS_adpater_R)
	alt_F<-paste0(CRS_adpater_F,test_control_SNPs$ALT_seq[i],CRS_adpater_R)

	ref_R<-paste0(CRS_adpater_R_rev_comp,test_control_SNPs$REF_seq[i],CRS_adpater_F_rev_comp)
	alt_R<-paste0(CRS_adpater_R_rev_comp,test_control_SNPs$ALT_seq[i],CRS_adpater_F_rev_comp)
	return(list(REF_fwd=ref_F,ALT_fwd=alt_F,REF_rev=ref_R,ALT_rev=alt_R))
})

names(library_adapters)<-test_control_SNPs$rsid

# carry out some final checks
library_DNA<-DNAStringSet(unlist(library_adapters))
homo<-homopolymerFinder(library_DNA)

homo_TF<-sapply(1:length(homo), function(i) {
	any(width(homo[[i]]) > max_hp_len)
}) 

#mutate the homopolyers final if there are any homopolyers this dont overlap with the adapters this is exit
has_homo<-which(homo_TF)

#due to unlisting the library this converts the row number back to a list row and slot format
convert_2_list_ind<-t(as.data.frame(sapply(has_homo, function(x)  { 
	resid=(x-1)/4
	row=floor(resid)
	slot=(resid-row)/.25
	return(c(row,slot)+1)
})))

colnames(convert_2_list_ind) <- c("library_row","library_slot") 

########################################
library_alter<-lapply(has_homo, function(ind) {

	homos<-homo[[ind]]
	homos<-homos[width(homos) > max_hp_len]
	keep=T
	library_seq<-as.character(library_DNA[ind])
	end<-nchar(library_seq) - nchar(CRS_adpater_R)
	genome_seq<-str_sub(library_seq,nchar(CRS_adpater_R)+1,end)
	oligo_name<-names(library_DNA)[ind]
	oligo_name<-str_split(oligo_name,"(\\.|_)")

	rsid<-unlist(oligo_name)[1]
	ref_alt<-unlist(oligo_name)[2]
	fwd_rev<-unlist(oligo_name)[3]
	test_cont_ind<-which(test_control_SNPs$rsid == rsid)
	all_var_ind<-which(all_variants$rsid == rsid)
	if ( length(all_var_ind) == 1) {
		orig_var<-all_variants[all_var_ind,]
		col_suff<-ifelse (is.na(orig_var$REF_mod),"seq","mod")
		coly<-paste(ref_alt,col_suff,sep="_")
		orig_seq<-unlist(orig_var[,..coly])
		if ( nchar(orig_seq) != nchar(genome_seq) ) stop ("error processing homopolymers in adapter ligated oligoes for SNP ", oligo_name)
		if ( orig_seq != genome_seq )  warning ("error processing homopolymers in adapter ligated oligoes for SNP ", oligo_name ," seqences between orginal sequence and adapter ligated sequence differ")
		} else all_var_ind<-NA
	for ( i in 1:length(homos)) {
		polybase<-mcols(homos)$base
		low_ol_TF<-start(homos)[i] < 16 & end(homos)[i] > 16
		high_ol_TF<-start(homos)[i] < end & end(homos)[i] > end
		if ( ! (low_ol_TF | high_ol_TF)) {
			keep=F
			return(list(library_seq,test_cont_ind,all_var_ind,ref_alt,keep,rsid))
		}
		#extract the homoplymer that overlaps with the genomic seq 
		ifelse(low_ol_TF, hp_pos<-unname(c(16,end(homos)[i])),hp_pos<-unname(c(start(homos)[i],end)))

		hom_width<-hp_pos[2]-hp_pos[1] # homopolymer width
		if (hom_width > 7) {
			keep=F
			return(list(library_seq,test_cont_ind,all_var_ind,ref_alt,keep,rsid))
		}
		ifelse(high_ol_TF, shift<-ceiling(hom_width/2),shift<-floor(hom_width/2))
		base_to_change=hp_pos[1]+shift
		str_sub(library_seq,base_to_change,base_to_change)<-swap_base_list[[polybase]]
	}
	return(list(library_seq,test_cont_ind,all_var_ind,ref_alt,keep,rsid))
})

update_seq<-sapply(library_alter,function(row) row[[1]])
test_cont_ind<-sapply(library_alter,function(row)	row[[2]])
all_var_ind<-sapply(library_alter,function(row)	row[[3]])
ref_alt_mod<-sapply(library_alter,function(row) row[[4]])
keep=sapply(library_alter,function(row) row[[5]])
rsid=sapply(library_alter,function(row) row[[6]])

test_TF=! is.na(all_var_ind)
test_cont=ifelse(test_TF, "test","cont")

#due to unlisting the library this converts the row number back to a list row and slot format
convert_2_list_ind<-t(as.data.frame(sapply(has_homo, function(x)  { 
	resid=(x-1)/4
	row=floor(resid)
	slot=(resid-row)/.25
	return(c(row,slot)+1)
})))

final_update<-data.frame(rsid=rsid,all_var_index=all_var_ind,test_cont_ind=test_cont_ind,test_cont=test_cont,allele=ref_alt_mod,modified_seq=update_seq, library_row=convert_2_list_ind[,1] , library_slot=convert_2_list_ind[,2],keep=keep )

message("	performing final checks")
# 
for ( snp in unique(final_update$rsid)) { 
	#print(snp)
	changed_seqs<-filter(final_update, rsid == {{snp}} )
	all_var_ind=unique(changed_seqs$all_var_index)
	test_cont_ind<-unique(changed_seqs$test_cont_ind)
	lib_rsid=names(library_adapters[unique(changed_seqs$library_row)])
	test_cont_rsid=test_control_SNPs$rsid[test_cont_ind]
	if  ( ! lib_rsid == snp & snp == test_cont_rsid ) stop ( "all_variants vs oligo library rsid mismatch")
	oligo_names<-names(library_adapters[[unique(changed_seqs$library_row)]][changed_seqs$library_slot])
	fwd_rev<-sapply(oligo_names,function(name) unlist(str_split(name,"_"))[2])
	row=unique(changed_seqs$library_row)

	if( ! nrow(changed_seqs) %% 2 == 0 ) stop ( "Edited and odd number of oliogs for ", snp, " must change both ref and alt or ") 

	if (all(changed_seqs$keep)) {
		if ( ! (sum ( fwd_rev == "fwd") == 2 | sum ( fwd_rev == "rev") == 2 )) stop ("final update error")	
		for ( i in 1:nrow(changed_seqs)) {
			slot=changed_seqs$library_slot[i]
			new_seq<-changed_seqs$modified_seq[i]
			ref_alt=changed_seqs$allele[i]
			col=paste0(ref_alt,"_mod")
			unalt_col=paste0(ref_alt,"_seq")

			old_seq<-ifelse(is.na(all_var_ind),test_control_SNPs[test_cont_ind,..unalt_col],all_variants[all_var_ind,..unalt_col])
			short_seq<-str_sub(new_seq,16,nchar(new_seq)-15)
			if ( ! nchar(old_seq) == nchar(short_seq)) stop ( snp, " original and replacment sequence have different length")
			test_control_SNPs[test_cont_ind,(unalt_col) := short_seq ] 
			test_control_SNPs$use_modified[test_cont_ind]<- "yes"
			library_adapters[[row]][slot]<-new_seq
			if ( ! is.na(all_var_ind))	all_variants[all_var_ind, col] <- short_seq 
		}
	} else {
		if ( ! (sum ( changed_seqs$allele == "REF") == 2 | sum ( changed_seqs$allele == "ALT") == 2 )) stop ("final update error")
	}
}

## deal with sequences that can being kept 
update_all_var_ind<-unique(final_update$all_var_ind[final_update$keep == T])
update_all_var_ind<-update_all_var_ind[ ! is.na(update_all_var_ind)]
if ( length(update_all_var_ind) > 0) { 
	all_variants$comments<-append_to_column(all_variants$comments,"Sequence modified near adapter due to creation of homopolyer",update_all_var_ind)
}

## deal with the sequences that need to be removed in the all_variants table
# update comments and filter in the all variants table and remove sequences from the library table
update_all_var_ind<-unique(final_update$all_var_ind[final_update$keep == F])
update_all_var_ind<-update_all_var_ind[ ! is.na(update_all_var_ind)]
update_test_cont_ind<-unique(final_update$test_cont_ind[final_update$keep == F])
update_library_ind<-unique(final_update$library_row [final_update$keep == F])

if ( ! (all(all_variants$rsid[update_all_var_ind] %in% names(library_adapters[update_library_ind])) & all(all_variants$rsid[update_all_var_ind] %in% test_control_SNPs$rsid[update_test_cont_ind]))) stop ( "all_variants vs oligo library mismatch")

if ( length(update_library_ind) > 0) { 
	#all_variants$filtered[update_all_var_ind] <- "yes" 
	all_variants$comments<-append_to_column(all_variants$comments,"homopolymer, cant be mutated within constraints",update_all_var_ind)
	library_adapters<-library_adapters[-update_library_ind]
	test_control_SNPs<-test_control_SNPs[-update_test_cont_ind,]
	filtered_snps<-rbind(filtered_snps,all_variants[update_all_var_ind,colnames(all_variants)[ colnames(all_variants) %in% colnames(filtered_snps) ]])
	all_variants=all_variants[-update_all_var_ind,]
}

# preform final check on library 
# make sure ref and alt are the same outside snp 
# check that those unaltered can be blasted back to the correct location 
message( "checking CRS sequence consistency")
snp_OK<-lapply(1:length(library_adapters), function(i) { 
	snp_OK=T
	snp_seq=NA
	rsid<-names(library_adapters)[i]
	#skip controls
	if ( ! rsid %in% controls$rsid ) {
		pc<-round(i/length(library_adapters),2)*100
		if ( pc %% 1 == 0 ) {
			cat("working - ", pc ," % complete \r")
			flush.console()
		}
		#print(rsid)
		row=filter(all_variants, rsid == {{rsid}})
		all_var_ind<-which(all_variants$rsid == rsid)
		ref_allele=unlist(row %>% dplyr::select(REF))
		alt_allele=unlist(row %>% dplyr::select(ALT))

		#next bit will extract the sequence around snps in each ref/alt fwd/rev format by leaving out the bit at the edge which may have been revised by in only the fwd or rev orientation. These may otherwise cause an error as the adapter ligation may introduce a homopolymer 
		ref_fwd_s<-str_sub(library_adapters[[i]]$REF_fwd,nchar(CRS_adpater_F)+1+max_hp_len,nchar(library_adapters[[i]]$REF_fwd)-(nchar(CRS_adpater_R)+max_hp_len))
		alt_fwd_s<-str_sub(library_adapters[[i]]$ALT_fwd,nchar(CRS_adpater_F)+1+max_hp_len,nchar(library_adapters[[i]]$ALT_fwd)-(nchar(CRS_adpater_R)+max_hp_len))
		ref_rev_s<-str_sub(library_adapters[[i]]$REF_rev,nchar(CRS_adpater_R)+1+max_hp_len,nchar(library_adapters[[i]]$REF_rev)-(nchar(CRS_adpater_F)+max_hp_len))
		alt_rev_s<-str_sub(library_adapters[[i]]$ALT_rev,nchar(CRS_adpater_R)+1+max_hp_len,nchar(library_adapters[[i]]$ALT_rev)-(nchar(CRS_adpater_F)+max_hp_len))

		if ( ! ref_fwd_s == ref_rev_s ) snp_OK=F #stop ( "ref sequences for ", rsid, " dont match")  
		if ( ! alt_fwd_s == alt_rev_s )snp_OK=F # stop ( "alt sequences for ", rsid, " dont match") 

		ref_fwd<-str_sub(library_adapters[[i]]$REF_fwd,nchar(CRS_adpater_F)+1,nchar(library_adapters[[i]]$REF_fwd)-nchar(CRS_adpater_R))
		alt_fwd<-str_sub(library_adapters[[i]]$ALT_fwd,nchar(CRS_adpater_F)+1,nchar(library_adapters[[i]]$ALT_fwd)-nchar(CRS_adpater_R))

		if ( ! nchar(ref_fwd) - nchar(ref_allele) == insert_size ) snp_OK=F # stop ( "length ref seq is not correct for ", rsid)
		if ( ! nchar(alt_fwd) - nchar(alt_allele) == insert_size ) snp_OK=F # stop ( "length alt seq is not correct for ", rsid)

		ref<-str_sub(ref_fwd,downstream+1,downstream + nchar(ref_allele))
		alt<-str_sub(alt_fwd,downstream+1,downstream + nchar(alt_allele))
		if ( ! (ref == ref_allele & alt == alt_allele)) snp_OK=F

		ref_no_SNP<-paste0(str_sub(ref_fwd,1,downstream),str_sub(ref_fwd,downstream+ 1 +nchar(ref),nchar(ref_fwd)))
		alt_no_SNP<-paste0(str_sub(alt_fwd,1,downstream),str_sub(alt_fwd,downstream+ 1 +nchar(alt),nchar(alt_fwd)))
		if ( ! ref_no_SNP == alt_no_SNP ) snp_OK=F #  stop ( rsid, " ref and alt sequncenes dont match")
		# find the ref seq in the ref genome 

		if ( is.na(row$REF_mod) & is.na(row$ALT_mod)) {
			snp_seq=ref_fwd_s
		}
	}
	return(list(snp_OK,snp_seq,all_var_ind))
})
cat("\n")

seq_OK<-sapply(1:length(snp_OK), function(i) snp_OK[[i]][[1]][1])
ref_seq<-sapply(1:length(snp_OK), function(i) snp_OK[[i]][[2]][1])
var_ind<-sapply(1:length(snp_OK), function(i) snp_OK[[i]][[3]][1])

bad_all_vars_ind<-var_ind[!seq_OK]
bad_lib_ind<-which(!seq_OK)

if ( length(bad_all_vars_ind) > 0 ) stop ( "some library sequences contain errors") 

message( "mapping the SNPs sequences back to the genome")
# test the postiion of the ref seqs 
ref_seq_ind<- which(! is.na(ref_seq))

ref_pos_check<-sapply(1:length(ref_seq_ind), function(i) {
	map_OK=T
	ind<-ref_seq_ind[i]
	pc<-round(i/length(ref_seq_ind),2)*100
	if ( pc %% 1 == 0 ) {
		cat("working - ", pc ," % complete \r")
		flush.console()
	}
	rsid<-names(library_adapters)[ind]
	#skip controls
	row=filter(all_variants, rsid == {{rsid}})
	chr=unlist(row$chr)
	chr=paste0("chr",chr)
	pos=as.numeric(unlist(row$pos))
	ref_fwd<-ref_seq[ind]
	snp_match<-vmatchPattern(ref_fwd,ref_genome[chr])
	if( length(snp_match) > 1 ) {
		#snp_OK=F # stop ( "mapped ", rsid , "to more than one loci")
		warning ( rsid, " maps to multiple locations")
	}
	if ( ! any ( sapply(start(snp_match), function(start) start + (downstream - max_hp_len) == pos))){
		map_OK=F #stop ( rsid, "seq map to the wrong loci")
		message( rsid, " maps to the wrong loci")
	}		
	return(map_OK)
})


if ( all(ref_pos_check)) { 
	message("	all snps look good")
} else stop ( "some variants dont map to the correct location or map twice")

adapter_F_up<-unique(sapply(1:length(library_adapters), function(i) {
unique(str_sub(library_adapters[[i]]$REF_fwd,1,nchar(CRS_adpater_F)),str_sub(library_adapters[[i]]$ALT_fwd,1,nchar(CRS_adpater_F)))
}))
if ( ! (length(adapter_F_up) == 1 & adapter_F_up == CRS_adpater_F)) stop ( "adapters have errors ") 
 
adapter_F_down<-unique(sapply(1:length(library_adapters), function(i) {
ref_len<-nchar(library_adapters[[i]]$REF_fwd)
alt_len<-nchar(library_adapters[[i]]$ALT_fwd)
unique(str_sub(library_adapters[[i]]$REF_fwd,ref_len-(nchar(CRS_adpater_R)-1),ref_len),str_sub(library_adapters[[i]]$ALT_ref,alt_len-(nchar(CRS_adpater_R)-1),alt_len))
}))
if ( ! (length(adapter_F_down) == 1 & adapter_F_down == CRS_adpater_R )) stop ( "adapters have errors ") 
 
adapter_R_up<-unique(sapply(1:length(library_adapters), function(i) {
unique(str_sub(library_adapters[[i]]$REF_rev,1,nchar(CRS_adpater_R)),str_sub(library_adapters[[i]]$ALT_rev,1,nchar(CRS_adpater_R)))
}))
if ( ! (length(adapter_R_up) == 1 & adapter_R_up == as.character(reverseComplement(DNAString(CRS_adpater_R))))) stop ( "adapters have errors ") 

adapter_R_down<-unique(sapply(1:length(library_adapters), function(i) {
ref_len<-nchar(library_adapters[[i]]$REF_rev)
alt_len<-nchar(library_adapters[[i]]$ALT_rev)
unique(str_sub(library_adapters[[i]]$REF_rev,ref_len-(nchar(CRS_adpater_F)-1),ref_len),str_sub(library_adapters[[i]]$ALT_rev,alt_len-(nchar(CRS_adpater_F)-1),alt_len))
}))
if ( ! length(adapter_R_down) == 1 & adapter_R_down == as.character(reverseComplement(DNAString(CRS_adpater_F )))) stop ( "adapters have errors ")

library_out<-unlist(library_adapters)
library_out<-as.data.frame(library_out)

#library_out
library_out$name <- rownames(library_out)
library_out<-library_out[,c(2,1)]
colnames(library_out)<- c("name","sequence")

if (any(duplicated(library_out$sequence))) stop ("library contains duplicated sequences")
message("library complete")

#other_qc_filt<-all_variants$rsid[! is.na(all_variants$filtered)]
not_filt<-all_variants$rsid

prefilt<-rbindlist(loci_snps_r2)
all_pre_P_filt<-unique(prefilt$rsid)

#do a little characterisation 
u_p_snp<-unique(p_val_lead_snp)

message("	generating summary")
numbers_remaining<-sapply(names(loci_snps_r2), function(lead) {
	pre_snps<-loci_snps_r2[[lead]]
	if ( ! lead %in%  all_variants$rsid ) {
		print(lead)
		message( lead, " is not in all_variants$rsid" )
		post_snps<-all_variants[all_variants$orig_lead == lead,]
		p_val_lead<-post_snps$rsid[which.min(post_snps$P_value)]
	} else {
		p_val_lead<-all_variants$P_lead[all_variants$rsid == lead]
		post_snps<-all_variants[all_variants$P_lead == p_val_lead,]
	}
	
	n_pre<-nrow(pre_snps)
	n_post<-nrow(post_snps)

	#post_p<-snps[snps %in% all_post_p_filt]
	post_qc<-post_snps$rsid[post_snps$rsid %in% not_filt]
	n_post_qc<-length(post_qc)

	min_old_p<-min(pre_snps$P_value)
	max_old_p<-max(pre_snps$P_value)
	min_new_p<-min(post_snps$P_value)
	max_new_p<-max(post_snps$P_value)

	return(c(n_pre,n_post,n_post_qc,min_old_p,max_old_p,min_new_p,max_new_p,p_val_lead))
})

orig<-sapply(1:length(numbers_remaining) , function(i) numbers_remaining[[i]][1])
post_p<-sapply(1:length(numbers_remaining) , function(i) numbers_remaining[[i]][2])
post_qc<-sapply(1:length(numbers_remaining) , function(i)numbers_remaining[[i]][3])

min_old_p<-sapply(1:length(numbers_remaining) , function(i) numbers_remaining[[i]][4])
max_old_p<-sapply(1:length(numbers_remaining) , function(i) numbers_remaining[[i]][5])
min_new_p<-sapply(1:length(numbers_remaining) , function(i) numbers_remaining[[i]][6])
max_new_p<-sapply(1:length(numbers_remaining) , function(i) numbers_remaining[[i]][7])

p_val_lead<-sapply(1:length(numbers_remaining) , function(i) numbers_remaining[[i]][8])

numbers_remining<-data.frame(lead_snp=names(loci_snps_r2), p_val_lead = p_val_lead,original_n = orig, post_p_filtering = post_p, post_qc_filtering=post_qc,old_max_p=max_old_p,old_min_p=min_old_p, new_max_p=max_new_p)

summary_cols<- c("total_sequences","n_SNPs","n_test_SNPs","n_contol_SNPs","n_loci")
summary_res<-c(nrow(library_out),length(library_adapters),sum(test_control_SNPs$test_control == "test"),sum(test_control_SNPs$test_control == "control"),length(lead_snps2))

summary_filt<-as.data.frame(table(unlist(filtered_snps$comments)),stringsAsFactors=F)
colnames(summary_filt)<-c("filter","count")

summary_df<-matrix(summary_res,nrow=1)
summary_df<-as.data.frame(summary_df)
colnames(summary_df)<- summary_cols

message("outputing library design")

write.table(library_out,file.path(outdir,paste0("MPRA_library_sequences",suffix,".txt")),sep="\t",row.names=F, quote=F)
write.table(summary_df,file.path(outdir,paste0("MPRA_summary",suffix,".txt")),sep="\t",row.names=F, quote=F)
fwrite(all_variants,file.path(outdir,paste0("MPRA_library_all_variants",suffix,".txt")),sep="\t",row.names=F, quote=F) # doesnt include those filtered for p-val or if indel
fwrite(filtered_snps,file.path(outdir,paste0("MPRA_filtered_SNPs",suffix,".txt")), sep="\t",row.names=F, quote=F)
write.table(prefilt,file.path(outdir,paste0("MPRA_prefilter_sequences",suffix,".txt")), sep="\t",row.names=F, quote=F)
write.table(test_control_SNPs,file.path(outdir,paste0("MPRA_library_included_variants_and_cotrols",suffix,".txt")),sep="\t",row.names=F, quote=F)
write.table(numbers_remining,file.path(outdir,paste0("MPRA_summary_stats",suffix,".txt")),sep="\t",row.names=F, quote=F)
write.table(summary_filt,file.path(outdir,paste0("MPRA_summary_filtered",suffix,".txt")),sep="\t",row.names=F, quote=F)

message( "finished library design")

if ( output_figures == T ) { 
	message ("Part 4")
	message ( "	outputing figures" ) 
	oldmar<-par()$mar
	filt_reason<-summary_filt$filter
	filt_reason<-c("retained",filt_reason)
	palette = brewer.pal(length(filt_reason),'Paired')
	names(palette)<-filt_reason

	for ( snp in lead_snps)  {
	  old_snps<-filter(prefilt,  lead == snp)
	  n_vars_old=nrow(old_snps)
	  
	  new_snps<-filter(all_variants, pre_filt_lead == {{snp}})
	  n_vars_new=nrow(new_snps)

	  n_filtered=n_vars_old-n_vars_new
	  #test all remove vars are lower than the threshold
	  filt_snps<-filter(old_snps, ! rsid %in% new_snps$rsid)
	  if ( nrow( filt_snps) > 0 & nrow(new_snps) > 0) { 
		chrom=unique(old_snps$chr)
		old_snps$group=NA
		old_snps$col=NA
		old_snps$group[ old_snps$rsid %in% new_snps$rsid] <- 16
		old_snps$group[ ! old_snps$rsid %in% new_snps$rsid] <- 2
		filt_rsid<-filt_snps$rsid
		old_snps$col[ old_snps$rsid %in% new_snps$rsid] <- palette["retained"] # retained
		resons<-filtered_snps[ rsid %in% filt_rsid , ]
		
		for ( reason in filt_reason) {
			filt_reasin_rsid<-resons$rsid[resons$comment==reason]
			old_snps$col[old_snps$rsid %in% filt_reasin_rsid] <-palette[reason]
		}
		old_snps$col[ is.na(old_snps$col)] <- "yellow"
		
		pdf(file.path(outdir,paste0(snp,"_filtered_plot.pdf")))
		  par(mar=c(5, 4, 4, 12), xpd=TRUE)
		  plot(old_snps$pos,-log(old_snps$P_value,10), pch=old_snps$group,col=old_snps$col,
			xlab=paste("Position Chr",chrom),
			ylab="-log10 Pvalue",
			cex=1.5,
			main=paste(snp,paste("filtered=",n_filtered),paste("kept=",n_vars_new),sep="\n")
			)
		legend("topright",inset=c(-0.6, 0),
			   legend= c(names(palette)[palette %in% unique(old_snps$col) ]),
			   fill = c(palette[palette %in% unique(old_snps$col)])
			   )
		legend("topleft",legend=c("retained","filtered"),pch=c(16,2))
		dev.off()
		}
	}
}

message ( "script complete") 

#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################