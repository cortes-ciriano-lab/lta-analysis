# LTA detection
# To sum up, we used the following criteria to detect LTA: 
#   
# * Biallelic _TP53_ disruption. 
# * Connections between chromosome 17p and other chromosome arms with CGRs.
# * Presence of an SV in the _TP53_ gene body, or a translocation breakpoint or a CGR in the region between start of _TP53_ and the centromere.
# * Loss of the 17p terminal segment harbouring _TP53_: a  LOH in the 17p terminal segment, minimum of 75% LOH in the genomic region from the start of the chromosome up to and including _TP53_. 
# * Alternatively: LOH in >33% of the terminal segment of 17p together with an amplification (minimum 1kbp of 9 copies or more)
# * Absence of full chromosome LOH and chromosome arm LOH
# * Removing cases where the chromosome (arm) shows extensive LOH (>90%) and exists in a single copy number state (>90%). This removes aneuploidies but does not discard extensive LOH due to rearrangements.
# * Removing all cases with full chromosome arm LOH (>99.5%)

## USAGE:
# source("/nfs/research/icortes/belzen/src/TCGA_WGD_analysis/lta_detection.conf")
# dataset_selection_label="TCGA-SARC"
# or 
# flag_run_single_sample=T
# target_sample="TCGA-SARC_TCGA-PC-A5DM_TCGA-PC-A5DM-01A-11D-A813-36"
# source("/nfs/research/icortes/belzen/src/lta_detection.R")

#!/usr/bin/env bash
# SINGULARITYENV_R_MAX_VSIZE=80Gb
# singularity exec --bind /nfs/research/icortes/ /nfs/research/icortes/belzen/src/structural_variation_202405_amd.sif R -e "source('/nfs/research/icortes/belzen/src/TCGA_WGD_analysis/lta_detection.conf');flag_run_single_sample=T;target_sample='TCGA-SARC_TCGA-PC-A5DM_TCGA-PC-A5DM-01A-11D-A813-36';source('/nfs/research/icortes/belzen/src/lta_detection.R')"

## OS reproduce
#source('/nfs/research/icortes/belzen/src/TCGA_WGD_analysis/osteos.lta_detection.conf');dataset_selection_label='osteos';source('/nfs/research/icortes/belzen/src/lta_detection.R')


suppressPackageStartupMessages({
  library(GenomicRanges)
  library(VariantAnnotation)
  library(StructuralVariantAnnotation)
  library(rtracklayer)
  library(tidyverse, quietly=TRUE)
  library(stringr, quietly=TRUE)
  library(stringi)
  library(dplyr)
})

# Settings ----
chrom_order=c(paste0("chr",1:22),"chrX")
options(expressions= 500000)
max_svs_display=1000
sv_id_col="bp_name"


if(!exists("chromothripsis_clusters_path")) {
  #local
  resources_dir = "/Users/belzen/resources/"
  metadata_dir="/Users/belzen/data/metadata/"
  results_dir="/Users/belzen//results/LTA-analysis/multi_tsg/"
  
  dataset_selection_label="osteos"

if(dataset_selection_label=="osteos") {
  #osteos
  driver_files_dir="/Users/belzen/data/ebi/osteosarcoma_analysis/driver_files/"
  driver_germline_files_dir=paste0("~/data/ebi/osteosarcoma_analysis/driver_germline_files/")
  sv_filtered_files_dir=paste0("~/data/sv_filtered_files/")
  
  cn_files_dir="/Users/belzen/data/cn_files/"
  chromothripsis_clusters_path="/Users/belzen/data/ebi/osteosarcoma_analysis/shatterseek/chromothripsis_clusters.tsv"
  
  dataset_selection_label="osteos"
  
  sample_table_path=paste0(metadata_dir,"OS_master_clinical_data_v2.csv")
  driver_masterfile_path = paste0(metadata_dir,"PerDonor_master_withBurdens_mainDrivers_28_Oct2023.csv")
  suptable1_path=paste0(metadata_dir,"OS_master_table_s1.tsv")
  
}  else{
  #tcga
  #/nfs/research/icortes/belzen/src/TCGA_WGD_analysis/lta_detection.conf
  dataset_selection_label="TCGA-SARC"
  
  flag_run_single_sample=T
  target_sample="TCGA-SARC_TCGA-PC-A5DM_TCGA-PC-A5DM-01A-11D-A813-36"
  if(flag_run_single_sample && exists("target_sample")) {
    dataset_selection_label=target_sample
  }
  
  #todo make into config local later
  filesystem="/Users/belzen/data/TCGA/"
  
  wgd_score_dir = "/Users/belzen/results/TCGA_WGD_analysis/v1/"
  cn_files_dir =  paste0(filesystem,"hmf/hmf_output_files/cn_files/")
  driver_files_dir=paste0(filesystem,"hmf/hmf_output_files/driver_files/")
  driver_germline_files_dir=paste0(filesystem,"hmf/hmf_output_files/driver_germline_files/")
  sv_filtered_files_dir =paste0(filesystem,"hmf/hmf_output_files/sv_filtered_files/")
  
  chromothripsis_clusters_path = paste0(filesystem,"hmf/data_freeze_20240617.chromothripsis_clusters.tsv")
  cohort_path="~/data/metadata/TCGA.data_freeze_20240617.purity_table.tsv"
}
  
}

if(!exists("flag_run_single_sample")) { flag_run_single_sample=F }

if(!exists("dataset_selection_label")) {
  print("Warning, no dataset_selection_label, setting default")
  dataset_selection_label="osteos"
}

print(paste0("Running: ",dataset_selection_label))


# Paths ----

# templates 
chromosome_bands_path = paste0(resources_dir,"chromosome_bands.gz")
gtf_path= paste0(resources_dir,"gencode.v38.annotation.gtf.gz")
#gtf_path = stri_replace_all_fixed(gtf_path_template,names(map_template_vars), map_template_vars,vectorize=F)

#input

sv_filtered_path_template = "${sv_filtered_files_dir}/${patient_basename}.purple.sv.filtered.vcf.gz"
cn_path_template = "${cn_files_dir}/${patient_basename}.purple.cnv.somatic.tsv"
drivers_path_template = paste0("${driver_files_dir}/${patient_basename}.purple.driver.catalog.tsv")
drivers_germline_path_template = paste0("${driver_germline_files_dir}/${patient_basename}.purple.driver.catalog.germline.tsv")

#output
gene_cn_sv_disruptions_path_template = paste0("${processed_output_dir}/${patient_basename}.gene_cn_sv_disruptions.tsv")
clusters_instability_region_path_template= paste0("${processed_output_dir}/${patient_basename}.clusters_instability_region.tsv")
svs_connecting_tsg_instability_path_template= paste0("${processed_output_dir}/${patient_basename}.svs_connecting_tsg_instability.tsv")
connected_regions_tsg_instability_path_template= paste0("${processed_output_dir}/${patient_basename}.connected_regions_tsg_instability.tsv")
cohort_call_lta_path_template= paste0("${processed_output_dir}/${patient_basename}.lta_detection_overview.tsv")
call_lta_annot_path_template = paste0("${processed_output_dir}/${patient_basename}.lta_detection_overview.all_tsgs.tsv")
oncogene_amp_annot_path_template = paste0("${processed_output_dir}/${patient_basename}.oncogene_amp.annot.tsv")

reconplot_path_template = paste0("${plot_dir}/${patient_basename}.${cgr_id}.reconplot.pdf")

#path names
plot_dir=paste0(results_dir,"plots/")

map_template_vars=c('${cn_files_dir}'=cn_files_dir,
                    '${sv_filtered_files_dir}'=sv_filtered_files_dir,
                    '${driver_files_dir}'=driver_files_dir,
                    '${driver_germline_files_dir}'=driver_germline_files_dir,
                    "${plot_dir}"=plot_dir,
                    "${processed_output_dir}"=results_dir)


#output
if(flag_run_single_sample && exists("target_sample")) {
  map_template_vars_cohort = c(map_template_vars,"${patient_basename}"=target_sample)
  
}  else {
  map_template_vars_cohort = c(map_template_vars,"${patient_basename}"=dataset_selection_label)
}

gene_cn_sv_disruptions_path = stri_replace_all_fixed(gene_cn_sv_disruptions_path_template,names(map_template_vars_cohort), map_template_vars_cohort,vectorize=F)
clusters_instability_region_path = stri_replace_all_fixed(clusters_instability_region_path_template,names(map_template_vars_cohort), map_template_vars_cohort,vectorize=F)
svs_connecting_tsg_instability_path = stri_replace_all_fixed(svs_connecting_tsg_instability_path_template,names(map_template_vars_cohort), map_template_vars_cohort,vectorize=F)
connected_regions_tsg_instability_path = stri_replace_all_fixed(connected_regions_tsg_instability_path_template,names(map_template_vars_cohort), map_template_vars_cohort,vectorize=F)
cohort_call_lta_path= stri_replace_all_fixed(cohort_call_lta_path_template,names(map_template_vars_cohort), map_template_vars_cohort,vectorize=F)
call_lta_annot_path = stri_replace_all_fixed(call_lta_annot_path_template,names(map_template_vars_cohort), map_template_vars_cohort,vectorize=F)
oncogene_amp_annot_path = stri_replace_all_fixed(oncogene_amp_annot_path_template,names(map_template_vars_cohort), map_template_vars_cohort,vectorize=F)
# Functions ----
get_chr_arms = function(chromosome_bands_df,split_giestain=T,autosomes= c(paste("chr",1:22,sep=""),"chrX","chrY")) {
  chromosome_bands_df = chromosome_bands_df %>% 
    dplyr::filter(seqnames %in% autosomes) %>% 
    dplyr::mutate(chr_arm = ifelse(grepl("p",cytoband),paste0(seqnames,"p"),paste0(seqnames,"q")))
  
  if(split_giestain==T){
    chromosome_bands_df = chromosome_bands_df %>% 
      dplyr::mutate(chr_arm = ifelse(grepl("gneg",gieStain)|grepl("gpos",gieStain),chr_arm,paste0(chr_arm,"_",gieStain))) 
  }
  chr_arms = chromosome_bands_df %>% dplyr::group_by(seqnames,chr_arm) %>% 
    dplyr::summarize(start = min(start)+1,
                     end = max(end)) %>% as.data.frame() %>% ungroup()
  
  return(chr_arms %>% as.data.frame())
}


overlap_fraction_gr = function(set1,set2,ignore_strand=TRUE){
  overlaps_intersect = pintersect(set1, set2,ignore.strand=ignore_strand)
  overlap_fraction = width(overlaps_intersect) / width(set1)
  return(overlap_fraction)
}
get_reciprocal_overlap_pairs = function(set1,set2,reciprocal_overlap=0.5,svtype_matching=TRUE,ignore_strand=TRUE,flag_allow_overlaps_self=F){  
  if(is.null(names(set1))) {
    names(set1)=paste0("set1_",1:length(set1))
  }
  if(is.null(names(set2))) {
    names(set2)=paste0("set2_",1:length(set2))
  }
  # make pairs, so note that this makes the lists larger
  overlap_pairs = findOverlaps(set1,set2,type="any",minoverlap = 1,ignore.strand=ignore_strand)
  
  if(svtype_matching==TRUE){
    svtype_match = set1[overlap_pairs@from]$svtype==set2[overlap_pairs@to]$svtype
    overlap_pairs = overlap_pairs[svtype_match==TRUE]
  }
  
  #remove if not at least overlaps => reciprocal overlap
  overlap_fraction = overlap_fraction_gr(set1[overlap_pairs@from], set2[overlap_pairs@to],ignore_strand) #width of intersect/width of first
  overlap_pairs = overlap_pairs[overlap_fraction > reciprocal_overlap]
  overlap_fraction = overlap_fraction_gr(set2[overlap_pairs@to], set1[overlap_pairs@from],ignore_strand)
  overlap_pairs = overlap_pairs[overlap_fraction > reciprocal_overlap]
  
  #length of intersect divided by length of first in function 
  overlap_fraction_metadata = overlap_fraction_gr(set1[overlap_pairs@from],set2[overlap_pairs@to],ignore_strand) 
  overlap_fraction_metadata_2 = overlap_fraction_gr(set2[overlap_pairs@to],set1[overlap_pairs@from],ignore_strand)
  overlap_pairs_df = as.data.frame(overlap_pairs)
  colnames(overlap_pairs_df) = c("from","to")
  
  overlap_pairs_df$set1 = names(set1[overlap_pairs@from])
  if("svtype" %in% names(mcols(set1))) {
    overlap_pairs_df$set1_svtype = set1[overlap_pairs@from]$svtype
  }
  
  overlap_pairs_df$set2 = names(set2[overlap_pairs@to])
  if("svtype" %in% names(mcols(set2))) {
    overlap_pairs_df$set2_svtype = set2[overlap_pairs@to]$svtype
  }
  overlap_pairs_df$overlap_set1_set2 = overlap_fraction_metadata
  overlap_pairs_df$overlap_set2_set1 = overlap_fraction_metadata_2
  
  #remove overlaps with self
  if(flag_allow_overlaps_self ==F) {
    overlap_pairs_df = overlap_pairs_df[overlap_pairs_df$set1!=overlap_pairs_df$set2,]
  }
  
  return(overlap_pairs_df)
}

get_svtype <- function(gr) {
  # CTX if seq names are the same
  ## return as complex if unpartnered
  
  return_gr = gr
  unpartnered_gr=gr[!gr$partner %in% names(gr)]
  gr = gr[gr$partner %in% names(gr)]
  
  svtype = ifelse(seqnames(gr) != seqnames(gr[gr$partner]), "CTX", # inter-chromosomosal
                  ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS", # TODO: improve classification of complex events
                         ifelse(strand(gr) == strand(gr[gr$partner]), "INV",
                                ifelse(xor(start(gr) < start(gr[gr$partner]), strand(gr) == "-"), "DEL",
                                       "DUP"))))
  
  if(length(unpartnered_gr)>0){
    return_gr[names(unpartnered_gr)]$svtype="complex"
  }
  return_gr[names(gr)]$svtype=svtype
  return(  return_gr$svtype)
}

rbind_no_colmatch = function (df1,df2)  {
  if(nrow(df1)==0) {
    return(df2)
  }
  cols_1 = names(df1)
  cols_2 = names(df2)
  col_diff_1 = cols_2[!cols_2 %in% cols_1]
  col_diff_2 = cols_1[!cols_1 %in% cols_2]
  df1[,col_diff_1]=NA
  df2[,col_diff_2]=NA
  
  join_df=rbind(df1,df2)
  return(join_df)
}

get_cn_abs_call = function(segments_df,cn_abs_homloss=0.5,cn_abs_amp=8,cn_abs_euploid_thresh=0.95,cn_abs_euploid=2,cn_value_col="copyNumber",CN_round_to_digits=1) {
  if(!"shift_to_ploidy" %in% names(segments_df)){
    warning(paste0("shift_to_ploidy missing, assuming ",cn_abs_euploid))
    segments_df$shift_to_ploidy=NA
  }
  warning(paste0("CN calling uses >= and round to digits ",CN_round_to_digits))
  warning(paste0("cn_abs_amp: ",cn_abs_amp,"; cn_abs_euploid_thresh: ",cn_abs_euploid_thresh,"; cn_abs_homloss: ",cn_abs_homloss))
  
  segments_df = segments_df %>% dplyr::mutate(
    !!sym(cn_value_col):=round(!!sym(cn_value_col),CN_round_to_digits),
    cn_abs_target_euploid = ifelse(!is.na(shift_to_ploidy),readr::parse_number(shift_to_ploidy),cn_abs_euploid),
    call= ifelse(!!sym(cn_value_col) <= cn_abs_homloss,"homloss",
                 ifelse(!!sym(cn_value_col) >= cn_abs_amp,"amp",
                        ifelse(!!sym(cn_value_col) <= (cn_abs_target_euploid-cn_abs_euploid_thresh),"loss",
                               ifelse(!!sym(cn_value_col) >= (cn_abs_target_euploid+cn_abs_euploid_thresh),"gain","0")))),
    #cr_l2fc=log2(copyNumber/cn_abs_target_euploid), #that doesnt make sense
  )
  return(segments_df)
}
get_reciprocal_overlap_pairs_start_end = function(svs,properties,reciprocal_overlap=0,svtype_matching=F,ignore_strand=T){
  ## starting bp 
  svs_start = svs
  end(svs_start)=start(svs_start)
  
  ## ending bp 
  svs_end = svs
  start(svs_end)=end(svs_end)
  
  start_overlaps = get_reciprocal_overlap_pairs(svs_start,properties,reciprocal_overlap = reciprocal_overlap,svtype_matching = svtype_matching,ignore_strand=ignore_strand)
  if(nrow(start_overlaps)>0){
    start_overlaps$sv_breakpoint_orientation="start"
  }
  
  end_overlaps = get_reciprocal_overlap_pairs(svs_end,properties,reciprocal_overlap = reciprocal_overlap,svtype_matching = svtype_matching,ignore_strand=ignore_strand)
  if(nrow(end_overlaps)>0){
    end_overlaps$sv_breakpoint_orientation="end"
  }
  overlaps= rbind(start_overlaps,end_overlaps)
  
  return(overlaps)
}
# Resources -----
## chr arms ----

chromosome_bands_df = read.table(chromosome_bands_path,header=T,sep="\t",comment.char = "")
chromosome_bands_df = chromosome_bands_df %>% dplyr::rename("seqnames"="X.chrom","start"="chromStart","end"="chromEnd","cytoband"="name") 
chromosome_bands_df$cytoband = paste0(chromosome_bands_df$seqnames,chromosome_bands_df$cytoband)
chr_arms = get_chr_arms(chromosome_bands_df,autosomes = chrom_order,split_giestain = F)
chr_arms = GRanges(chr_arms)
names(chr_arms) = chr_arms$chr_arm
chr_arms_df = as.data.frame(chr_arms)
chr_arms_df = chr_arms_df %>% dplyr::mutate(chr_arm_start=start,chr_arm_end=end) 

gtf <- rtracklayer::import(gtf_path)
gene_properties = gtf[gtf$type=="gene"]
names(gene_properties) = gene_properties$gene_id
gene_properties = gene_properties[gene_properties@seqnames %in% chrom_order]
gene_properties_df = as.data.frame(gene_properties)
#gene_properties_df$gene_coord = paste0(gene_properties_df$seqnames,":",gene_properties_df$start,"-",gene_properties_df$end)
gene_properties_df = gene_properties_df %>% dplyr::mutate(gene_start = start, gene_end=end)

# Cohort  ----

##for osteos ----
if(dataset_selection_label=="osteos") {
  sample_table = read.csv(sample_table_path,header=T) 
  sample_table$basename=sample_table$name
  sample_table = sample_table %>% filter(qc_final!="blacklist" & !is.na(basename))
  
  suptable1 = read.table(suptable1_path,sep = "\t", header=T) 
  suptable1 = suptable1 %>% dplyr::mutate(Region_Rcode=Region,
                                          Region=str_replace(Region,"R","tumour-"),
                                          Region=ifelse(Region==".",NA,Region))
  sample_table = sample_table %>% dplyr::mutate(Region_Rcode=str_replace(Region,"tumour-","R"),
                                                Region_Rcode=str_replace(Region_Rcode,"metastasis-","R"),
                                                Region_Rcode=ifelse(is.na(Region_Rcode),".",Region_Rcode))
  
  sample_table = sample_table %>% left_join(suptable1[,c("final_identifier","Region_Rcode","LTA")],by=c("final_id"="final_identifier","Region_Rcode"))
  
  driver_masterfile = read.csv(driver_masterfile_path)
  sample_table = sample_table %>% mutate(selected_sample = names2 %in% driver_masterfile$names2)
  sample_table = sample_table %>% left_join(driver_masterfile %>% select(final_id,Region,CGR_17p,contains("TP53"),contains("OncoAmp")),by=c("final_id","Region")) 
  
  sample_table %>% filter(selected_sample) %>% nrow() ==  sample_table  %>% select(donor_id) %>% unique() %>% nrow()
  cohort = sample_table %>% filter(selected_sample)
  
  cohort$cohort_id=cohort$Cohort
  driver_masterfile = driver_masterfile %>% left_join(cohort[,c("names2","basename")]) 

}

##for tcga  ----
if(grepl("TCGA",dataset_selection_label)) {
cohort=read.table(cohort_path,header = T)
cohort$basename=cohort$sample_id
cohort = cohort %>% rowwise() %>% mutate(cohort_id=unlist(strsplit(sample_id, "_"))[1])

cohort = cohort %>% filter(qc=="PASS" )
cohort = cohort %>% filter(!((snv_n_chromosomes <= 15 | sv_n_chromosomes <= 15) & (snv_n < 100 & sv_n < 10)))

if(flag_run_single_sample==F & exists("dataset_selection_label")) {
  cohort = cohort %>% filter(cohort_id==dataset_selection_label)
  print(paste0("RUNNING: ",dataset_selection_label))
}
print("Number of samples")
print(cohort %>% nrow())

}

if(flag_run_single_sample && exists("target_sample")) {
  cohort = cohort %>% filter(basename==target_sample)
  if(cohort %>% nrow() != 1 ) { print("ERROR: not single sample"); q()}
  print(paste0("RUNNING: ",cohort$basename))
  
}

if(cohort %>% nrow() == 0 ) {
  print("No cohort, quitting")
  q()
}

# Load CGRs  ----
chromothripsis_clusters=read.table(chromothripsis_clusters_path,sep = "\t",header=T)

##for osteos ----
if(dataset_selection_label=="osteos") {
  #for os: subset to selected
cohort = cohort %>% mutate(region=ifelse(is.na(Region),".",Region))
chromothripsis_clusters=chromothripsis_clusters %>% merge(cohort[,c("donor_id","region","basename")])
}

#for tcga
if(grepl("TCGA",dataset_selection_label)) {
  chromothripsis_clusters$basename = chromothripsis_clusters$sample_id
  chromothripsis_clusters=chromothripsis_clusters %>% merge(cohort[,c("cohort_id","basename")])
}

#for both
chromothripsis_clusters = chromothripsis_clusters %>% mutate(cluster_id=paste0("cluster_",basename,"_",chrom))
chromothripsis_clusters = chromothripsis_clusters %>% mutate(multichromosomal=grepl("_",chrom_all))


chromothripsis_clusters_gr = GRanges(chromothripsis_clusters)
names(chromothripsis_clusters_gr)=chromothripsis_clusters_gr$cluster_id


if(chromothripsis_clusters %>% nrow() == 0) {
  print("No chromothripsis_clusters, quitting")
q()
}

# Load svs, cns, drivers ----

svs_df = data.frame()
drivers = data.frame()
drivers_germline = data.frame()
cohort_segments_df = data.frame()

for(target_sample in cohort$basename) {
  map_template_vars_patient=c(map_template_vars,'${patient_basename}'=target_sample)
  
  #sv_wgd_score_path = stri_replace_all_fixed(sv_wgd_score_path_template,names(map_template_vars_patient), map_template_vars_patient,vectorize=F)
  sv_filtered_path = stri_replace_all_fixed(sv_filtered_path_template,names(map_template_vars_patient), map_template_vars_patient,vectorize=F)
  
  if(length(Sys.glob(sv_filtered_path))==1){
    gridss_vcf = readVcf(sv_filtered_path, "hg38")
    svs_gr = breakpointRanges(gridss_vcf,nominalPosition=T,inferMissingBreakends=F)
    elementMetadata(svs_gr)["svtype"] = get_svtype(svs_gr)
    svs_sample = as.data.frame(svs_gr)
    if(svs_sample %>% nrow() == 0) {
      print(paste0("Missing SVs:",target_sample))
      next()
    }
    svs_sample = svs_sample %>% mutate(!!sym(sv_id_col) := rownames(.)) 
    
    #important to use nominal position and remove unpartnered breakends
    
    ## Annotate SVs with EBI pipeline output 
    info_vcf = as.data.frame(info(gridss_vcf)) %>% mutate(!!sym(sv_id_col) := rownames(.)) 
    join_cols = names(info_vcf)[!names(info_vcf) %in% names(svs_sample)]
    
    selected_metadata_cols = join_cols[grepl("PURPLE",join_cols)]
    
    svs_sample = svs_sample %>% left_join(info_vcf[,unique(c(sv_id_col,selected_metadata_cols))],by=sv_id_col)
    
    
    #get first item of purple AF , CN and CN change = current breakpoint
    #otherwise also cannot save
    #to indicate the difference use lowercase letters
    
    svs_sample = svs_sample %>% rowwise() %>% dplyr::mutate(
      purple_af_start=ifelse(length(PURPLE_AF)>0,PURPLE_AF[[1]] %>% as.character(),NA),
      purple_cn_start=PURPLE_CN[[1]] %>% as.character(),
      purple_cn_change_start=PURPLE_CN_CHANGE[[1]] %>% as.character(),
      purple_af_end=ifelse(length(PURPLE_AF)>0,PURPLE_AF[2][] %>% as.character(),NA),
      purple_cn_end=PURPLE_CN[2][] %>% as.character(),
      purple_cn_change_end=PURPLE_CN_CHANGE[2][] %>% as.character()) %>% as.data.frame()
    #dplyr::select(-PURPLE_AF,-PURPLE_CN,-PURPLE_CN_CHANGE)
    
    # svs_sample = read.table(sv_wgd_score_path,sep = "\t",header = T)
     svs_sample$basename=target_sample
     svs_df = rbind(svs_df,svs_sample)  
    
  } else {
    print(paste0("Missing SVs:",target_sample))
    next()
  }
  
  cn_path = stri_replace_all_fixed(cn_path_template,names(map_template_vars_patient), map_template_vars_patient,vectorize=F)
  if(Sys.glob(cn_path) %>% length() == 0) { 
    print(target_sample)
    next() 
  }
  segments_df = read.table(cn_path,sep="\t",header=T)
  segments_df = segments_df %>% dplyr::rename(seqnames=chromosome)
  segments_df$basename=target_sample
  segments_df$segment_id = paste0("seg_",1:nrow(segments_df))
  cohort_segments_df=rbind(cohort_segments_df,segments_df)
  
  drivers_path = stri_replace_all_fixed(drivers_path_template,names(map_template_vars_patient), map_template_vars_patient,vectorize=F)
  if(Sys.glob(drivers_path) %>% length() == 0) { next() }
  drivers_sample = read.table(drivers_path,sep="\t",header=T)
  if(nrow(drivers_sample)==0) {next()}
  drivers_sample$basename=target_sample
  drivers=rbind(drivers,drivers_sample)
 
  drivers_germline_path = stri_replace_all_fixed(drivers_germline_path_template,names(map_template_vars_patient), map_template_vars_patient,vectorize=F)
  if(Sys.glob(drivers_germline_path) %>% length() == 0) { next() }
  drivers_germline_sample = read.table(drivers_germline_path,sep="\t",header=T)
  if(nrow(drivers_germline_sample)==0) {next()}
  drivers_germline_sample$basename=target_sample
  drivers_germline=rbind(drivers_germline,drivers_germline_sample)
  
}

svs_df = svs_df %>% mutate(sv_id=paste0(basename,"_",bp_name),partner_sv_id=paste0(basename,"_",partner),sv_event=paste0(basename,"_",event))
svs_gr=GRanges(svs_df)
names(svs_gr)=svs_df$sv_id

cohort_segments_df = cohort_segments_df %>% mutate(cn_seg_id=paste0(basename,"_",segment_id))
cn_seg_gr = GRanges(cohort_segments_df)
names(cn_seg_gr) = cohort_segments_df$cn_seg_id

drivers = rbind(drivers,drivers_germline) %>% unique()
drivers = drivers %>% filter(driverLikelihood>=0.5) 
drivers = drivers %>% filter(chromosome != "chrY")
drivers = drivers %>% dplyr::mutate(gene_name=gene) %>% left_join(gene_properties_df %>% select(gene_name,gene_id,gene_start,gene_end)) 
drivers = drivers %>% dplyr::mutate(chr_arm = ifelse(grepl("p",chromosomeBand),paste0(chromosome,"p"),paste0(chromosome,"q")))
drivers = drivers %>% left_join(chr_arms_df %>% select(chr_arm,chr_arm_start,chr_arm_end))
drivers = drivers %>% mutate(sample_gene = paste0(basename,"_",gene_name))

#collapse multiple
driver_summary = drivers %>% group_by(basename,gene_name) %>% 
  summarize(hmf_biallelic=any(biallelic=="true"),likelihoodMethod=toString(likelihoodMethod),hmf_driver=toString(driver)) %>% ungroup() %>% as.data.frame()


#for initial screening consider all 
if(!exists("lta_tsg_disruption_lst")) {
  # driver_freq = drivers %>% filter(category=="TSG" & biallelic=="true") %>% select(basename,gene_name) %>% unique() %>% group_by(gene_name) %>% count() %>% arrange(-n) 
  # sample_cnt = cohort$basename %>% unique() %>% length()
  # lta_tsg_disruption_lst = filter(driver_freq,n>sample_cnt*0.1)$gene_name
 tsg_disruption_lst = filter(drivers,category=="TSG" & biallelic=="true")$gene_name %>% unique()
} else {
  tsg_disruption_lst=lta_tsg_disruption_lst
}

#

# TSG disruption ----

# biallelic disruption
# + any SV evidence
# overlapping SV or CGR
# LOH (minor CN <0.5)  and downstream SV 
## => any downstream SV or  more specific: CTX or CGR

#overlap gene with CN for LOH and homz loss
gene_cn_overlap = get_reciprocal_overlap_pairs(gene_properties[gene_properties$gene_name %in% tsg_disruption_lst],cn_seg_gr,reciprocal_overlap = 0,svtype_matching = F,ignore_strand = T)
gene_cn_overlap = gene_cn_overlap  %>% dplyr::rename(gene_id=set1,cn_seg_id=set2) 
gene_cn_overlap = gene_cn_overlap %>% left_join(cohort_segments_df %>% dplyr::select(seqnames,
                                                                                     start,end,
                                                                                     basename,cn_seg_id,
                                                                                     copyNumber,minorAlleleCopyNumber))
gene_cn_overlap = gene_cn_overlap  %>% left_join(gene_properties_df[,c("gene_id","gene_name")] %>% unique()) 
gene_cn_overlap = gene_cn_overlap %>% mutate(sample_gene = paste0(basename,"_",gene_name))
gene_cn_overlap = gene_cn_overlap %>% get_cn_abs_call(cn_abs_euploid_thresh=0.85,cn_abs_amp = 9)

gene_cn = gene_cn_overlap %>% group_by(seqnames,basename,gene_name,sample_gene) %>% summarize(weighted_minorCN = sum(overlap_set1_set2*minorAlleleCopyNumber),
                                                                                   weighted_CN= sum(overlap_set1_set2*copyNumber),
                                                                                   seg_cnt = length(unique(cn_seg_id)))
gene_cn_wide = gene_cn_overlap %>% group_by(basename,gene_name,call,sample_gene) %>% 
  summarize(frac = sum(overlap_set1_set2,na.rm = T),
            weighted_CN= sum(overlap_set1_set2*copyNumber,na.rm = T)/frac)

gene_cn_wide = gene_cn_wide %>% pivot_wider(names_from = "call",values_from = c("weighted_CN","frac"))

#overlap gene with sv
#check if gene is knocked out by sv
multiplicity_gene_region_sv.min=0.9
biallelic_sv_bp_tumor_af.min=0.9
sv_bp_tumor_af.min=0.5

loh_minorCN.max=0.5

gene_sv_bp_overlap = get_reciprocal_overlap_pairs(svs_gr,gene_properties[gene_properties$gene_name %in% tsg_disruption_lst],reciprocal_overlap = 0,svtype_matching = F,ignore_strand = T)
if(gene_sv_bp_overlap %>% nrow() > 0 ) {
  
gene_sv_bp_overlap = gene_sv_bp_overlap  %>% dplyr::rename(sv_id=set1,gene_id=set2)  %>% select(-to,-from)
gene_sv_bp_overlap = gene_sv_bp_overlap %>% left_join(svs_df %>% select(sv_id,basename,purple_af_start,PURPLE_JCN,sv_event))  %>% mutate(tumor_af=purple_af_start)
gene_sv_bp_overlap = gene_sv_bp_overlap %>% left_join(gene_properties_df[,c("gene_id","gene_name")] %>% unique()) %>% mutate(sample_gene = paste0(basename,"_",gene_name))

gene_sv_bp_max_af = gene_sv_bp_overlap %>% 
  group_by(sample_gene,basename,gene_name) %>% 
  summarize(max_tumor_af = max(tumor_af,na.rm = T),
            max_multiplicity = max(PURPLE_JCN,na.rm = T),
            sv_bp_cnt = length(unique(sv_id))
  ) 

gene_sv_event = gene_sv_bp_overlap %>% select(sample_gene,basename,gene_name,sv_event,PURPLE_JCN) %>% unique() %>%
  group_by(sample_gene,basename,gene_name) %>% 
  summarize(sum_sv_event_multiplicity = sum(PURPLE_JCN,na.rm = T),
            sv_event_cnt = length(unique(sv_event))
  ) 


gene_sv_count_clonal = gene_sv_bp_overlap %>% 
  filter(tumor_af>= biallelic_sv_bp_tumor_af.min | PURPLE_JCN >= multiplicity_gene_region_sv.min) %>%
  group_by(sample_gene,basename,gene_name) %>% 
  summarize(sv_event_clonal_cnt=length(unique(sv_event))
            )
            

#sum multiplicity of events 
} 
#check for arm level disruptions
#for P arm: gene overlapping or downsteam so gene start to chr arm end; for Q arm chr arm start to gene end
#chr arm from driver file or mappiing if want to check hom_disruptions

map_genes_chr_arms = get_reciprocal_overlap_pairs(gene_properties[gene_properties$gene_name %in% tsg_disruption_lst],chr_arms,reciprocal_overlap = 0,svtype_matching = F,ignore_strand = T)
map_genes_chr_arms = map_genes_chr_arms  %>% dplyr::rename(gene_id=set1,chr_arm=set2) %>% select(-from,-to,-overlap_set1_set2,-overlap_set2_set1)
map_genes_chr_arms = map_genes_chr_arms  %>% left_join(gene_properties_df[,c("gene_id","gene_name","gene_start","gene_end")] %>% unique()) 
map_genes_chr_arms = map_genes_chr_arms %>% left_join(chr_arms_df %>% select(seqnames,chr_arm,chr_arm_start,chr_arm_end))

region_gene_sv_disruption_chr_arm = map_genes_chr_arms %>% mutate(start = ifelse(grepl("p",chr_arm),gene_start,chr_arm_start),
                              end = ifelse(grepl("p",chr_arm),chr_arm_end,gene_end))

region_gene_sv_disruption_chr_arm_gr = GRanges(region_gene_sv_disruption_chr_arm)
names(region_gene_sv_disruption_chr_arm_gr) = region_gene_sv_disruption_chr_arm_gr$gene_id

#overlap with SV > enable filtering later for CTX and AF
#overlap with cgr 
gene_region_cgr_overlap = get_reciprocal_overlap_pairs(region_gene_sv_disruption_chr_arm_gr,chromothripsis_clusters_gr,reciprocal_overlap = 0,svtype_matching = F,ignore_strand = T)
gene_region_cgr_overlap = gene_region_cgr_overlap  %>% dplyr::rename(gene_id=set1,cluster_id=set2)  %>% select(-to,-from)
gene_region_cgr_overlap = gene_region_cgr_overlap %>% left_join(chromothripsis_clusters %>% select(cluster_id,basename)) 
gene_region_cgr_overlap = gene_region_cgr_overlap %>% left_join(gene_properties_df[,c("gene_id","gene_name")] %>% unique())  %>% mutate(sample_gene = paste0(basename,"_",gene_name))


gene_region_sv_bp_overlap = get_reciprocal_overlap_pairs(svs_gr,region_gene_sv_disruption_chr_arm_gr,reciprocal_overlap = 0,svtype_matching = F,ignore_strand = T)
gene_region_sv_bp_overlap = gene_region_sv_bp_overlap %>% dplyr::rename(sv_id=set1,gene_id=set2,svtype=set1_svtype)  %>% select(-to,-from)
gene_region_sv_bp_overlap = gene_region_sv_bp_overlap %>% left_join(svs_df %>% select(sv_id,basename,purple_af_start,PURPLE_JCN,sv_event))  %>% mutate(tumor_af=purple_af_start)
gene_region_sv_bp_overlap = gene_region_sv_bp_overlap %>% left_join(gene_properties_df[,c("gene_id","gene_name")] %>% unique())  %>% mutate(sample_gene = paste0(basename,"_",gene_name))



#add info on chromosomes
#to assess specificity of LOH
region_tsg_chr_arm_terminal = map_genes_chr_arms %>% mutate(start = ifelse(grepl("p",chr_arm),chr_arm_start,gene_start),
                                                            end = ifelse(grepl("p",chr_arm),gene_end,chr_arm_end))
region_tsg_chr_arm_terminal_gr = GRanges(region_tsg_chr_arm_terminal)
names(region_tsg_chr_arm_terminal_gr) = region_tsg_chr_arm_terminal$gene_id

region_tsg_chr_arm_terminal_cn_overlap = get_reciprocal_overlap_pairs(region_tsg_chr_arm_terminal_gr,cn_seg_gr,reciprocal_overlap = 0,svtype_matching = F,ignore_strand = T)
region_tsg_chr_arm_terminal_cn_overlap = region_tsg_chr_arm_terminal_cn_overlap  %>% dplyr::rename(gene_id=set1,cn_seg_id=set2) 
region_tsg_chr_arm_terminal_cn_overlap = region_tsg_chr_arm_terminal_cn_overlap %>% 
  left_join(cohort_segments_df %>%   dplyr::mutate(cn_seg_width=end-start) %>%
              dplyr::select(cn_seg_id,basename,copyNumber,minorAlleleCopyNumber,cn_seg_width)) %>%
  left_join(region_tsg_chr_arm_terminal %>% dplyr::mutate(region_width=end-start) %>% dplyr::select(gene_id,gene_name,region_width))
region_tsg_chr_arm_terminal_cn_overlap = region_tsg_chr_arm_terminal_cn_overlap %>% mutate(sample_gene = paste0(basename,"_",gene_name))

cn_tsg_chr_arm_terminal = region_tsg_chr_arm_terminal_cn_overlap %>% 
  group_by(basename,gene_name,sample_gene,region_width,gene_id) %>% 
  filter(minorAlleleCopyNumber<=loh_minorCN.max) %>% 
  summarize(seg_loh = sum(overlap_set2_set1*cn_seg_width),
            chr_arm_terminal_seg_cnt = n()) %>%
  dplyr::mutate(chr_arm_terminal_frac_loh = seg_loh / region_width) %>% as.data.frame() %>%
  dplyr::select(basename,gene_name,sample_gene,chr_arm_terminal_frac_loh)


cohort_segments_df = cohort_segments_df %>% dplyr::mutate(copyNumberInt=round(copyNumber))

chr_arms_cn_overlap = get_reciprocal_overlap_pairs(chr_arms,cn_seg_gr,reciprocal_overlap = 0,svtype_matching = F,ignore_strand = T)
chr_arms_cn_overlap = chr_arms_cn_overlap  %>% dplyr::rename(chr_arm=set1,cn_seg_id=set2) 
chr_arms_cn_overlap = chr_arms_cn_overlap %>% 
  left_join(cohort_segments_df %>%   dplyr::mutate(cn_seg_width=end-start) %>%
              dplyr::select(seqnames,cn_seg_id,basename,copyNumber,copyNumberInt,minorAlleleCopyNumber,cn_seg_width,
                            start,end)) %>%
  left_join(chr_arms_df %>% dplyr::mutate(chr_arm_width=chr_arm_end-chr_arm_start)%>% dplyr::select(chr_arm,chr_arm_width,chr_arm_start,chr_arm_end))

#todo remove chr_arm_seg_cnt because misleading is LOH seg only
cn_chr_arm = chr_arms_cn_overlap %>% 
  filter(minorAlleleCopyNumber<=loh_minorCN.max) %>% 
  group_by(basename,chr_arm,chr_arm_width) %>% 
  summarize(seg_loh = sum(overlap_set2_set1*cn_seg_width),
            chr_arm_seg_cnt = n()) %>%
  dplyr::mutate(chr_arm_frac_loh = seg_loh / chr_arm_width) %>% as.data.frame() %>%
  dplyr::select(-seg_loh,-chr_arm_width)

chromosomes_df = chr_arms_df %>%  group_by(seqnames) %>% summarize(start = min(start),end=max(end),.groups="drop")
chromosomes_df$chrom_width=chromosomes_df$end-chromosomes_df$start

cn_full_chrom = chr_arms_cn_overlap %>% left_join(chromosomes_df %>% select(seqnames,chrom_width)) %>% 
  filter(minorAlleleCopyNumber<=loh_minorCN.max) %>% 
  group_by(basename,seqnames,chrom_width) %>% 
  summarize(seg_loh = sum(overlap_set2_set1*cn_seg_width),
            chrom_frac_seg_cnt=n()) %>%
  dplyr::mutate(chrom_frac_loh = seg_loh / chrom_width) %>% as.data.frame() %>%
  dplyr::select(-seg_loh,-chrom_width)

cn_chr_arm_amp = chr_arms_cn_overlap %>% 
  filter(copyNumber>=9) %>% 
  group_by(basename,chr_arm) %>% 
  summarize(seg_amp9_kbp = sum(overlap_set2_set1*cn_seg_width)/1e3) %>% as.data.frame() 

# per sample, get chr arm frac most common CN integer state, to check for full arm loss and LOH
#only take part within chr arm region
chr_arm_cn_states = chr_arms_cn_overlap %>% 
  filter(copyNumberInt<=8 & copyNumberInt>=0) %>%
  rowwise() %>% 
  dplyr::mutate(seg_start=max(chr_arm_start,start),seg_end = min(chr_arm_end,end),seg_width=seg_end-seg_start,
                sample_chr_arm=paste0(basename,"_",chr_arm)) %>%
  group_by(basename,chr_arm,seqnames,sample_chr_arm,copyNumberInt,chr_arm_width) %>% summarize(
    bp=sum(seg_width),.groups = "drop") %>% dplyr::mutate(frac=bp/chr_arm_width) 

select_largest_frac_chr_arm = chr_arm_cn_states[chr_arm_cn_states$frac == ave(chr_arm_cn_states$frac,chr_arm_cn_states$sample_chr_arm,FUN=max),]
select_largest_frac_chr_arm = select_largest_frac_chr_arm %>% 
  dplyr::rename(chr_arm_largest_frac_cn_state_int=copyNumberInt,
                chr_arm_largest_frac_cn_state=frac) %>%
  dplyr::select(-sample_chr_arm,-chr_arm_width,-bp,-seqnames)

#only take part within chr arm region still needed to prevent frac>1 for segs in both p and q
chrom_cn_states = chr_arms_cn_overlap %>% 
  filter(copyNumberInt<=8 & copyNumberInt>=0) %>%
  left_join(chromosomes_df %>% select(seqnames,chrom_width)) %>% 
  rowwise() %>% 
  dplyr::mutate(seg_start=max(chr_arm_start,start),seg_end = min(chr_arm_end,end),seg_width=seg_end-seg_start,
                sample_chrom=paste0(basename,"_",seqnames)) %>%
  group_by(basename,seqnames,sample_chrom,copyNumberInt,chrom_width) %>% summarize(
    bp=sum(seg_width),.groups="drop") %>% dplyr::mutate(frac=bp/chrom_width) 

select_largest_frac_chrom = chrom_cn_states[chrom_cn_states$frac == ave(chrom_cn_states$frac,chrom_cn_states$sample_chrom,FUN=max),]
select_largest_frac_chrom = select_largest_frac_chrom %>% 
  dplyr::rename(chrom_largest_frac_cn_state_int=copyNumberInt,
                chrom_largest_frac_cn_state=frac) %>%
  dplyr::select(-sample_chrom,-chrom_width,-bp)


if(exists("gene_sv_bp_max_af")) {
  gene_cn_sv_disruptions = gene_cn %>% ungroup() %>%
    left_join(gene_sv_bp_max_af) %>% 
    left_join(gene_sv_count_clonal) %>%
    left_join(gene_sv_event)
} else {
  #otherwise add empty cols
  gene_cn_sv_disruptions = gene_cn %>% ungroup() 
  gene_cn_sv_disruptions$max_tumor_af=NA
  gene_cn_sv_disruptions$max_multiplicity=NA
  gene_cn_sv_disruptions$sum_sv_event_multiplicity=NA
  gene_cn_sv_disruptions$sv_event_cnt=NA
  gene_cn_sv_disruptions$sv_event_clonal_cnt=NA
  
  
}

# additional criteria for evidence of dicentric chromosomes:
# 1) loss of terminal chromosomal fraction extending from the TSG: tsg_chr_arm_terminal_frac_loh >= 0.75
# 2) not full chromosome (arm) LOH, but local loss due to sv with breakpoint at the TSG side of centromere: 
# filtering removes simple focal deletions and cases where gene is KO by single bp in haploid chromosome arm
# does not discard LOH due to rerrangements due to requirement of same copy number state

gene_cn_sv_disruptions = gene_cn_sv_disruptions %>%
  left_join(driver_summary) %>%
   dplyr::mutate( hmf_biallelic =ifelse(is.na(hmf_biallelic),F,hmf_biallelic),
  gene_sv_bp= ifelse(!is.na(max_tumor_af), (max_tumor_af >= sv_bp_tumor_af.min | (max_multiplicity/weighted_CN) >= sv_bp_tumor_af.min), F),
  gene_loh=weighted_minorCN<=loh_minorCN.max,
  gene_sv_bp_biallelic = ifelse(!is.na(sum_sv_event_multiplicity), ( (gene_loh & max_tumor_af >= biallelic_sv_bp_tumor_af.min) | 
                                  (sum_sv_event_multiplicity/weighted_CN) >= biallelic_sv_bp_tumor_af.min ), F),
  
  gene_region_cgr = sample_gene %in% gene_region_cgr_overlap$sample_gene,
  gene_region_ctx = sample_gene %in% filter(gene_region_sv_bp_overlap,svtype=="CTX")$sample_gene) %>%
  left_join(map_genes_chr_arms) %>% 
  left_join(cn_tsg_chr_arm_terminal,by=c("basename", "gene_name", "sample_gene")) %>% 
  left_join(cn_chr_arm,by=c("basename", "chr_arm")) %>%
  left_join(cn_full_chrom,by=c("basename", "seqnames")) %>%
  left_join(cn_chr_arm_amp,by=c("basename", "chr_arm")) %>%
  left_join(select_largest_frac_chrom,by=c("basename", "seqnames")) %>%
  left_join(select_largest_frac_chr_arm,by=c("basename", "chr_arm"))

gene_cn_sv_disruptions = gene_cn_sv_disruptions %>%
  dplyr::mutate(seg_amp9_kbp=ifelse(is.na(seg_amp9_kbp),0,seg_amp9_kbp),
                chr_arm_largest_frac_cn_state=ifelse(is.na(chr_arm_largest_frac_cn_state),0,chr_arm_largest_frac_cn_state),
                chrom_largest_frac_cn_state=ifelse(is.na(chrom_largest_frac_cn_state),0,chrom_largest_frac_cn_state)
  ) %>%
  dplyr::mutate(
    gene_knockout =  hmf_biallelic | (!is.na(hmf_driver) & gene_sv_bp) | gene_sv_bp_biallelic | weighted_CN<=loh_minorCN.max,
    gene_knockout_sv_component = gene_loh & gene_knockout,
    gene_knockout_dicentric = ((gene_sv_bp_biallelic | (gene_knockout & gene_loh & (gene_region_cgr | gene_region_ctx | gene_sv_bp) )) & 
                                 (chr_arm_terminal_frac_loh >= 0.75 | (chr_arm_terminal_frac_loh>=0.33 & seg_amp9_kbp>1e3)) &
                                 ( ! (chr_arm_frac_loh >0.9 & chr_arm_largest_frac_cn_state>0.9 ) & chr_arm_frac_loh<0.995 ) &
                                 ( ! (chrom_frac_loh >0.9 & chrom_largest_frac_cn_state >0.9) )
    

    )
  )

#notes
#gene_knockout_sv_component_strict is obsolete 
# no longer filter on multiplicity for dicentric v3 gene_region_ctx = sample_gene %in% filter(gene_region_sv_bp_overlap,svtype=="CTX"&PURPLE_JCN>=multiplicity_gene_region_sv.min)$sample_gene,

if(!exists("lta_tsg_disruption_lst")) {
  tsg_biallelic_ko_freq = gene_cn_sv_disruptions %>% filter(gene_knockout_sv_component) %>% select(basename,gene_name) %>% unique() %>% group_by(gene_name) %>% count() %>% arrange(-n) 
  sample_cnt = cohort$basename %>% unique() %>% length()
  print(tsg_biallelic_ko_freq)
  lta_tsg_disruption_lst = filter(tsg_biallelic_ko_freq,n>sample_cnt*0.10)$gene_name
}
  

# Regions of instability / oncogene amp ----
# for OS didnt require oncogene but CGR "region of instability"  
# make regions of all CGRs then filter on region of interest: oncogene, multichromosomal etc

instability_region = chromothripsis_clusters  %>% dplyr::rename(cluster_start=start,cluster_end=end)

#can be multiple chr arms
map_clusters_chr_arms = get_reciprocal_overlap_pairs(chromothripsis_clusters_gr,chr_arms,reciprocal_overlap = 0,svtype_matching = F,ignore_strand = T)
map_clusters_chr_arms = map_clusters_chr_arms %>% dplyr::rename(cluster_id=set1,chr_arm=set2)
instability_region = instability_region %>% left_join(map_clusters_chr_arms[,c("cluster_id","chr_arm")] %>% unique())
instability_region = instability_region %>% left_join(chr_arms_df %>% select(chr_arm,chr_arm_start,chr_arm_end))
instability_region$region_id=paste0(instability_region$cluster_id,"_",instability_region$chr_arm)

#from drivers, get amplified oncogenes
# oncogene to chr arm
# cgr in chrom arm
# for region take chrom arm  OR CGR + margin

drivers_amplified = drivers %>% filter(driver=="AMP")  %>% select(basename,gene_name,maxCopyNumber,chr_arm) %>% unique()

instability_region_oncogene_amp_long = instability_region %>% merge(drivers_amplified,by=c("basename","chr_arm"))
instability_region_oncogene_amp = instability_region_oncogene_amp_long %>% group_by(basename,chr_arm,cluster_id) %>% 
  summarize(onco_gene_amp_lst = toString(unique(paste0(gene_name," (",round(maxCopyNumber,2),") "))))

cgr_multichrom_oncogene_amp = instability_region_oncogene_amp_long %>% group_by(basename,chrom_all) %>% 
  summarize(cgr_multichrom_onco_gene_amp_lst = toString(unique(paste0(gene_name," (",round(maxCopyNumber,2),") "))))

instability_region = instability_region %>% left_join(instability_region_oncogene_amp) %>% left_join(cgr_multichrom_oncogene_amp) %>% 
  mutate(onco_amp_in_chr_arm = !is.na(onco_gene_amp_lst),
         onco_amp_in_cgr_multichrom = !is.na(cgr_multichrom_onco_gene_amp_lst))


# Connections: Link tsg region to onco region ----
#get disrupted tsgs, look up chromosome arm of genem append chrom arm location and make that into range 
#for question of disruptive SV can use for p arm: gene start to chrom end (towards centomere), and for q arm chrom start to gene end

if(dataset_selection_label=="osteos") {
  tsg_region = gene_cn_sv_disruptions %>% filter(gene_name %in% lta_tsg_disruption_lst) %>% as.data.frame()
} else {
  tsg_region = gene_cn_sv_disruptions %>% filter(gene_name %in% lta_tsg_disruption_lst & gene_knockout_sv_component) %>%  as.data.frame()
}

tsg_region$region_id = paste0(tsg_region$basename,"_",tsg_region$gene_name)
tsg_region$region_type="tsg_region"

coord_start_col="chr_arm_start"
coord_end_col="chr_arm_end"
tsg_region_id_col="region_id"
tsg_region_gr = tsg_region %>% dplyr::rename(start=!!sym(coord_start_col),end=!!sym(coord_end_col)) %>% GRanges()
names(tsg_region_gr) = tsg_region[,tsg_region_id_col]


instability_region$region_type="instability_region"
coord_start_col="chr_arm_start"
coord_end_col="chr_arm_end"
instability_region_id_col="region_id"

instability_region_gr = instability_region %>% dplyr::rename(start=!!sym(coord_start_col),end=!!sym(coord_end_col)) %>% GRanges()
names(instability_region_gr) = instability_region[,instability_region_id_col]


## add regions together
#map svs inside and see if start/end connect tsg and onco 
region_cols=c("region_id","region_type","chr_arm","basename")
region_collection = rbind_no_colmatch(
   tsg_region %>% select(all_of(region_cols)),
   instability_region %>% select(all_of(region_cols)))

regions_gr = c(instability_region_gr[,region_cols],tsg_region_gr[,region_cols])


map_svs_regions = get_reciprocal_overlap_pairs(svs_gr[svs_gr$svtype=="CTX"],regions_gr,reciprocal_overlap = 0,svtype_matching = F,ignore_strand = T)
map_svs_regions = map_svs_regions %>% dplyr::rename(sv_id=set1,region_id=set2,svtype=set1_svtype) %>% select(-to,-from,-overlap_set1_set2,-overlap_set2_set1)
map_svs_regions = map_svs_regions %>% left_join(region_collection %>% select(region_id,basename,region_type) %>% unique())
#for the basename match
map_svs_regions = map_svs_regions %>% merge(svs_df %>% select(sv_id,basename,partner_sv_id))

#nb: chromosome arms etc are from perspective of instability region
map_svs_regions = map_svs_regions %>% left_join(instability_region %>% select(basename,region_id,chr_arm))

#annotation tsg chr arm for plots
svs_connecting_tsg_instability = map_svs_regions %>% filter(region_type=="instability_region") %>% 
  filter(partner_sv_id %in% filter(map_svs_regions,region_type=="tsg_region")$sv_id) %>%
  left_join( filter(map_svs_regions,region_type=="tsg_region") %>% select(region_id,sv_id) %>% unique() %>%
           dplyr::rename(tsg_region_id=region_id),by=c("partner_sv_id"="sv_id")) %>%
  left_join(tsg_region %>% dplyr::rename(tsg_chr_arm=chr_arm) %>% 
              select(region_id,tsg_chr_arm) %>% unique(),
            by=c("tsg_region_id"="region_id")) %>% 
  left_join(svs_df %>% select(sv_id,basename,purple_af_start,PURPLE_JCN)) 

#+ annotation about disrupted tsg 
#+ instability region
selected_connected_tsg_regions = svs_connecting_tsg_instability %>%
#  select(basename,region_id,chr_arm,multichromosomal,onco_amp_in_chr_arm,onco_gene_amp_lst,tsg_region_id,tsg_chr_arm) %>% unique() %>%
  group_by(basename,region_id,chr_arm,tsg_region_id,tsg_chr_arm) %>%
  summarize(connecting_svs = length(unique(sv_id)),
            ctx_max_tumor_af = max(purple_af_start,na.rm = T),
            ctx_max_multiplicity = max(PURPLE_JCN,na.rm = T),.groups = "drop") %>% 
  left_join(instability_region %>% select(region_id,chrom_all,classification,cluster_id,multichromosomal,onco_amp_in_chr_arm,onco_gene_amp_lst,onco_amp_in_cgr_multichrom,cgr_multichrom_onco_gene_amp_lst,chrom) ) %>%
  dplyr::rename_with(.fn=function(x){paste0("cgr_",x)},.cols=c("multichromosomal","chrom_all","classification")) %>%
  left_join(tsg_region %>% dplyr::mutate(chrom=seqnames) %>%
              select(region_id,chrom,gene_name,weighted_CN,weighted_minorCN,max_tumor_af,max_multiplicity,contains("hmf"),
                     gene_loh,gene_sv_bp_biallelic,gene_region_cgr,gene_region_ctx,gene_region_ctx,gene_sv_bp,
                     #gene_knockout_dicentric_prev,gene_knockout_dicentric2,gene_knockout_sv_component_strict
                     gene_knockout_sv_component,gene_knockout_dicentric,
                     ) %>% unique() %>%
              dplyr::rename_with(.fn=function(x){paste0("tsg_",x)}),
            by=c("tsg_region_id"))  %>%
  as.data.frame()

#for "multichromosomal CGR"  uses "chrom_all" field so can be 1 chrom + the chromosome of the tsg itself...

assess_multichromosomal = chromothripsis_clusters %>% select(chrom_all,basename,chrom) %>% dplyr::rename(cgr_chrom_all=chrom_all) %>% 
  merge(selected_connected_tsg_regions[,c("basename","cgr_chrom_all","tsg_chrom","tsg_region_id")] %>% unique()) %>%
  anti_join(selected_connected_tsg_regions %>% select(basename,tsg_region_id,tsg_chrom,cgr_chrom_all) %>% 
              dplyr::rename(chrom=tsg_chrom))

assess_multichromosomal = assess_multichromosomal %>% group_by(basename,cgr_chrom_all,tsg_chrom) %>% summarize(chrom_cnt=length(unique(chrom)),
                                                                                                     chrom_lst_excl_tsg_region=toString(sort(unique(factor(chrom,levels=chrom_order)))))
selected_connected_tsg_regions = selected_connected_tsg_regions %>% left_join(assess_multichromosomal) %>% 
  mutate(cgr_multichromosomal_excl_tsg_region = !is.na(chrom_cnt) & chrom_cnt>1)

#LTA multichromosomal if CGR is multichromosomal but also if tsg connected to two different clusters
#from perspective of tsg region: connected instability regions
count_connected_tsg_regions = selected_connected_tsg_regions %>% group_by(basename,tsg_region_id) %>% 
  summarize(instability_region_cnt = length(unique(cluster_id)),instability_region_chrom_cnt = length(unique(chrom)))

selected_connected_tsg_regions = selected_connected_tsg_regions %>% left_join(count_connected_tsg_regions) %>% mutate(tsg_to_multiple_chrom=ifelse(!is.na(instability_region_chrom_cnt),instability_region_chrom_cnt>1,F))
#selected_connected_tsg_regions %>% filter(basename=="gel-100k_215000400-tumour-1_LP3000428-DNA_A01")


# Call LTA ----
#extend to amp in other cgr arms for multichromosomal ltas 
#onco_amp_in_cgr_multichrom includes direct onco_amp_in_chr_arm
#kept meaning of lta_onco and lta_multichrom the same for legacy purposes

gene_cn_sv_disruptions.call_lta = gene_cn_sv_disruptions %>% 
  dplyr::mutate(tsg_chrom=seqnames,tsg_gene_name=gene_name) %>% 
  select(basename,tsg_gene_name,tsg_chrom) %>% 
  left_join(selected_connected_tsg_regions) %>% 
  group_by(basename,tsg_gene_name,tsg_chrom,tsg_region_id) %>% 
  summarize( 
    lta_any = any(!is.na(region_id)), #instability region from the selected table
    lta = any(!is.na(tsg_gene_knockout_sv_component) & tsg_gene_knockout_sv_component), #the knockout sv component is from connected table, so know that lta is detected
    lta_multichrom =  any( ((!is.na(tsg_to_multiple_chrom) & tsg_to_multiple_chrom) | (!is.na(cgr_multichromosomal_excl_tsg_region) & cgr_multichromosomal_excl_tsg_region) ) &
                             !is.na(tsg_gene_knockout_sv_component) & tsg_gene_knockout_sv_component),
    lta_onco =  any(!is.na(onco_amp_in_chr_arm) & onco_amp_in_chr_arm & !is.na(tsg_gene_knockout_sv_component) & tsg_gene_knockout_sv_component),
    lta_onco_direct =  lta_onco,
    lta_onco_incl_indirect =  any(!is.na(onco_amp_in_cgr_multichrom) & onco_amp_in_cgr_multichrom & !is.na(tsg_gene_knockout_sv_component) & tsg_gene_knockout_sv_component),
    lta_onco_multichrom =  any( !is.na(onco_amp_in_chr_arm) & onco_amp_in_chr_arm &  
                                  ((!is.na(tsg_to_multiple_chrom) & tsg_to_multiple_chrom) | (!is.na(cgr_multichromosomal_excl_tsg_region) & cgr_multichromosomal_excl_tsg_region) ) &
                                  !is.na(tsg_gene_knockout_sv_component) & tsg_gene_knockout_sv_component),
    lta_onco_multichrom_incl_indirect =  any( !is.na(onco_amp_in_cgr_multichrom) & onco_amp_in_cgr_multichrom &  
                                  ((!is.na(tsg_to_multiple_chrom) & tsg_to_multiple_chrom) | (!is.na(cgr_multichromosomal_excl_tsg_region) & cgr_multichromosomal_excl_tsg_region) ) &
                                  !is.na(tsg_gene_knockout_sv_component) & tsg_gene_knockout_sv_component),
    .groups="drop"
  )

#knockout of gene regardless of lta
gene_cn_sv_disruptions.call_lta = gene_cn_sv_disruptions.call_lta %>% 
  left_join(gene_cn_sv_disruptions %>% select(basename,gene_name,weighted_CN,weighted_minorCN,max_tumor_af,max_multiplicity,contains("hmf"),
                                              gene_loh,gene_sv_bp_biallelic,gene_region_cgr,gene_region_ctx,gene_sv_bp,
                                              chr_arm_terminal_frac_loh,chr_arm_frac_loh,chrom_frac_loh,seg_amp9_kbp,chr_arm_largest_frac_cn_state,chrom_largest_frac_cn_state,
                                              #gene_knockout_sv_component_strict,gene_knockout_dicentric_prev,gene_knockout_dicentric2,
                                              gene_knockout,gene_knockout_sv_component,gene_knockout_dicentric
                                              ) %>% unique() %>%
              dplyr::rename_with(.cols=-basename,.fn=function(x){paste0("tsg_",x)}) )

gene_cn_sv_disruptions.call_lta = gene_cn_sv_disruptions.call_lta %>% dplyr::mutate(connected_to_tp53_tsg_region = 
                                                                                      tsg_region_id %in% filter(selected_connected_tsg_regions,region_id %in% filter(selected_connected_tsg_regions,tsg_gene_name=="TP53")$region_id)$tsg_region_id)

gene_cn_sv_disruptions.call_lta = gene_cn_sv_disruptions.call_lta %>% mutate( across(.cols=c(contains("knockout"),contains("lta"),connected_to_tp53_tsg_region), ~replace_na(.x, F)) )


cohort_call_lta = cohort  %>% 
  select(cohort_id,basename,contains("final_id"),contains("TP53"),contains("LTA"),contains("subtype_short")) %>%
  dplyr::mutate(sample_analysed=basename %in% svs_df$basename) %>%
  left_join(gene_cn_sv_disruptions.call_lta %>% filter(tsg_gene_name=="TP53") %>% select(basename,tsg_gene_name,contains("lta"),contains("knockout")))

cohort_call_lta = cohort_call_lta %>% mutate( across(.cols=c(contains("knockout"),contains("lta")), ~replace_na(.x, F)) )


selected_connected_tsg_regions_annot = selected_connected_tsg_regions %>% left_join(gene_cn_sv_disruptions.call_lta %>% select(basename,tsg_region_id,contains("lta"),connected_to_tp53_tsg_region))


oncogene_amp_annot = drivers_amplified %>% 
  left_join(instability_region_oncogene_amp_long %>% select(basename,chr_arm,gene_name,cluster_id,region_id),by=c("basename", "gene_name", "chr_arm")) %>%
  left_join(selected_connected_tsg_regions_annot %>%
              select(region_id,cgr_chrom_all,cgr_classification,tsg_region_id,tsg_gene_name,tsg_chr_arm,tsg_gene_knockout_sv_component,tsg_gene_knockout_dicentric,
                     #tsg_gene_knockout_sv_component_strict,tsg_gene_knockout_dicentric_prev,tsg_gene_knockout_dicentric2,
                     tsg_to_multiple_chrom,cgr_multichromosomal_excl_tsg_region,
                     contains("lta"),connected_to_tp53_tsg_region),relationship = "many-to-many")
#TODO add indirect amp here too but is tricky 
#selected_connected_tsg_regions_annot %>% filter(tsg_gene_name=="TP53") %>% filter(basename %in% filter(cohort_call_lta,tsg_gene_name=="TP53" & lta_onco_multichrom & !lta_onco_direct)$basename) %>% select(basename,cgr_chrom_all) %>% unique()

# oncogene_amp_annot %>% filter(basename %in% filter(cohort_call_lta,tsg_gene_name=="TP53" & lta_onco_multichrom & !lta_onco_direct)$basename) %>%
#   select(gene_name,cgr_chrom_all,chrom_all,lta_onco_multichrom)

oncogene_amp_annot = oncogene_amp_annot %>% mutate(
  across(.cols=c(contains("knockout"),contains("lta"),connected_to_tp53_tsg_region), ~replace_na(.x, F))
)


#look at instability regions connected to multiple tsg disrupted regions indicating same event
#selected_connected_tsg_regions_annot %>% filter(lta) %>% group_by(region_id) %>% summarize(toString(unique(tsg_region_id)))


#gene disruptions but then ltas only
call_lta_annot = gene_cn_sv_disruptions.call_lta %>% filter(lta_any)

#Export -----

write.table(gene_cn_sv_disruptions.call_lta,gene_cn_sv_disruptions_path,sep = "\t",row.names = F,col.names = T)
write.table(instability_region,clusters_instability_region_path,sep = "\t",row.names = F,col.names = T)
write.table(svs_connecting_tsg_instability,svs_connecting_tsg_instability_path,sep = "\t",row.names = F,col.names = T)
write.table(selected_connected_tsg_regions_annot,connected_regions_tsg_instability_path,sep = "\t",row.names = F,col.names = T)
write.table(cohort_call_lta,cohort_call_lta_path,sep = "\t",row.names = F,col.names = T)
write.table(call_lta_annot,call_lta_annot_path,sep = "\t",row.names = F,col.names = T)
write.table(oncogene_amp_annot,oncogene_amp_annot_path,sep = "\t",row.names = F,col.names = T)

if(!exists("skip_plots")) { skip_plots=F }
if(skip_plots) {
   q()
}

# Reconplots for multi TSG ----
library(ReConPlot)

if(dataset_selection_label=="osteos") {
  plot_cases = call_lta_annot %>% filter(lta_any)
  plot_cases = call_lta_annot %>% filter(lta)
  
} else {
  plot_cases = call_lta_annot %>% filter(lta_onco&tsg_gene_knockout_dicentric)
  #plot_cases = call_lta_annot %>% filter(lta_onco_multichrom)
}

plot_cases = plot_cases  %>% arrange(desc(tsg_gene_name)) %>% arrange(-tsg_gene_knockout_dicentric,-lta_onco)

for(target_tsg_region_id in plot_cases$tsg_region_id) {
  sample=filter(call_lta_annot,tsg_region_id==target_tsg_region_id)
  target_sample=sample$basename
  plot_path_label=paste0("lta")
  if("LTA" %in% names(sample)) {
    plot_path_label=paste0(plot_path_label,".truth_set=",sample$LTA)
  }
  plot_path_label=paste0(plot_path_label,".detected=",sample$lta)
  if(sample$lta) {
    plot_path_label=paste0(plot_path_label,".onco_multichrom=",sample$lta_onco_multichrom,".multichrom=",sample$lta_multichrom,".onco=",sample$lta_onco)
  }

  if("tsg_gene_knockout_dicentric" %in% names(sample)) {
    plot_path_label=paste0(plot_path_label,".tsg_gene_knockout_dicentric=",sample$tsg_gene_knockout_dicentric)
  }
  #if("tsg_gene_knockout_sv_component_strict" %in% names(sample)) {
  #   plot_path_label=paste0(plot_path_label,".tsg_gene_knockout_sv_component_strict=",sample$tsg_gene_knockout_sv_component_strict)
  # }
  #  # if("connected_to_tp53_tsg_region" %in% names(sample)) {
   #   plot_path_label=paste0(plot_path_label,".connected_to_tp53_tsg_region=",sample$connected_to_tp53_tsg_region)
   # }
  
  map_template_vars_patient=c(map_template_vars,'${patient_basename}'=target_tsg_region_id)
  
#chr arm region 
#include connections found
#all svs of cgrs of the multichrom constructs 

selected_instability_region = instability_region %>% filter(region_id %in% filter(selected_connected_tsg_regions,tsg_region_id==target_tsg_region_id)$region_id)
selected_connecting_svs = svs_connecting_tsg_instability %>% 
  filter(tsg_region_id==target_tsg_region_id) %>%
  filter(region_id %in% selected_instability_region$region_id) %>%
  left_join(selected_instability_region %>% select(region_id,chrom_all))

map_template_vars_cgr = c(map_template_vars_patient,'${cgr_id}'=plot_path_label)
reconplot_path = stri_replace_all_fixed(reconplot_path_template,names(map_template_vars_cgr), map_template_vars_cgr,vectorize=F)

plot_chrom = str_count(selected_instability_region$chrom_all,"_") %>% max() + 1
plot_width = ifelse(plot_chrom<4, 35, 35+(2.5*(plot_chrom-3)))

if(length(Sys.glob(reconplot_path))==1) {next()}

pdf(reconplot_path, width = plot_width/2.2, height = 8/2.2,pointsize = 8)


#plot the entire thing and plot per multichrom group (chrom_all)
multichrom_group_lst = c("all",selected_instability_region$chrom_all %>% unique())
multichrom_group="all"
for(multichrom_group in multichrom_group_lst) {
  if(multichrom_group=="all" ) { 
    plot_selected_instability_region = selected_instability_region 
  } else {
    plot_selected_instability_region = selected_instability_region %>% filter(chrom_all==multichrom_group)
  }

#expand to multichrom 
selected_cgrs = chromothripsis_clusters %>% merge(plot_selected_instability_region[,c("basename","chrom_all")])
selected_cgrs = selected_cgrs %>% left_join(map_clusters_chr_arms[,c("cluster_id","chr_arm")])

map_svs_target_cgrs = get_reciprocal_overlap_pairs(svs_gr[svs_gr$basename == target_sample],chromothripsis_clusters_gr[chromothripsis_clusters_gr$cluster_id %in% selected_cgrs$cluster_id,],reciprocal_overlap = 0,svtype_matching = F,ignore_strand = T)
map_svs_target_cgrs = map_svs_target_cgrs %>% dplyr::rename(sv_id=set1,cluster_id=set2,svtype=set1_svtype) %>% select(-to,-from,-overlap_set1_set2,-overlap_set2_set1)

plot_selected_connecting_svs = selected_connecting_svs %>% filter(chrom_all %in% plot_selected_instability_region$chrom_all)

plot_svs = c(plot_selected_connecting_svs$sv_id,
             plot_selected_connecting_svs$partner_sv_id,
             map_svs_target_cgrs$sv_id) %>% unique()

#plot_gene_lst=c("TP53")
plot_gene_lst=lta_tsg_disruption_lst

plot_genes = instability_region_oncogene_amp_long %>% merge(plot_selected_instability_region[,c("basename","chr_arm")])
if(nrow(plot_genes)>0) {
  plot_gene_lst = c(plot_gene_lst,plot_genes$gene_name) 
}
 
sv_data = svs_df %>% 
  filter(sv_id %in% plot_svs) %>% 
  select(basename,bp_name,start,seqnames,strand) %>% 
  dplyr::rename(chr1=seqnames,pos1=start) %>%
  left_join(svs_df %>% select(basename,partner,start,seqnames,strand) %>% 
              dplyr::rename(chr2=seqnames,pos2=start),by=c("bp_name"="partner","basename")) %>%
  dplyr::mutate(strands=paste0(strand.x,strand.y)) %>% select(-strand.x,-strand.y) %>%
  select(chr1,pos1,chr2,pos2,strands) %>% as.data.frame()

  sv_data = sv_data %>% filter(chr1 %in% chrom_order & chr2 %in% chrom_order)

  plot_region_instability_only = c(plot_selected_connecting_svs$chr_arm,plot_selected_connecting_svs$tsg_chr_arm,plot_selected_instability_region$chr_arm) %>% unique()
  
  plot_region_instability_only = chr_arms_df %>%
    filter(chr_arm %in% plot_region_instability_only) %>% arrange(seqnames) %>%
    select(seqnames,start,end) %>% dplyr::mutate(chr1=as.character(seqnames),chr=chr1,start=1,end=500e6) %>% as.data.frame() %>% unique()
  sv_data_instability_only = sv_data %>% filter(chr1 %in% plot_region_instability_only$seqnames)
  
  
  cn_data = cohort_segments_df %>% 
    dplyr::rename(chr=seqnames) %>% 
    filter(basename==target_sample & chr %in% sv_data$chr1) %>%
    mutate(copyNumber=round(copyNumber),minorAlleleCopyNumber=round(minorAlleleCopyNumber)) %>%
    select(chr,start,end,copyNumber,minorAlleleCopyNumber) %>% as.data.frame()
  
  plot_region_arms = c(plot_selected_connecting_svs$chr_arm,plot_selected_connecting_svs$tsg_chr_arm,selected_cgrs$chr_arm)
  cgrs_chr_arm_coord = chr_arms_df %>%
    filter(chr_arm %in% plot_region_arms) %>%
    dplyr::mutate(chr1=seqnames,chr=seqnames) %>%
    group_by(chr1,chr) %>% summarize(start=min(start),end=max(end)) %>% ungroup() %>% as.data.frame()
  
  #otherwise error
  cgrs_chr_arm_coord$chr = as.character(cgrs_chr_arm_coord$chr)
  cgrs_chr_arm_coord$chr1 = as.character(cgrs_chr_arm_coord$chr1)

  #request from jose full chromosomes
  full_chrom = cgrs_chr_arm_coord 
  full_chrom$start=1
  full_chrom$end=500e6
  
  if(nrow(sv_data) > max_svs_display) {
    print(paste0("Too many SVs to display, skipping full plot", nrow(sv_data)))
  } else {
  
  p = ReConPlot(genes = plot_gene_lst,
                sv_data,
                cn_data,
                chr_selection=full_chrom,
                legend_SV_types=T,
                pos_SVtype_description=1000000,
                scale_separation_SV_type_labels=1/23,
                title=paste0(multichrom_group," ",target_sample),
                max.cn = 12,
                label_interchr_SV = TRUE,
  )
  print(p)
  }
  
  if(full_chrom %>% filter(!chr %in% plot_region_instability_only$chr) %>% nrow() == 0 ) {next()}
  if(nrow(sv_data_instability_only) > max_svs_display) {
    print(paste0("Too many SVs to display, skipping instability only plot", nrow(sv_data_instability_only)))
    print("Last resort: only TSG arm")
    plot_tsg_region_only = c(plot_selected_connecting_svs$tsg_chr_arm) %>% unique()
    
    plot_tsg_region_only = chr_arms_df %>%
      filter(chr_arm %in% plot_tsg_region_only) %>% arrange(seqnames) %>%
      select(seqnames,start,end) %>% dplyr::mutate(chr1=as.character(seqnames),chr=chr1,start=1,end=500e6) %>% as.data.frame() %>% unique()
    sv_data_tsg_region_only = sv_data %>% filter(chr1 %in% plot_tsg_region_only$seqnames)
    
    if(nrow(sv_data_tsg_region_only) <= max_svs_display) {
      print(paste0("Last resort: only TSG arm"))
      
    p = ReConPlot(genes = plot_gene_lst,
                  sv_data_tsg_region_only,
                  cn_data %>% filter(chr %in% plot_tsg_region_only$seqnames),
                  chr_selection=plot_tsg_region_only,
                  legend_SV_types=T,
                  pos_SVtype_description=1000000,
                  scale_separation_SV_type_labels=1/23,
                  title=paste0(multichrom_group," ",target_sample," - plot tsg region only"),
                  max.cn = 12,
                  label_interchr_SV = TRUE,
    )
    print(p)
    }
  } else {
  p = ReConPlot(genes = plot_gene_lst,
                sv_data_instability_only,
                cn_data %>% filter(chr %in% plot_region_instability_only$seqnames),
                chr_selection=plot_region_instability_only,
                legend_SV_types=T,
                pos_SVtype_description=1000000,
                scale_separation_SV_type_labels=1/23,
                title=paste0(multichrom_group," ",target_sample," - plot chrom instability only"),
                max.cn = 12,
                label_interchr_SV = TRUE,
  )
  print(p)
  }
  
  if(full_chrom %>% filter(!chr %in% plot_region_instability_only$chr) %>% nrow() == 0 ) {next()}
  if(nrow(sv_data_instability_only) > max_svs_display) {
    print(paste0("Too many SVs to display, skipping instability only plot", nrow(sv_data)))
  } else {
    p = ReConPlot(genes = plot_gene_lst,
                  sv_data_instability_only,
                  cn_data %>% filter(chr %in% plot_region_instability_only$seqnames),
                  chr_selection=plot_region_instability_only,
                  legend_SV_types=T,
                  pos_SVtype_description=1000000,
                  scale_separation_SV_type_labels=1/23,
                  title=paste0(multichrom_group," ",target_sample," - plot chrom instability only"),
                  max.cn = 12,
                  label_interchr_SV = TRUE,
    )
    print(p)
  }
}

dev.off()

}

