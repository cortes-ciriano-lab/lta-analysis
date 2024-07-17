# LTA detection

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

#gliobblastoma lung breast => disruption of TSGs connected to oncogene amps. esoph 

suppressPackageStartupMessages({
  library(GenomicRanges, quietly=TRUE)
  library(rtracklayer, quietly=TRUE)
  library(tidyverse, quietly=TRUE)
  library(stringr, quietly=TRUE)
  library(stringi)
  library(dplyr)
})

# Settings ----
chrom_order=c(paste0("chr",1:22),"chrX")
options(expressions= 500000)
max_svs_display=1000

if(!exists("chromothripsis_clusters_path")) {
  #local
  resources_dir = "/Users/belzen/resources/"
  metadata_dir="/Users/belzen/data/metadata/"
  results_dir="/Users/belzen//results/LTA-analysis/"
  
  # resources_dir = "/nfs/research/icortes/belzen/resources/"
  # results_dir="/nfs/research/icortes/belzen/results/lta_detection/"
  # metadata_dir="/nfs/research/icortes/belzen/data/metadata/"
  # 
  if(FALSE){
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
  
  chromothripsis_clusters_path = paste0(filesystem,"hmf/data_freeze_20240512.chromothripsis_clusters.tsv")
  cohort_path="~/data/metadata/TCGA.data_freeze_20240512.purity_table.tsv"
  
  # chromothripsis_clusters_path = "/nfs/research/icortes/DATA/TCGA_WGS/hmf/data_freeze_20240512/shatterseek/chromothripsis_clusters.tsv"
  # cohort_path="/nfs/research/icortes/DATA/TCGA_WGS/hmf/data_freeze_20240512/purity_table.tsv"
  # 
  }

if(FALSE) {
  #osteos
  wgd_score_dir="~/results/osteos_gel-100k/"
  wgd_score_dir="~/results/osteos_all/"
  driver_files_dir="/Users/belzen/data/ebi/osteosarcoma_analysis/driver_files/"
  driver_germline_files_dir=paste0("~/data/ebi/osteosarcoma_analysis/driver_germline_files/")
  
  cn_files_dir="/Users/belzen/data/cn_files/"
  chromothripsis_clusters_path="/Users/belzen/data/ebi/osteosarcoma_analysis/shatterseek/chromothripsis_clusters.tsv"
  
  dataset_selection_label="osteos"
  
  sample_table_path=paste0(metadata_dir,"OS_master_clinical_data_v2.csv")
  driver_masterfile_path = paste0(metadata_dir,"PerDonor_master_withBurdens_mainDrivers_28_Oct2023.csv")
  suptable1_path=paste0(metadata_dir,"OS_master_table_s1.tsv")
  
  #os driver files
  #/nfs/research/icortes/DATA/osteosarcoma/analysis/output_files/driver_files
  
  #shatterseek output
  #/nfs/research/icortes/DATA/osteosarcoma/analysis/shatterseek
  
  #svs from 
  #/nfs/research/icortes/DATA/osteosarcoma/analysis/output_files/sv_filtered_files
  
  #instead of 05 could be run on 
  #/nfs/research/icortes/DATA/TCGA_WGS/hmf/data_freeze_20240617/purity_table.tsv
  #todo proper sbatch launcer for  sbatch --mem=50Gb --time=35:00:00 run_wgd_score.sh 
  
  
}
}

if(!exists("flag_run_single_sample")) { flag_run_single_sample=F }

if(!exists("dataset_selection_label")) {
  print("Warning, no dataset_selection_label, setting default")
  dataset_selection_label="TCGA-SARC"
}

print(paste0("Running: ",dataset_selection_label))


# Paths ----

# templates 
chromosome_bands_path = paste0(resources_dir,"chromosome_bands.gz")
gtf_path= paste0(resources_dir,"gencode.v38.annotation.gtf.gz")
#gtf_path = stri_replace_all_fixed(gtf_path_template,names(map_template_vars), map_template_vars,vectorize=F)

#input

sv_wgd_score_path_template = paste0("${wgd_score_dir}/${patient_basename}.sv.filtered.wgd_score.tsv")
cn_path_template = "${cn_files_dir}/${patient_basename}.purple.cnv.somatic.tsv"
drivers_path_template = paste0("${driver_files_dir}/${patient_basename}.purple.driver.catalog.tsv")
drivers_germline_path_template = paste0("${driver_germline_files_dir}/${patient_basename}.purple.driver.catalog.germline.tsv")

#output
gene_cn_sv_disruptions_path_template = paste0("${processed_output_dir}/${patient_basename}.gene_cn_sv_disruptions.tsv")
clusters_instability_region_path_template= paste0("${processed_output_dir}/${patient_basename}.clusters_instability_region.tsv")
svs_connecting_tsg_instability_path_template= paste0("${processed_output_dir}/${patient_basename}.svs_connecting_tsg_instability.tsv")
connected_regions_tsg_instability_path_template= paste0("${processed_output_dir}/${patient_basename}.connected_regions_tsg_instability.tsv")
cohort_call_lta_path_template= paste0("${processed_output_dir}/${patient_basename}.lta_detection_overview.tsv")

reconplot_path_template = paste0("${plot_dir}/${patient_basename}.${cgr_id}.reconplot.pdf")

#path names
plot_dir=paste0(results_dir,"plots/")

map_template_vars=c('${cn_files_dir}'=cn_files_dir,
                    '${wgd_score_dir}'=wgd_score_dir,
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

# Functions ----
get_chr_arms = function(chromosome_bands_df,split_giestain=T,autosomes= c(paste("chr",1:22,sep=""),"chrX","chrY")) {
  ## note; 1 based coordinates??? 
  
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
chr_arms = get_chr_arms(chromosome_bands_df,split_giestain = F)
chr_arms = GRanges(chr_arms)
names(chr_arms) = chr_arms$chr_arm
chr_arms_df = as.data.frame(chr_arms)
chr_arms_df = chr_arms_df %>% dplyr::mutate(chr_arm_start=start,chr_arm_end=end) 

gtf <- rtracklayer::import(gtf_path)
gene_properties = gtf[gtf$type=="gene"]
names(gene_properties) = gene_properties$gene_id
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

cohort = cohort %>% filter(purity>=0.3 & qc=="PASS" )

if(flag_run_single_sample==F & exists("dataset_selection_label")) {
  cohort = cohort %>% filter(cohort_id==dataset_selection_label)
  print(paste0("RUNNING: ",dataset_selection_label))
}
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

sv_wgd_score = data.frame()
drivers = data.frame()
drivers_germline = data.frame()
cohort_segments_df = data.frame()

for(target_sample in cohort$basename) {
  map_template_vars_patient=c(map_template_vars,'${patient_basename}'=target_sample)
  
  sv_wgd_score_path = stri_replace_all_fixed(sv_wgd_score_path_template,names(map_template_vars_patient), map_template_vars_patient,vectorize=F)

  if(length(Sys.glob(sv_wgd_score_path))==1){
    sv_wgd_score_sample = read.table(sv_wgd_score_path,sep = "\t",header = T)
    sv_wgd_score_sample$basename=target_sample
    sv_wgd_score = rbind(sv_wgd_score,sv_wgd_score_sample)  
    
  } else {
    print(target_sample)
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

cohort %>% filter(basename %in% cohort_segments_df$basename)

svs_df = sv_wgd_score
svs_df = svs_df %>% mutate(sv_id=paste0(basename,"_",bp_name),partner_sv_id=paste0(basename,"_",partner))
svs_gr=GRanges(svs_df)
names(svs_gr)=svs_df$sv_id

cohort_segments_df = cohort_segments_df %>% mutate(cn_seg_id=paste0(basename,"_",segment_id))
cn_seg_gr = GRanges(cohort_segments_df)
names(cn_seg_gr) = cohort_segments_df$cn_seg_id

drivers = rbind(drivers,drivers_germline) %>% unique()
drivers = drivers %>% dplyr::mutate(gene_name=gene) %>% left_join(gene_properties_df %>% select(gene_name,gene_id,gene_start,gene_end)) 
drivers = drivers %>% dplyr::mutate(chr_arm = ifelse(grepl("p",chromosomeBand),paste0(chromosome,"p"),paste0(chromosome,"q")))
drivers = drivers %>% left_join(chr_arms_df %>% select(chr_arm,chr_arm_start,chr_arm_end))
drivers = drivers %>% mutate(sample_gene = paste0(basename,"_",gene_name))

#collapse multiple
driver_summary = drivers %>% group_by(basename,gene_name) %>% 
  summarize(hmf_biallelic=any(biallelic=="true"),likelihoodMethod=toString(likelihoodMethod),hmf_driver=toString(driver)) %>% ungroup() %>% as.data.frame()


    
# TSG disruption ----

# biallelic disruption
# + any SV evidence
# overlapping SV or CGR
# LOH (minor CN <0.5)  and downstream SV 
## => any downstream SV or  more specific: CTX or CGR

#todo replace by driver file TSGs
lta_tsg_disruption_lst = c("TP53")

#overlap gene with CN for LOH and homz loss
gene_cn_overlap = get_reciprocal_overlap_pairs(gene_properties[gene_properties$gene_name %in% lta_tsg_disruption_lst],cn_seg_gr,reciprocal_overlap = 0,svtype_matching = F,ignore_strand = T)
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

#overlap gene with cgr => not used
# gene_cgr_overlap = get_reciprocal_overlap_pairs(gene_properties[gene_properties$gene_name %in% lta_tsg_disruption_lst],chromothripsis_clusters_gr,reciprocal_overlap = 0,svtype_matching = F,ignore_strand = T)
# gene_cgr_overlap = gene_cgr_overlap  %>% dplyr::rename(gene_id=set1,cluster_id=set2)  %>% select(-to,-from)
# gene_cgr_overlap = gene_cgr_overlap %>% left_join(chromothripsis_clusters %>% select(cluster_id,basename)) 

#overlap gene with sv
gene_sv_bp_overlap = get_reciprocal_overlap_pairs_start_end(svs_gr,gene_properties[gene_properties$gene_name %in% lta_tsg_disruption_lst],reciprocal_overlap = 0,svtype_matching = F,ignore_strand = T)
if(gene_sv_bp_overlap %>% nrow() > 0 ) {
  
gene_sv_bp_overlap = gene_sv_bp_overlap  %>% dplyr::rename(sv_id=set1,gene_id=set2)  %>% select(-to,-from)
gene_sv_bp_overlap = gene_sv_bp_overlap %>% left_join(svs_df %>% select(sv_id,basename,purple_af_start,purple_af_end,PURPLE_JCN)) %>% mutate(tumor_af = ifelse(sv_breakpoint_orientation=="start",purple_af_start,purple_af_end))
gene_sv_bp_overlap = gene_sv_bp_overlap %>% left_join(gene_properties_df[,c("gene_id","gene_name")] %>% unique()) %>% mutate(sample_gene = paste0(basename,"_",gene_name))

gene_sv_bp_max_af = gene_sv_bp_overlap %>% 
  group_by(sample_gene,basename,gene_name) %>% 
  summarize(max_tumor_af = max(tumor_af,na.rm = T),
            max_multiplicity = max(PURPLE_JCN,na.rm = T))
} 
#check for arm level disruptions
#for P arm: gene overlapping or downsteam so gene start to chr arm end; for Q arm chr arm start to gene end
#chr arm from driver file or mappiing if want to check hom_disruptions

map_genes_chr_arms = get_reciprocal_overlap_pairs(gene_properties[gene_properties$gene_name %in% lta_tsg_disruption_lst],chr_arms,reciprocal_overlap = 0,svtype_matching = F,ignore_strand = T)
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


gene_region_sv_bp_overlap = get_reciprocal_overlap_pairs_start_end(svs_gr,region_gene_sv_disruption_chr_arm_gr,reciprocal_overlap = 0,svtype_matching = F,ignore_strand = T)
gene_region_sv_bp_overlap = gene_region_sv_bp_overlap %>% dplyr::rename(sv_id=set1,gene_id=set2,svtype=set1_svtype)  %>% select(-to,-from)
gene_region_sv_bp_overlap = gene_region_sv_bp_overlap %>% left_join(svs_df %>% select(sv_id,basename,purple_af_start,purple_af_end,PURPLE_JCN)) %>% mutate(tumor_af = ifelse(sv_breakpoint_orientation=="start",purple_af_start,purple_af_end))
gene_region_sv_bp_overlap = gene_region_sv_bp_overlap %>% left_join(gene_properties_df[,c("gene_id","gene_name")] %>% unique())  %>% mutate(sample_gene = paste0(basename,"_",gene_name))

#check if gene is knocked out by sv
multiplicity_gene_region_sv.min=0.9
biallelic_sv_bp_tumor_af.min=0.9
sv_bp_tumor_af.min=0.5

loh_minorCN.max=0.5

#TODO: if not empty ...
if(exists("gene_sv_bp_max_af")) {
  gene_cn_sv_disruptions = gene_cn %>% ungroup() %>%
    left_join(gene_sv_bp_max_af)
} else {
  #otherwise add empty cols
  gene_cn_sv_disruptions = gene_cn %>% ungroup() 
  gene_cn_sv_disruptions$max_tumor_af=NA
  gene_cn_sv_disruptions$max_multiplicity=NA
  
}

gene_cn_sv_disruptions = gene_cn_sv_disruptions %>%
  left_join(driver_summary) %>%
   dplyr::mutate(
     hmf_biallelic =ifelse(is.na(hmf_biallelic),F,T),
  gene_sv_bp= ifelse(!is.na(max_tumor_af), (max_tumor_af >= sv_bp_tumor_af.min | (max_multiplicity/weighted_CN) >= sv_bp_tumor_af.min), F),
  gene_loh=weighted_minorCN<=loh_minorCN.max,
  gene_sv_bp_biallelic = ifelse(!is.na(max_tumor_af), gene_loh & (max_tumor_af >= biallelic_sv_bp_tumor_af.min | (max_multiplicity/weighted_CN) >= biallelic_sv_bp_tumor_af.min), F),
  gene_region_cgr = sample_gene %in% gene_region_cgr_overlap$sample_gene,
 # gene_region_sv_bp = sample_gene %in% filter(gene_region_sv_bp_overlap,PURPLE_JCN>=multiplicity_gene_region_sv.min)$sample_gene,
  gene_region_ctx = sample_gene %in% filter(gene_region_sv_bp_overlap,svtype=="CTX"&PURPLE_JCN>=multiplicity_gene_region_sv.min)$sample_gene) %>%
  dplyr::mutate(
    gene_knockout_sv_component = gene_loh & (hmf_biallelic | (!is.na(hmf_driver) & gene_sv_bp) | gene_sv_bp_biallelic),
    #gene_knockout_sv_component_strict = (gene_sv_bp_biallelic | (hmf_biallelic & (gene_loh & (gene_region_cgr | gene_region_ctx))))
    gene_knockout_sv_component_strict = gene_knockout_sv_component & (gene_region_cgr | gene_region_ctx | gene_sv_bp)
  ) %>%
  left_join(map_genes_chr_arms) #for coordinates
    

if(FALSE) {
  #gene_cn_sv_disruptions_bk=gene_cn_sv_disruptions
  missed_disruptions = gene_cn_sv_disruptions %>% filter(basename %in% 
                                                           filter(assess_lta,LTA=="Yes" &  lta==F & has_tsg_disruption==F)$basename)
  
  missed_disruptions %>% as.data.frame()
  
  
  
  missed_disruptions %>% filter(!basename %in% gene_region_sv_bp_overlap$basename)
  missed_disruptions %>% filter(!basename %in% filter(gene_region_sv_bp_overlap,svtype=="CTX")$basename)
  
  gene_region_sv_bp_overlap %>% filter(svtype=="CTX" & basename %in% missed_disruptions$basename)
  
  missed_disruptions
  gene_cn_sv_disruptions %>% filter(basename %in% missed_disruptions$basename) %>% filter(gene_knockout_sv_component==T)
  gene_cn_sv_disruptions %>% filter(basename %in% missed_disruptions$basename) %>% filter(gene_knockout_sv_component==F)
  
  missed_disruptions = gene_cn_sv_disruptions %>% filter(gene_knockout_sv_component==F) %>% merge(cohort[,c("basename","LTA","final_id","TP53status","TP53")]) %>% 
    filter(LTA=="Yes")
  
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

#from drivers, get amplified oncogenes
# oncogene to chr arm
# cgr in chrom arm
# for region take chrom arm  OR CGR + margin


instability_region_oncogene_amp_long = instability_region %>% merge(drivers %>% filter(driver=="AMP")  %>% select(basename,gene_name,maxCopyNumber,chr_arm),
                                                              by=c("basename","chr_arm"))
instability_region_oncogene_amp = instability_region_oncogene_amp_long %>% group_by(basename,chr_arm,cluster_id) %>% 
  summarize(onco_gene_amp_lst = toString(unique(paste0(gene_name," (",round(maxCopyNumber,2),") "))))


instability_region = instability_region %>% left_join(instability_region_oncogene_amp) %>% mutate(onco_amp_in_chr_arm = !is.na(onco_gene_amp_lst))


# Connections: Link tsg region to onco region ----
#get disrupted tsgs, look up chromosome arm of genem append chrom arm location and make that into range 
#for question of disruptive SV can use for p arm: gene start to chrom end (towards centomere), and for q arm chrom start to gene end

tsg_region = gene_cn_sv_disruptions %>% filter(gene_knockout_sv_component) %>% as.data.frame()

tsg_region$region_id = paste0(tsg_region$basename,"_",tsg_region$gene_name)
tsg_region$region_type="tsg_region"

coord_start_col="chr_arm_start"
coord_end_col="chr_arm_end"
tsg_region_id_col="region_id"
tsg_region_gr = tsg_region %>% dplyr::rename(start=!!sym(coord_start_col),end=!!sym(coord_end_col)) %>% GRanges()
names(tsg_region_gr) = tsg_region[,tsg_region_id_col]


instability_region$region_type="instability_region"
instability_region$region_id=paste0(instability_region$cluster_id,"_",instability_region$chr_arm)

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

#multi to multi because of cgrs spanning chroms > region id with 
#instability_region %>% filter(region_id=="cluster_gel-100k_215000258-tumour-1_LP3000758-DNA_F04_chr9")

#todo refactor to carry less on the svs annotated df 
#nb: chromosome arms etc are from perspective of instability region
map_svs_regions = map_svs_regions %>% left_join(instability_region %>% select(basename,region_id,chr_arm))
                                                 #,multichromosomal,onco_amp_in_chr_arm,onco_gene_amp_lst))
# map_svs_regions[is.na(map_svs_regions$onco_amp_in_chr_arm),]$onco_amp_in_chr_arm=F
# map_svs_regions[is.na(map_svs_regions$multichromosomal),]$multichromosomal=F

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
  left_join(instability_region %>% select(region_id,chrom_all,classification,cluster_id,multichromosomal,onco_amp_in_chr_arm,onco_gene_amp_lst) ) %>%
  dplyr::rename_with(.fn=function(x){paste0("cgr_",x)},.cols=c("multichromosomal","chrom_all","classification")) %>%
  left_join(tsg_region %>% 
              select(region_id,gene_name,weighted_CN,weighted_minorCN,max_tumor_af,max_multiplicity,contains("hmf"),
                     gene_loh,gene_sv_bp_biallelic,gene_region_cgr,gene_region_ctx,gene_knockout_sv_component,gene_knockout_sv_component_strict) %>% unique() %>%
              dplyr::rename_with(.fn=function(x){paste0("tsg_",x)}),
            by=c("tsg_region_id"))  %>%
  as.data.frame()


if(FALSE) {
  #todocheck multichromosomal if they all show up in the connection table
  #check attributes
  selected_connected_tsg_regions %>% filter(basename=="TCGA-SARC_TCGA-PC-A5DM_TCGA-PC-A5DM-01A-11D-A813-36")
  svs_connecting_tsg_instability %>% filter(basename=="TCGA-SARC_TCGA-PC-A5DM_TCGA-PC-A5DM-01A-11D-A813-36")
  
}

# Call LTA ----
#multichromosomal if CGR is multichromosomal but also if tsg connected to two different clusters

count_connected_tsg_regions = selected_connected_tsg_regions %>% group_by(basename,tsg_region_id) %>% 
  summarize(instability_region_cnt = length(unique(cluster_id)))

cohort_call_lta = cohort  %>% 
  select(cohort_id,basename,contains("final_id"),contains("TP53"),contains("LTA"),contains("subtype_short")) %>%
  mutate(has_tsg_disruption = basename %in% tsg_region$basename,
         has_tsg_disruption_strict = basename %in% filter(tsg_region,gene_knockout_sv_component_strict)$basename,
         has_instability_region = basename %in% instability_region$basename,
         lta = basename %in% selected_connected_tsg_regions$basename,
         lta_multichrom = basename %in% c(filter(selected_connected_tsg_regions,cgr_multichromosomal)$basename,filter(count_connected_tsg_regions,instability_region_cnt>1)$basename),
         lta_onco = basename %in% filter(selected_connected_tsg_regions,onco_amp_in_chr_arm)$basename,
         lta_onco_multichrom = basename %in% filter(selected_connected_tsg_regions,onco_amp_in_chr_arm&cgr_multichromosomal)$basename,
  )

#Export -----

write.table(gene_cn_sv_disruptions %>% left_join(cohort_call_lta),gene_cn_sv_disruptions_path,sep = "\t",row.names = F,col.names = T)
write.table(instability_region %>% left_join(cohort_call_lta),clusters_instability_region_path,sep = "\t",row.names = F,col.names = T)
write.table(svs_connecting_tsg_instability %>% left_join(cohort_call_lta),svs_connecting_tsg_instability_path,sep = "\t",row.names = F,col.names = T)
write.table(selected_connected_tsg_regions %>% left_join(cohort_call_lta),connected_regions_tsg_instability_path,sep = "\t",row.names = F,col.names = T)
write.table(cohort_call_lta,cohort_call_lta_path,sep = "\t",row.names = F,col.names = T)

# Reconplots ----
library(ReConPlot)

# target_sample="TCGA-SARC_TCGA-DX-A8BK_TCGA-DX-A8BK-01A-11D-A811-36"
# target_sample = "TCGA-SARC_TCGA-DX-A6YR_TCGA-DX-A6YR-01A-33D-A811-36"
# target_sample="TCGA-SARC_TCGA-DX-AB2Q_TCGA-DX-AB2Q-01A-11D-A812-36"
# target_sample="TCGA-SARC_TCGA-PC-A5DM_TCGA-PC-A5DM-01A-11D-A813-36"
#for(target_sample in unique(selected_connected_tsg_regions$basename))



#plot_path_label=paste0("lta_multichrom")
#for(target_sample in filter(cohort_call_lta,LTA=="No" & lta_multichrom)$basename) {
  #how to viz missed??
for(target_sample in filter(cohort_call_lta,lta)$basename) {
  sample=filter(cohort_call_lta,basename==target_sample)
  plot_path_label=paste0("lta")
  if("LTA" %in% names(sample)) {
    plot_path_label=paste0(plot_path_label,".truth_set=",sample$LTA)
  }
  plot_path_label=paste0(plot_path_label,".detected=",sample$lta)
  if(sample$lta) {
    plot_path_label=paste0(plot_path_label,".onco_multichrom=",sample$lta_onco_multichrom,".multichrom=",sample$lta_multichrom,".onco=",sample$lta_onco)
  }
  
  map_template_vars_patient=c(map_template_vars,'${patient_basename}'=target_sample)
  
#chr arm region 
#include connections found
#all svs of cgrs of the multichrom constructs 

selected_instability_region = instability_region %>% filter(region_id %in% filter(selected_connected_tsg_regions,basename==target_sample)$region_id)
selected_connecting_svs = svs_connecting_tsg_instability %>% filter(region_id %in% selected_instability_region$region_id) %>%
  left_join(selected_instability_region %>% select(region_id,chrom_all))

map_template_vars_cgr = c(map_template_vars_patient,'${cgr_id}'=plot_path_label)
reconplot_path = stri_replace_all_fixed(reconplot_path_template,names(map_template_vars_cgr), map_template_vars_cgr,vectorize=F)

plot_chrom = str_count(selected_instability_region$chrom_all,"_") %>% max() + 1
plot_width = ifelse(plot_chrom<4, 35, 35+(2.5*(plot_chrom-3)))

if(length(Sys.glob(reconplot_path))==1) {next()}

pdf(reconplot_path, width = plot_width/2.2, height = 8/2.2,pointsize = 8)


#plot the entire thing and plot per multichrom group (chrom_all)
multichrom_group_lst = c("all",selected_instability_region$chrom_all %>% unique())
multichrom_group="chr3_chr5_chr7_chr8_chr9_chr11_chr13_chr19_chrX"
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

plot_gene_lst=c("TP53")
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

  if(nrow(sv_data) > max_svs_display) {
    print(paste0("Too many SVs to display, skipping plot", nrow(sv_data)))
    next()
  }
  
  cn_data = cohort_segments_df %>% 
    dplyr::rename(chr=seqnames) %>% 
    filter(basename==target_sample & chr %in% sv_data$chr1) %>%
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
  
  plot_region_instability_only = c(plot_selected_connecting_svs$chr_arm,plot_selected_connecting_svs$tsg_chr_arm,plot_selected_instability_region$chr_arm) %>% unique()
  
  plot_region_instability_only = chr_arms_df %>%
    filter(chr_arm %in% plot_region_instability_only) %>% arrange(seqnames) %>%
    select(seqnames,start,end) %>% dplyr::mutate(chr1=as.character(seqnames),chr=chr1,start=1,end=500e6) %>% as.data.frame() %>% unique()
  
  if(full_chrom %>% filter(!chr %in% plot_region_instability_only$chr) %>% nrow() == 0 ) {next()}
  p = ReConPlot(genes = plot_gene_lst,
                sv_data %>% filter(chr1 %in% plot_region_instability_only$seqnames),
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

dev.off()

}

