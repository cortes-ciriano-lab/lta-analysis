#analyse results
#load config osteos
suppressPackageStartupMessages({
  library(GenomicRanges, quietly=TRUE)
  library(rtracklayer, quietly=TRUE)
  library(tidyverse, quietly=TRUE)
  library(stringr, quietly=TRUE)
  library(stringi)
  library(dplyr)
})
# Functions ----

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

# paths ----
if(!exists("results_dir")) {
  metadata_dir="/Users/belzen/data/metadata/"
  results_dir="/Users/belzen//results/LTA-analysis/multi_tsg/"
  sample_table_path=paste0(metadata_dir,"OS_master_clinical_data_v2.csv")
  driver_masterfile_path = paste0(metadata_dir,"PerDonor_master_withBurdens_mainDrivers_28_Oct2023.csv")
  suptable1_path=paste0(metadata_dir,"OS_master_table_s1.tsv")
  
  
}

#templates
gene_cn_sv_disruptions_path_template = paste0("${processed_output_dir}/${patient_basename}.gene_cn_sv_disruptions.tsv")
clusters_instability_region_path_template= paste0("${processed_output_dir}/${patient_basename}.clusters_instability_region.tsv")
svs_connecting_tsg_instability_path_template= paste0("${processed_output_dir}/${patient_basename}.svs_connecting_tsg_instability.tsv")
connected_regions_tsg_instability_path_template= paste0("${processed_output_dir}/${patient_basename}.connected_regions_tsg_instability.tsv")
cohort_call_lta_path_template= paste0("${processed_output_dir}/${patient_basename}.lta_detection_overview.tsv")
call_lta_annot_path_template = paste0("${processed_output_dir}/${patient_basename}.lta_detection_overview.all_tsgs.tsv")
oncogene_amp_annot_path_template = paste0("${processed_output_dir}/${patient_basename}.oncogene_amp.annot.tsv")

#path names
plot_dir=paste0(results_dir,"plots/")

map_template_vars=c(
                   "${plot_dir}"=plot_dir,
                    "${processed_output_dir}"=results_dir)


#Cohort ----
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
  

# Compare to OS manual ----

dataset_selection_label = "osteos"
map_template_vars_cohort = c(map_template_vars,"${patient_basename}"=dataset_selection_label)
cohort_call_lta_path= stri_replace_all_fixed(cohort_call_lta_path_template,names(map_template_vars_cohort), map_template_vars_cohort,vectorize=F)

cohort_call_lta=read.table(cohort_call_lta_path,sep = "\t",header = T)

#agreement?
#subset to HGOS 
assess_lta = cohort_call_lta  %>% filter(subtype_short=="HGOS")
assess_lta = assess_lta %>% left_join(driver_masterfile %>% select(basename,contains("OncoAmp"))) 

assess_lta %>% nrow()
assess_lta %>% filter(LTA=="Yes") %>% nrow()

assess_lta %>% filter(LTA=="Yes" ) %>% 
  filter(lta & tsg_gene_knockout_dicentric) %>% nrow()



#oncogene amp
#aim 45/68 but can only get either 47 with any amp, or amp+tra, or 44/42 CGR assoc
if(FALSE) {
  assess_lta %>% filter(LTA=="Yes" &OncoAmp_status=="OncoAmp") %>% nrow()
  assess_lta %>% filter(LTA=="Yes" &OncoAmp_TRA_status=="OncoAmp+TRA") %>% nrow()
  
  #CGR in onco amp arm
  assess_lta %>% filter(LTA=="Yes" & OncoAmp_allCGR_status=="OncoAmp+CGR") %>% nrow()
  
  #entry point for my analysis is "connected to CGR" 
  assess_lta %>% filter(LTA=="Yes" & OncoAmp_CGR_w_TRA_status=="OncoAmp+CGRassoc_TRA") %>% nrow()
  
  assess_lta %>% filter(LTA=="Yes" & OncoAmp_CGR_w_TRA_status=="OncoAmp+CGRassoc_TRA") %>% 
    filter(lta_any) %>% nrow()
  
  assess_lta %>% filter(LTA=="Yes" & OncoAmp_CGR_w_TRA_status=="OncoAmp+CGRassoc_TRA") %>% 
    filter(lta) %>% nrow() #these should be linked
  
  assess_lta %>% filter(LTA=="Yes" & OncoAmp_CGR_w_TRA_status=="OncoAmp+CGRassoc_TRA") %>% 
    filter(lta & !lta_onco_incl_indirect) #missing 2 with indirect approach, 5 with direct connections 
  
  #onco amps not associated to CGR in my analysis
  #no shatterseek CGR there
  oncogene_amp_annot %>% filter(basename %in% filter(assess_lta,LTA=="Yes" & OncoAmp_CGR_w_TRA_status=="OncoAmp+CGRassoc_TRA" & lta & !lta_onco_incl_indirect)$basename)  %>% 
    select(basename,gene_name,chr_arm,region_id)
  chromothripsis_clusters %>% filter(basename %in% filter(assess_lta,LTA=="Yes" & OncoAmp_CGR_w_TRA_status=="OncoAmp+CGRassoc_TRA" & lta & !lta_onco_incl_indirect)$basename)  %>% 
    select(basename,chrom_all) %>% unique()
  
  
  assess_lta %>% filter(LTA=="Yes" & OncoAmp_allCGR_status=="OncoAmp+CGR") %>% 
    filter(lta & !lta_onco_incl_indirect) #missing 3 with indirect approach, 5 with direct connections 
  
  oncogene_amp_annot %>% filter(basename %in% filter(assess_lta,LTA=="Yes" & OncoAmp_allCGR_status=="OncoAmp+CGR" & lta & !lta_onco_incl_indirect)$basename)  %>% 
    select(basename,gene_name,chr_arm,region_id)
  chromothripsis_clusters %>% filter(basename %in% filter(assess_lta,LTA=="Yes" & OncoAmp_allCGR_status=="OncoAmp+CGR" & lta & !lta_onco_incl_indirect)$basename)  %>% 
    select(basename,chrom_all) %>% unique()
  #6              gel-100k_215001755-tumour-6_LP3001501-DNA_E02      EGFR   chr7p cluster_gel-100k_215001755-tumour-6_LP3001501-DNA_E02_chr7_chr7p
  #is linked to cgr, but also in Solanges table no CTX found
  
  
  assess_lta %>% filter(LTA=="Yes" & OncoAmp_CGR_w_TRA_status=="OncoAmp+CGRassoc_TRA") %>% 
    filter(lta & !lta_onco_direct) #missing 5 with direct connections 
  
}


assess_lta %>% filter(LTA=="Yes" & OncoAmp_CGR_w_TRA_status=="OncoAmp+CGRassoc_TRA") %>% 
  filter(lta & lta_onco_direct&tsg_gene_knockout_dicentric) %>% nrow()


assess_lta %>% filter(LTA=="Yes") %>% 
  filter(lta & lta_onco_direct&tsg_gene_knockout_dicentric) %>% nrow()

#keep the 2/3rds of the 49 high confidence 29 have onco amplified directly on partner 


#found any LTA, including samples without sufficient evidence tp53 knockout by sv
assess_lta %>% filter(LTA=="Yes"  & lta_any==T ) %>% nrow()

assess_lta %>% filter(LTA=="Yes"  & lta==T ) %>% nrow() 

assess_lta %>% filter(LTA=="Yes"  & lta==T & tsg_gene_knockout_dicentric) %>% nrow()
assess_lta %>% filter(LTA=="Yes"  & lta==T & tsg_gene_knockout_dicentric) %>% nrow()

assess_lta %>% filter(LTA=="Yes"  & lta==T & lta_onco_direct & tsg_gene_knockout_dicentric) %>% nrow()
assess_lta %>% filter(LTA=="Yes"  & lta==T & !lta_onco_direct  & lta_onco_incl_indirect & tsg_gene_knockout_dicentric) %>% nrow()

assess_lta %>% filter(LTA=="Yes"  & lta==T & lta_multichrom & tsg_gene_knockout_dicentric) %>% nrow()

assess_lta %>% filter(LTA=="Yes"  & lta==T & (lta_onco_multichrom) ) %>% nrow()
assess_lta %>% filter(LTA=="Yes"  & lta==T & (lta_onco_multichrom) & tsg_gene_knockout_dicentric) %>% nrow()

# would be filtered out of TCGA analysis because not oncogene
assess_lta %>% filter(LTA=="Yes" & lta==T & !(lta_onco_direct &tsg_gene_knockout_dicentric)) %>% nrow()

#new cases or false positives also outside of HGOS?
assess_lta %>% filter(LTA=="No" & lta & (tsg_gene_knockout_dicentric)) %>% select(subtype_short,basename,final_id,lta_onco,tsg_gene_knockout_dicentric)
cohort_call_lta %>% filter(LTA=="No" & lta & (tsg_gene_knockout_dicentric)) %>% select(subtype_short,basename,final_id,lta_onco,tsg_gene_knockout_dicentric)

#missed cases
assess_lta %>% filter(LTA=="Yes" &  lta==F)  %>% select(basename,final_id,lta_any,tsg_gene_knockout)


if(FALSE) {
  #additional true positive
  gene_cn_sv_disruptions %>% filter(gene_name=="TP53" & basename %in% 
                                      filter(assess_lta,LTA=="Yes"  & lta==T & !tsg_gene_knockout_dicentric)$basename) %>% 
    filter(gene_knockout_dicentric2 &gene_knockout_dicentric==F) %>% as.data.frame()
  
 
  #missing true positive
  gene_cn_sv_disruptions %>% filter(gene_name=="TP53" & basename %in% 
                                      filter(assess_lta,LTA=="Yes"  & lta==T & tsg_gene_knockout_dicentric)$basename) %>% 
    filter(gene_knockout_dicentric2==F) %>% as.data.frame()
  
  
  #additional false positive
  gene_cn_sv_disruptions %>% filter(gene_name=="TP53" & basename %in% 
                                      filter(assess_lta,LTA=="No"  & lta==T & !tsg_gene_knockout_dicentric)$basename) %>% 
    filter(gene_knockout_dicentric2 &gene_knockout_dicentric==F) %>% as.data.frame()
 
 
  #checkfor FP with other genes 
  gene_cn_sv_disruptions.call_lta %>% filter(tsg_gene_name=="CDKN2A" & basename %in% 
                                      filter(assess_lta,LTA=="No")$basename) %>% 
    filter(lta & tsg_gene_knockout_dicentric2) %>% as.data.frame()
  #think this one is acceptable mdacc_MDACC36_MDA-SMF-12-33_RB1 also gel-100k_215000645-tumour-1_LP3000544-DNA_G10_CDKN2A

  verified_misses = c("mdacc_MDACC41_MDA-SMF-12-9",
                      "gel-100k_215001602-tumour-1_LP3001058-DNA_A05","gel-100k_215001679-tumour-2_LP3000850-DNA_F01",
                      "gel-100k_215003146-tumour-1_LP3001502-DNA_C05","gel-100k_215003563-tumour-1_LP3001609-DNA_D04",                    
                      "kidsfirst_PT_870WE3XB-tumour-1_H_VC-UC0027-GMKF-40-UC0027-01A-01E","kidsfirst_PT_X8TWRDEH-tumour-1_H_VC-GMKF-40-PASTRH-08A-01D",
                      "pcawg_DO52642_f83fc777-5416-c3e9-e040-11ac0d482c8e")

  gene_cn_sv_disruptions.call_lta %>% 
    filter(tsg_gene_name=="TP53") %>%
    filter(basename %in% filter(cohort_call_lta,LTA=="Yes"&subtype_short=="HGOS")$basename) %>% 
    filter(!basename %in% verified_misses) %>%
    filter(lta) %>%
    # filter(lta_onco_multichrom)  %>% 
    # filter(lta_onco_amp)  %>% 
     mutate(tsg_seg_amp9_kbp=ifelse(is.na(tsg_seg_amp9_kbp),0,tsg_seg_amp9_kbp)) %>%
    
    filter(tsg_gene_knockout_dicentric==F & tsg_gene_knockout_dicentric2) %>%
    #filter( ! (((tsg_gene_sv_bp_biallelic) | (tsg_gene_knockout & tsg_gene_loh & (tsg_gene_region_cgr | tsg_gene_region_ctx_notaf | tsg_gene_sv_bp))) & 
    #            ((tsg_chr_arm_terminal_frac_loh >= 0.75 | tsg_seg_amp9_kbp>1e3) & tsg_chr_arm_frac_loh <=0.9 ) & tsg_chrom_frac_loh<=0.9)) %>% #dicentric chromosome is likely
    
    #filter( !(tsg_gene_region_cgr | tsg_gene_region_ctx_notaf | tsg_gene_sv_bp)) %>%
    #filter( ! (((tsg_chr_arm_terminal_frac_loh >= 0.75 | tsg_seg_amp9_kbp>1e3))) ) %>%
    as.data.frame() %>% select(basename,tsg_chr_arm_frac_loh,tsg_chrom_frac_loh)
    arrange(-tsg_gene_knockout_sv_component_strict)
    group_by(tsg_gene_name,connected_to_tp53_tsg_region) %>% summarize(n=length(unique(basename)),samples=toString(unique(basename))) %>% arrange(-n) %>% as.data.frame()
  #group_by(cohort_id) %>% summarize(n=length(unique(basename)),samples=toString(unique(basename)),tsgs=toString(unique(tsg_gene_name))) %>% arrange(-n)

    #good  LSAMP still false
  gene_cn_sv_disruptions %>% filter(basename=="kidsfirst_PT_JMJ6ZHPM-tumour-1_H_VC-UC0029-GMKF-40-UC0029-01A-01E"& gene_name=="TP53") %>% as.data.frame()
}

