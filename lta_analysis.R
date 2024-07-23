#analyse results
#on server:
#start singularity session:
#singularity exec --bind /nfs/research/icortes/ /nfs/research/icortes/belzen/src/structural_variation_202405_amd.sif R -e "source('/nfs/research/icortes/belzen/src/TCGA_WGD_analysis/lta_detection.conf');source('/nfs/research/icortes/belzen/src/lta-analysis/lta_analysis.R')"
#default output for plots here
#/nfs/research/icortes/belzen/results/lta_detection/multi_tsg/plots/
#specify other plot dir that is writable

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

#read lta_detection.conf 
if(!exists("results_dir")) {
  tcga_cohort_path="~/data/metadata/TCGA.data_freeze_20240617.purity_table.tsv"
  results_dir="/Users/belzen//results/LTA-analysis/multi_tsg/"
  plot_dir=paste0(results_dir,"plots/")
  #results_dir="/nfs/research/icortes/belzen/results/lta_detection/multi_tsg/"
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

map_template_vars=c(
                   "${plot_dir}"=plot_dir,
                    "${processed_output_dir}"=results_dir)

# TCGA-----

tcga_cohort=read.table(tcga_cohort_path,header = T)
tcga_cohort$basename=tcga_cohort$sample_id
tcga_cohort = tcga_cohort %>% rowwise() %>% mutate(cohort_id=unlist(strsplit(sample_id, "_"))[1])

tcga_cohort = tcga_cohort %>% filter(qc=="PASS" )
tcga_cohort = tcga_cohort %>% filter(!((snv_n_chromosomes <= 15 | sv_n_chromosomes <= 15) & (snv_n < 100 & sv_n < 10)))
tcga_cohort %>% nrow()

cohort_lst = c("osteos","TCGA-SARC","TCGA-OV","TCGA-GBM")
cohort_lst = c(unique(tcga_cohort$cohort_id))
cohort_lst = c(unique(tcga_cohort$cohort_id),"osteos")
dataset_selection_label="TCGA-BRCA"

tcga_cohort_call_lta = data.frame()
tcga_call_lta_annot = data.frame()
tcga_oncogene_amp_annot = data.frame()
tcga_selected_connected_tsg_regions_annot = data.frame()
tcga_gene_cn_sv_disruptions = data.frame()

for(dataset_selection_label in cohort_lst) {
map_template_vars_cohort = c(map_template_vars,"${patient_basename}"=dataset_selection_label)

gene_cn_sv_disruptions_path = stri_replace_all_fixed(gene_cn_sv_disruptions_path_template,names(map_template_vars_cohort), map_template_vars_cohort,vectorize=F)
cohort_call_lta_path= stri_replace_all_fixed(cohort_call_lta_path_template,names(map_template_vars_cohort), map_template_vars_cohort,vectorize=F)
call_lta_annot_path = stri_replace_all_fixed(call_lta_annot_path_template,names(map_template_vars_cohort), map_template_vars_cohort,vectorize=F)
oncogene_amp_annot_path = stri_replace_all_fixed(oncogene_amp_annot_path_template,names(map_template_vars_cohort), map_template_vars_cohort,vectorize=F)
connected_regions_tsg_instability_path = stri_replace_all_fixed(connected_regions_tsg_instability_path_template,names(map_template_vars_cohort), map_template_vars_cohort,vectorize=F)


if(length(Sys.glob(cohort_call_lta_path))==0){ print(paste0("Missing ",dataset_selection_label)); next() } 
  
cohort_call_lta=read.table(cohort_call_lta_path,sep = "\t",header = T)
cohort_call_lta$cohort_id_original = cohort_call_lta$cohort_id
cohort_call_lta$cohort_id=dataset_selection_label
tcga_cohort_call_lta = rbind_no_colmatch(tcga_cohort_call_lta,cohort_call_lta)

call_lta_annot = read.table(call_lta_annot_path,sep = "\t",header = T)
if(call_lta_annot %>% nrow() > 0) {
call_lta_annot$cohort_id=dataset_selection_label
tcga_call_lta_annot = rbind_no_colmatch(tcga_call_lta_annot,call_lta_annot)
}

oncogene_amp_annot = read.table(oncogene_amp_annot_path,sep = "\t",header = T)
if(oncogene_amp_annot %>% nrow() > 0) {
  oncogene_amp_annot$cohort_id=dataset_selection_label
tcga_oncogene_amp_annot = rbind_no_colmatch(tcga_oncogene_amp_annot,oncogene_amp_annot)
}

selected_connected_tsg_regions_annot = read.table(connected_regions_tsg_instability_path,sep = "\t",header= T)
if(selected_connected_tsg_regions_annot %>% nrow() > 0) {
selected_connected_tsg_regions_annot$cohort_id=dataset_selection_label
tcga_selected_connected_tsg_regions_annot = rbind_no_colmatch(tcga_selected_connected_tsg_regions_annot,selected_connected_tsg_regions_annot)
}
gene_cn_sv_disruptions = read.table(gene_cn_sv_disruptions_path,sep = "\t",header = T)
gene_cn_sv_disruptions$cohort_id=dataset_selection_label
tcga_gene_cn_sv_disruptions = rbind_no_colmatch(tcga_gene_cn_sv_disruptions,gene_cn_sv_disruptions)

}


#Quantification ----

#remove non hgos
if("subtype_short" %in% names(tcga_cohort_call_lta)) {
tcga_cohort_call_lta = tcga_cohort_call_lta %>% filter(subtype_short=="HGOS"|is.na(subtype_short))
#ostoes is 49 of manual set and 8 additional cases
} else {
  tcga_cohort_call_lta$subtype_short=NA
  
}

#cohort overview, dicentric LTA per gene
tcga_gene_cn_sv_disruptions %>%  
  filter(basename %in% tcga_cohort_call_lta$basename) %>%
  filter(lta&tsg_gene_knockout_dicentric) %>%
  group_by(cohort_id,tsg_gene_name) %>% summarize(n=length(unique(basename)),samples=toString(unique(basename))) %>% arrange(-n)

#TSGs most commonly affected and associated with oncogene amps
tcga_gene_cn_sv_disruptions %>%  
  filter(cohort_id!="osteos") %>%
  filter(basename %in% tcga_cohort_call_lta$basename) %>%
  filter(lta&tsg_gene_knockout_dicentric&lta_onco) %>%
  group_by(tsg_gene_name) %>% summarize(n=length(unique(basename)),samples=toString(unique(basename))) %>% arrange(-n)


#additional cases have not got amp so this is subset of manual set
tcga_gene_cn_sv_disruptions %>%  
  filter(basename %in% tcga_cohort_call_lta$basename) %>%
  filter(lta_onco&tsg_gene_knockout_dicentric) %>% 
  #group_by(cohort_id) %>% summarize(n=length(unique(basename)),samples=toString(unique(basename))) %>% arrange(-n)
  group_by(cohort_id,tsg_gene_name) %>% summarize(n=length(unique(basename)),samples=toString(unique(basename))) %>% arrange(-n)


## LTA with other TSGs -----
#In HGOS, the 17 samples that had other TSGs affected by LTA consisted for the 
# vast majority of samples with TP53 LTA detected (n= 12)
#indicating that LTA can also disrupt additional TSGs 
#in 7 of 12 cases we could find connections to the TP53 LTA event

tcga_gene_cn_sv_disruptions = tcga_gene_cn_sv_disruptions %>% 
  dplyr::mutate(sample_has_tp53_lta = basename %in% filter(tcga_cohort_call_lta,lta)$basename)

lta_no_tp53_breakdown = tcga_gene_cn_sv_disruptions %>% 
  filter(basename %in% filter(tcga_cohort_call_lta,subtype_short=="HGOS")$basename) %>%
  filter(tsg_gene_name!="TP53" & lta & tsg_gene_knockout_dicentric) %>%
  #select(basename) %>% unique() %>% nrow()
  group_by(basename,sample_has_tp53_lta) %>% summarize(any_tp53_conn = any(connected_to_tp53_tsg_region)) %>%
  group_by(sample_has_tp53_lta,any_tp53_conn) %>% summarize(sample_cnt=length(unique(basename))) %>% arrange(-sample_cnt)

print(lta_no_tp53_breakdown)
print(tcga_gene_cn_sv_disruptions %>% filter(connected_to_tp53_tsg_region==F) %>%
        filter(basename %in% filter(tcga_cohort_call_lta,subtype_short=="HGOS")$basename) %>%
        filter(tsg_gene_name!="TP53" & lta & tsg_gene_knockout_dicentric) %>%
        filter(sample_has_tp53_lta==F))

#Hence, looking at other TSGs, we only focused on cases without TP53 LTA
#for split plot:
tsg_lta_no_tp53_breakdown = tcga_gene_cn_sv_disruptions %>% 
  #filter(basename %in% filter(tcga_cohort_call_lta,subtype_short=="HGOS")$basename) %>% #check OS first
  filter(sample_has_tp53_lta==F) %>% #cohort call only contains TP53 LTAs
  filter(lta & tsg_gene_knockout_dicentric) %>% 
  group_by(cohort_id,tsg_gene_name) %>% summarize(sample_cnt=length(unique(basename))) %>% arrange(-sample_cnt)

print(tsg_lta_no_tp53_breakdown)

##applying the same to TCGA cohort 
#oncogene amplifications :
#any evidence of tp53 lta or only high conf? => selected not any to look at independant LTA events
tcga_oncogene_amp_annot = tcga_oncogene_amp_annot %>% 
  dplyr::mutate(sample_has_tp53_lta = basename %in% filter(tcga_cohort_call_lta,lta)$basename)

tsg_with_onco_amp_lta_no_tp53_breakdown = tcga_oncogene_amp_annot  %>% 
  #  filter(basename %in% filter(tcga_cohort_call_lta,subtype_short=="HGOS")$basename) %>% #check OS first
  filter(sample_has_tp53_lta==F) %>% #cohort call only contains TP53 LTAs
  filter(lta & tsg_gene_knockout_dicentric) %>% 
  group_by(tsg_gene_name) %>% summarize(sample_cnt=length(unique(basename))) %>% arrange(-sample_cnt)

print(tsg_with_onco_amp_lta_no_tp53_breakdown)

# LTA associated with oncogene amplification ----

#remove HGOS and samples without processing
tcga_oncogene_amp_annot = tcga_oncogene_amp_annot %>% filter(basename %in% filter(tcga_cohort_call_lta,sample_analysed)$basename)
tcga_oncogene_amp_annot = tcga_oncogene_amp_annot %>% dplyr::mutate(lta_dicentric = lta & tsg_gene_knockout_dicentric)
  

#just to check 
tcga_oncogene_amp_annot %>% filter(lta_dicentric & !lta_onco)


tcga_oncogene_amp_annot = tcga_oncogene_amp_annot %>% mutate(amp_cgr = !is.na(region_id))

#for the overall plot regardless of TSG disrupted
amp_samples = tcga_oncogene_amp_annot %>%  group_by(cohort_id,basename) %>% 
  summarize(onco_amp = length(unique(gene_name)),
            amp_cgr=any(amp_cgr),
            lta_onco_amp=any(lta_onco),
            lta_onco_multichrom=any(lta_onco_multichrom),
            lta_onco_dicentric=any(lta_dicentric),
            lta_onco_multichrom_dicentric=any(lta_onco_multichrom&lta_dicentric))


amp_samples = amp_samples %>% dplyr::mutate(label = ifelse(lta_onco_multichrom_dicentric,"Amp_LTA-multi",
                                             ifelse(lta_onco_dicentric,"Amp_LTA-single",
                                             ifelse(amp_cgr,"Amp_CGR","Amp_Other"))),
                                            label_simple = ifelse(lta_onco_dicentric,"Amp_LTA", 
                                                                  ifelse(amp_cgr,"Amp_CGR","Amp_Other")))

amp_samples_overall = amp_samples

## Plot oncogene amplifications with LTA component -----
#tcga_cohort_size = tcga_cohort %>% group_by(cohort_id) %>% summarize(Cohort_Size=n())
#analysed samples:
tcga_cohort_size = tcga_cohort_call_lta %>% filter(sample_analysed) %>% group_by(cohort_id) %>% summarize(Cohort_Size=n(),sample_cnt=Cohort_Size)
cohort_order=c("osteos",tcga_cohort_size %>% arrange(-Cohort_Size) %>% dplyr::select(cohort_id) %>% flatten_chr()) %>% unique()

label_order=c("Amp_LTA","Amp_LTA-multi","Amp_LTA-single","Amp_CGR","Amp_Other")
label_colors=c("Amp_LTA"="purple3","Amp_LTA-multi"="purple3","Amp_LTA-single"="hotpink","Amp_CGR"="darkolivegreen3","Amp_Other"="skyblue2")

pdf(paste0(plot_dir,"LTA_oncogene_amp.pdf"),height=12,width=10)

p = ggplot(amp_samples_overall) + 
  geom_col(data=tcga_cohort_size,aes(y=factor(cohort_id,levels=rev(cohort_order)),x=sample_cnt), fill="grey",alpha=0.5,color="black", linewidth=.1) +
  ggtitle("# samples with oncogene amplifications, connected to LTA event disrupting any TSG") +
  scale_fill_manual(values = label_colors) + 
  theme_bw() + ylab("") + guides(fill=guide_legend(title=""))

print(p +   geom_bar(aes(y=factor(cohort_id,levels=rev(cohort_order)),fill=factor(label,levels=names(label_colors))),color="black", linewidth=.1,position = position_stack(reverse = TRUE)))

p = p +   geom_bar(aes(y=factor(cohort_id,levels=rev(cohort_order)),fill=factor(label_simple,levels=names(label_colors))),color="black", linewidth=.1,position = position_stack(reverse = TRUE))
print(p)

saveRDS(p, file=paste0(plot_dir,"LTA_oncogene_amp.any_tsg.rds"))


target_tsg="TP53"
for(target_tsg in c("TP53",tsg_with_onco_amp_lta_no_tp53_breakdown$tsg_gene_name)) {
#tp53 only
amp_samples_single_tsg = tcga_oncogene_amp_annot %>% 
  filter(tsg_gene_name==target_tsg | is.na(tsg_gene_name)) %>%
  group_by(cohort_id,tsg_gene_name,basename,sample_has_tp53_lta) %>% 
  summarize(lta_onco_dicentric=any(lta_dicentric),
            lta_onco_multichrom_dicentric=any(lta_onco_multichrom&lta_dicentric))
#for genes other than tp53, only count if sample has no tp53 LTA
amp_samples_single_tsg = amp_samples_single_tsg %>% dplyr::mutate(assign_lta_event = tsg_gene_name=="TP53" | sample_has_tp53_lta==F)

amp_samples_single_tsg = amp_samples_overall %>% select(basename,onco_amp,amp_cgr) %>% 
  left_join(amp_samples_single_tsg) 
amp_samples_single_tsg = amp_samples_single_tsg %>% dplyr::mutate(
  assign_lta_event = ifelse(is.na(assign_lta_event),F,assign_lta_event),
  label = ifelse(lta_onco_multichrom_dicentric&assign_lta_event,"Amp_LTA-multi", 
                 ifelse(lta_onco_dicentric&assign_lta_event,"Amp_LTA-single", 
                        ifelse(amp_cgr,"Amp_CGR","Amp_Other"))),
  label_simple = ifelse(lta_onco_dicentric&assign_lta_event,"Amp_LTA", 
                               ifelse(amp_cgr,"Amp_CGR","Amp_Other"))
)


p = ggplot(amp_samples_single_tsg) + 
  geom_col(data=tcga_cohort_size,aes(y=factor(cohort_id,levels=rev(cohort_order)),x=sample_cnt), fill="grey",alpha=0.5,color="black", linewidth=.1) +
  ggtitle(paste0("# samples with oncogene amplifications, connected to LTA event disrupting ",target_tsg)) +
  scale_fill_manual(values = label_colors) +
  theme_bw() + ylab("") + guides(fill=guide_legend(title=""))

p = p + geom_bar(aes(y=factor(cohort_id,levels=rev(cohort_order)),fill=factor(label_simple,levels=names(label_colors))),color="black", linewidth=.1,,position = position_stack(reverse = TRUE)) 
print(p)
#print(p + geom_bar(aes(y=factor(cohort_id,levels=rev(cohort_order)),fill=factor(label,levels=names(label_colors))),color="black", linewidth=.1,,position = position_stack(reverse = TRUE)) )
if(target_tsg=="TP53") {
  amp_samples_tp53=amp_samples_single_tsg
  plot_amp_samples_tp53=p
}


}

#all non tp53 together
#for genes other than tp53, only count if sample has no tp53 LTA

amp_samples_no_tp53 = tcga_oncogene_amp_annot %>% 
  filter(tsg_gene_name!="TP53" & sample_has_tp53_lta==F) %>%
  group_by(cohort_id,basename) %>% 
  summarize(lta_onco_dicentric=any(lta_dicentric),
            lta_onco_multichrom_dicentric=any(lta_onco_multichrom&lta_dicentric))
amp_samples_no_tp53 = amp_samples_overall %>% select(basename,onco_amp,amp_cgr) %>% 
  left_join(amp_samples_no_tp53) 
amp_samples_no_tp53 = amp_samples_no_tp53 %>% dplyr::mutate(
  lta_onco_dicentric = ifelse(is.na(lta_onco_dicentric),F,lta_onco_dicentric),
  lta_onco_multichrom_dicentric=ifelse(is.na(lta_onco_multichrom_dicentric),F,lta_onco_multichrom_dicentric),
  label = ifelse(lta_onco_multichrom_dicentric,"Amp_LTA-multi", 
                 ifelse(lta_onco_dicentric,"Amp_LTA-single", 
                        ifelse(amp_cgr,"Amp_CGR","Amp_Other"))),
  label_simple = ifelse(lta_onco_dicentric,"Amp_LTA", 
                        ifelse(amp_cgr,"Amp_CGR","Amp_Other"))
)


p = ggplot(amp_samples_no_tp53) + 
  geom_col(data=tcga_cohort_size,aes(y=factor(cohort_id,levels=rev(cohort_order)),x=sample_cnt), fill="grey",alpha=0.5,color="black", linewidth=.1) +
  ggtitle(paste0("# samples with oncogene amplifications, with LTA involving other TSGs (excl TP53 LTA)")) +
  scale_fill_manual(values = label_colors) +
  theme_bw() + ylab("") + guides(fill=guide_legend(title=""))

p = p + geom_bar(aes(y=factor(cohort_id,levels=rev(cohort_order)),fill=factor(label_simple,levels=names(label_colors))),color="black", linewidth=.1,,position = position_stack(reverse = TRUE)) 
print(p)
#print(p + geom_bar(aes(y=factor(cohort_id,levels=rev(cohort_order)),fill=factor(label,levels=names(label_colors))),color="black", linewidth=.1,,position = position_stack(reverse = TRUE)) )

plot_amp_samples_no_tp53=p

dev.off()

pdf(paste0(plot_dir,"LTA_oncogene_amp.summary.landscape.pdf"),height=5,width=15)

print(plot_amp_samples_tp53  + coord_flip() +  scale_y_discrete(limits=rev) + theme(axis.text.x = element_text(angle=90)))
print(plot_amp_samples_no_tp53  + coord_flip() +  scale_y_discrete(limits=rev) + theme(axis.text.x = element_text(angle=90)))

dev.off()

saveRDS(plot_amp_samples_tp53, file=paste0(plot_dir,"LTA_oncogene_amp.tp53.rds"))
saveRDS(plot_amp_samples_no_tp53, file=paste0(plot_dir,"LTA_oncogene_amp.no_tp53.rds"))



# Extended plot show all LTAs ----
#TODO move this up
#LTA overall captures all TSG disruptions by LTA 
#remove HGOS and samples without processing
tcga_gene_cn_sv_disruptions = tcga_gene_cn_sv_disruptions %>% filter(basename %in% filter(tcga_cohort_call_lta,sample_analysed)$basename)
tcga_gene_cn_sv_disruptions = tcga_gene_cn_sv_disruptions %>% dplyr::mutate(lta_dicentric = lta & tsg_gene_knockout_dicentric)

#equiv of amp_cgr for tsgs is tsg_gene_knockout_sv_component

#for the overall plot regardless of TSG disrupted
tsg_disruption_samples= tcga_gene_cn_sv_disruptions %>%  
  #filter(basename %in% filter(tcga_cohort_call_lta,subtype_short=="HGOS")$basename) %>% #check OS first
  # filter(tsg_gene_name=="TP53") %>%
  filter(tsg_gene_knockout_sv_component) %>% 
  group_by(cohort_id,basename) %>% 
  summarize(tsg_sv_ko = length(unique(tsg_gene_name)),
            lta_dicentric=any(lta_dicentric),
            lta_multichrom_dicentric=any(lta_dicentric&lta_multichrom)) %>% 
  dplyr::mutate(label = ifelse(lta_dicentric,"LTA","TSG_koSV"),label_simple = label)


tsg_per_ct =  tsg_disruption_samples %>% group_by(cohort_id) %>% summarize(TSG_koSV=length(unique(basename)), #should be identical to n(),
                                                                                     LTA=sum(lta_dicentric),
                                                                                     LTA_multi=sum(lta_multichrom_dicentric))

#tp53 only version
tsg_disruption_samples_tp53 = tcga_gene_cn_sv_disruptions %>%  
  #filter(basename %in% filter(tcga_cohort_call_lta,subtype_short=="HGOS")$basename) %>% #check OS first
  filter(tsg_gene_name=="TP53") %>%
  filter(tsg_gene_knockout_sv_component) %>% 
  group_by(cohort_id,basename) %>% 
  summarize(tsg_sv_ko = length(unique(tsg_gene_name)),
            lta_dicentric=any(lta_dicentric),
            lta_multichrom_dicentric=any(lta_dicentric&lta_multichrom)) %>% 
  dplyr::mutate(label = ifelse(lta_dicentric,"LTA","TSG_koSV"),label_simple = label)

#other tsgs no tp53  => excluding tumors with tp53 
tsg_disruption_samples_no_tp53= tcga_gene_cn_sv_disruptions %>%  
  #filter(basename %in% filter(tcga_cohort_call_lta,subtype_short=="HGOS")$basename) %>% #check OS first
  filter(tsg_gene_name!="TP53" & sample_has_tp53_lta==F) %>%
  filter(lta_dicentric) %>%
  group_by(cohort_id,basename) %>% 
  summarize(lta_dicentric=any(lta_dicentric),
            lta_multichrom_dicentric=any(lta_dicentric&lta_multichrom)) %>% 
  dplyr::mutate(label = ifelse(lta_dicentric,"LTA",NA),label_simple = label)



tsg_per_ct_tp53 =  tsg_disruption_samples_tp53 %>% group_by(cohort_id) %>% summarize(TSG_koSV=length(unique(basename)), #should be identical to n(),
                                                                           LTA=sum(lta_dicentric),
                                                                           LTA_multi=sum(lta_multichrom_dicentric))

tsg_per_ct_no_tp53 =  tsg_disruption_samples_no_tp53 %>% group_by(cohort_id) %>% summarize(
                                                                                           LTA=sum(lta_dicentric),
                                                                                           LTA_multi=sum(lta_multichrom_dicentric))



amp_per_ct =  amp_samples %>% group_by(cohort_id) %>% summarize(OncoAmp=length(unique(basename)), #should be identical to n(),
                                                                          OncoAmp_LTA=sum(lta_onco_dicentric,na.rm = T),
                                                                          OncoAmp_LTA_multi=sum(lta_onco_multichrom_dicentric,na.rm = T))


amp_per_ct_tp53 =  amp_samples_tp53 %>% group_by(cohort_id) %>% summarize(OncoAmp=length(unique(basename)), #should be identical to n(),
                                                                OncoAmp_LTA=sum(lta_onco_dicentric,na.rm = T),
                                                                OncoAmp_LTA_multi=sum(lta_onco_multichrom_dicentric,na.rm = T))


amp_per_ct_no_tp53 =  amp_samples_no_tp53 %>% group_by(cohort_id) %>% summarize(OncoAmp=length(unique(basename)), #should be identical to n(),
                                                                          OncoAmp_LTA=sum(lta_onco_dicentric,na.rm = T),
                                                                          OncoAmp_LTA_multi=sum(lta_onco_multichrom_dicentric,na.rm = T))


lta_prevalence_tp53 = tcga_cohort_size  %>% select(-sample_cnt) %>% left_join(amp_per_ct_tp53) %>% left_join(tsg_per_ct_tp53) %>% dplyr::mutate( across(.cols=everything(), ~replace_na(.x, F)) )
lta_prevalence = tcga_cohort_size  %>% select(-sample_cnt) %>% left_join(amp_per_ct) %>% left_join(tsg_per_ct) %>% dplyr::mutate( across(.cols=everything(), ~replace_na(.x, F)) )
lta_prevalence_no_tp53 = tcga_cohort_size  %>% select(-sample_cnt) %>% left_join(amp_per_ct_no_tp53) %>% left_join(tsg_per_ct_no_tp53) %>% dplyr::mutate( across(.cols=everything(), ~replace_na(.x, F)) )


label_order=c("Cohort_Size","TSG_koSV","OncoAmp","LTA","OncoAmp_LTA")
prevalence_label_colors=c("Cohort_Size"="grey","TSG_koSV"="skyblue2","OncoAmp"="red3","LTA"="purple3","OncoAmp_LTA"="hotpink")
prevalence_label_colors=c("Cohort_Size"="grey","LTA"="purple3","OncoAmp"="red3","OncoAmp_LTA"="hotpink")


pdf(paste0(plot_dir,"LTA_prevalence_breakdown.pdf"),height=20,width=10)

p_lta_tp53 = ggplot(lta_prevalence_tp53 %>% 
         select(-contains("multi"),-TSG_koSV) %>%
         pivot_longer(-cohort_id,values_to = "cnt",names_to="attr") %>%
         dplyr::mutate(attr = factor(attr,levels=names(prevalence_label_colors)))
         ) +
  geom_col(aes(x=cnt,y=factor(attr,levels=names(prevalence_label_colors) %>% rev()),fill=attr),color="black", linewidth=.1,position = position_stack(reverse = TRUE)) +
  facet_grid(rows=vars(factor(cohort_id,levels=cohort_order))) + 
  theme_bw() + ylab("") + xlab("# samples") +
  scale_fill_manual(values=prevalence_label_colors) +
  ggtitle("LTA prevalence - events disrupting TP53")
print(p_lta_tp53)

p_lta_no_tp53 = ggplot(lta_prevalence_no_tp53 %>% 
             select(-contains("multi")) %>%
             pivot_longer(-cohort_id,values_to = "cnt",names_to="attr") %>%
             dplyr::mutate(attr = factor(attr,levels=names(prevalence_label_colors)))
) +
  geom_col(aes(x=cnt,y=factor(attr,levels=names(prevalence_label_colors) %>% rev()),fill=attr),color="black", linewidth=.1,position = position_stack(reverse = TRUE)) +
  facet_grid(rows=vars(factor(cohort_id,levels=cohort_order))) + 
  theme_bw() + ylab("") + xlab("# samples") +
  scale_fill_manual(values=prevalence_label_colors) +
  ggtitle("LTA prevalence - excluding tumors with TP53 LTA")
print(p_lta_no_tp53)

p = ggplot(lta_prevalence %>% 
         select(-contains("multi"),-TSG_koSV) %>%
         pivot_longer(-cohort_id,values_to = "cnt",names_to="attr") %>%
         dplyr::mutate(attr = factor(attr,levels=names(prevalence_label_colors)))
) +
  geom_col(aes(x=cnt,y=factor(attr,levels=names(prevalence_label_colors) %>% rev()),fill=attr),color="black", linewidth=.1,position = position_stack(reverse = TRUE)) +
  facet_grid(rows=vars(factor(cohort_id,levels=cohort_order))) + 
  theme_bw() + ylab("") + xlab("# samples") +
  scale_fill_manual(values=prevalence_label_colors) +
  ggtitle("LTA prevalence - events disrupting any TSG")
print(p)

saveRDS(p, file=paste0(plot_dir,"LTA_prevalence_breakdown.any_tsgs.rds"))
saveRDS(p_lta_tp53, file=paste0(plot_dir,"LTA_prevalence_breakdown.tp53.rds"))
saveRDS(p_lta_no_tp53, file=paste0(plot_dir,"LTA_prevalence_breakdown.no_tp53.rds"))

dev.off()

pdf(paste0(plot_dir,"LTA_prevalence_breakdown.landscape.pdf"),height=5,width=20)

print(p_lta_tp53  + coord_flip() +  facet_grid(cols=vars(factor(cohort_id,levels=cohort_order))) + scale_y_discrete(limits=rev) + theme(axis.text.x = element_text(angle=90)))
print(p_lta_no_tp53  + coord_flip() + facet_grid(cols=vars(factor(cohort_id,levels=cohort_order))) + scale_y_discrete(limits=rev) + theme(axis.text.x = element_text(angle=90)))

dev.off()



#rate of LTA------

print("In all cases, LTA detected after application of the dicentric filter")
print("TP53 LTA")

lta_prevalence_tp53_export = lta_prevalence_tp53 %>% select(-TSG_koSV,-contains("multi")) 

tcga_colsums_tp53 = lta_prevalence_tp53_export %>% filter(cohort_id!="osteos") %>% select(-cohort_id) %>% colSums() %>% t() %>% as.data.frame()
tcga_colsums_tp53$cohort_id="TCGA_overall"
lta_prevalence_tp53_export = rbind(lta_prevalence_tp53_export,tcga_colsums_tp53)
lta_prevalence_tp53_export = lta_prevalence_tp53_export %>% mutate(LTA_frac = LTA/Cohort_Size, 
                                           LTA_OncoAmp_frac = OncoAmp_LTA/Cohort_Size) %>% arrange(-LTA_OncoAmp_frac)
lta_prevalence_tp53_export[is.na(lta_prevalence_tp53_export)]=0

print(lta_prevalence_tp53_export)
print(lta_prevalence_tp53_export %>% filter(cohort_id=="TCGA_overall"))

print("other TSG with LTA, excluding tumors with TP53 LTA")


lta_prevalence_no_tp53_export = lta_prevalence_no_tp53 %>% select(-contains("multi")) 
tcga_colsums_no_tp53 = lta_prevalence_no_tp53_export %>% filter(cohort_id!="osteos") %>% select(-cohort_id) %>% colSums() %>% t() %>% as.data.frame()
tcga_colsums_no_tp53$cohort_id="TCGA_overall"
lta_prevalence_no_tp53_export = rbind(lta_prevalence_no_tp53_export,tcga_colsums_no_tp53)
lta_prevalence_no_tp53_export = lta_prevalence_no_tp53_export %>% mutate(LTA_frac = LTA/Cohort_Size, 
                                                                   LTA_OncoAmp_frac = OncoAmp_LTA/Cohort_Size) %>% arrange(-LTA_OncoAmp_frac)

lta_prevalence_no_tp53_export[is.na(lta_prevalence_no_tp53_export)]=0
print(lta_prevalence_no_tp53_export)

print(lta_prevalence_no_tp53_export %>% filter(cohort_id=="TCGA_overall"))


print("any type of LTA ")


lta_prevalence_export = lta_prevalence %>% select(-TSG_koSV,-contains("multi")) 
tcga_colsums = lta_prevalence_export %>% filter(cohort_id!="osteos") %>% select(-cohort_id) %>% colSums() %>% t() %>% as.data.frame()
tcga_colsums$cohort_id="TCGA_overall"
lta_prevalence_export = rbind(lta_prevalence_export,tcga_colsums)
lta_prevalence_export = lta_prevalence_export %>% mutate(LTA_frac = LTA/Cohort_Size, 
                                                                         LTA_OncoAmp_frac = OncoAmp_LTA/Cohort_Size) %>% arrange(-LTA_OncoAmp_frac)

lta_prevalence_export[is.na(lta_prevalence_export)]=0
print(lta_prevalence_export)
print(lta_prevalence_export %>% filter(cohort_id=="TCGA_overall"))


write.table(lta_prevalence_export,paste0(plot_dir,"LTA_prevalence.any_tsgs.tsv"),sep="\t",row.names = F,col.names = T)
write.table(lta_prevalence_no_tp53_export, paste0(plot_dir,"LTA_prevalence.no_tp53.tsv"),sep="\t",row.names = F,col.names = T)
write.table(lta_prevalence_tp53_export, paste0(plot_dir,"LTA_prevalence.tp53.tsv"),sep="\t",row.names = F,col.names = T)

