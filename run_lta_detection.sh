#! /usr/bin/sh

singularity_img="/nfs/research/icortes/belzen/src/structural_variation_202405_amd.sif"
singularity_bind_dir="/nfs/research/icortes/"
script_dir="/nfs/research/icortes/belzen/src/"
script="/nfs/research/icortes/belzen/src/lta_detection.R"
resources="--mem=150Gb --time=35:00:00"

#run_pattern="TCGA-PRAD_*"
#run_pattern=$1

cohort_table="/nfs/research/icortes/DATA/TCGA_WGS/hmf/data_freeze_20240617/purity_table.tsv"
config_cohort="/nfs/research/icortes/belzen/src/TCGA_WGD_analysis/lta_detection.conf"

#run_cohort=/nfs/research/icortes/belzen/src/run_shatterseek_cohort.lst
#cat ${cohort_table} |cut -f 1 -d "_" | uniq  > $run_cohort
#lst=`echo ${run_cohort}`
#echo $lst

lst=`cat ${cohort_table} |cut -f 1 -d "_" | uniq`

for ct in $lst; do

    #add  underscore to prevent overgrapping eg for TCGA-OV
    #run_pattern="${ct}_*"
#    singularity_settings="source('${config_cohort}');dataset_selection_label='${ct}';skip_plots=T;"
    singularity_settings="source('${config_cohort}');dataset_selection_label='${ct}';"
    
    runscript="${script_dir}run_lta_detection.${ct}.sh"

     echo "#!/usr/bin/env bash" > ${runscript}
echo "SINGULARITYENV_R_MAX_VSIZE=120Gb" >> ${runscript}
     echo "singularity exec --bind ${singularity_bind_dir} ${singularity_img} R -e \"${singularity_settings}source('${script}')\"" >> ${runscript}

    sbatch ${resources} ${runscript}


done


