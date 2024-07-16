#!/usr/bin/env bash
SINGULARITYENV_R_MAX_VSIZE=150Gb
singularity exec --bind /nfs/research/icortes/ /nfs/research/icortes/belzen/src/structural_variation_202405_amd.sif R -e "source('/nfs/research/icortes/belzen/src/TCGA_WGD_analysis/osteos.lta_detection.conf');dataset_selection_label='osteos';source('/nfs/research/icortes/belzen/src/lta_detection.R')"

