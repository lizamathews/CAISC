setwd("/data/YoungLab_Bioinfo/DENDROsimulation/infercnvsubclustersHMM3_4fold_addingchr59")
rm(list=ls())
library(infercnv)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="expressionmatrix.txt",
                                    annotations_file="annoationCell.txt",
                                    delim="\t",
                                    gene_order_file="annoationGene.txt",
                                    ref_group_names=c("Healthy")) 
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="./", 
                             cluster_by_groups=TRUE,
                             analysis_mode = "subclusters",
                             k_obs_groups = 4,
                             denoise=TRUE,
                             HMM=TRUE,
			     HMM_type = "i3",
                             num_threads=6,
                             no_plot=FALSE
                             )
save(infercnv_obj, file="DENDROsimulation_Result.RData")

