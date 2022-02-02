# LUSC_integrated_network

The source code of the paper "Analyzing Integrated Network of Methylation and
Gene Profile in Lung Squamous Cell Carcinoma" by Yusri Dwi Heryanto, Kotoe Katayama, and Seiya Imoto.

The steps to run the analyses:
1. Using R, run the DEGs_data_preparation_and_DMCs_ELMER_pipeline.R. It will download all LUSC TCGA gene expresion, methylation, and clinical data.
Here, we preprocessing and standardize DEGs data, perform DMCs network reconstruction using ELMER, and save the result to: 
     - File TCGA-LUSC_Tumor.tsv (DEGs list) 
     - File TCGA-LUSC_gene_location.csv (DEGs location and strand)
     - Directory ./Data for all TCGA methylation, clinical, and gene expresion data downloaded by ELMER
     - Directory ./Result for all ELMER result (DMCs network, motif, and master regulator TF)
2. Using Julia, run DEGs_network_reconstruction.jl to reconstruct DEGs network. The output is a file output_TCGA-LUSC_Tumor.tsv
3. Using Julia, run network_integration_centrality_measurement.jl to integrate the DMCs-DEGs network, perform centrality measurement, and community detection analysis. The outputs were saved in directory "graph visualize". The output files:
     - TCGA-LUSC_integrated_network_nodes.csv (nodes list)
     - TCGA-LUSC_integrated_network_edges.csv (edges list)
4. Using R, run enrichment_survival_analysis_and_plotting.R to perform GSEA, survival analysis, motif and master regulator TFs identification for each community.
