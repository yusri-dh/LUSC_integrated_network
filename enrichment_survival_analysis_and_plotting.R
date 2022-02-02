# initialization -----------------------------------------------------------

library(tidyverse)
library(gprofiler2)
library(TCGAbiolinks)
library(ggplot2)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(ELMER)
library(pheatmap)
library(RColorBrewer)
library(grid)
working_dir = "/home/yusri/Documents/project/LUSC_integrated_network/"
setwd(working_dir)
df = read_csv("graph_visualize/TCGA-LUSC_methylation_tumor_gene_loc_giant.csv")
df_network = read_csv("graph_visualize/TCGA-LUSC_methylation_tumor_gene_giant.csv")

# functional gene  enrichment -----------------------------------------------------------
# Count the number of member of each community and select the communities that have > 200 member
count_membership = dplyr::count(df,leiden_membership)
selected_membership = count_membership$leiden_membership[count_membership$n > 200]
cumsum(count_membership[selected_membership,]$n)/sum(count_membership$n)

query_gene_member <- function(df,member){
  res <-df %>% 
    filter(leiden_membership == member) %>%
    filter(type == "gene") %>%
    select(id,label)
  
  return(res)
}

query_gene_list_ = list()
query_gene_list_label = list()
for (i in selected_membership) {
  print(i)
  query <- query_gene_member(df,i)
  identifier = paste("Community",as.character(i),sep="_")
  query_gene_list_[[identifier]] = query$id
  query_gene_list_label[[identifier]] = query$label
}

query_gene_list = query_gene_list_[order(lengths(query_gene_list_),decreasing=TRUE)]
query_gene_list_label = query_gene_list_label[order(lengths(query_gene_list_),decreasing=TRUE)]

all_gene_label <- (df %>% 
                     filter(type == "gene"))$label
all_gene_id <- (df %>% 
                  filter(type == "gene"))$id

# functional enrichment using g:profiler
gostres <- gost(query = query_gene_list,organism="hsapiens",sources = c("GO","KEGG","REAC","WP"),as_short_link = FALSE)

# find the 3 lowest p value functional class for every data sources
low_p_value_df  <- gostres$result %>% 
  group_by(query,source) %>%
  slice_min(order_by=p_value,n=3) %>%
  select(query,source,term_name,term_size,intersection_size,p_value)

# Save all significant result
write_csv(gostres$result,"./paper/Supplemental/significant_gsea.csv")

# Save the 3 lowest p value functional class
write_csv(low_p_value_df,"./paper/Supplemental/low_p_value_enrichment.csv")


# centrality ranking -----------------------------------------------------------------

high_gene_ranking <- df %>% 
  slice_max(order_by=betweenness,n=20) %>%
  pull(label)

# DEGs data preparation---------------------------------------------------------
query.exp <- GDCquery(project = "TCGA-LUSC", 
                      legacy = FALSE,
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - FPKM-UQ",
                      sample.type = c("Primary Tumor","Solid Tissue Normal"))

clinical_patient_cancer <- GDCquery_clinic("TCGA-LUSC","clinical")

cancer.exp <- GDCprepare(query = query.exp,directory="Data/LUSC/Raw/RNA/")
dataPrep <- TCGAanalyze_Preprocessing(object = cancer.exp, cor.cut = 0.6)
dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfoHT,
                                      method = "gcContent")
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)
# Which samples are Primary Tumor
dataSmTP <- TCGAquery_SampleTypes(getResults(query.exp,cols="cases"),"TP") 
# which samples are solid tissue normal
dataSmNT <- TCGAquery_SampleTypes(getResults(query.exp,cols="cases"),"NT")
dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,dataSmNT],
                            mat2 = dataFilt[,dataSmTP],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")



dataFiltDE <- subset(dataFilt, subset = rownames(dataFilt) %in% rownames(dataDEGs))

# survival analysis of betweenness centrality-----------------------------------
selected_column = "betweenness"
lower_limit = quantile(df[[selected_column]],c(0.1))
upper_limit = quantile(df[[selected_column]],c(0.9))
selected_surv_gene_low_rank<-df %>%
  filter(type=="gene") %>%
  filter(!!as.symbol(selected_column)<=lower_limit) %>%
  arrange(desc(!!as.symbol(selected_column)))

selected_surv_gene_high_rank<-df %>%
  filter(type=="gene") %>%
  filter(!!as.symbol(selected_column)>=upper_limit)  %>%
  arrange(desc(!!as.symbol(selected_column)))



surv_an_low <-TCGAanalyze_SurvivalKM(clinical_patient_cancer,
                                     dataFilt,
                                     Genelist = selected_surv_gene_low_rank$id,
                                     Survresult = FALSE,
)

surv_an_high <-TCGAanalyze_SurvivalKM(clinical_patient_cancer,
                                      dataFilt,
                                      Genelist = (selected_surv_gene_high_rank)$id,
                                      Survresult = FALSE,
)

surv_an_tf <-TCGAanalyze_SurvivalKM(clinical_patient_cancer,
                                      dataFilt,
                                      Genelist = gconvert(c("tp63", "klf5", "sox2"))$target,
                                      Survresult = TRUE,
)
#write_csv(surv_an_low,"./paper/Supplemental/survival_analysis_low_betweenness.csv")
#write_csv(surv_an_high,"./paper/Supplemental/survival_analysis_high_betweenness.csv")

# connection analysis --------------------------------------------
df_network_selected <- df_network %>% 
  filter(link_type == "gene-gene") %>%
  filter(Target_member <= 10)
count_connection <- function(df_network,member){
  df_network_member <- df_network %>% 
    filter(Source_member == member) %>%
    select(Target_member,Source_member) %>% 
    dplyr::count(Target_member)
  return(df_network_member)
}

test <- count_connection(df_network_selected)

connection_matrix <- matrix(0,10,10)
for (member in 1:10) {
  df_member <- count_connection(df_network_selected,member)
  for (i in df_member$Target_member){
    connection_matrix[member,i] = df_member %>%
      filter(Target_member ==i) %>%
      pull()
  }
  
}

diag(connection_matrix) <- 0
ratio_connection_matrix <- connection_matrix/rowSums(connection_matrix)
#ratio_connection_matrix <- (ratio_connection_matrix + t(ratio_connection_matrix))/2

colnames(ratio_connection_matrix) <- seq(1,10,by=1)
rownames(ratio_connection_matrix) <- seq(1,10,by=1)
Cxy_heatmap <- pheatmap(ratio_connection_matrix,
                        display_numbers = TRUE,
                        number_color = "black", 
                        fontsize_number = 8,
                        colorRampPalette(rev(brewer.pal(n =3, name =
                                                          "RdYlBu")))(101),
                        cluster_cols=F,
                        cluster_rows=F,
                        angle_col = 0,
                        legend_breaks = c(0.1, 0.2, 0.3, 0.4, 0.5,max(ratio_connection_matrix)),
                        legend_labels = c("0.1", "0.2", "0.3", "0.4", "0.5", "Ratio Cx(y)"),
                        gaps_row = seq(1,10,by=1)
)
# heatmap --------------------------------------------------------------------
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")

Cxy_heatmap <- pheatmap(ratio_connection_matrix,
         display_numbers = TRUE,
         number_color = "black", 
         fontsize_number = 8,
         colorRampPalette(rev(brewer.pal(n =3, name =
                                           "RdYlBu")))(101),
         cluster_cols=F,
         cluster_rows=F,
         angle_col = 0,
         legend_breaks = c(0.1, 0.2, 0.3, 0.4, 0.5,max(ratio_connection_matrix)),
         legend_labels = c("0.1", "0.2", "0.3", "0.4", "0.5", "Ratio Cx(y)"),
         gaps_row = seq(1,10,by=1)
         )

setHook("grid.newpage", NULL, "replace")
grid.text("Target Community", x=0.4,y=-0.02, gp=gpar(fontsize=16))
grid.text("Source Community", x=0.9,y=0.3, rot=-90, gp=gpar(fontsize=16))

# methylation analysis preparation --------------------------------------------

query_methylation_member <- function(df,df_network,member,link.type="all"){
  if (link.type=="all" ) {
    df_network_methylation = df_network %>% 
      filter(link_type != "gene-gene") %>%
      filter(Source_member == member) %>%
      select(Target,Source)
  }else{
    df_network_methylation = df_network %>% 
      filter(link_type == link.type) %>%
      filter(Source_member == member) %>%
      select(Target,Source)
  }
  
  
  return(df_network_methylation)
}

# find the methylated genes and probes for each community
query_probe_list_all = list()
query_probe_list_hypermethylation = list()
query_probe_list_hypomethylation = list()

query_dmg_list_all = list()
query_dmg_list_hypermethylation = list()
query_dmg_list_hypomethylation = list()

for (i in selected_membership) {
  print(i)
  identifier = paste("Community",as.character(i),sep="_")
  query_probe_list_all[[identifier]] = query_methylation_member(df,df_network,i,"all")  %>% 
    select(Target) %>% 
    distinct() %>%
    pull()
  query_probe_list_hypermethylation[[identifier]] = query_methylation_member(df,df_network,i,"hypermethylation")  %>% 
    select(Target) %>% 
    distinct()  %>%
    pull()
  query_probe_list_hypomethylation[[identifier]] = query_methylation_member(df,df_network,i,"hypomethylation")  %>% 
    select(Target) %>% 
    distinct()  %>%
    pull()
  query_dmg_list_all[[identifier]] = query_methylation_member(df,df_network,i,"all")  %>% 
    select(Source) %>% 
    distinct()  %>%
    pull()
  query_dmg_list_hypermethylation[[identifier]] = query_methylation_member(df,df_network,i,"hypermethylation")  %>% 
    select(Source) %>% 
    distinct()  %>%
    pull()
  query_dmg_list_hypomethylation[[identifier]] = query_methylation_member(df,df_network,i,"hypomethylation")  %>% 
    select(Source) %>% 
    distinct()  %>%
    pull()
}

df_dmg_hypermethyl = sapply(query_dmg_list_hypermethylation, length) %>%
  as_tibble() %>%
  rownames_to_column() %>%
  dplyr::rename(Community = rowname, "Hypermethylated genes" = value )
  
df_dmg_hypomethyl = sapply(query_dmg_list_hypomethylation, length)%>%
  as_tibble() %>%
  rownames_to_column() %>%
  dplyr::rename(Community = rowname, "Hypomethylated genes" = value )

df_probe_hypermethyl = sapply(query_probe_list_hypermethylation, length)%>%
  as_tibble() %>%
  rownames_to_column() %>%
  dplyr::rename(Community = rowname, "Hypermethylated probes" = value )

                                     
df_probe_hypomethyl = sapply(query_probe_list_hypomethylation, length)%>%
  as_tibble() %>%
  rownames_to_column() %>%
  dplyr::rename(Community = rowname, "Hypomethylated probes" = value )

                  


df_methylation_community <- Reduce(function(x, y) merge(x, y), list(df_dmg_hypermethyl, 
                                                                              df_dmg_hypomethyl, 
                                                                              df_probe_hypermethyl,
                                                                              df_probe_hypomethyl)) %>%
  mutate(Community = as.integer(Community)) %>%
  filter(Community <= 10) %>%
  mutate(Community = as.factor(Community)) %>%
  gather(measurement_type, n,"Hypermethylated genes":"Hypomethylated probes",factor_key=TRUE)%>%
  arrange(Community)

write_csv(df_methylation_community,"df_methylation_community.csv")
fig <- ggplot(df_methylation_community,aes(x = Community, y = n, fill = measurement_type))
fig + 
  geom_bar(position="dodge", stat="identity") +
  theme_classic()+
  theme(legend.position = c(1, 1),
        legend.justification = c(1,1),
        legend.direction = "horizontal",
        legend.title = element_blank()) +
  ylab("Number of genes/probes")
ggsave("./paper/Figure/methylation_community.png", dpi = 600)

# find enriched motif and master regulator TF for each communtiy----------------
mae <- get(load("./Result/LUSC/LUSC_mae_hg38.rda"))

get_enriched_motif_community <- function(data, probe_list, member,direction){
  if (direction == "hypo") {
    result_dir = "Result/LUSC/TN_Tumor_vs_Normal/hypo"
    sig.diff <- read.csv(file.path(result_dir,"getMethdiff.hypo.probes.significant.csv"))
    pair <- read.csv(file.path(result_dir,"getPair.hypo.pairs.significant.csv"))
  } else if (direction == "hyper"){
    result_dir = "Result/LUSC/TN_Tumor_vs_Normal/hyper"
    sig.diff <- read.csv(file.path(result_dir,"getMethdiff.hyper.probes.significant.csv"))
    pair <- read.csv(file.path(result_dir,"getPair.hyper.pairs.significant.csv"))
  }
  
  # Identify enriched motif for significantly hypo/hypermethylated probes which 
  # have putative target genes.
  label <- paste(direction,member,sep="_")
  enriched.motif <- ELMER::get.enriched.motif(data = data,
                                              probes = probe_list[[member]], 
                                              dir.out = file.path(result_dir,"community"), 
                                              label = label,
  )
  
  TF <- ELMER::get.TFs(data = data, 
                       group.col = "definition",
                       group1 =  "Primary solid Tumor",
                       group2 = "Solid Tissue Normal",
                       mode = "unsupervised",
                       enriched.motif = enriched.motif,
                       dir.out = file.path(result_dir,"community"), 
                       cores = 1, 
                       label = label)
} 



