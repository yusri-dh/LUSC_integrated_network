
library(MultiAssayExperiment)
library(ELMER)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(EDASeq)
library(edgeR)
library(dplyr)
working_dir = "/home/yusri/Documents/project/LUSC_integrated_network/"

setwd(working_dir)
genome <- "hg38"
cancer_name = "LUSC"
project_name = "TCGA-LUSC"
data_path = paste("Data/",cancer_name)
if(!dir.exists(file.path("Data",cancer_name))){
  getTCGA(disease = cancer_name,
          Meth = TRUE,
          RNA = TRUE,
          Clinic = TRUE,
          basedir="./Data",
          genome=genome)
}

ELMER::TCGA.pipe(cancer_name,
          cores = 10,
          mode = "unsupervised",
          Data = data_path,
          analysis = c("createMAE","distal.probes","diffMeth","pair","motif","TF.search"),
          group.col = "definition",
          group1 = "Primary solid Tumor",
          group2 = "Solid Tissue Normal",
          diff.dir = "hyper",
)
print("hypermethylation analysis finished")
ELMER::TCGA.pipe(cancer_name,
          cores = 10,
          mode = "unsupervised",
          Data = data_path,
          group.col = "definition",
          analysis = c("createMAE","distal.probes","diffMeth","pair","motif","TF.search"),
          group1 = "Primary solid Tumor",
          group2 = "Solid Tissue Normal",
          diff.dir = "hypo",
)
print("hypomethylation analysis finished")



query.gene.exp <- GDCquery(project = project_name, 
                      legacy = FALSE,
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - FPKM-UQ",
                      sample.type = c("Primary Tumor","Solid Tissue Normal"))




gene.exp <- GDCprepare(query = query.gene.exp,
                       save = TRUE,
                       directory = file.path("Data/LUSC","Raw","RNA"),
                       save.filename = "lusc_gene.rda",
                       summarizedExperiment = TRUE)


# Gene expression data preparation
dataPrep <- TCGAanalyze_Preprocessing(object = gene.exp, cor.cut = 0.6)
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
logdataFiltDE <- log1p(dataFiltDE)

dataFiltDETumor<- logdataFiltDE[,colnames(dataFilt) %in% dataSmTP]
dataFiltDENormal<- logdataFiltDE[,colnames(dataFilt) %in% dataSmNT]

normlogdataFiltDETumor <- t(scale(t(dataFiltDETumor)))
normlogdataFiltDENormal <- t(scale(t(dataFiltDENormal)))

normlogdataFiltDETumor <- normlogdataFiltDETumor[complete.cases(normlogdataFiltDETumor), ] 
normlogdataFiltDENormal <- normlogdataFiltDENormal[complete.cases(normlogdataFiltDENormal), ] 


write.table(normlogdataFiltDETumor,paste(project_name,"Tumor.tsv",sep="_"),sep=" ",col.names=NA)
write.table(normlogdataFiltDENormal,paste(project_name,"Normal.tsv",sep="_"),sep=" ",col.names=NA)
write.csv(as.data.frame(rowRanges(gene.exp))[c("ensembl_gene_id","external_gene_name","seqnames","start","end","strand")],paste(project_name,"gene_location.csv",sep="_"))
