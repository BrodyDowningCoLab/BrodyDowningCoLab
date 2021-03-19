###########################
## Update h2o
if ("package:h2o" %in% search()) { detach("package:h2o", unload=TRUE) }
if ("h2o" %in% rownames(installed.packages())) { remove.packages("h2o") }
pkgs <- c("RCurl","jsonlite")
for (pkg in pkgs) {
  if (! (pkg %in% rownames(installed.packages()))) { install.packages(pkg) }
}
install.packages("h2o", type="source", 
                 repos="http://h2o-release.s3.amazonaws.com/h2o/rel-zermelo/4/R")
###########################
library("h2o")
library("bigrquery")
library("tidyr")
library("dplyr")


## Set up Google Bigquery project name
billing <- "tcga-masked"

## Downlod clinicalData for masked CNV data
clinicalQuery <- paste(
  "SELECT
     cd.case_gdc_id, cd.project_short_name, cd.age_at_diagnosis, cd.race, cd.gender, cd.ethnicity
   FROM
   `isb-cgc.TCGA_bioclin_v0.Clinical` cd")
# clinicalData_masked is used to extract metadata, i.e. cancer_type, for masked CNV dataset
clinicalData_masked <- query_exec(clinicalQuery, project,max_pages=Inf, use_legacy_sql=FALSE)

# clinicalData is used to extract gender and cancer type data for unmasked CNVs
#clinicalData<-readRDS("clinicalData.rds")

##### Upstream/Downstream analysis
## Make the dataframe to store data
cnv_auc <- data.frame(matrix(ncol=6, nrow=20)) # Store auc for each CNV
rocPlot <- data.frame() # Store tpr/fpr values for ROC Plots
metCount <- data.frame() # Store methylation probe count data

# Declare DNA segment (from )
geneName <- cnvList[8] 
chr <- paste0('chr',sapply(strsplit(geneName,"_"),`[`,1))
start <- as.numeric(sapply(strsplit(geneName,"_"),`[`,2))
end <- as.numeric(sapply(strsplit(geneName,"_"),`[`,3))
#plotName <- paste0(chr,":",start,"-",end)

# (1) Random sampling of the same size up/downstream + within CNV
# (2) percentage of CNV length instead of 1MB
cnv_len <- end-start

################
## Plot methylation probe distribution
# Extract probe IDs within 2MB of CNVs
billing <- "tcga-masked"
cpg_query <- paste("SELECT
                 cn.CpG_probe_id AS probe_ID,
                 cn.chromosome AS Chr,
                 cn.position AS Start_Pos,
                 FROM
                 `isb-cgc.platform_reference.GDC_hg38_methylation_annotation` cn
                 WHERE   		         	
                 (cn.chromosome = 'chr20' AND			
                   cn.position>",lb, "AND 	
                   cn.position<",ub,")")

tb <- bq_project_query(billing, cpg_query)
cpg_spread <- bq_table_download(tb)
cpg_spread$cpg_name <- geneName
metCount <- rbind(metCount,cpg_spread)

# Plot density of cpg probes
library(ggplot2)
ggplot(cpg_spread) + 
  geom_histogram(aes(x=Start_Pos),binwidth=cnv_len)+theme_bw()+
  annotate("rect", xmin = start, xmax = end, ymin = 0, ymax = Inf,
           alpha = .2)+annotate(geom="text", x=start, y=50, label=plotName)+
  xlab('Chromosome Location')+ggtitle(paste("Methylation Probe Count within 2Mb Window of",plotName))
################

## Extract methylation data for a specific gene or for a specific CNV segment
query <- paste("WITH PID as (SELECT
                 cn.CpG_probe_id AS probe_ID,
                 cn.chromosome AS Chr,
                 cn.position AS Start_Pos,
                 FROM
                 `isb-cgc.platform_reference.GDC_hg38_methylation_annotation` cn
                 WHERE   		         	
                 (cn.chromosome = 'chr4' AND			
                   cn.position>",start-cnv_len, "AND 	
                   cn.position<",start,"))
               
               SELECT
               met.probe_ID, met.beta_value, met.case_gdc_id, met.sample_barcode
               FROM `isb-cgc.TCGA_hg38_data_v0.DNA_Methylation_chr4` met
               JOIN PID ON (PID.probe_ID = met.probe_ID)")

# Get data from Google Bigquery
tb <- bq_project_query(billing, query)
gene <- bq_table_download(tb,page_size=1e6)

##### Store all cnv-met data for modeling
gene <- gene[!gene$probe_ID=="",]
gene$Tissue <- clinicalData_masked[match(gene$case_gdc_id,
                        clinicalData_masked$case_gdc_id),"project_short_name"]
gene <- gene[!is.na(gene$Tissue),]
gene$sampleType <- sapply(strsplit(sapply(strsplit(gene$sample_barcode,"-"), 
                                          `[`, 4),"[aA-zZ]+"),`[`,1)
gene$chr_loc <- "4_9459504_9477699:2 CNVs downstream"
head(gene)
met_dat <- rbind(met_dat,gene)
saveRDS(met_dat,"metDat.rds")

cnvLen_auc <- data.frame(matrix(ncol=7, nrow=0)) # Store auc for each CNV
names(cnvLen_auc) <- c("cnv","loc","1","2","3","4","5")
cnv_rocPlot <- data.frame()

glm_auc <- data.frame(matrix(ncol=7, nrow=20)) # Store auc for each CNV
names(glm_auc) <- c("cnv","loc","1","2","3","4","5")
glm_rocPlot <- data.frame()
### TODO
# For each CNV: (1) CNV (2) CNV+1 upstream (3) CNV+2 upstream (4) CNV+1 downstream
# (5) CNV+2 downstream (6) CNV+1 both (7) CNV+2 both
# Modeling: GBM vs. GLM

cnvList <- unique(sapply(strsplit(met_dat$chr_loc,':'),`[`,1))

h2o.init(nthreads=6,min_mem_size = "4g")
# Pick one CNV and its surrounding regions


for (n in 2:length(cnvList)){
j=nrow(cnvLen_auc)+1

cnv <- cnvList[n]
cnvLen_auc[j,1] <- cnv
cnvLen_auc[j,2] <- "clean"

glm_auc[j,1] <- cnv
glm_auc[j,2] <- "clean"
# Extract CNV-only data to find baseline AUC and number of predictors to use
clean <- met_dat[met_dat$chr_loc==cnv,]
clean_wide <- clean %>% pivot_wider(names_from=probe_ID,values_from=beta_value,
                                    values_fn=list(beta_value=mean))
numProbe <- ncol(clean_wide)-5
clean_wide <- clean_wide[clean_wide$sampleType=="01",]
clean_wide$Tissue <- ifelse(clean_wide$Tissue=="TCGA-COAD","COAD","Normal")
response <- "Tissue"
predictors <- setdiff(names(clean_wide), c(response,"case_gdc_id",
                        "sample_barcode","aliquot_id","sampleType","chr_loc"))
if (length(predictors)>50){
  predictors <- predictors[sample(length(predictors),50)]
  numProbe <- 50
}
clean_wide[[response]] <- as.factor(clean_wide[[response]]) 


train.hex<- as.h2o(clean_wide, destination_frame = "train.hex")

plotNum <- sample(3:7,1)
for (i in 3:7){
  
  model <- h2o.gbm(x=predictors,
                  y=response,
                  training_frame = train.hex,
                  nfolds=10,
                  ntrees=100)
  
  model1 <- h2o.glm(x=predictors,
                   y=response,
                   training_frame = train.hex,
                   nfolds=10)
  # Record the AUC
  cnvLen_auc[j,i] <- h2o.auc(model, train=FALSE, xval=TRUE)
  glm_auc[j,i] <- h2o.auc(model1, xval=TRUE)
  
  if (i==plotNum){
    tpr <- h2o.performance(model,xval=T)@metrics$thresholds_and_metric_scores["tpr"]
    fpr <- h2o.performance(model,xval=T)@metrics$thresholds_and_metric_scores["fpr"]
    tmp <- cbind(tpr,fpr)
    tmp$cnv <- cnv
    tmp$loc <- "clean"
    cnv_rocPlot <- rbind(cnv_rocPlot,tmp)
    
    tpr1 <- h2o.performance(model1,xval=T)@metrics$thresholds_and_metric_scores["tpr"]
    fpr1 <- h2o.performance(model1,xval=T)@metrics$thresholds_and_metric_scores["fpr"]
    tmp1 <- cbind(tpr1,fpr1)
    tmp1$cnv <- cnv
    tmp1$loc <- "clean"
    glm_rocPlot <- rbind(glm_rocPlot,tmp1)
  }
}


full <- met_dat[grep(cnv,met_dat$chr_loc),]
cnv1u <- full[grep("1 CNV upstream",full$chr_loc),]
cnv2u <- full[grep("2 CNVs upstream",full$chr_loc),]
cnv1d <- full[grep("1 CNV downstream",full$chr_loc),]
cnv2d <- full[grep("2 CNVs downstream",full$chr_loc),]

full2 <- rbind(clean,cnv1u)
if (nrow(full2)==nrow(clean)){
  full2 <- data.frame()
}
full3 <- rbind(clean,cnv2u)
if (nrow(full3)==nrow(clean)){
  full3 <- data.frame()
}
full4 <- rbind(clean,cnv1d)
if (nrow(full4)==nrow(clean)){
  full4 <- data.frame()
}
full5 <- rbind(clean,cnv2d)
if (nrow(full5)==nrow(clean)){
  full5 <- data.frame()
}
full6 <- rbind(clean,cnv1u,cnv1d)
if (nrow(full6)==nrow(clean)){
  full6 <- data.frame()
}

dfName <- c("full","1 CNV upstream","2 CNVs upstream","1 CNV downstream","2 CNVs downstream",
            "CNV+1 CNV upstream","CNV+2 CNVs upstream","CNV+1 CNV downstream","CNV+2 CNVs downstream",
            "CNV+1 CNV up/downstream")
toRun <- list(full,cnv1u,cnv2u,cnv1d,cnv2d,
           full2,full3,full4,full5,full6)
# Transform data from long to wide format (one probe/column)

for (i in 1:length(toRun)){
  j=j+1
  current_name <- dfName[i]
  current <- toRun[[i]]
  
  if (nrow(current)==0){
    next
  }
  
  cnvLen_auc[j,"cnv"] <- cnv
  cnvLen_auc[j,"loc"] <- current_name
  glm_auc[j,"cnv"] <- cnv
  glm_auc[j,"loc"] <- current_name
  
  
  current_wide <- current %>% pivot_wider(names_from=probe_ID,values_from=beta_value,
                                          values_fn=list(beta_value=mean))
  current_wide <- current_wide[current_wide$sampleType=="01",]
  current_wide$Tissue <- ifelse(current_wide$Tissue=="TCGA-COAD","COAD","Normal")
  response <- "Tissue"
  predictors_current <- setdiff(names(current_wide),c(response,"case_gdc_id",
                          "sample_barcode","aliquot_id","sampleType","chr_loc"))
  if (length(predictors_current) > numProbe){
    predictors_current <- predictors_current[sample(length(predictors_current),numProbe)]
  } 
  current_wide[[response]] <- as.factor(current_wide[[response]])

  train1.hex<- as.h2o(current_wide, destination_frame = "train1.hex")
  plotNum <- sample(3:7,1)
  
  for (k in 3:7){
    
    model <- h2o.gbm(x=predictors_current,
                     y=response,
                     training_frame = train1.hex,
                     nfolds=10,
                      ntrees=100)
    
    model1 <- h2o.glm(x=predictors_current,
                      y=response,
                      training_frame = train1.hex,
                      nfold=10)
    # Record the AUC
    cnvLen_auc[j,k] <- h2o.auc(model, train=FALSE, xval=TRUE)
    glm_auc[j,k] <- h2o.auc(model1, xval=TRUE)
    
    if (k==plotNum){
      tpr <- h2o.performance(model,xval=T)@metrics$thresholds_and_metric_scores["tpr"]
      fpr <- h2o.performance(model,xval=T)@metrics$thresholds_and_metric_scores["fpr"]
      tmp <- cbind(tpr,fpr)
      tmp$cnv <- cnv
      tmp$loc <- current_name
      cnv_rocPlot <- rbind(cnv_rocPlot,tmp)
      
      tpr1 <- h2o.performance(model1,xval=T)@metrics$thresholds_and_metric_scores["tpr"]
      fpr1 <- h2o.performance(model1,xval=T)@metrics$thresholds_and_metric_scores["fpr"]
      tmp1 <- cbind(tpr1,fpr1)
      tmp1$cnv <- cnv
      tmp1$loc <- current_name
      glm_rocPlot <- rbind(glm_rocPlot,tmp1)
    }
  }
}

}

h2o.shutdown(prompt = F)



################# 
## t-test
#################
cnvList <- unique(combAUC$cnv)
combAUC$mean_auc <- round(rowMeans(combAUC[,3:7]),5)
sig_gbm <- as.data.frame(matrix(ncol=10))
sig_glm <- as.data.frame(matrix(ncol=10))

names(sig_gbm) <- c(names(combAUC))
names(sig_glm) <- c(names(combAUC))
sig_gbm$significance <- " "
sig_glm$significance <- " "

for (i in 1:length(cnvList)){
  curr_cnv <- cnvList[i]
  current_gbm <- combAUC %>% filter(cnv==curr_cnv, model=="GBM")
  current_glm <- combAUC %>% filter(cnv==curr_cnv, model=="GLM")
  baseline_gbm <- combAUC %>% filter(cnv==curr_cnv, loc=="clean", model=="GBM")
  baseline_glm <- combAUC %>% filter(cnv==curr_cnv, loc=="clean", model=="GLM")
  for (j in 1:nrow(current_gbm)){
    current_gbm$p_value[j] <- round(t.test(baseline_gbm[,c(3:7)],current_gbm[j,c(3:7)])$p.value,5)
    current_glm$p_value[j] <- round(t.test(baseline_glm[,c(3:7)],current_glm[j,c(3:7)])$p.value,5)
    current_gbm$significance <- ifelse(current_gbm$p_value<0.05 & 
                                         current_gbm$mean_auc>baseline_gbm$mean_auc,"+",
                                       current_gbm$significance)
    current_gbm$significance <- ifelse(current_gbm$p_value<0.05 & 
                                         current_gbm$mean_auc<baseline_gbm$mean_auc,"-",
                                       current_gbm$significance)
    current_glm$significance <- ifelse(current_glm$p_value<0.05 & 
                                         current_glm$mean_auc>baseline_glm$mean_auc,"+",
                                       current_glm$significance)
    current_glm$significance <- ifelse(current_glm$p_value<0.05 & 
                                         current_glm$mean_auc<baseline_glm$mean_auc,"-",
                                       current_glm$significance)
  }
  sig_gbm <- rbind(sig_gbm,current_gbm)
  sig_glm <- rbind(sig_glm,current_glm)
}


##########
## Differentiate COAD tumor tissue vs normal tissue
# coad_dat <- filter(met_dat,Tissue=="COAD")
# coad_dat$Tissue <- ifelse(coad_dat$sampleType==11,"Normal",coad_dat$Tissue)
# coad_dat1 <- coad_dat[,c(random_bv,which(colnames(coad_dat)=="Tissue"))]
# #coad_dat1 <- coad_dat
# response <- "Tissue"
# predictors_coad <- setdiff(names(met_dat1), c(response,"case_gdc_id",
#                            "sample_barcode","aliquot_id","sampleType"))
# coad_dat1[[response]] <- as.factor(coad_dat1[[response]])  
##########

# Start h2o
h2o.init(nthreads=6,min_mem_size = "4g")

# CHANGE FOR EVERY CNV
j=1
cnv_auc[j,1] <- geneName
# Pick random number to decide which iteration to save tpr,fpr
plotNum <- sample(2:6,1)
for (i in 2:6){
  
  ## COAD vs. other cancer types
  # Load data into h2o
  train1.hex<- as.h2o(clean_wide, destination_frame = "train1.hex")
  model1 <- h2o.gbm(x=predictors,
                   y=response,
                   training_frame = train1.hex,
                   nfolds=10,
                   ntrees=100)
  # Record the AUC
  h2o.auc(model1, train=FALSE, xval=TRUE)
  
  # Load data into h2o
  model2 <- h2o.glm(x=predictors,
                    y=response,
                    training_frame = train1.hex)
  # Record the AUC
  h2o.auc(model2)
  cnv_auc[j,i] <- h2o.auc(model1, train=FALSE, xval=TRUE)
 
  train_full.hex <- as.h2o(full2,destination_frame = "train_full.hex")
  model_full <- h2o.gbm(x=predictors_full2,
                    y=response,
                    training_frame = train_full.hex,
                    nfolds=10,ntrees=100)
  # Record the AUC
  h2o.auc(model_full,xval=TRUE)
  ## COAD normal vs. tumor tissue
  # train_coad.hex<- as.h2o(coad_dat1, destination_frame = "train_coad.hex")        
  # model_coad <- h2o.gbm(x=predictors_coad,
  #                   y=response,
  #                   training_frame = train_coad.hex,
  #                   nfolds=10,
  #                   ntrees=100)
  # # Record the AUC
  # coad_auc[j,i] <- h2o.auc(model_coad, train=FALSE, xval=TRUE)
  
  if (i==plotNum){
    tpr <- h2o.performance(model1,xval=T)@metrics$thresholds_and_metric_scores["tpr"]
    fpr <- h2o.performance(model1,xval=T)@metrics$thresholds_and_metric_scores["fpr"]
    tmp <- cbind(tpr,fpr)
    tmp$chrLoc <- geneName
    rocPlot <- rbind(rocPlot,tmp)

    # tpr_coad <- h2o.performance(model_coad,xval=T)@metrics$thresholds_and_metric_scores["tpr"]
    # fpr_coad <- h2o.performance(model_coad,xval=T)@metrics$thresholds_and_metric_scores["fpr"]
    # tmp1 <- cbind(tpr_coad,fpr_coad)
    # tmp1$chrLoc <-  geneName
    # coad_rocPlot <- rbind(coad_rocPlot,tmp1)
  }
  gc()
  h2o:::.h2o.garbageCollect()
  h2o:::.h2o.garbageCollect()
  h2o:::.h2o.garbageCollect()
}

saveRDS(cnv_auc,"cnv_auc.rds")
saveRDS(rocPlot,"cnv_rocPlot.rds")

h2o.shutdown()


# Plot upstream/downstream regions of a CNV
ggplot(cnvROC1,aes(fpr,tpr,colour=chrLoc))+geom_line() + 
 geom_segment(aes(x=0,y=0,xend = 1, yend = 1),linetype = 2,col='grey')+
  theme_bw()


saveRDS(gene_auc,file="cnv_up_downstream_auc.rds")
saveRDS(cnvROC,file="cnv_up_downstream_rocPlot.rds")
saveRDS(coad_auc,file="coadOnly_cnv_up_downstream_auc.rds")
saveRDS(coad_rocPlot,file="coad_Only_cnv_up_downstream_rocPlot.rds")
h2o.shutdown(prompt=FALSE)

#############
## Processing for DMR-overlapped CNV
# cnv_query <- paste("#standardSQL
#      with topCNVs as
#                      (SELECT
#                      Chr,
#                      Start_Pos,
#                      End_Pos,
#                      COUNT(*) AS total
#                      FROM
#                      `tcga-unmasked.cnvMet.cnv_DMR` cn
#                      GROUP BY
#                      Chr,
#                      Start_Pos,
#                      End_Pos
#                      ORDER BY
#                      total DESC
#                      limit 50)
#         SELECT
#                      Aliquot_ID AS Aliquot_ID,
#                      cn.Chr AS Chr,
#                      cn.Start_Pos AS Start_Pos,
#                      cn.End_Pos AS End_Pos,
#                      cn.Num_Probes AS Num_Probes,
#                      avg(Segment_Mean) AS Segment_Mean,
#                      cn.Cancer_Type AS Cancer_Type
#                      FROM
#                      `tcga-unmasked.TCGA_unmasked_copyNum.TCGA_CNV_unmasked` cn
#                      JOIN
#                      topCNVs top ON (top.Chr=cn.Chr and top.Start_Pos=cn.Start_Pos and top.End_Pos=cn.End_Pos)
#                      group by Aliquot_ID,cn.Chr,cn.Start_Pos,cn.End_Pos,Num_Probes,Cancer_Type")
# 
# cnv_gene <- query_exec(cnv_query,project,max_pages=Inf, use_legacy_sql=FALSE)
# cnv_gene1 <- unite(cnv_gene, "chromosome_loc", c("Chr","Start_Pos","End_Pos"))
# cnv_gene_wide <- cnv_gene1 %>% pivot_wider(names_from = chromosome_loc,values_from=Segment_Mean)
# 
# 
# # Isolate COAD samples
# cnv_gene_wide$Cancer_Type <- ifelse(cnv_gene_wide$Cancer_Type=="TCGA-COAD","COAD","Normal")
# response <- "Cancer_Type"
# predictors <- setdiff(names(cnv_gene_wide), c(response, "Aliquot_ID","Num_Probes"))
# cnv_gene_wide[[response]] <- as.factor(cnv_gene_wide[[response]])  
# 
# # Start h2o
# h2o.init(nthreads=4,min_mem_size = "4g")
# 
# # Load data into h2o
# train.hex<- as.h2o(cnv_gene_wide, destination_frame = "train.hex")        
# 
# model <- h2o.gbm(x=predictors,
#                  y=response,
#                  training_frame = train.hex,
#                  nfolds=10,
#                  ntrees=100)
# h2o.auc(model,xval=TRUE)
# tpr <- h2o.performance(model,xval=T)@metrics$thresholds_and_metric_scores["tpr"]
# fpr <- h2o.performance(model,xval=T)@metrics$thresholds_and_metric_scores["fpr"]
# tmp <- cbind(tpr,fpr)
# ggplot(tmp,aes(fpr,tpr))+geom_line() + 
#   geom_segment(aes(x=0,y=0,xend = 1, yend = 1),
#                linetype = 2,col='grey')+theme_bw()+
#   ggtitle(geneIn)
# h2o.shutdown(prompt=FALSE)
