if (!requireNamespace("bigrquery", quietly = TRUE))
    install.packages("bigrquery")
if (!requireNamespace("h2o", quietly = TRUE))
    install.packages("h2o")
if (!requireNamespace("dplyr", quietly = TRUE))
    install.packages("dplyr")
if (!requireNamespace("tidyr", quietly = TRUE))
    install.packages("tidyr")
if (!requireNamespace("ggplot2", quietly = TRUE))
    install.packages("ggplot2")
if (!requireNamespace("openssl", quietly = TRUE))
    install.packages("openssl")
if (!requireNamespace("data.table", quietly = TRUE))
  install.packages("data.table")
    
#
#
# This is the Google cloud project ID.  Change it to *your* project ID, or it won't work
# When this program makes the call to Bigquery to retrieve the TCGA data, it will prompt you to 
# login to your google account.

project<- "symbolic-rope-182022"

# Google charges for Bigquery access, but gives you 1TB of queries free per month.
# The query here will consume about 250 MB, or 1/4000 th of your monthly quota.



library(dplyr)
library(h2o)
library(bigrquery)
library(ggplot2)
library(tidyr)
library(openssl)
library(data.table)

###Get Methylation Data ###
BigQuerySQL_start_end <- paste(
"SELECT bs.case_barcode, methyl.sample_barcode, bs.sample_type, methyl.probe_id, methyl.beta_value, annotate.CpG_island_coord, annotate.chromosome, annotate.position, annotate.CGI_feature_type, clin.gender, clin.ethnicity
FROM
  `isb-cgc.TCGA_bioclin_v0.Biospecimen` bs
JOIN
  `isb-cgc.TCGA_hg38_data_v0.DNA_Methylation_chr1` methyl
ON
  bs.sample_barcode = methyl.sample_barcode
JOIN
  `isb-cgc.TCGA_bioclin_v0.Clinical` clin
ON
  bs.case_barcode = clin.case_barcode
JOIN
  `isb-cgc.platform_reference.GDC_hg38_methylation_annotation` annotate
ON 
  methyl.probe_id = annotate.CpG_probe_id
WHERE
  bs.project_short_name ='TCGA-COAD'"                         
)

#Get the data
TCGA_DNA_Methylation <- bq_project_query(project, BigQuerySQL_start_end)
data <- bq_table_download(TCGA_DNA_Methylation, page_size = 100000000)
system.time(
  fwrite(data, "/home/data/Shared/TCGA/TCGA_COAD_Methylation_Chr1.csv", sep = ",", col.names = TRUE)
  )



###Get CNV Data ###

BigQuerySQL_start_end_cnv <- paste(
  "SELECT bs.case_barcode, cnv.sample_barcode, bs.sample_type, cnv.chromosome, cnv.start_pos, cnv.end_pos, cnv.segment_mean
FROM
  `isb-cgc.TCGA_bioclin_v0.Biospecimen` bs
JOIN
  `isb-cgc.TCGA_hg38_data_v0.Copy_Number_Segment_Masked_r14` cnv
ON
  bs.sample_barcode = cnv.sample_barcode
WHERE
  bs.project_short_name ='TCGA-COAD'"                         
)


#Get the data
TCGA_DNA_CNV <- bq_project_query(project, BigQuerySQL_start_end_cnv)
data <- bq_table_download(TCGA_DNA_CNV, page_size = 100000000)
system.time(
  fwrite(data, "/home/data/Shared/TCGA/TCGA_COAD_CNV.csv", sep = ",", col.names = TRUE)
)

rm(data)