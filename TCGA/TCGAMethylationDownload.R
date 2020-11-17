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
    
#
#
# This is the Google cloud project ID.  Change it to *your* project ID, or it won't work
# When this program makes the call to Bigquery to retrieve the TCGA data, it will prompt you to 
# login to your google account.

project<- "symbolic-rope-182022"

# Google charges for Bigquery access, but gives you 1TB of queries free per month.
# The query here will consume about 250 MB, or 1/4000 th of your monthly quota.



library("dplyr")
library("h2o")
library("bigrquery")
library("ggplot2")
library("tidyr")
library("openssl")

BigQuerySQL_start_end <- paste(
"SELECT methyl.sample_barcode, beta_value",
"FROM",
  "`isb-cgc.TCGA_hg38_data_v0.DNA_Methylation_chr21` methyl",
"JOIN",
  "`isb-cgc.TCGA_bioclin_v0.Biospecimen` bs",
"ON",
  "methyl.sample_barcode = bs.sample_barcode",
"WHERE",
  "bs.project_short_name ='TCGA-COAD'"                         
)


#Get the data
TCGA_DNA_Methylation <- bq_project_query(project, BigQuerySQL_start_end)
data <- bq_table_download(TCGA_DNA_Methylation)