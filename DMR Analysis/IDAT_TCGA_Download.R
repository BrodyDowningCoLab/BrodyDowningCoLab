# Using sesame  http://bioconductor.org/packages/sesame/
# Please cite 10.1093/nar/gky691 and doi: 10.1093/nar/gkt090.
library(TCGAbiolinks)

# bars=file.list
# 
# for(i in 1:length(file.list)){
#   #retrieve the sample ID: separate by "_" to get sample ID, then separate by "-" to get sample type, then substring to get rid of letter at the end
#   bars[i] = unlist(strsplit(file.list[i], "_"))[2]
# }

#-------------------------------------------------------
# Example to idat files from TCGA projects
#-------------------------------------------------------
projects <- TCGAbiolinks:::getGDCprojects()$project_id
projects <- "TCGA-COAD"
match.file.cases.all <- NULL
for(proj in projects){
  print(proj)
  query <- GDCquery(project = proj,
                    data.category = "Raw microarray data",
                    data.type = "Raw intensities", 
                    experimental.strategy = "Methylation array",
                    legacy = TRUE,
                    file.type = ".idat",
                    platform = "Illumina Human Methylation 450")
  match.file.cases <- getResults(query,cols=c("cases","file_name"))
  match.file.cases$project <- proj
  match.file.cases.all <- rbind(match.file.cases.all,match.file.cases)
  tryCatch(GDCdownload(query, method = "api", files.per.chunk = 20),
           error = function(e) GDCdownload(query, method = "client"))
}
# This will create a map between idat file name, cases (barcode) and project
readr::write_tsv(match.file.cases.all, file =  "idat_filename_case.txt")
# code to move all files to local folder
for(file in dir(".",pattern = ".idat", recursive = T)){
  TCGAbiolinks:::move(file,basename(file))
}
