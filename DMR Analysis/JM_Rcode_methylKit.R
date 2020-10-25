# Code for methylKit analysis in R (based off of ):
# 
# In Bash:
# module load R/3.6.0
# 
# In R:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
# Just need the line below
BiocManager::install("methylKit")

file.list=list("/dfs3/pub/jmorival/Zaragoza/RRBS_21020/SortedSAM/P02_sort.bam", "/dfs3/pub/jmorival/Zaragoza/RRBS_21020/SortedSAM/P01_sort.bam")
moyobj=processBismarkAln(file.list,sample.id=list("P1","C1"), assembly="hg19", treatment=c(1,0), read.context="CpG",mincov=5)
fMeth=unite(fMyObj,destrand=TRUE)

###############################################################
# Calculate methylation difference (DMRs)
# fmyDiff=calculateDiffMeth(fMeth)
# two groups detected:
#  will calculate methylation difference as the difference of
# treatment (group: 1) - control (group: 0)


# Find DMRs with 30% meth diff
fmyDiff30p=getMethylDiff(fmyDiff,difference=30,qvalue=0.01)
dd=as.data.frame(imyDiff30p)
write.table(dd,"/dfs3/pub/jmorival/Zaragoza/DMRs/iPSC_baseDMR.txt",sep="\t",quote=FALSE,row.names=FALSE)

# gene annotation
gene.obj=readTranscriptFeatures("/dfs3/pub/jmorival/Zaragoza/DMRs/refseq.hg19.bed.txt")
annotateWithGeneParts(as(imyDiff30p,"GRanges"),gene.obj)

# CpG annotation
cpg.obj=readFeatureFlank("/dfs3/pub/jmorival/Zaragoza/DMRs/CpGIslands.hg19.bed.txt",feature.flank.name=c("CpGi","shores"))
annotateWithFeatureFlank(as(imyDiff30p,"GRanges"),cpg.obj$CpGi,cpg.obj$shores,feature.name="CpGi",flank.name="shores")

# Find closest gene
idiffAnn=annotateWithGeneParts(as(imyDiff30p,"GRanges"),gene.obj)
head(getAssociationWithTSS(idiffAnn))
write.table(getAssociationWithTSS(idiffAnn),"/dfs3/pub/jmorival/Zaragoza/DMRs/iPSC_closestGene_hg19_5x.txt",sep="\t",quote=FALSE,row.names=FALSE)

# regional tiles
fTiles=tileMethylCounts(fMyObj,win.size=1000,step.size=1000,cov.bases=5)

# reassign treatments
getTreatment(iMyObj)=c(0,0,1,1,0,1,0,0,0,0)