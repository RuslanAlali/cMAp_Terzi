#cMAP
setwd("/home/ros/Desktop/cMAP_data")
library(DESeq2)
library(pheatmap)
library(dplyr)

#----------------------------

# #The aim is to compare sham with Nx(kidney damaged)
# the condition is tgfa(m) is better than tgfa(p) when the kidney is damaged



# DESeq2

#load data
RNAseq=read.table("DataCount.txt",header =  T,row.names = 1)

#prepare metadata
m_data=colnames(RNAseq)
meta_data <- as.vector(factor(grepl("sham",m_data,ignore.case = T),labels =  c("Nx","Sham")))
meta_data <- rbind( meta_data, as.vector(factor(grepl("EGF",m_data,ignore.case = T),labels =  c("Tgfa","Egf"))))
meta_data <- rbind( meta_data, as.vector(factor(grepl("_p_",m_data,ignore.case = T),labels =  c("mut","parent"))))
colnames(meta_data)=m_data
row.names(meta_data)=c("condition","GFactor","Genotype")
meta_data=t(meta_data)
meta_data<-as.data.frame(meta_data)
remove(m_data)   

#----------------------

m_data=colnames(RNAseq)
egf_RNA= RNAseq[grepl("egf",m_data,ignore.case = T)] 
tgfa_RNA= RNAseq[grepl("tgfa",m_data,ignore.case = T)] 

egf_dds=DESeqDataSetFromMatrix(egf_RNA, colData=meta_data[meta_data$GFactor=="Egf",], design = ~Genotype)
egf_dds<- vst(egf_dds, blind = T)
v_assay<-assay(egf_dds)
v_assay<-cor(v_assay)
pheatmap(v_assay,annotation =  meta_data)

tgfa_dds=DESeqDataSetFromMatrix(tgfa_RNA, colData=meta_data[meta_data$GFactor!="Egf",], design = ~Genotype)
tgfa_dds<- vst(tgfa_dds, blind = T)
v_assay<-assay(tgfa_dds)
v_assay<-cor(v_assay)
pheatmap(v_assay,annotation =  meta_data)



#################### EGF(p)
m_data=colnames(egf_RNA)
EGF_Nx=egf_RNA[grepl("Nx",m_data,ignore.case = T)] 
EGF_Nx=DESeqDataSetFromMatrix(EGF_Nx, colData=meta_data[meta_data$GFactor=="Egf"& meta_data$condition=="Nx",], design = ~Genotype)
EGF_Nx<- vst(EGF_Nx, blind = T)
v_assay<-assay(EGF_Nx)
v_assay<-cor(v_assay)
pheatmap(v_assay,annotation =  meta_data)

m_data=colnames(egf_RNA)
EGF_sham=egf_RNA[grepl("sham",m_data,ignore.case = T)] 
EGF_sham=DESeqDataSetFromMatrix(EGF_sham, colData=meta_data[meta_data$GFactor=="Egf"& meta_data$condition!="Nx",], design = ~Genotype)
EGF_sham<- vst(EGF_sham, blind = T)
v_assay<-assay(EGF_sham)
v_assay<-cor(v_assay)
pheatmap(v_assay,annotation =  meta_data)


plotPCA(EGF_Nx,intgroup="Genotype")
plotPCA(EGF_sham,intgroup="Genotype")


#################### Tgfa(p)
m_data=colnames(tgfa_RNA)
Tgfa_Nx=tgfa_RNA[grepl("Nx",m_data,ignore.case = T)] 
Tgfa_Nx=DESeqDataSetFromMatrix(Tgfa_Nx, colData=meta_data[meta_data$GFactor!="Egf"& meta_data$condition=="Nx",], design = ~Genotype)
Tgfa_Nx<- vst(Tgfa_Nx, blind = T)
v_assay<-assay(Tgfa_Nx)
v_assay<-cor(v_assay)
pheatmap(v_assay,annotation =  meta_data)

m_data=colnames(tgfa_RNA)
Tgfa_sham=tgfa_RNA[grepl("sham",m_data,ignore.case = T)] 
Tgfa_sham=DESeqDataSetFromMatrix(Tgfa_sham, colData=meta_data[meta_data$GFactor!="Egf"& meta_data$condition!="Nx",], design = ~Genotype)
Tgfa_sham<- vst(Tgfa_sham, blind = T)
v_assay<-assay(Tgfa_sham)
v_assay<-cor(v_assay)
pheatmap(v_assay,annotation =  meta_data)


plotPCA(Tgfa_Nx,intgroup="Genotype")
plotPCA(Tgfa_sham,intgroup="Genotype")


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@
#@       run DESeq2 wrong
#@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
m_data=colnames(egf_RNA)

EGF_Nx=egf_RNA[grepl("Nx",m_data,ignore.case = T)] 
EGF_Nx=DESeqDataSetFromMatrix(EGF_Nx, colData=meta_data[meta_data$GFactor=="Egf"& meta_data$condition=="Nx",], design = ~Genotype)

EGF_sham=egf_RNA[grepl("sham",m_data,ignore.case = T)] 
EGF_sham=DESeqDataSetFromMatrix(EGF_sham, colData=meta_data[meta_data$GFactor=="Egf"& meta_data$condition!="Nx",], design = ~Genotype)

m_data=colnames(tgfa_RNA)

Tgfa_Nx=tgfa_RNA[grepl("Nx",m_data,ignore.case = T)] 
Tgfa_Nx=DESeqDataSetFromMatrix(Tgfa_Nx, colData=meta_data[meta_data$GFactor!="Egf"& meta_data$condition=="Nx",], design = ~Genotype)

Tgfa_sham=tgfa_RNA[grepl("sham",m_data,ignore.case = T)] 
Tgfa_sham=DESeqDataSetFromMatrix(Tgfa_sham, colData=meta_data[meta_data$GFactor!="Egf"& meta_data$condition!="Nx",], design = ~Genotype)


EGF_Nx<-DESeq(EGF_Nx)
Tgfa_Nx<-DESeq(Tgfa_Nx)
EGF_sham<-DESeq(EGF_sham)
Tgfa_sham<-DESeq(Tgfa_sham)

plotDispEsts(EGF_Nx)
plotDispEsts(EGF_sham)
plotDispEsts(Tgfa_Nx)
plotDispEsts(Tgfa_sham)


resEGF_Nx<-results(EGF_Nx, alpha=0.05)
resTgfa_Nx<-results(Tgfa_Nx, alpha=0.05)
resEGF_sham<-results(EGF_sham, alpha=0.05)
resTgfa_sham<-results(Tgfa_sham, alpha=0.05)


resEGF_Nx1<-lfcShrink(EGF_Nx,contrast = c("Genotype","parent","mut"),res = resEGF_Nx)
resTgfa_Nx1<-lfcShrink(Tgfa_Nx,contrast = c("Genotype","parent","mut"),res = resTgfa_Nx)
resEGF_sham1<-lfcShrink(EGF_sham,contrast = c("Genotype","parent","mut"),res = resEGF_sham)
resTgfa_sham1<-lfcShrink(Tgfa_sham,contrast = c("Genotype","parent","mut"),res = resTgfa_sham)

summary(resEGF_Nx1)
summary(resTgfa_Nx1)
summary(resEGF_sham1)
summary(resTgfa_sham1)


temp=resEGF_Nx1
colnames(temp)=paste0("EGF_Nx_", colnames( resEGF_Nx1))
total_res=temp
temp=resTgfa_Nx1
colnames(temp)=paste0("Tgfa_Nx_", colnames( resTgfa_Nx1))
total_res= cbind(total_res,temp)
temp=resEGF_sham1
colnames(temp)=paste0("EGF_sham_", colnames( resEGF_sham1))
total_res= cbind(total_res,temp)
temp=resTgfa_sham1
colnames(temp)=paste0("Tgfa_sham_", colnames( resTgfa_sham1))
total_res= cbind(total_res,temp)

write.csv(total_res,file = "total_res_mouse.csv")


#---- try two factors
Tgfa_=tgfa_RNA
Tgfa_=DESeqDataSetFromMatrix(Tgfa_, colData=meta_data[meta_data$GFactor!="Egf",], design = ~condition+Genotype)
Tgfa_<-DESeq(Tgfa_)
resTgfa_<-results(Tgfa_, alpha=0.05)
resTgfa_1<-lfcShrink(Tgfa_,contrast = c("Genotype","parent","mut"),res = resTgfa_)
summary(resTgfa_1)

write.csv(resTgfa_1,file = "resTgfa_1_mouse.csv")


#------ read data

resTgfa_1<-read.csv("resTgfa_1_mouse.csv",row.names = 1)
total_res<-read.csv("total_res_mouse.csv",row.names = 1)

temp_tgfa<-resTgfa_1[which(resTgfa_1$padj<0.05),]  
temp_total<-total_res[which(total_res$Tgfa_Nx_padj<0.05),]  

library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(2, "Pastel2")
set1 <- rownames(temp_tgfa)
set2 <- rownames(temp_total)

venn.diagram(
  x = list(set1, set2),
  category.names = c("Tgfa_all" , "Tgfa_Nx"),
  filename = 'venn_diagramm.tiff',
  output=TRUE,
  fill = myCol[1:2]
  
)


#######################3
#Convert Ensembl gene names into Gene Symbols using annotables package.
#
######################

#install.packages("devtools")
#devtools::install_github("stephenturner/annotables")
#install.packages("tibble",dependencies = TRUE, repos = "http://cran.us.r-project.org")

library(annotables)
library(tibble)
library(dplyr)

wt_res<-data.frame(temp_total)%>%
  rownames_to_column(var="ensgene")
wt_res<-  left_join(x=wt_res, y=grcm38[,c("ensgene","symbol","description")],by="ensgene")
wt_res<- wt_res[!is.na(wt_res$symbol),]
wt_res<- wt_res[!duplicated(wt_res$symbol),]
View(wt_res)

res_sig<-subset(wt_res,Tgfa_Nx_padj<0.05)
res_sig<-res_sig%>%arrange(Tgfa_Nx_padj)
View(res_sig)


