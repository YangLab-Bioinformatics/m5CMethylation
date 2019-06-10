rm(list = ls())
setwd("/pnas/yangyg_group/hanyn/2-patient-bs/Tissue_all/Diff")
mc   <-read.table(file="2-Diff_m5C_all.txt",header=T,stringsAsFactors=F)
mc1  <-read.table(file="2-Diff_m5C.txt",header=T,stringsAsFactors=F)
mcup  <-read.table(file="2-Diff_m5C_up_rmintron.bed",header=T,stringsAsFactors=F)
colnames(mcup)<-c("chrom","start","end","strand","site" ,"Gene" ,"diff","pvalue","FDR","pair_num","pair_up","pair_dn","mean")
cor  <-read.table(file="3-cor_positive_737.txt",header=T,stringsAsFactors=F)

rpkm <-read.table(file="/asnas/yangyg_group/hanyn/Tissue-RNA/RNA-Seq-CP/NGT-20180627/DEG/DESeq2_N29T36/2-N29T36-rpkm.txt",header=T,stringsAsFactors=F)
expu <-read.table(file="/asnas/yangyg_group/hanyn/Tissue-RNA/RNA-Seq-CP/NGT-20180627/DEG/DESeq2_N29T36/DESeq2/TvsN-up-fc1.2-padj05.txt",header=F,stringsAsFactors=F)[,c(1,3,7)]
colnames(expu) <- c("Gene","log2FC","FDR")

name <-read.table(file="/pnas/yangyg_group/hanyn/data/reference/Locus_tag.txt",header=T,stringsAsFactors=F)
onco <- read.csv(file="../Diff_3036_v0/uniprot-oncogene.csv",header=T,stringsAsFactors=F)
tcga <-read.table(file="../Diff_3036_v0/TCGA_SurviveCandidate.txt",header=F,stringsAsFactors=F)
library(gdata)
marker <-read.xls("marker_all.xlsx",sheet="Sheet1",stringsAsFactors=F,header=F)
nor=29
tum=36
sam <- 65
sample <- c("N0156","N0453","N0466",
            "N0474","N0531","N0555","N0565","N0897","N1037","N1063",
            "N1067","N1495","N1664","N1691","N2224","N2297","N3111",
            "N3138","N3740","N3880","P2347","P2688","P2759","P2804",
            "P3143","P3451","P3467",        "P3522","P3635","T0156",
            "T0453","T0466","T0474","T0531","T0555","T0565","T0982",
            "T1037","T1067","T1495","T1664","T2224","T2226","T2271",
            "T2297","T2347","T2644","T2759","T2804","T3138","T3143",
            "T3196","T3282","T3357","T3371","T3403","T3451","T3515",
            "T3560","T3614","T3635","T3660","T3740","T3750","T3880")

##################################################################################### 1. diff m5C     #####
rm(list = ls())
setwd("/pnas/yangyg_group/hanyn/2-patient-bs/Tissue_all/")
site <-read.table(file="summary_m5C_gene.txt",header=F,stringsAsFactors=F)
setwd("/pnas/yangyg_group/hanyn/2-patient-bs/Tissue_all/Diff")  
nor <- 29
tum <- 36
site$covn <- apply(site[,c(seq(6,nor*3+4,by=3))],1,function(x){sum(!is.na(x)&x>=10)})
site$covt <- apply(site[,c(seq(6+nor*3,4+nor*3+tum*3,by=3))],1,function(x){sum(!is.na(x)&x>=10)})
site0 <- site[(site$covn+site$covt)>=30&site$covn>=10&site$covt>=10,]
site1 <- site0[,c(1:4,seq(7,4+nor*3+tum*3,by=3),4+nor*3+tum*3+6)]
site1$meann <- apply(site1[,5:(4+nor)],1,function(x){mean(x[!is.na(x)])})
site1$meant <- apply(site1[,(5+nor):(4+nor+tum)],1,function(x){mean(x[!is.na(x)])})
site1$diff <- site1$meant-site1$meann
site2 <- site1[abs(site1$diff)>=0.05,]
site2$pvalue <- as.numeric(apply(site2[,5:(4+nor+tum)],1,function(x){wilcox.test(x[1:nor],x[(1+nor):(nor+tum)],paired=F)$p.value}))
library(qvalue)
site2$qvalue <- qvalue(site2$pvalue)$qvalue
site3 <- site2[site2$pvalue < 0.05,]
head(site3)

colnames(site3) <- c("chrom","start","end","strand",sample,"Gene","meann","meant","diff","pvalue","FDR")
head(site3[site3$Gene=="ENSG00000143321",],n=10)
dim(site)  #1425509
dim(site0) #1005357
dim(site2) #  48636
dim(site3) #  12085   ??????
length(unique(site3$Gene))# 
length(unique(site3$Gene[site3$diff>0]))# 2979
length(unique(site3$Gene[site3$diff<0]))#  915
nrow(unique(site3[site3$diff>0,1:4]))# 8346
nrow(unique(site3[site3$diff<0,1:4]))# 2123 
## paired 
head(site3)
site3$site <- paste(site3$chrom,site3$start,sep="_")
pair <-read.table(file="1-pair_555.txt",header=T,stringsAsFactors=F)
dim(site3)# 
dim(pair) # 
sitep <- merge(site3,pair,by.x="site",by.y="site",all.x=F,all.y=F, suffixes = c(".x",".y"))
dim(sitep)# 
length(unique(sitep$Gene))# 
length(unique(sitep$Gene[sitep$meant-sitep$meann>0]))#  
length(unique(sitep$Gene[sitep$meant-sitep$meann<0]))# 
nrow(unique(sitep[sitep$meant-sitep$meann>0,1:5]))#  
nrow(unique(sitep[sitep$meant-sitep$meann<0,1:5]))# 1154 site 504 gene
head(sitep[sitep$Gene=="ENSG00000143321",],n=10) 
write.table(sitep,"2-Diff_m5C_all.txt",quote=F,sep="\t",row.names = F)
write.table(unique(sitep[,c(1,(5+nor+tum+1),(5+nor+tum+4):ncol(sitep))]),"2-Diff_m5C.txt",quote=F,sep="\t",row.names = F)
write.table(unique(sitep[sitep$diff>0,c(2:5,1,(5+nor+tum+1),(5+nor+tum+4):ncol(sitep))]),"2-Diff_m5C_up.bed",quote=F,sep="\t",row.names = F,col.names=F)
write.table(unique(sitep[sitep$diff<0,c(71)]),"2-Diff_m5C_down_genelist.txt",quote=F,sep="\t",row.names = F,col.names=F)

##################################################################################### 1.1 heatmap     ###### 
setwd("/pnas/yangyg_group/hanyn/2-patient-bs/Tissue_all/Diff")
mc   <-read.table(file="2-Diff_m5C_all.txt",header=T,stringsAsFactors=F)
sitep <- mc[,-1]
## heatmap
diff <- unique(sitep[,c(1:(4+nor+tum),(8+nor+tum))])
diff1 <- diff[order(diff$diff,decreasing=T),5:(4+nor+tum)]
# pdf("FigureS1G_diff-m5C-heatmap.pdf",width=6,height=7)
library(pheatmap)
annotation_col = data.frame(Sample =c(rep("Normal",nor),rep("Tumor",tum)))
rownames(annotation_col)=colnames(diff1)
pheatmap(as.matrix(diff1), 
         color = colorRampPalette(c("blue","white","red"))(n = 500),
         cellwidth = 4,cellheight =0.06,
         show_rownames=F,border_color="gray",font_size=15,show_colnames=F,
         cluster_cols=F,scale="row",cluster_rows=F,breaks= unique(c(seq(-2,2, length=500))),
         legend = T,annotation_col=annotation_col,
         annotation_colors=list(Sample = c(Normal = "darkgreen", Tumor="orange"))[1],
         gaps_row = nrow(diff[diff$diff>0,]),gaps_col=nor
)
dev.off()
library(ggplot2)
diff2 <- data.frame(diff=sort(diff$diff),fill=0)
diff2$fill[diff2$diff>0] <- "Hyper-methylated m5C"
diff2$fill[diff2$diff<0] <- "Hypo-methylated m5C"
sum(diff2$diff>0)#4126
sum(diff2$diff<0)#1154
diff2$perc <- ((1:5280)/5280-(1154/5280))*100
pdf("FigureS1_diff-m5C-bar_sort-0805-hyn.pdf",width=7,height=3)
ggplot(data=diff2)+geom_bar(aes(x=perc,y=diff,fill=fill),stat="identity")+theme_classic()+
  scale_fill_manual(values=c("red","blue"))+labs(x="Percentage (%)",y="m5C level difference",fill="")+
  annotate("text",x=-20,y=-0.3,label="21.86%",vjust=0,hjust=0,size=3)+
  annotate("text",x=63,y=0.35,label="78.14%",vjust=0,hjust=0,size=3)+
  theme_bw()+
  theme(panel.grid.minor.x=element_blank()) +
  theme(panel.grid.major.x=element_blank())+
  theme(panel.grid.minor.y=element_blank()) +
  theme(panel.grid.major.y=element_blank())
dev.off()
## 2019 rebuttal heatmap
setwd("/pnas/yangyg_group/hanyn/2-patient-bs/Tissue_all/Diff")
mc   <-read.table(file="2-Diff_m5C_all.txt",header=T,stringsAsFactors=F)
gene   <-read.table(file="pathway/Figure1a_pre_gene.txt",header=T,stringsAsFactors=F)
head(mc)
pre <-merge(mc,gene,by.x="Gene",by.y="gene.y",all.x=F,all.y=F, suffixes = c(".x",".y"))
#285

length(unique(gene$gene)) # one gene without ensemble gene id

length(unique(pre$Gene))# 
pre1 = pre[pre$diff>0,]
nrow(pre1)# 
length(unique(pre1$Gene))# 
nrow(pre1[,7:71])
nrow(unique(pre1[,7:71]))
pdf("pre3.pdf",width=6,height=7)
library(pheatmap)
#annotation_col = data.frame(Sample =c(rep("Normal",nor),rep("Tumor",tum)))
#rownames(annotation_col)=colnames(diff1)
pheatmap(as.matrix(unique(pre1[,7:71])), 
         color = colorRampPalette(c("blue","white","red"))(n = 500),
   #      cellwidth = 4,cellheight =0.06,
         show_rownames=F,border_color="gray",font_size=15,show_colnames=F,
         cluster_cols=F,scale="row",cluster_rows=F,
        #breaks= unique(c(seq(-2,2, length=500))),
        # legend = T#,#annotation_col=annotation_col,
        # annotation_colors=list(Sample = c(Normal = "darkgreen", Tumor="orange"))[1],
         gaps_col=29
)
dev.off()


##
if (!require("devtools")) {
  install.packages("devtools", dependencies = TRUE)
  library(devtools)
}
install_github("raivokolde/pheatmap")
pheatmap(..., na_col = "grey", ...)
##
####################################################################################  2. ten pathways heatmap boxplot #####
setwd("/pnas/yangyg_group/hanyn/2-patient-bs/Tissue_all/Diff/pathway")
mc   <-read.table(file="../2-Diff_m5C_all.txt",header=T,stringsAsFactors=F)
name <-read.table(file="/pnas/yangyg_group/hanyn/data/reference/Locus_tag.txt",header=T,stringsAsFactors=F)
path <-read.table(file="path_0730.txt",header=T,stringsAsFactors=F)
mcup <-read.table(file="../2-Diff_m5C_up_rmintron.bed",header=T,stringsAsFactors=F)
colnames(mcup)<-c("chrom","start","end","strand","site" ,"Gene" ,"diff","pvalue","FDR","pair_num","pair_up","pair_dn","mean")
head(mc)
head(name)
head(path)
head(mcup)
sam <- 65
mer1 <-merge(mcup,name,by.x="Gene",by.y="gene",all.x=F,all.y=F, suffixes = c(".x",".y"))
mer2 <-merge(mer1,path,by.x="locus_tag",by.y="Gene",all.x=F,all.y=F, suffixes = c(".x",".y"))
mer3 <-unique(merge(mer2[,c(7,1,2,15)],mc[,c(1,(9+sam),(10+sam),(12+sam),(13+sam),6:(5+sam))],by.x="site",by.y="site",all.x=F,all.y=F, suffixes = c(".x",".y")))
length(unique(mer3$locus_tag))
#write.csv(mer3,"path_53genes.csv",quote=F,row.names = F)

library(pheatmap)
heat<- read.csv("path_53genes_filtered.csv",stringsAsFactors=F)  
nor=29
tum=36
sam <- 65
rownames(heat) <- heat$locus_tag
heat2 <- heat[order(heat$Path),]
heat3 <- heat2[,c(9:(8+sam))]
(num_path <- as.numeric(tapply(heat2$Path,heat2$Path,length)))
#4  5  9 14 11  7
gapsrow <- c(4,9,18,32,43)##!!!!!!!!!!!!!!!!!!!!!  gai
annotation_col = data.frame(Sample =c(rep("Normal",nor),rep("Tumor",tum)))
rownames(annotation_col)=colnames(heat3)
annotation_row = data.frame(Path   =c(rep("JAK/Stat",  num_path[1]),
                                      rep("VEGF",    num_path[2]),
                                      rep("ERBB",  num_path[3]),
                                      rep("PI3K/AKT",      num_path[4]),
                                      rep("EMT",       num_path[5]),
                                      rep("ERK/MAPK",  num_path[6])
                                      ))
rownames(annotation_row)=rownames(heat3)
library(pheatmap)
#pdf("Figure1A_path5_0807_hyn.pdf",width=16,height=10)
pheatmap(as.matrix(heat3[,c(1:29,30,31,32,65,33,34,64,35,36,63,37,38,62,39,40,61,41,51,42,60,43,44,
                            59,45,50,46,58,47,48,57,49,52,56,55,53,54)]), 
         color = colorRampPalette(c("darkblue","white","red"))(n = 500),
         cellwidth = 10,cellheight =9,
         show_rownames=T,border_color="gray",font_size=15,show_colnames=F,
         cluster_cols=F,scale="row",cluster_rows=F,gaps_col = nor,gaps_row=gapsrow,
         legend = T,
         annotation_col=annotation_col,
         annotation_colors=list(Sample = c(Normal = "darkblue", Tumor="yellow")),
         annotation_row=annotation_row,
         breaks= unique(c(seq(-3,3, length=500)))
)
dev.off() 

### test genes 
heat$medien_diff <- apply(heat[,-1:-8],1,function(x){median(x[30:65],na.rm=T)-median(x[1:29],na.rm=T)})
heat9<- heat[,c(4:8,ncol(heat))]
heat9["AKT2",]
heat9["LAMC2",]
heat9["PDPK1",]
heat_md <- heat[heat$medien_diff>0.1,]
write.csv(heat_md,"Figure1B_mediandiff0.1.csv",quote=F,row.names = F)
heat_md[order(heat_md$Path),c(4:6,ncol(heat_md))]
heat[heat$medien_diff>0.1 & heat$Path=="5EMT",c(4:8,74)]
###
#######  boxplot 
mygene <- c("IL6","TYK2",
            "TBX1","PTK2",
            "PRKD1","PRKCZ","CHAD","RPS6KB2","ZDHHC8",
           "DVL2","HDGF","GRB2")
heat2[mygene,]$pvalue  ## check the star number
# [1] 6.450000e-05 8.584725e-03 3.580000e-05 1.443009e-02 2.040000e-05
# [6] 7.230000e-06 4.098510e-04 2.560000e-06 4.120000e-09 9.050000e-06
# [11] 3.630000e-07 6.270000e-08

#### 20190510###
mean(heat2[mygene,]$diff) # 0.1598568
write.table(heat2[mygene,],file="figure1Btable.txt",sep="\t",quote=FALSE);
######
onco <- heat2[mygene,c(2,9:(8+sam))]
gene_order <- mygene
oncos <- onco[gene_order,]
oncos$locus_tag=factor(oncos$locus_tag, levels=gene_order)
library(reshape2) 
onco1 <- melt(oncos) 
gnum <- length(mygene)
onco1$type <- c(rep("Normal",nor*gnum),rep("Tumor",tum*gnum))
library(ggplot2)
onco1$x <- c(rep(1:gnum-0.2,nor),rep(1:gnum+0.2,tum))
max<-tapply(onco1$value,onco1$locus_tag,function(x){max(x,na.rm=T)})
pdf("Figure1B_boxplot5_0808-hyn.pdf",width=10,height=3)
ggplot(onco1)+
  geom_boxplot(aes(x=locus_tag,y=value,fill=type),position=position_dodge(width=0.8),outlier.shape = NA)+
  geom_jitter(aes(x=x,y=value),size=0.4,width=0.07)+
  theme_classic()+ylim(c(0,1.2))+
  scale_fill_manual(values=c("blue","red"))+labs(x="",y="m5C level",fill="")+
 # theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))+
 annotate("text",x=1:gnum,y=max+0.06,label=c("***","**","***","*",rep("***",8)))+
 annotate("segment", x = 1:gnum-0.2, xend = 1:gnum+0.2, y = max+0.055, yend = max+0.055)+
 annotate("segment", x = c(1,3,5,7,10,11)-0.2, xend = c(2,4,6,9,10,12)+0.2, y = 1.07, yend = 1.07)+
 annotate("text",x = c(1.5,3.5,5.5,8,10,11.5),y=1.14,label=c("JAK/STAT","VEGF","ERBB","PI3K/AKT","EMT","ERK/MAPK"))+
  scale_y_continuous(breaks=seq(0,1,by=0.2))
dev.off()
##################################################################################### 4. early late  ########
rm(list = ls())
setwd("/pnas/yangyg_group/hanyn/2-patient-bs/Tissue_all/Diff")
mc   <-read.table(file="2-Diff_m5C_all.txt",header=T,stringsAsFactors=F)
sitep <- mc[(mc$diff)>0,-1]
#### diff #### 
tumo <- data.frame(
  sitep$chrom,sitep$start,sitep$end,sitep$strand,    
  sitep$T0466,  sitep$T0555,   sitep$T2226,  sitep$T2347,  sitep$T2804,  sitep$T3196,
  sitep$T3282,  sitep$T3403,  sitep$T3451,  sitep$T3560,  sitep$T3614,  sitep$T3660,
  sitep$T0531,  sitep$T0565,  sitep$T2224,  sitep$T2297,  sitep$T3515,  sitep$T3740,
  sitep$T3880,  sitep$T0453,  sitep$T0474,  sitep$T2271,  sitep$T2644,  sitep$T3138,
  sitep$T3357,  sitep$T3371,  sitep$T3635,  sitep$T0156,  sitep$T0982,  sitep$T1037,
  sitep$T1067,  sitep$T1495,  sitep$T1664,  sitep$T2759,  sitep$T3143,  sitep$T3750,  sitep$Gene)
earl=12
late=24
## diff 
tumo$meann <- apply(tumo[,5:(4+earl)],1,function(x) {mean(x[!is.na(x)])})
tumo$meant <- apply(tumo[,(5+earl):(4+earl+late)],1,function(x){mean(x[!is.na(x)])})
tumo2 <- tumo[!is.na(tumo$meann) & (tumo$meant-tumo$meann)>=0.05,]## 
tumo2$pvalue <- as.numeric(apply(tumo2[,5:(4+earl+late)],1,function(x){wilcox.test(x[1:earl],x[(1+earl):(earl+late)],paired=F)$p.value}))
tumo3 <- tumo2[tumo2$pvalue<0.05,]
#write.table(tumo3,"Stage/Stage_diff_all.txt",quote=F,sep="\t",row.names = F)#208

nrow(unique(tumo3[,1:2]))
length(unique(tumo3$sitep.Gene))
#### count
head(tumo3)
head(tumo3[tumo3$sitep.Gene=="ENSG00000143321",],n=10)# yes
length(unique(tumo3$sitep.Gene)) 
length(unique(tumo3$sitep.Gene[tumo3$meant-tumo3$meann>0])) 
length(unique(tumo3$sitep.Gene[tumo3$meant-tumo3$meann<0])) 
nrow(unique(tumo3[(tumo3$meant-tumo3$meann)>0,1:4]))#  
nrow(unique(tumo3[(tumo3$meant-tumo3$meann)<0,1:4]))#  

#### heatmap all #### 
matrix <- tumo3[order((tumo3$meant-tumo3$meann),decreasing=T),c(5:(4+earl+late))]
gaps_row <- nrow(unique(tumo3[(tumo3$meant-tumo3$meann)>0,1:4]))
dev.off()
pdf("FigureS1J_stage-heatmap.pdf",width=7,height=7)
library(pheatmap)
annotation_col = data.frame(Stage =c(rep("NMIBC",earl),rep("MIBC",late)))
rownames(annotation_col)=colnames(matrix)
pheatmap(as.matrix(unique(matrix[,c(12,1:5,7,8,9,10,11,6,13:15,36,16:34,35)])), 
 #        color = colorRampPalette(c("lightgreen","white","red"))(n = 500),
        color = c(
         colorRampPalette(colors=c("navy","lavender"))(n=50),
         colorRampPalette(colors=c("lavender","white","lightpink"))(n=100)[2:99],
          colorRampPalette(colors=c("lightpink","firebrick4"))(n=50)),
         cellwidth = 6,cellheight =6*36/275,
         show_rownames=F,border_color="gray",font_size=15,show_colnames=F,
         cluster_cols=F,scale="row",cluster_rows=F,breaks= unique(c(seq(-2,2, length=200))),
         legend = T,annotation_col=annotation_col,
         annotation_colors=list(Stage = c(NMIBC = "darkgreen", MIBC="orange")),
         gaps_col = earl,gaps_row = gaps_row
)
dev.off()
#### heatmap part 
stage <-read.table(file="./Stage/Stage_filtered.txt",header=F,stringsAsFactors=F)
name <-read.table(file="/pnas/yangyg_group/hanyn/data/reference/Locus_tag.txt",header=T,stringsAsFactors=F)
mer1<-merge(tumo3,name,by.x="sitep.Gene",by.y="gene",all.x=F,all.y=F, suffixes = c(".x",".y"))
mer2<-merge(mer1,stage,by.x="locus_tag",by.y="V1",all.x=F,all.y=F, suffixes = c(".x",".y"))
mer5 <- mer2[c(-2,-4,-9),]

library(pheatmap)
matrix <- as.matrix(mer5[,c(7:(6+earl+late))])
rownames(matrix) <- mer5$locus_tag
gene_sort <- names(sort(apply(matrix[,c((earl+1):(earl+late))],1,function(x){median(x,na.rm=T)}),decreasing=T))
pdf("Figure1C_stage_heatmap5.pdf",width=7,height=5)
annotation_col = data.frame(Stage =c(rep("NMIBC",earl),rep("MIBC",late)))
rownames(annotation_col)=colnames(matrix)
pheatmap(matrix[gene_sort[c(-9,-11,-13:-15)],c(12,1:5,7,8,9,10,11,6,13:15,36,16:34,35)], 
#         color = colorRampPalette(c("darkblue","white","red"))(n = 500),
         color = c(
           colorRampPalette(colors=c("navy","lavender"))(n=50),
           colorRampPalette(colors=c("lavender","white","lightpink"))(n=100)[2:99],
           colorRampPalette(colors=c("lightpink","firebrick4"))(n=50)),
         cellwidth = 9,cellheight =9,
         show_rownames=T,border_color="gray",font_size=15,show_colnames=F,
         cluster_cols=F,scale="row",cluster_rows=F,gaps_col = earl,
         annotation_col=annotation_col,
         annotation_colors=list(Stage = c(NMIBC = "darkgreen", MIBC="orange")),
         legend = T,breaks= unique(c(seq(-2,2, length=200)))
)
dev.off()
#### boxplot part
onco <- mer5   ### ?????? ????????????????????matrix
rownames(onco) <- onco$locus_tag
onco[gene_sort,c("locus_tag","pvalue")]#####   check if match with star in the plot
onco[gene_sort,c("pvalue")]
# locus_tag      pvalue
# AKT2        AKT2 0.009241761
# ZDHHC8    ZDHHC8 0.004112815
# GRB2        GRB2 0.017183517
# HDGF        HDGF 0.024543212
# GIPC1      GIPC1 0.026431146
# SRRT        SRRT 0.030406307
# MMP15      MMP15 0.017248213
# MMP19      MMP19 0.047464417
# RCN1        RCN1 0.008422988
# LRP5        LRP5 0.049557536
# LAMC2      LAMC2 0.015686111


oncoa <- onco[gene_sort[c(-9,-11,-13:-15)],c(1,7:(6+earl+late))]
#### 20190510
# mean(onco[gene_sort[c(-9,-11,-13:-15)],]$meant - onco[gene_sort[c(-9,-11,-13:-15)],]$meann)
# 0.102123
#write.table(onco[gene_sort[c(-9,-11,-13:-15)],],file="Figure1cd_table.txt",sep="\t",quote=FALSE);


####
oncoa$locus_tag <- factor(oncoa$locus_tag,levels=gene_sort[c(-9,-11,-13:-15)])
library(reshape2)
onco1 <- melt(oncoa)
gnum <- 10
onco1$type <- factor(c(rep("NMIBC",earl*gnum),rep("MIBC",late*gnum)),levels=c("NMIBC","MIBC"))
library(ggplot2)
onco1$x <- c(rep(1:gnum-0.2,earl),rep(1:gnum+0.2,late))
max<-tapply(onco1$value,onco1$locus_tag,function(x){max(x,na.rm=T)})
pdf("Figure1D-stage_boxplot4.pdf",width=8,height=5)
ggplot(onco1)+
  geom_boxplot(aes(x=locus_tag,y=value,fill=type),position=position_dodge(width=0.8),outlier.shape = NA)+
  geom_jitter(aes(x=x,y=value),size=0.4,width=0.07)+
  theme_classic()+ylim(c(0,1.05))+
#  scale_fill_manual(values=c("blue","red"))+labs(x="",y="m5C level",fill="")+
  scale_fill_manual(values=c("navy","firebrick4"))+labs(x="",y="m5C level",fill="")+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))+
  annotate("text",x=1:gnum,y=max+0.03,label=c("**","**",rep("*",6),"**",rep("*",4),rep("**",2))[c(-9,-11,-13:-15)])+
  annotate("segment", x = 1:gnum-0.2, xend = 1:gnum+0.2, y = max+0.025, yend = max+0.025)
dev.off()
################################################################################ 5. lymph stage  ########
rm(list = ls())
setwd("/pnas/yangyg_group/hanyn/2-patient-bs/Tissue_all/Diff")
mc   <-read.table(file="2-Diff_m5C_all.txt",header=T,stringsAsFactors=F)
sitep <- mc[(mc$diff)>0,-1]
#### diff #####
tumo <- data.frame(
  sitep$chrom,sitep$start,sitep$end,sitep$strand,    
  sitep$T0531,  sitep$T0555,  sitep$T0565,  sitep$T1037,  sitep$T1067,  sitep$T2226,
  sitep$T2297,  sitep$T2347,  sitep$T2759,  sitep$T2804,  sitep$T3138,  sitep$T3282,
  sitep$T3403,  sitep$T3451,  sitep$T3515,  sitep$T3560,  sitep$T3614,  sitep$T3635,
  sitep$T3660,  sitep$T3740,  sitep$T3880,  sitep$T0453,  sitep$T0466,  sitep$T1664,
  sitep$T2224,  sitep$T2644,  sitep$T3143,  sitep$T3196,  sitep$T0156,  sitep$T0474,
  sitep$T0982,  sitep$T1495,  sitep$T2271,  sitep$T3371,  sitep$T3750,  sitep$T3357,  sitep$Gene)
earl=21
late=15
tumo$meann <- apply(tumo[,5:(4+earl)],1,function(x) {mean(x[!is.na(x)])})
tumo$meant <- apply(tumo[,(5+earl):(4+earl+late)],1,function(x){mean(x[!is.na(x)])})
tumo2 <- tumo[!is.na(tumo$meann) & (tumo$meant-tumo$meann)>=0.05,]
tumo2$pvalue <- as.numeric(apply(tumo2[,5:(4+earl+late)],1,function(x){wilcox.test(x[1:earl],x[(1+earl):(earl+late)],paired=F)$p.value}))
tumo3 <- tumo2[tumo2$pvalue<0.05,]
#write.table(tumo3,"Stage/Lymph_diff_all.txt",quote=F,sep="\t",row.names = F)



cd /pnas/yangyg_group/hanyn/2-patient-bs/Tissue_all/Diff/Stage
cut -f 41 Stage_diff_all.txt|sort -u |wc -l
cut -f 1,2 Stage_diff_all.txt|sort -u |wc -l

cut -f 41 Lymph_diff_all.txt|sort -u |wc -l
cut -f 1,2 Lymph_diff_all.txt|sort -u |wc -l

#### count ####
head(tumo3)
head(tumo3[tumo3$sitep.Gene=="ENSG00000143321",],n=10)# no
length(unique(tumo3$sitep.Gene)) 
length(unique(tumo3$sitep.Gene[tumo3$meant-tumo3$meann>0])) 
length(unique(tumo3$sitep.Gene[tumo3$meant-tumo3$meann<0])) 
nrow(unique(tumo3[(tumo3$meant-tumo3$meann)>0,1:4]))#  
nrow(unique(tumo3[(tumo3$meant-tumo3$meann)<0,1:4]))#  

#### heatmap all ####
matrix <- tumo3[order(tumo3$meant-tumo3$meann,decreasing=T),c(5:(4+earl+late))]
gaps_row <- nrow(unique(tumo3[(tumo3$meant-tumo3$meann)>0,1:4]))
pdf("FigureS1J_lymph-heatmap.pdf",width=7,height=7)
library(pheatmap)
annotation_col = data.frame("LM" =c(rep("NonLM",earl),rep("LM",late)))
rownames(annotation_col)=colnames(matrix)
pheatmap(as.matrix(unique(matrix)[c(1:109,186:189,110:119,190:195,120:129,180:185,130:139,140:149,206:209,150:159,200:205,160:169,196:199,170:179),
                                 c(1,3,5,8,9,10,2,11,20,12,
                             13,6,14:16,4,17,21,19,18,7,
                             24:29,22,31:36,23,30)]), 
         color = colorRampPalette(c("cyan","white","red"))(n = 500),
         cellwidth = 6,cellheight =6*36/213,
         show_rownames=F,border_color="gray",font_size=15,show_colnames=F,
         cluster_cols=F,scale="row",cluster_rows=F,breaks= unique(c(seq(-2,2, length=500))),
         legend = T,annotation_col=annotation_col,
         annotation_colors=list("LM" = c(NonLM = "burlywood1", LM="brown"))[1],
         gaps_col = earl,gaps_row = gaps_row
)
dev.off()
#### heatmap part ####
#mcup  <-read.table(file="2-Diff_m5C_up_rmintron.bed",header=T,stringsAsFactors=F)
#colnames(mcup)<-c("chrom","start","end","strand","site" ,"Gene" ,"diff","pvalue","FDR","pair_num","pair_up","pair_dn","mean")
name  <-read.table(file="/pnas/yangyg_group/hanyn/data/reference/Locus_tag.txt",header=T,stringsAsFactors=F)
tumo3$site <- paste(tumo3$sitep.chrom,tumo3$sitep.start,sep="_")
tumo4  <-merge(tumo3,name,by.x="sitep.Gene",by.y="gene",all.x=F,all.y=F, suffixes = c(".x",".y"))
#tumo41 <-merge(tumo4,mcup,by.x="site",by.y="site",all.x=F,all.y=F, suffixes = c(".x",".y"))

tumo5 <-tumo4[(tumo4$meant-tumo4$meann)>0,]
mer5  <- tumo5[order(tumo5$pvalue)[1:50],c(46,1,42:45,6:41)]

onco <- read.csv(file="uniprot-oncogene.csv",header=T,stringsAsFactors=F)
tumo6 <-merge(mer5,onco,by.x="locus_tag",by.y="Gene.names...primary..",all.x=F,all.y=F, suffixes = c(".x",".y"))

stage <-read.table(file="./Stage/lymph_filtered.txt",header=F,stringsAsFactors=F)
tumo7<-merge(mer5,stage,by.x="locus_tag",by.y="V1",all.x=F,all.y=F, suffixes = c(".x",".y"))

#tumo9<- tumo5[tumo5$locus_tag=="ATP8B1"|tumo5$locus_tag=="NETO2"|tumo5$locus_tag=="PSMB9",c(46,1,42:45,6:41)]
tumo9<- tumo5[tumo5$locus_tag=="UBE2Q2"|tumo5$locus_tag=="NETO2"|tumo5$locus_tag=="PSMB9",c(46,1,42:45,6:41)]


tumo8 <- unique(rbind(tumo6[,1:42],tumo7[,1:42],tumo9[,1:42]))[c(-1,-4),]
matrix <- as.matrix(tumo8[,c(7:(6+earl+late))])
rownames(matrix) <- tumo8$locus_tag
gene_sort <- names(sort(apply(matrix[,c((earl+1):(earl+late))],1,function(x){median(x,na.rm=T)}),decreasing=T))
annotation_col = data.frame("LM" =c(rep("NonLM",earl),rep("LM",late)))
rownames(annotation_col)=colnames(matrix)
library(pheatmap)
pdf("Figure1E_lymph_heatmap5_20181101_1.pdf",width=7,height=5)
pheatmap(matrix[c(9,4,2,7,3,6,5,1,8,10),c(1:2,4:21,3,30,22,32,35,23,24,28,25,26,27,31,33,34,36,29)], 
         color = colorRampPalette(c("cyan","white","red"))(n = 500),
         cellwidth = 9,cellheight =9,
         show_rownames=T,border_color="gray",font_size=15,show_colnames=F,
         cluster_cols=F,scale="row",cluster_rows=F,gaps_col = earl,
         annotation_col=annotation_col,
         annotation_colors=list("LM" = c(NonLM = "burlywood1", LM="brown"))[1],
         legend = T,breaks= unique(c(seq(-2,2, length=500)))
)
dev.off()
#### boxplot part ####
### ?????? ????????????????????matrix
rownames(tumo8) <- tumo8$locus_tag
oncoa <- tumo8[gene_sort,c(1,7:(6+earl+late))]
tumo8[gene_sort,c(5)]  ### check star 
#### 20190510 
# mean(tumo8[gene_sort,]$meant-tumo8[gene_sort,]$meann)
#write.table(tumo8[gene_sort,],file="Figure1LM.txt",sep="\t",quote=FALSE);
oncoa$locus_tag <- factor(oncoa$locus_tag,levels=gene_sort[1:10])
library(reshape2)
onco1 <- melt(oncoa)
gnum <- 10
onco1$type <- factor(c(rep("NonLM",earl*gnum),rep("LM",late*gnum)),levels=c("NonLM","LM"))
library(ggplot2)
onco1$x <- c(rep(1:gnum-0.2,earl),rep(1:gnum+0.2,late))
max<-tapply(onco1$value,onco1$locus_tag,function(x){max(x,na.rm=T)})
pdf("Figure1F_lymph_boxplot5.pdf",width=6,height=5)
ggplot(onco1)+
  geom_boxplot(aes(x=locus_tag,y=value,fill=type),position=position_dodge(width=0.8),outlier.shape = NA)+
  geom_jitter(aes(x=x,y=value),size=0.4,width=0.07)+
  theme_classic()+ylim(c(0,0.8))+
  scale_fill_manual(values=c("cyan","red"))+labs(x="",y="m5C level",fill="")+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))+
  annotate("text",x=1:gnum,y=max+0.03,label=c("*","**",rep("*",3),"**","*","**","*","**"))+
  annotate("segment", x = 1:gnum-0.2, xend = 1:gnum+0.2, y = max+0.025, yend = max+0.025)
dev.off()
################################ Figure 2 ########
##################################################################################### 1. ovlp gene expr #######
rm(list = ls())
setwd("/pnas/yangyg_group/hanyn/2-patient-bs/Tissue_all/Diff")
mc1  <-read.table(file="2-Diff_m5C.txt",header=T,stringsAsFactors=F)
exp <-read.csv(file="/asnas/yangyg_group/hanyn/Tissue-RNA/RNA-Seq-CP/NGT-20180627/DEG/DESeq2_N29T36/DESeq2/TvsN.csv",header=T,stringsAsFactors=F)
exp1 <- unique(exp[abs(exp$log2FoldChange)>=log2(1.2)&exp$padj<0.05,])
ovlp  <-merge(mc1,exp1,by.x="Gene",by.y="X",all.x=F,all.y=F, suffixes = c(".x",".y")) #7
dim(ovlp)
length(unique(ovlp$site))#
length(unique(ovlp$Gene))#
## point quadrant  ##
ovlp$group <- "NA"
ovlp$group[ovlp$diff>0 & ovlp$log2FoldChange>0] <-"uu"
ovlp$group[ovlp$diff>0 & ovlp$log2FoldChange<0] <-"ud" 
ovlp$group[ovlp$diff<0 & ovlp$log2FoldChange>0] <-"du"
ovlp$group[ovlp$diff<0 & ovlp$log2FoldChange<0] <-"dd" 
length(unique(ovlp$site[ovlp$group=="uu"]))#  
length(unique(ovlp$site[ovlp$group=="ud"]))#  
length(unique(ovlp$site[ovlp$group=="du"]))#  
length(unique(ovlp$site[ovlp$group=="dd"]))#  
length(unique(ovlp$Gene[ovlp$group=="uu"]))#  
length(unique(ovlp$Gene[ovlp$group=="ud"]))#  
length(unique(ovlp$Gene[ovlp$group=="du"]))# 
length(unique(ovlp$Gene[ovlp$group=="dd"]))#
write.table(ovlp,"Figure2B-quantile.txt",quote=F,sep="\t",row.names = F)

pdf("Figure2A-point-quadrant1-0807-hyn.pdf",width=5,height=4)
library(ggplot2)
ggplot()+
  geom_point(data=ovlp,aes(x=diff,y=log2FoldChange,colour=group),size=1)+
  scale_colour_manual(values=c("blue", "green","yellow", "red"))+
   annotate("text", x = 0.3, y = 3.5,  label = "Hyper-up
            600",size=5)+
   annotate("text", x = 0.3, y = -3 ,  label = "Hyper-down
            244",size=5)+
   annotate("text", x = -0.3, y = 3.5, label = "Hypo-up
            111",size=5)+
   annotate("text", x = -0.3, y = -3,  label = "Hypo-down
            91",size=5)+
  # annotate("text", x = 0.3, y = 3,    label = "Gene=475",size=5)+
  # annotate("text", x = 0.3, y = -3.5, label = "Gene=182",size=5)+
  # annotate("text", x = -0.3, y = 3,   label = "Gene=92",size=5)+
  # annotate("text", x = -0.3, y = -3.5,label = "Gene=72",size=5)+
   annotate("segment", x = -0.5, y = 0,xend = 0.5, yend = 0,color="gray90")+
   annotate("segment", x = 0, y = -4,xend = 0, yend = 4,color="gray90")+
  labs(x="m5C level difference (Tumor-Normal)",y="Relatiive gene expression (log2 FoldChange)",color="",size=9)+
  theme_bw()+
  theme(panel.grid.minor.x=element_blank()) +
  theme(panel.grid.major.x=element_blank())+
  theme(panel.grid.minor.y=element_blank()) +
  theme(panel.grid.major.y=element_blank())+
  scale_x_continuous(breaks=c(-0.5,0,0.5),limits=c(-0.5,0.5))+
  scale_y_continuous(breaks=c(-3,0,3),limits=c(-4,4))+
  theme(legend.position="none")
dev.off()
write.table(unique(ovlp$Gene[ovlp$group=="uu"]),"Figure2B-genelist_upup_0807.txt",quote=F,sep="\t",row.names = F)


##################################################################################### 2. function
##################################################################################### 2. ipa function #####
setwd("/pnas/yangyg_group/hanyn/2-patient-bs/Tissue_all/Diff")
mc1  <-read.table(file="2-Diff_m5C.txt",header=T,stringsAsFactors=F)
write.table(unique(mc1[mc1$diff>0,2]),"FigureS1H_genelist_unique_0807.txt",quote=F,sep="\t",row.names = F)
###### function
rm(list = ls())
setwd("G:\\2-patient-bs\\2-Tissue_all\\Diff\\Pathway")
path <-read.table(file="pathway.txt",header=T,sep="\t",stringsAsFactors=F)
pvalue <- path$PValue
names(pvalue)<-path$Path
all_pv<-sort(unique(pvalue));
col_data_pvalue <- NULL;
for(i in 1:length(pvalue)){
  col_data_pvalue<-c(col_data_pvalue,which(all_pv==pvalue[i]));
}
pdf("FigureS1G_0806_hyn.pdf",width=5,height=3)
par(mar=c(4,16,2,1),lwd=3);
barplot(rev(pvalue),horiz=T,las=1,cex.names=1,space=0.5,
        col="red",border="red",xlab="-log10(pvalue)",xaxt="n");
axis(1,at=seq(from=0,to=10,by=1),lwd=1,tcl=-1)
dev.off()
##########

##################################################################################### 3. listA ###############
rm(list = ls())
setwd("/pnas/yangyg_group/hanyn/2-patient-bs/Tissue_all/Diff")
mcup  <-read.table(file="2-Diff_m5C_up_rmintron.bed",header=T,stringsAsFactors=F)
colnames(mcup)<-c("chrom","start","end","strand","site" ,"Gene" ,"diff","pvalue","FDR","pair_num","pair_up","pair_dn","mean")
expu <-read.table(file="/asnas/yangyg_group/hanyn/Tissue-RNA/RNA-Seq-CP/NGT-20180627/DEG/DESeq2_N29T36/DESeq2/TvsN-up-fc1.2-padj05.txt",header=F,sep="\t",stringsAsFactors=F)[,c(1,3,7)]
colnames(expu) <- c("Gene","log2(FC)","FDR")
cor  <-read.table(file="3-cor_positive_737.txt",header=T,stringsAsFactors=F)
head(mcup)
head(expu)
head(cor)
mcupexpu <-merge(mcup,expu,by.x="Gene",by.y="Gene",all.x=F,all.y=F, suffixes = c(".x",".y"))
list <-merge(mcupexpu,cor,by.x="site",by.y="site",all.x=F,all.y=F, suffixes = c(".x",".y"))
head(list)
list1 <- list[list$Gene.x==list$Gene.y,]
length(unique(mcupexpu$Gene))# 
length(unique(mcupexpu$site))#  
length(unique(list1$Gene.x))#  
length(unique(list1$site))#  
name <-read.table(file="/pnas/yangyg_group/hanyn/data/reference/Locus_tag.txt",header=T,stringsAsFactors=F)
mer1 <-merge(list1,name,by.x="Gene.x",by.y="gene",all.x=F,all.y=F, suffixes = c(".x",".y"))
mer2 <- mer1[mer1$cor.s>=0.2,]
length(unique(mer2$locus_tag))#  
write.table(unique(mer1[,c(1:2,ncol(mer1))]),file="4-list-uup_124g.txt",quote=FALSE,sep="\t",row.names = F)

##################################################################################### 4. cor & box ##########
rm(list = ls())
setwd("/pnas/yangyg_group/hanyn/2-patient-bs/Tissue_all/Diff")
mc   <-read.table(file="2-Diff_m5C_all.txt",header=T,stringsAsFactors=F)
mc_pair   <-read.table(file="1-paired_detail_all.txt",header=T,stringsAsFactors=F)
rpkm <-read.table(file="/asnas/yangyg_group/hanyn/Tissue-RNA/RNA-Seq-CP/NGT-20180627/DEG/DESeq2_N29T36/2-N29T36-rpkm.txt",header=T,stringsAsFactors=F)
expu <-read.table(file="/asnas/yangyg_group/hanyn/Tissue-RNA/RNA-Seq-CP/NGT-20180627/DEG/DESeq2_N29T36/DESeq2/TvsN-up-fc1.5-padj01.txt",header=F,stringsAsFactors=F)[,c(1,3,6)]
colnames(expu) <- c("Gene","log2FC","Pvalue")
cor  <-read.table(file="3-cor_positive_737.txt",header=T,stringsAsFactors=F)
name <-read.table(file="/pnas/yangyg_group/hanyn/data/reference/Locus_tag.txt",header=T,stringsAsFactors=F)

sam <- 65
nor <- 29
tum <- 36
mcrp <-merge(mc[(mc$diff)>0,c((6+sam),1,(10+sam),6:(5+sam))],rpkm,by.x="Gene",by.y="GeneID",all.x=F,all.y=F, suffixes = c(".x",".y"))

sample <- c("N0156","N0453","N0466",
            "N0474","N0531","N0555","N0565","N0897","N1037","N1063",
            "N1067","N1495","N1664","N1691","N2224","N2297","N3111",
            "N3138","N3740","N3880","P2347","P2688","P2759","P2804",
            "P3143","P3451","P3467",        "P3522","P3635","T0156",
            "T0453","T0466","T0474","T0531","T0555","T0565","T0982",
            "T1037","T1067","T1495","T1664","T2224","T2226","T2271",
            "T2297","T2347","T2644","T2759","T2804","T3138","T3143",
            "T3196","T3282","T3357","T3371","T3403","T3451","T3515",
            "T3560","T3614","T3635","T3660","T3740","T3750","T3880")
############################################   cor ##########
mycor <- function(myname,x1,x2,y1,y2){
    mygene <- name$gene[name$locus_tag==myname]
    mysite <- cor[!is.na(cor$Gene)&cor$Gene==mygene,"site"]
    mycor  <- paste("R=",format(cor$cor.s[!is.na(cor$Gene)&cor$Gene==mygene],digit=2),sep="")
    mycorp <- paste("p=",format(cor$cor.p.s[!is.na(cor$Gene)&cor$Gene==mygene],digit=3),sep="")
    info <- data.frame(
      m5c= as.numeric(mcrp[mcrp$Gene==mygene & mcrp$site==mysite,4:(3+sam)]),
      rpk= as.numeric(mcrp[mcrp$Gene==mygene & mcrp$site==mysite,(4+sam):(3+sam*2)])
    )
    library(ggplot2)
    ggplot(info,aes(x=m5c,y=log2(rpk)))+
      geom_point()+
      annotate("text",x=x2-0.2,y=y1+.5,label=mycor,vjust=0,hjust=0,size=4)+
      annotate("text",x=x2-0.2,y=y1,label=mycorp,vjust=0,hjust=0,size=4)+
      stat_smooth(method='lm',col="blue",size=1)+
      labs(x="m5C level",y="Expression level (Log2 RPKM)",title=myname)+
      theme_bw()+
      scale_x_continuous(expand = c((x2-x1)/40,(x2-x1)/40),limits=c(x1,x2),breaks=c(seq(0,1,by=0.1)))+
      scale_y_continuous(expand = c((y2-y1)/60,(y2-y1)/60),limits=c(y1,y2),breaks=c(-2:8))+
      theme(panel.grid.major.y=element_blank()) +
      theme(panel.grid.major.x=element_blank())+
      theme(panel.grid.minor.y=element_blank())+
      theme(panel.grid.minor.x=element_blank())
} 

p1 <- mycor("HDGF",0,0.6,4,8)
p3 <- mycor("SMG7",0,0.6,2,5)
p4 <- mycor("MLLT11",0,0.9,-1,7)
p5 <- mycor("PSMD3",-0.02,0.33,3,7)
p6 <- mycor("RAB11FIP4",0,0.62,-1,5)
### F11R matched two ensembl genes
mygene <- "ENSG00000158769" #
myname <- "F11R"
x1 <- 0
x2 <- 0.8
y1 <- 2
y2 <- 7
mysite <- cor[!is.na(cor$Gene)&cor$Gene==mygene,"site"]
mycor  <- paste("R=",format(cor$cor.s[!is.na(cor$Gene)&cor$Gene==mygene],digit=2),sep="")
mycorp <- paste("p=",format(cor$cor.p.s[!is.na(cor$Gene)&cor$Gene==mygene],digit=3),sep="")
info <- data.frame(
  m5c= as.numeric(mcrp[mcrp$Gene==mygene & mcrp$site==mysite,4:(3+sam)]),
  rpk= as.numeric(mcrp[mcrp$Gene==mygene & mcrp$site==mysite,(4+sam):(3+sam*2)])
)
library(ggplot2)
p2<-ggplot(info,aes(x=m5c,y=log2(rpk)))+
  geom_point()+
  annotate("text",x=x2-0.3,y=y1+.5,label=mycor,vjust=0,hjust=0,size=4)+
  annotate("text",x=x2-0.3,y=y1,label=mycorp,vjust=0,hjust=0,size=4)+
  stat_smooth(method='lm',col="blue",size=1)+
  labs(x="m5C level",y="Expression level (Log2 RPKM)",title=myname)+
  theme_bw()+
  scale_x_continuous(expand = c((x2-x1)/40,(x2-x1)/40),limits=c(x1,x2))+
  scale_y_continuous(expand = c((y2-y1)/60,(y2-y1)/60),limits=c(y1,y2))+
  theme(panel.grid.major.y=element_blank()) +
  theme(panel.grid.major.x=element_blank())+
  theme(panel.grid.minor.y=element_blank())+
  theme(panel.grid.minor.x=element_blank())
pdf("Figure2C_cor5_axis.pdf",width=10,height=6)
library("grid")
pushViewport(viewport(layout=grid.layout(2,3)))
vplayout<-function(x,y) viewport(layout.pos.row=x,layout.pos.col=y)
print(p1,vp=vplayout(1,1))
print(p2,vp=vplayout(1,2))
print(p3,vp=vplayout(1,3))
print(p4,vp=vplayout(2,1))
print(p5,vp=vplayout(2,2))
print(p6,vp=vplayout(2,3))
dev.off()

##################################################   cor remove ourliers##########
mycor <- function(myname,x1,x2,y1,y2){
  mygene <- name$gene[name$locus_tag==myname]
  mysite <- cor[!is.na(cor$Gene)&cor$Gene==mygene,"site"]
  mycor  <- paste("R=",format(cor$cor.s[!is.na(cor$Gene)&cor$Gene==mygene],digit=2),sep="")
  mycorp <- paste("p=",format(cor$cor.p.s[!is.na(cor$Gene)&cor$Gene==mygene],digit=3),sep="")
  info <- data.frame(
    m5c= as.numeric(mcrp[mcrp$Gene==mygene & mcrp$site==mysite,4:(3+sam)]),
    rpk= as.numeric(mcrp[mcrp$Gene==mygene & mcrp$site==mysite,(4+sam):(3+sam*2)])
  )
  library(ggplot2)
  x1 <- 
  
  ggplot(info,aes(x=m5c,y=log2(rpk)))+
    geom_point()+
    annotate("text",x=x2-0.2,y=y1+.5,label=mycor,vjust=0,hjust=0,size=4)+
    annotate("text",x=x2-0.2,y=y1,label=mycorp,vjust=0,hjust=0,size=4)+
    stat_smooth(method='lm',col="blue",size=1)+
    labs(x="m5C level",y="Expression level (Log2 RPKM)",title=myname)+
    theme_bw()+
    scale_x_continuous(expand = c((x2-x1)/40,(x2-x1)/40),limits=c(x1,x2),breaks=c(seq(0,1,by=0.2)))+
    scale_y_continuous(expand = c((y2-y1)/60,(y2-y1)/60),limits=c(y1,y2),breaks=c(-2:8))
} 

p1 <- mycor("HDGF",0,0.6,4,8)
p3 <- mycor("SMG7",0,0.6,2,5)
p4 <- mycor("MLLT11",0,0.9,-1,7)
p5 <- mycor("PSMD3",-0.02,0.23,3,7)
p6 <- mycor("RAB11FIP4",0,0.62,-1,5)
### F11R matched two ensembl genes
mygene <- "ENSG00000158769" #
myname <- "F11R"
x1 <- 0
x2 <- 0.8
y1 <- 2
y2 <- 7
mysite <- cor[!is.na(cor$Gene)&cor$Gene==mygene,"site"]
mycor  <- paste("R=",format(cor$cor.s[!is.na(cor$Gene)&cor$Gene==mygene],digit=2),sep="")
mycorp <- paste("p=",format(cor$cor.p.s[!is.na(cor$Gene)&cor$Gene==mygene],digit=3),sep="")
info <- data.frame(
  m5c= as.numeric(mcrp[mcrp$Gene==mygene & mcrp$site==mysite,4:(3+sam)]),
  rpk= as.numeric(mcrp[mcrp$Gene==mygene & mcrp$site==mysite,(4+sam):(3+sam*2)])
)
library(ggplot2)
p2<-ggplot(info,aes(x=m5c,y=log2(rpk)))+
  geom_point()+
  annotate("text",x=x2-0.3,y=y1+.5,label=mycor,vjust=0,hjust=0,size=4)+
  annotate("text",x=x2-0.3,y=y1,label=mycorp,vjust=0,hjust=0,size=4)+
  stat_smooth(method='lm',col="blue",size=1)+
  labs(x="m5C level",y="Expression level (Log2 RPKM)",title=myname)+
  theme_bw()+
  scale_x_continuous(expand = c((x2-x1)/40,(x2-x1)/40),limits=c(x1,x2))+
  scale_y_continuous(expand = c((y2-y1)/60,(y2-y1)/60),limits=c(y1,y2))
pdf("Figure2C_cor1.pdf",width=10,height=6)
library("grid")
pushViewport(viewport(layout=grid.layout(2,3)))
vplayout<-function(x,y) viewport(layout.pos.row=x,layout.pos.col=y)
print(p1,vp=vplayout(1,1))
print(p2,vp=vplayout(1,2))
print(p3,vp=vplayout(1,3))
print(p4,vp=vplayout(2,1))
print(p5,vp=vplayout(2,2))
print(p6,vp=vplayout(2,3))
dev.off()

######
############################################   m5c only paired samples ####
mym5c <- function(myname,py,y2){
    mygene <- name$gene[name$locus_tag==myname]
    mysite <- cor[!is.na(cor$Gene)&cor$Gene==mygene,"site"]
    mym5cp <- paste("p=",format(mc_pair$pair_p[mc_pair$Gene==mygene&mc_pair$site==mysite],digit=3),sep="")
    myup <- paste("Up:",format(mc$pair_up[mc$Gene==mygene&mc$site==mysite],digit=3),sep="")
    mydw <- paste("Down:",format(mc$pair_dn[mc$Gene==mygene&mc$site==mysite],digit=3),sep="")
    m5c <- data.frame(m5c=as.numeric(mc_pair[mc_pair$Gene==mygene&mc_pair$site==mysite,-1:-5]),sample=c(rep("Normal",22),rep("Tumor",22)))
    ggplot()+
      geom_boxplot(data=m5c,aes(x=sample,y=m5c,fill=sample),outlier.shape = NA)+
      geom_line(data=m5c,aes(x=sample,y=m5c,group=c(1:22,1:22)),linetype="dashed",size=0.6)+
      theme_bw()+theme(legend.position="none")+
      annotate("text",x=0.5,y=py, label=mym5cp,vjust=0,hjust=0,size=4)+
      annotate("text",x=0.5,y=py-0.05, label=myup,vjust=0,hjust=0,size=4)+
      annotate("text",x=0.5,y=py-0.1, label=mydw,vjust=0,hjust=0,size=4)+
      scale_fill_manual(values=c("aquamarine","salmon"))+labs(title=myname,x="",y="m5C level",fill="")+
      scale_y_continuous(expand = c(y2/30,y2/30),limits=c(0,y2))+
      theme(panel.grid.major.y=element_blank()) +
      theme(panel.grid.major.x=element_blank())+
      theme(panel.grid.minor.y=element_blank())+
      theme(panel.grid.minor.x=element_blank())
    
}
library(ggplot2)
p1 <- mym5c("HDGF",0.55,0.6)
p3 <- mym5c("SMG7",0.45,0.5)
p4 <- mym5c("MLLT11",0.55,0.6)
p5 <- mym5c("PSMD3",0.25,0.3)
p6 <- mym5c("RAB11FIP4",0.5,0.6)
##
### F11R matched two ensembl genes
mygene <- "ENSG00000158769" #F11R
myname <- "F11R"
py <- 0.75
y2 <- 0.8
mysite <- cor[!is.na(cor$Gene)&cor$Gene==mygene,"site"]
mym5cp <- paste("p=",format(mc_pair$pair_p[mc_pair$Gene==mygene&mc_pair$site==mysite],digit=3),sep="")
myup <- paste("Up:",format(mc$pair_up[mc$Gene==mygene&mc$site==mysite],digit=3),sep="")
mydw <- paste("Down:",format(mc$pair_dn[mc$Gene==mygene&mc$site==mysite],digit=3),sep="")
m5c <- data.frame(m5c=as.numeric(mc_pair[mc_pair$Gene==mygene&mc_pair$site==mysite,-1:-5]),sample=c(rep("Normal",22),rep("Tumor",22)))
p2 <- ggplot()+
  geom_boxplot(data=m5c,aes(x=sample,y=m5c,fill=sample),outlier.shape = NA)+
  geom_line(data=m5c,aes(x=sample,y=m5c,group=c(1:22,1:22)),linetype="dashed",size=0.6)+
  theme_bw()+theme(legend.position="none")+
  annotate("text",x=0.5,y=py, label=mym5cp,vjust=0,hjust=0,size=4)+
  annotate("text",x=0.5,y=py-0.05, label=myup,vjust=0,hjust=0,size=4)+
  annotate("text",x=0.5,y=py-0.1, label=mydw,vjust=0,hjust=0,size=4)+
  scale_fill_manual(values=c("aquamarine","salmon"))+labs(title=myname,x="",y="m5C level",fill="")+
  scale_y_continuous(expand = c(y2/30,y2/30),limits=c(0,y2))+
  theme(panel.grid.major.y=element_blank()) +
  theme(panel.grid.major.x=element_blank())+
  theme(panel.grid.minor.y=element_blank())+
  theme(panel.grid.minor.x=element_blank())


pdf("Figure2D_m5c4.pdf",width=10,height=6)
library("grid")
pushViewport(viewport(layout=grid.layout(2,3)))
vplayout<-function(x,y) viewport(layout.pos.row=x,layout.pos.col=y)
print(p1,vp=vplayout(1,1))
print(p2,vp=vplayout(1,2))
print(p3,vp=vplayout(1,3))
print(p4,vp=vplayout(2,1))
print(p5,vp=vplayout(2,2))
print(p6,vp=vplayout(2,3))
dev.off()

############################################   EXP ####
rpkmp <- rpkm[,c("GeneID",
                 "S0156N","N0453N","N0466N","N0474N","N0531N","N0555N","N0565N","N1037N","N1067N","N1495N",
                 "N1664N","N2224N","N2297N","P2347P","P2759P","P2804P","N3138N","P3143P","P3451P","P3635P",
                 "N3740N","N3880N",
                 "S0156T","S0453T","S0466T","T0474T","T0531T","T0555T","T0565T","T1037T","T1067T","T1495T",
                 "T1664T","T2224T","T2297T","T2347T","T2759T","T2804T","T3138T","T3143T","T3451T","T3635T",
                 "T3740T","T3880T")]
rpkmp$pair_up <- apply(rpkmp[,2:45],1,function(x){sum(x[23:44]-x[1:22]>0,na.rm=T)})
rpkmp$pair_dn <- apply(rpkmp[,2:45],1,function(x){sum(x[23:44]-x[1:22]<0,na.rm=T)})
rpkmp$pair_p  <- apply(rpkmp[,2:45],1,function(x){wilcox.test(x[23:44],x[1:22],paired=T)$p.value})
myexp <- function(myname,ylim1,ylim2){
  mygene <- name$gene[name$locus_tag==myname]
  myexpp <- paste("p=",format(rpkmp$pair_p[rpkmp$Gene==mygene],digit=3),sep="")
  myup <- paste("Up:",rpkmp$pair_up[rpkmp$GeneID==mygene],sep="")
  mydw <- paste("Down:",rpkmp$pair_dn[rpkmp$GeneID==mygene],sep="")
  myrpkm <- data.frame(rpkm=as.numeric(rpkm[rpkm$GeneID==mygene,-1]),sample=c(rep("Normal",nor),rep("Tumor",tum)))
  rownames(myrpkm) <- sample
  myrpkmp <- myrpkm[c("N0156","N0453","N0466","N0474","N0531","N0555","N0565","N1037","N1067","N1495","N1664",
                  "N2224","N2297","P2347","P2759","P2804","N3138","P3143","P3451","P3635","N3740","N3880",
                  "T0156","T0453","T0466","T0474","T0531","T0555","T0565","T1037","T1067","T1495","T1664",
                  "T2224","T2297","T2347","T2759","T2804","T3138","T3143","T3451","T3635","T3740","T3880"),]
  library(ggplot2)
  ggplot()+
    geom_boxplot(data=myrpkm,aes(x=sample,y=log2(rpkm),fill=sample),outlier.shape = NA)+
    geom_line(data=myrpkmp,aes(x=sample,y=log2(rpkm),group=c(1:22,1:22)),linetype="dashed",size=0.6)+
    theme_bw()+theme(legend.position="none")+
    annotate("text",x=0.5,y=ylim2*0.95, label=myexpp,vjust=0,hjust=0,size=4)+
    annotate("text",x=0.5,y=ylim2*0.9, label=myup,vjust=0,hjust=0,size=4)+
    annotate("text",x=0.5,y=ylim2*0.85, label=mydw,vjust=0,hjust=0,size=4)+
    scale_fill_manual(values=c("cyan","brown1"))+labs(title=myname,x="",y="Expression level (Log2 RPKM)",fill="")+
    scale_y_continuous(expand = c((ylim2-ylim1)/50,(ylim2-ylim1)/50),limits=c(ylim1,ylim2))+
    theme(panel.grid.major.y=element_blank()) +
    theme(panel.grid.major.x=element_blank())+
    theme(panel.grid.minor.y=element_blank())+
    theme(panel.grid.minor.x=element_blank())
  
}
library(ggplot2)
p1 <- myexp("HDGF",5,8)
p3 <- myexp("SMG7",2,5)
p4 <- myexp("MLLT11",-1,5)
p5 <- myexp("PSMD3",3,7)
p6 <- myexp("RAB11FIP4",-1,5)
### F11R matched two ensembl genes
mygene <- "ENSG00000158769" #F11R
myname <- "F11R"
ylim1 <- 2
ylim2 <- 8
myexpp <- paste("p=",format(rpkmp$pair_p[rpkmp$Gene==mygene],digit=3),sep="")
myup <- paste("Up:",rpkmp$pair_up[rpkmp$GeneID==mygene],sep="")
mydw <- paste("Down:",rpkmp$pair_dn[rpkmp$GeneID==mygene],sep="")
myrpkm <- data.frame(rpkm=as.numeric(rpkm[rpkm$GeneID==mygene,-1]),sample=c(rep("Normal",nor),rep("Tumor",tum)))
rownames(myrpkm) <- sample
myrpkmp <- myrpkm[c("N0156","N0453","N0466","N0474","N0531","N0555","N0565","N1037","N1067","N1495","N1664",
                    "N2224","N2297","P2347","P2759","P2804","N3138","P3143","P3451","P3635","N3740","N3880",
                    "T0156","T0453","T0466","T0474","T0531","T0555","T0565","T1037","T1067","T1495","T1664",
                    "T2224","T2297","T2347","T2759","T2804","T3138","T3143","T3451","T3635","T3740","T3880"),]
library(ggplot2)
p2 <- ggplot()+
  geom_boxplot(data=myrpkm,aes(x=sample,y=log2(rpkm),fill=sample),outlier.shape = NA)+
  geom_line(data=myrpkmp,aes(x=sample,y=log2(rpkm),group=c(1:22,1:22)),linetype="dashed",size=0.6)+
  theme_bw()+theme(legend.position="none")+
  annotate("text",x=0.5,y=ylim2*0.95, label=myexpp,vjust=0,hjust=0,size=4)+
  annotate("text",x=0.5,y=ylim2*0.9, label=myup,vjust=0,hjust=0,size=4)+
  annotate("text",x=0.5,y=ylim2*0.85, label=mydw,vjust=0,hjust=0,size=4)+
  scale_fill_manual(values=c("cyan","brown1"))+labs(title=myname,x="",y="Expression level (Log2 RPKM)",fill="")+
  scale_y_continuous(expand = c((ylim2-ylim1)/50,(ylim2-ylim1)/50),limits=c(ylim1,ylim2))+
  theme(panel.grid.major.y=element_blank()) +
  theme(panel.grid.major.x=element_blank())+
  theme(panel.grid.minor.y=element_blank())+
  theme(panel.grid.minor.x=element_blank())


pdf("Figure2E_exp.pdf",width=10,height=6)
library("grid")
pushViewport(viewport(layout=grid.layout(2,3)))
vplayout<-function(x,y) viewport(layout.pos.row=x,layout.pos.col=y)
print(p1,vp=vplayout(1,1))
print(p2,vp=vplayout(1,2))
print(p3,vp=vplayout(1,3))
print(p4,vp=vplayout(2,1))
print(p5,vp=vplayout(2,2))
print(p6,vp=vplayout(2,3))
dev.off()



########  file 5th ##########
rm(list = ls())
setwd("/pnas/yangyg_group/hanyn/2-patient-bs/Tissue_all/Diff")
mc   <-read.table(file="2-Diff_m5C_all.txt",header=T,stringsAsFactors=F)
name <-read.table(file="/pnas/yangyg_group/hanyn/data/reference/Locus_tag.txt",header=T,stringsAsFactors=F)
mymc <-merge(mc,name,by.x="Gene",by.y="gene",all.x=F,all.y=F, suffixes = c(".x",".y"))
listA <- data.frame(Gene=c("HDGF","F11R","MLLT11","PSMD3","RAB11FIP4"))
data <-merge(mymc,listA,by.x="locus_tag",by.y="Gene",all.x=F,all.y=F, suffixes = c(".x",".y"))
write.table(data[-5,],"Figure2_Gene_m5c.txt",quote=F,sep="\t",row.names = F)

rpkm <-read.table(file="/asnas/yangyg_group/hanyn/Tissue-RNA/RNA-Seq-CP/NGT-20180627/DEG/DESeq2_N29T36/2-N29T36-rpkm.txt",header=T,stringsAsFactors=F)
mymc <-merge(rpkm,name,by.x="GeneID",by.y="gene",all.x=F,all.y=F, suffixes = c(".x",".y"))
data <-merge(mymc,listA,by.x="locus_tag",by.y="Gene",all.x=F,all.y=F, suffixes = c(".x",".y"))
write.table(data[-2,],"Figure2_Gene_RPKM.txt",quote=F,sep="\t",row.names = F)

#################


############################ 20181030 ##################
# GO of hypo methylated genes
setwd("/Users/hanyanan/1-Project/BladderCancer/2-patient-bs/2-Tissue_all/Diff/")
kegg <-read.table(file="2-Diff_m5C_down_genelist_GO.txt",header=T,sep="\t",stringsAsFactors=F)
go1 <- kegg[kegg$FDR<=0.05,]
pvalue <- go1$PValue
pvalue <- (-log10(as.numeric(pvalue)))
names(pvalue)<-substr(go1$Term,12,100);
all_pv<-sort(unique(pvalue));
col_data_pvalue <- NULL;
for(i in 1:length(pvalue)){
  col_data_pvalue<-c(col_data_pvalue,which(all_pv==pvalue[i]));
}
pdf("2-hypom5c-gene-kegg.pdf",width=9,height=3)
par(mar=c(4,32,2,2),lwd=3);
barplot(rev(pvalue),horiz=T,las=1,cex.names=1,space=0.5,col="blue",border="blue",xlab="-log10(pvalue)",xaxt="n",main=go1$V3[1]);
axis(1,at=seq(from=0,to=20,by=1),lwd=2,tcl=-0.3)
dev.off()

(pvalue=-log10(c(1.60E-09,5.11E-08,3.21E-04,0.004598368,0.01250396)))
names(pvalue)=c("Ribosome","Translation","Oxidative phosphorylation","Endocytosis","Ubiquitin mediated proteolysis")
pdf("2-hypom5c-gene-kegg_2019??????1.pdf",width=7,height=3)
par(mar=c(4,20,2,2),lwd=3);
barplot(rev(pvalue),horiz=T,las=1,cex.names=1,space=0.5,col="blue",border="blue",xlab="-log10(pvalue)",xaxt="n",main=go1$V3[1]);
axis(1,at=seq(from=0,to=10,by=2),lwd=2,tcl=-0.3)
dev.off()

########################################################    20190510    Figure1 BDF overall change   ########################################################    
# tumor normal
setwd("/pnas/yangyg_group/hanyn/2-patient-bs/Tissue_all/Diff/pathway")
mc   <-read.table(file="../2-Diff_m5C_all.txt",header=T,stringsAsFactors=F)
nrow(unique(mc[,c(1:70,74)]))
head(unique(mc[,c(1:70,74)])$diff)
mean(unique(mc[,c(1:70,74)])$diff) #0.05279981
abs=abs(unique(mc[,c(1:70,74)])$diff)
mean(abs) # 0.09596022

mygene   <-read.table(file="figure1Btable.txt",header=T,stringsAsFactors=F)
data1 = data.frame(type=c(rep("Rep Genes",length(mygene$diff)),rep("Control",length(abs))),
                   Diff=c(mygene$diff,abs)
                    )


library(ggplot2)
pdf("rm-diff1.pdf",width=4,height=3)
ggplot(data1,aes(x=type,y=Diff,fill=type))+
  geom_boxplot()+geom_jitter()+theme_bw()
dev.off()
pdf("rm-diff1a.pdf",width=4,height=3)
ggplot(data1,aes(x=type,y=Diff,fill=type))+
  geom_boxplot()+theme_bw()
dev.off()

#### ##### #### #### ####    early late
rm(list = ls())
setwd("/pnas/yangyg_group/hanyn/2-patient-bs/Tissue_all/Diff")
mc   <-read.table(file="2-Diff_m5C_all.txt",header=T,stringsAsFactors=F)
sitep <- mc[,-1]
tumo <- data.frame(
  sitep$chrom,sitep$start,sitep$end,sitep$strand,    
  sitep$T0466,  sitep$T0555,   sitep$T2226,  sitep$T2347,  sitep$T2804,  sitep$T3196,
  sitep$T3282,  sitep$T3403,  sitep$T3451,  sitep$T3560,  sitep$T3614,  sitep$T3660,
  sitep$T0531,  sitep$T0565,  sitep$T2224,  sitep$T2297,  sitep$T3515,  sitep$T3740,
  sitep$T3880,  sitep$T0453,  sitep$T0474,  sitep$T2271,  sitep$T2644,  sitep$T3138,
  sitep$T3357,  sitep$T3371,  sitep$T3635,  sitep$T0156,  sitep$T0982,  sitep$T1037,
  sitep$T1067,  sitep$T1495,  sitep$T1664,  sitep$T2759,  sitep$T3143,  sitep$T3750,  sitep$Gene)
earl=12
late=24
tumo$meann <- apply(tumo[,5:(4+earl)],1,function(x) {mean(x[!is.na(x)])})
tumo$meant <- apply(tumo[,(5+earl):(4+earl+late)],1,function(x){mean(x[!is.na(x)])})
a = unique(tumo[,c(1:40,42:43)])
mean(abs(a$meant-a$meann),na.rm=T) # 0.04801195


mygene   <-read.table(file="Figure1cd_table.txt",header=T,stringsAsFactors=F)
data1 = data.frame(type=c(rep("Rep Genes",length(mygene$meant-mygene$meann)),rep("Control",length(a$meant))),
                   Diff=c(abs(mygene$meant-mygene$meann),abs(a$meant-a$meann))
)

library(ggplot2)
pdf("rm-diff2.pdf",width=4,height=3)
ggplot(data1,aes(x=type,y=Diff,fill=type))+
  geom_boxplot()+theme_bw()
dev.off()

# tumo2 <- tumo[!is.na(tumo$meann) & abs(tumo$meant-tumo$meann)>=0.05,]## only up
# tumo2$pvalue <- as.numeric(apply(tumo2[,5:(4+earl+late)],1,function(x){wilcox.test(x[1:earl],x[(1+earl):(earl+late)],paired=F)$p.value}))
# tumo3 <- tumo2[tumo2$pvalue<0.05,]
# head(tumo3)
# 
# a = unique(tumo3[,c(1:40,42:43)]) # 309
# mean(abs(a$meant-a$meann),na.rm=T) # 0.1222068

#### ##### #### #### ####    lymph
rm(list = ls())
setwd("/pnas/yangyg_group/hanyn/2-patient-bs/Tissue_all/Diff")
mc   <-read.table(file="2-Diff_m5C_all.txt",header=T,stringsAsFactors=F)
sitep <- mc[,-1]

tumo <- data.frame(
  sitep$chrom,sitep$start,sitep$end,sitep$strand,    
  sitep$T0531,  sitep$T0555,  sitep$T0565,  sitep$T1037,  sitep$T1067,  sitep$T2226,
  sitep$T2297,  sitep$T2347,  sitep$T2759,  sitep$T2804,  sitep$T3138,  sitep$T3282,
  sitep$T3403,  sitep$T3451,  sitep$T3515,  sitep$T3560,  sitep$T3614,  sitep$T3635,
  sitep$T3660,  sitep$T3740,  sitep$T3880,  sitep$T0453,  sitep$T0466,  sitep$T1664,
  sitep$T2224,  sitep$T2644,  sitep$T3143,  sitep$T3196,  sitep$T0156,  sitep$T0474,
  sitep$T0982,  sitep$T1495,  sitep$T2271,  sitep$T3371,  sitep$T3750,  sitep$T3357,  sitep$Gene)
earl=21
late=15
tumo$meann <- apply(tumo[,5:(4+earl)],1,function(x) {mean(x[!is.na(x)])})
tumo$meant <- apply(tumo[,(5+earl):(4+earl+late)],1,function(x){mean(x[!is.na(x)])})
a = unique(tumo[,c(1:40,42:43)])
nrow(a)
mean(abs(a$meant-a$meann),na.rm=T) # 0.0425948

mygene   <-read.table(file="Figure1LM.txt",header=T,stringsAsFactors=F)
data1 = data.frame(type=c(rep("Rep Genes",length(mygene$meant-mygene$meann)),rep("Control",length(a$meant))),
                   Diff=c(abs(mygene$meant-mygene$meann),abs(a$meant-a$meann))
)

library(ggplot2)
pdf("rm-diff_LM.pdf",width=4,height=3)
ggplot(data1,aes(x=type,y=Diff,fill=type))+
  geom_boxplot()+theme_bw()
dev.off()



# tumo2 <- tumo[!is.na(tumo$meann) & abs(tumo$meant-tumo$meann)>=0.05,]
# tumo2$pvalue <- as.numeric(apply(tumo2[,5:(4+earl+late)],1,function(x){wilcox.test(x[1:earl],x[(1+earl):(earl+late)],paired=F)$p.value}))
# tumo3 <- tumo2[tumo2$pvalue<0.05,]
# a = unique(tumo3[,c(1:40,42:43)])
# nrow(a)
# mean(abs(a$meant-a$meann),na.rm=T) # 0.1123028

########################################################    20190513    Figure1 BDF overall change   ########################################################    
setwd("/pnas/yangyg_group/hanyn/2-patient-bs/Tissue_all/Diff/pathway")
heat<- read.csv("path_53genes_filtered.csv",stringsAsFactors=F)  
meann = apply(heat[,9:37],1,function(x){mean(x,na.rm=T)})
names(meann)=heat$locus_tag
mygene<-read.table(file="figure1Btable.txt",header=T,stringsAsFactors=F)
mygen = meann[mygene$locus_tag]
other = meann[-which(names(meann) %in% mygene$locus_tag)]

min(mygen)
median(mygen)
mean(mygen)
max(mygen)

min(other)
median(other)
mean(other)
max(other)

#### early late  
rm(list = ls())
setwd("/pnas/yangyg_group/hanyn/2-patient-bs/Tissue_all/Diff/Stage")
stage <- read.table("Stage_diff_all.txt",stringsAsFactors=F,header=T)  

mygene   <-read.table(file="../Figure1cd_table.txt",header=T,stringsAsFactors=F)
mygen = mygene$meann
(other = stage$meann[-which(stage$sitep.Gene %in% mygene$sitep.Gene)])

min(mygen)
median(mygen)
mean(mygen)
max(mygen)

min(other)
median(other)
mean(other)
max(other)

#### lymph
rm(list = ls())
setwd("/pnas/yangyg_group/hanyn/2-patient-bs/Tissue_all/Diff/Stage")
lymph <- read.table("Lymph_diff_all.txt",stringsAsFactors=F,header=T)  
mygene   <-read.table(file="../Figure1LM.txt",header=T,stringsAsFactors=F)
mygen = mygene$meann
other = lymph$meann[-which(lymph$sitep.Gene %in% mygene$sitep.Gene)]

min(mygen)
median(mygen)
mean(mygen)
max(mygen)

min(other)
median(other)
mean(other)
max(other)



########################################################    20190514    Figure1 A all boxplot    ########################################################    
rm(list = ls())
setwd("/pnas/yangyg_group/hanyn/2-patient-bs/Tissue_all/Diff/pathway")
heat<- read.csv("path_53genes_filtered.csv",stringsAsFactors=F)  
mygene <- c("IL6","TYK2",
            "TBX1","PTK2",
            "PRKD1","PRKCZ","CHAD","RPS6KB2","ZDHHC8",
            "DVL2","HDGF","GRB2")
heat1 = heat[-which(heat$locus_tag %in% mygene),]

sam=65
nor=29
tum=36

library(reshape2) 
library(ggplot2)


onco <- heat1[1:12,c(2,9:(8+sam))]
onco$locus_tag=factor(onco$locus_tag, levels=onco$locus_tag)
onco1 <- melt(onco) 
gnum <- 12
onco1$type <- c(rep("Normal",nor*gnum),rep("Tumor",tum*gnum))
onco1$x <- c(rep(1:gnum-0.2,nor),rep(1:gnum+0.2,tum))
p1=ggplot(onco1)+
  geom_boxplot(aes(x=locus_tag,y=value,fill=type),position=position_dodge(width=0.8),outlier.shape = NA)+
  geom_jitter(aes(x=x,y=value),size=0.4,width=0.07)+
  theme_classic()+ylim(c(0,1))+
  scale_fill_manual(values=c("blue","red"))+labs(x="",y="m5C level",fill="")

onco <- heat1[13:24,c(2,9:(8+sam))]
onco$locus_tag=factor(onco$locus_tag, levels=onco$locus_tag)
onco1 <- melt(onco) 
gnum <- 12
onco1$type <- c(rep("Normal",nor*gnum),rep("Tumor",tum*gnum))
onco1$x <- c(rep(1:gnum-0.2,nor),rep(1:gnum+0.2,tum))
p2=ggplot(onco1)+
  geom_boxplot(aes(x=locus_tag,y=value,fill=type),position=position_dodge(width=0.8),outlier.shape = NA)+
  geom_jitter(aes(x=x,y=value),size=0.4,width=0.07)+
  theme_classic()+ylim(c(0,1))+
  scale_fill_manual(values=c("blue","red"))+labs(x="",y="m5C level",fill="")


onco <- heat1[25:36,c(2,9:(8+sam))]
onco$locus_tag=factor(onco$locus_tag, levels=onco$locus_tag)
onco1 <- melt(onco) 
gnum <- 12
onco1$type <- c(rep("Normal",nor*gnum),rep("Tumor",tum*gnum))
onco1$x <- c(rep(1:gnum-0.2,nor),rep(1:gnum+0.2,tum))
p3=ggplot(onco1)+
  geom_boxplot(aes(x=locus_tag,y=value,fill=type),position=position_dodge(width=0.8),outlier.shape = NA)+
  geom_jitter(aes(x=x,y=value),size=0.4,width=0.07)+
  theme_classic()+ylim(c(0,1))+
  scale_fill_manual(values=c("blue","red"))+labs(x="",y="m5C level",fill="")

pdf("rm_boxplot_rebuttal3.pdf",width=10,height=3)
print(p1)
print(p2)
print(p3)
dev.off()

pdf("rm_boxplot_rebuttal4.pdf",width=3,height=3)
onco <- heat1[37:38,c(2,9:(8+sam))]
onco$locus_tag=factor(onco$locus_tag, levels=onco$locus_tag)
onco1 <- melt(onco) 
gnum <- 2
onco1$type <- c(rep("Normal",nor*gnum),rep("Tumor",tum*gnum))
onco1$x <- c(rep(1:gnum-0.2,nor),rep(1:gnum+0.2,tum))
ggplot(onco1)+
  geom_boxplot(aes(x=locus_tag,y=value,fill=type),position=position_dodge(width=0.8),outlier.shape = NA)+
  geom_jitter(aes(x=x,y=value),size=0.4,width=0.07)+
  theme_classic()+ylim(c(0,1))+
  scale_fill_manual(values=c("blue","red"))+labs(x="",y="m5C level",fill="")
dev.off()
########################################################    20190514    Figure1 BDF overall change -2  ########################################################    
##### tumor normal
#
rm(list = ls())
setwd("/pnas/yangyg_group/hanyn/2-patient-bs/Tissue_all/Diff")
mc   <-read.table(file="2-Diff_m5C_all.txt",header=T,stringsAsFactors=F)
#mygene
mygene<-read.table(file="pathway/figure1Btable.txt",header=T,stringsAsFactors=F)
mygene$meann = apply(mygene[,9:37],1,function(x){mean(x,na.rm=T)})
mygene$meant = apply(mygene[,38:73],1,function(x){mean(x,na.rm=T)})
median(mygene$meann)#0.1409934
median(mygene$diff) #0.1337754
median(mygene$meant/mygene$meann,na.rm=T)#2.162284
#others
ipa   <-read.table(file="ipa_gene.txt",header=T,stringsAsFactors=F)
ipa1 = ipa$gene.y[-which(ipa$gene %in% rownames(mygene))]
ipamc <- mc[which(mc$Gene %in% ipa1),]
median(ipamc$meann)#0.05584013
median(abs(ipamc$meant-ipamc$meann))#0.0837713
median(ipamc$meant/ipamc$meann)# 2.561351

a=length(ipamc$diff)##1290
sum(abs(ipamc$diff)<0.2)/a

##### early late  
rm(list = ls())
setwd("/pnas/yangyg_group/hanyn/2-patient-bs/Tissue_all/Diff/Stage")
stage <- read.table("Stage_diff_all.txt",stringsAsFactors=F,header=T)  

mygene   <-read.table(file="../Figure1cd_table.txt",header=T,stringsAsFactors=F)
median(mygene$meann)#0.1389167
median(abs(mygene$meant-mygene$meann))#0.07410688
median(mygene$meant/mygene$meann)#1.620684

other = stage[-which(stage$sitep.Gene %in% mygene$sitep.Gene),]
median(other$meann)# 0.079875
median(abs(other$meant-other$meann))#  0.08856566
median(other$meant/other$meann)#  2.020081

other$diff=other$meant-other$meann
#
sum(other$diff<0.2)/length(other$diff)

##### lymph
rm(list = ls())
setwd("/pnas/yangyg_group/hanyn/2-patient-bs/Tissue_all/Diff/Stage")
lymph <- read.table("Lymph_diff_all.txt",stringsAsFactors=F,header=T)  

mygene   <-read.table(file="../Figure1LM.txt",header=T,stringsAsFactors=F)
median(mygene$meann)#  
median(abs(mygene$meant-mygene$meann))#  
median(mygene$meant/mygene$meann)#

other = lymph[-which(lymph$sitep.Gene %in% mygene$sitep.Gene),]
median(other$meann)#   
median(abs(other$meant-other$meann))#   
median(other$meant/other$meann)#

other$diff=other$meant-other$meann
length(other$diff)#
sum(other$diff<0.2)/length(other$diff)

sum(other$diff<0.1)/length(other$diff)


###############################   20190525      #########
rm(list = ls())
setwd("/pnas/yangyg_group/hanyn/2-patient-bs/Tissue_all/Diff")
mc   <-read.table(file="2-Diff_m5C_all.txt",header=T,stringsAsFactors=F)

length(unique(mc$Gene[mc$diff>0]))
length(unique(mc$Gene[mc$diff<0]))
length(unique(mc$Gene[mc$diff==0]))

###############################   20190527      #########
rm(list = ls())
setwd("/Users/hanyanan/1-Project/BladderCancer/2-patient-bs/2-Tissue_all/Diff")

down <- read.csv(file="down_IPA.csv",header=F,stringsAsFactors=F)
pvalue <- down$V2
names(pvalue)<-down$V1

pdf("down_ipa2.pdf",width=6,height=4)
par(mar=c(5,17,2,2),lwd=3);
barplot(rev(pvalue),horiz=T,las=1,cex.names=1,space=0.5,col="red",border="red",xlab="-log10(pvalue)",xaxt="n",main="Hypomethylated gene pathway",xlim=c(0,8));
axis(1,at=c(0,4,8),lwd=3,tcl=-1)
dev.off()


down <- read.table(file="Stage/ipa_stage_filtered.txt",header=F,sep="\t",stringsAsFactors=F)
pvalue <- down$V2
names(pvalue)<-down$V1
pdf("ipa_stage.pdf",width=6,height=4)
par(mar=c(5,17,2,2),lwd=3);
barplot(rev(pvalue),horiz=T,las=1,cex.names=1,space=0.5,col="red",border="red",xlab="-log10(pvalue)",xaxt="n",main="Stage IPA pathway",xlim=c(0,3));
axis(1,at=c(0,1,2,3),lwd=3,tcl=-1)
dev.off()


###################################### rebuttal upload table
rm(list = ls())
setwd("/Users/hanyanan/1-Project/BladderCancer/2-patient-bs/2-Tissue_all/Diff/")
mc   <-read.table(file="2-Diff_m5C_all.txt",header=T,stringsAsFactors=F)
site=unique(mc[,c(2,4:70,72:75)])
nrow(site[site$diff>0,])
nrow(site)
head(site)
str(site)
site$key=substr(site$chrom,4,8)
site$key[site$key=="X"]=50
site$key[site$key=="MT"]=60
site1=site[order(as.numeric(site$key),as.numeric(site$end)),]
head(site1)
nrow(site1)
site1[1:30,1:3]
ncol(site1)
write.table(site1[,1:72],"Diff_m5C_upload_5280.txt",quote=F,sep="\t",row.names = F)



setwd("/Users/hanyanan/1-Project/BladderCancer/2-patient-bs/2-Tissue_all/Diff/")
mc   <-read.table(file="2-Diff_m5C_all.txt",header=T,stringsAsFactors=F)
load =unique(mc[,c(2,4:75)])


nrow(load)
nrow(mc)
nrow(unique(mc))
length(unique(mc$Gene[mc$diff>0]))
length(unique(mc$Gene[mc$diff<0]))











