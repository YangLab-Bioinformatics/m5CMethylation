rm(list = ls())
setwd("/pnas/yangyg_group/hanyn/2-patient-bs/Tissue_all/")
site0 <-read.table(file="summary_m5C_gene.txt",header=F,stringsAsFactors=F)
setwd("/pnas/yangyg_group/hanyn/2-patient-bs/Tissue_all/Diff/")

colnames(site0) <- c("chrom","start","end","strand",
                     "met.N0156","cov.N0156","lev.N0156",
                     "met.N0453","cov.N0453","lev.N0453",
                     "met.N0466","cov.N0466","lev.N0466",
                     "met.N0474","cov.N0474","lev.N0474",
                     "met.N0531","cov.N0531","lev.N0531",
                     "met.N0555","cov.N0555","lev.N0555",
                     "met.N0565","cov.N0565","lev.N0565",
                     "met.N0897","cov.N0897","lev.N0897",
                     "met.N1037","cov.N1037","lev.N1037",
                     "met.N1063","cov.N1063","lev.N1063",
                     "met.N1067","cov.N1067","lev.N1067",
                     "met.N1495","cov.N1495","lev.N1495",
                     "met.N1664","cov.N1664","lev.N1664",
                     "met.N1691","cov.N1691","lev.N1691", 
                     "met.N2224","cov.N2224","lev.N2224",
                     "met.N2297","cov.N2297","lev.N2297",
                     "met.N3111","cov.N3111","lev.N3111",
                     "met.N3138","cov.N3138","lev.N3138",
                     "met.N3740","cov.N3740","lev.N3740",
                     "met.N3880","cov.N3880","lev.N3880",
                     "met.P2347","cov.P2347","lev.P2347",
                     "met.P2688","cov.P2688","lev.P2688",
                     "met.P2759","cov.P2759","lev.P2759",
                     "met.P2804","cov.P2804","lev.P2804",
                     "met.P3143","cov.P3143","lev.P3143",
                     "met.P3451","cov.P3451","lev.P3451",
                     "met.P3467","cov.P3467","lev.P3467",
                   #  "met.P3515","cov.P3515","lev.P3515",
                     "met.P3522","cov.P3522","lev.P3522",
                     "met.P3635","cov.P3635","lev.P3635",
                     
                     "met.T0156","cov.T0156","lev.T0156",
                     "met.T0453","cov.T0453","lev.T0453",
                     "met.T0466","cov.T0466","lev.T0466",
                     "met.T0474","cov.T0474","lev.T0474",
                     "met.T0531","cov.T0531","lev.T0531",
                     "met.T0555","cov.T0555","lev.T0555",
                     "met.T0565","cov.T0565","lev.T0565",
                     "met.T0982","cov.T0982","lev.T0982",
                     "met.T1037","cov.T1037","lev.T1037",
                     "met.T1067","cov.T1067","lev.T1067",
                     "met.T1495","cov.T1495","lev.T1495",
                     "met.T1664","cov.T1664","lev.T1664",
                     "met.T2224","cov.T2224","lev.T2224",
                     "met.T2226","cov.T2226","lev.T2226",
                     "met.T2271","cov.T2271","lev.T2271",
                     "met.T2297","cov.T2297","lev.T2297",
                     "met.T2347","cov.T2347","lev.T2347",
                     "met.T2644","cov.T2644","lev.T2644",
                     "met.T2759","cov.T2759","lev.T2759",
                     "met.T2804","cov.T2804","lev.T2804",
                     "met.T3138","cov.T3138","lev.T3138",
                     "met.T3143","cov.T3143","lev.T3143",
                     "met.T3196","cov.T3196","lev.T3196",
                     "met.T3282","cov.T3282","lev.T3282",
                     "met.T3357","cov.T3357","lev.T3357",
                     "met.T3371","cov.T3371","lev.T3371",
                     "met.T3403","cov.T3403","lev.T3403",
                     "met.T3451","cov.T3451","lev.T3451",
                     "met.T3515","cov.T3515","lev.T3515",
                     "met.T3560","cov.T3560","lev.T3560",
                     "met.T3614","cov.T3614","lev.T3614",
                     "met.T3635","cov.T3635","lev.T3635",
                     "met.T3660","cov.T3660","lev.T3660",
                     "met.T3740","cov.T3740","lev.T3740",
                     "met.T3750","cov.T3750","lev.T3750",
                     "met.T3880","cov.T3880","lev.T3880","chromg","startg","endg","strandg","Trans","Gene")
 
site <- data.frame(
  site0$chrom,site0$start,site0$end,site0$strand,    
  site0$lev.T0156-site0$lev.N0156,
  site0$lev.T0453-site0$lev.N0453,
  site0$lev.T0466-site0$lev.N0466,
  site0$lev.T0474-site0$lev.N0474,
  site0$lev.T0531-site0$lev.N0531,
  site0$lev.T0555-site0$lev.N0555,
  site0$lev.T0565-site0$lev.N0565,
  site0$lev.T1037-site0$lev.N1037,
  site0$lev.T1067-site0$lev.N1067,
  site0$lev.T1495-site0$lev.N1495,
  site0$lev.T1664-site0$lev.N1664,
  site0$lev.T2224-site0$lev.N2224,
  site0$lev.T2297-site0$lev.N2297,
  site0$lev.T2347-site0$lev.P2347,
  site0$lev.T2759-site0$lev.P2759,
  site0$lev.T2804-site0$lev.P2804,
  site0$lev.T3138-site0$lev.N3138,
  site0$lev.T3143-site0$lev.P3143,
  site0$lev.T3451-site0$lev.P3451,
  site0$lev.T3635-site0$lev.P3635,
  site0$lev.T3740-site0$lev.N3740,
  site0$lev.T3880-site0$lev.N3880,
  site0$Gene)

## diff 
pair=22
site$pair_num <- apply(site[,5:(4+pair)],1,function(x){sum(!is.na(x))})
site$pair_up  <- apply(site[,5:(4+pair)],1,function(x){sum(!is.na(x) & x>0)})
site$pair_dn  <- apply(site[,5:(4+pair)],1,function(x){sum(!is.na(x) & x<0)})
site$mean     <- apply(site[,5:(4+pair)],1,function(x){mean(x,na.rm=T)})
head(site)

site1 <- site[site$pair_num>=5,]
site2 <- site1[site1$pair_up>=5 | site1$pair_dn>=5,]
site3 <- site2[abs(site2$mean)>=0.05,]
site3$site <- paste(site3$site0.chrom,site3$site0.start,sep="_")
length(unique(site3$site))
head(unique(site3[,c(10+pair,(6+pair):(9+pair))]))
write.table(unique(site3[,c(10+pair,(6+pair):(9+pair))]),"1-pair_555.txt",quote=F,sep="\t",row.names = F)

### ¿¨ pairµÄpvalue
# site3$pvalue <- as.numeric(apply(site3[,(5+pair):(4+pair*3)],1,function(x){wilcox.test(x[1:pair],x[(1+pair):(pair*2)],paired=F)$p.value}))
# library(qvalue)
# site3$qvalue <- qvalue(site3$pvalue)$qvalues
# site4 <- site3[abs(site3$qvalue)< 0.05,]
# site4$site <- paste(site4$site0.chrom,site4$site0.start,sep="_")
# write.table(unique(site4[,c(12+pair*3,(6+pair*3):(9+pair*3))]),"pair_22_555.txt",quote=F,sep="\t",row.names = F)
# dim(site4)#317409



#################### detailed table 
setwd("/pnas/yangyg_group/hanyn/2-patient-bs/Tissue_all/Diff/")
mc   <-read.table(file="2-Diff_m5C_all.txt",header=T,stringsAsFactors=F)

pair <- mc[,c("Gene","site",
                 "N0156", "N0453", "N0466", "N0474", "N0531", "N0555", "N0565", "N1037", "N1067", "N1495",
                 "N1664", "N2224", "N2297", "P2347", "P2759", "P2804", "N3138", "P3143", "P3451", "P3635",
                 "N3740", "N3880",
                 "T0156", "T0453", "T0466", "T0474", "T0531", "T0555", "T0565", "T1037", "T1067", "T1495",
                 "T1664", "T2224", "T2297", "T2347", "T2759", "T2804", "T3138", "T3143", "T3451", "T3635",
                 "T3740", "T3880")]
sam_p <- 22
pair$up <- apply(pair[,-1:-2],1,function(x){sum((x[(1+sam_p):(sam_p*2)]-x[1:sam_p])>0,na.rm=T)})
pair$dn <- apply(pair[,-1:-2],1,function(x){sum((x[(1+sam_p):(sam_p*2)]-x[1:sam_p])<0,na.rm=T)})
pair$pair_p <- apply(pair[,-1:-2],1,function(x){wilcox.test(x[1:sam_p],x[(sam_p+1):(sam_p*2)],paired=T)$p.value})
pair1 <- pair[,c(1,2,47:49,3:46)]
write.table(pair1,"1-paired_detail_all.txt",quote=F,sep="\t",row.names = F)
# 
# pair1[pair1$Gene=="ENSG00000143321",1:5]
# pair1[pair1$Gene=="ENSG00000116698",1:5]
# pair1[pair1$Gene=="ENSG00000213190",1:5]
# pair1[pair1$Gene=="ENSG00000158769",1:5]
# pair1[pair1$Gene=="ENSG00000108344",1:5]
# pair1[pair1$Gene=="ENSG00000131242",1:5]
# 
# HDGF <- pair1[pair1$Gene=="ENSG00000143321",-1:-5]
# HDGF[,23:44]-HDGF[,1:22]
# 
# x <- as.numeric(HDGF)
# 


              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              