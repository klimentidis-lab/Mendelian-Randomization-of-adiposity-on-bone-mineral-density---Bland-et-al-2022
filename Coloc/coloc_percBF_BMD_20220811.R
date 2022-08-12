
## Amit Arora
## May 17, 2022

library(coloc)

### read data
dat<-read.table(file="C xxx /MRBase/Harmonized_Perc_BF_TB_BMD_20220516.txt",T)
dat<-dat[order(dat$CHR,dat$BP),] 

## create list of snps and corresponding genes
SNP<-c("rs13204965","rs3743347")
GENE<-c("RSPO3","IQCH")

colnames(dat)[1]<-"SNP"
tmp<-cbind(SNP,GENE)
tmp2<-merge(tmp,dat[,c(1,30,31)],by="SNP")


### pull all snps 500kn up and downstream of these snps
### do it for 1st row

myGWAS5=dat[(dat$CHR==tmp2[1,3]),]
upstream=myGWAS5[(myGWAS5$BP>=tmp2[1,4]) & (myGWAS5$BP<(tmp2[1,4]+500000)),]
downstream=myGWAS5[(myGWAS5$BP>(tmp2[1,4]-500000)) & (myGWAS5$BP<tmp2[1,4]),]
#middle=myGWAS5[(myGWAS5$BP>=39304985) & (myGWAS5$BP<=39323226),]

dat3=rbind(downstream,upstream)

#N for each henotype
#bf=450557
#bmd=66132
#bmd=66615

# rRun coloc

dat3$varbeta.bf<-((dat3$se.exposure)^2)
dat3$varbeta.bmd<-((dat3$se.outcome)^2)
my.res <- coloc.abf(dataset1=list(beta=dat3$beta.exposure, varbeta=dat3$varbeta.bf, N=450557,type="quant"),
                    dataset2=list(beta=dat3$beta.outcome, varbeta=dat3$varbeta.bmd, N=66132,type="quant"),
                    MAF=dat3$eaf.exposure)

## min and max base pair position for each chunk
min_bp<-min(dat3$BP)
max_bp<-max(dat3$BP)

## combine the results
out<-cbind(as.data.frame(t(my.res$summary)),min_bp,max_bp)

rm(dat3)
rm(my.res)
rm(myGWAS5)
rm(upstream)
rm(downstream)
rm(min_bp)
rm(max_bp)

######## do it for next snp


for (i in 2: nrow(tmp2)){
  
  myGWAS5=dat[(dat$CHR==tmp2[i,3]),]
  upstream=myGWAS5[(myGWAS5$BP>=tmp2[i,4]) & (myGWAS5$BP<(tmp2[i,4]+500000)),]
  downstream=myGWAS5[(myGWAS5$BP>(tmp2[i,4]-500000)) & (myGWAS5$BP<tmp2[i,4]),]
  #middle=myGWAS5[(myGWAS5$BP>=39304985) & (myGWAS5$BP<=39323226),]
  
  dat3=rbind(downstream,upstream)
  
  dat3$varbeta.bf<-((dat3$se.exposure)^2)
  dat3$varbeta.bmd<-((dat3$se.outcome)^2)
  my.res <- coloc.abf(dataset1=list(beta=dat3$beta.exposure, varbeta=dat3$varbeta.bf, N=450557,type="quant"),
                      dataset2=list(beta=dat3$beta.outcome, varbeta=dat3$varbeta.bmd, N=66615,type="quant"),
                      MAF=dat3$eaf.exposure)
  min_bp<-min(dat3$BP)
  max_bp<-max(dat3$BP)
  A<-cbind(as.data.frame(t(my.res$summary)),min_bp,max_bp)
  out<-rbind(out,A) ## add everything to out df
  print(i)
}

out2<-cbind(tmp2,out)

write.table(out2,file="tbl_2Snps_coloc.results_se2_500kb_perbf_TB-BMD_20220517.txt",row.names=F,quote=F)

######################