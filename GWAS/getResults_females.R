setwd('/xdisk/yann/mig2020/rsgrps/yann/Victoria/GWAS_BodyFatPer_SexCombo/Females')
dat2 <- read.table('bolt_460K_selfRepWhite.PerBF_Females_Gwas.bgenv3.stats.gz', header=T)

dim(dat2)
#[1] 19971372       14

#looks like there are rare snps for some reason
dat2<-dat2[!(dat2$A1FREQ<0.001),]
dat2<-dat2[!(dat2$A1FREQ>0.999),]

dim(dat2)
#[1] 16831149       14

table(dat2$CHR[which(dat2$P_BOLT_LMM_INF<5e-8)])
#   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
#1993 2198  914  747 1357 2250  472 1008 1378  442 1566 1412  514  314 1672 1620 
#17   18   19   20   21   22   23 
#931 1181  479  961   45  512   50

save(dat2, file='FullGWASresults_PerBF_Females_GWAS.rda')

###############################################################################################
#get top hits #################################################################################
dat3<-dat2[order(dat2$P_BOLT_LMM_INF), ]
#head(psortdata)
#length(psortdata$p[which(psortdata$p<0.00000005)])
datasig<-dat3[which(dat3$P_BOLT_LMM_INF<5e-8),] ## for 5X10-8 we have about 71 SNPs to go through and end up with 508 regions in 'mat'
#datasig<-datasig[1-86,]                     ##### there are 1500 2Mb regions in the genome, perhaps limit it to 1Mb regions (3000 of those)
#head(datasig)
#length(datasig$p)
#summary(datasig$chromosome)
datasig$SNP<-as.character(datasig$SNP)
datasig$CHR<-as.character(datasig$CHR)
datasig$BP<-as.numeric(datasig$BP)
nrow(datasig)
#[1] 24016

###########################################################################################
#GET 'INDEPENDENT' SNPs
SNP<-rep(NA, 4000)
Chrom<-rep(NA, 4000)
kb<-rep(NA, 4000)
p<-rep(NA, 4000)
mat<-data.frame(SNP, Chrom, kb, p)
#head(mat)
mat[1,1]<-datasig$SNP[1]
mat[1,2]<-datasig$CHR[1]
mat[1,3]<-datasig$BP[1]
mat[1,3]<-as.numeric(mat[1,3])
mat[1,4]<-datasig$P_BOLT_LMM_INF[1]

j=2
for(i in 2:length(datasig$P_BOLT_LMM_INF)){
  enter<-c()
  if (!datasig$CHR[i]%in%mat[,2][which(mat[,2]!='NA')]) {
    mat[j,1]<-datasig$SNP[i]
    mat[j,2]<-datasig$CHR[i]
    mat[j,3]<-datasig$BP[i]
    mat[j,4]<-datasig$P_BOLT_LMM_INF[i]
    j<-j+1
  }
  else
    
    for (k in 1:length(which(mat[,2]!='NA'))) {
      if (mat[k,2]==datasig$CHR[i]) {
        if (datasig$BP[i] < (mat[k,3] - 500000) | datasig$BP[i] > (mat[k,3] + 500000)) {
          enter[k]<-1 }} 
      if (mat[k,2]==datasig$CHR[i]){
        if (datasig$BP[i] > (mat[k,3] - 500000) && datasig$BP[i] < (mat[k,3] + 500000)) {
          enter[k]<-0}}
    }
  
  if(!0%in%enter && 1%in%enter) {
    mat[j,1]<-datasig$SNP[i]
    mat[j,2]<-datasig$CHR[i]
    mat[j,3]<-datasig$BP[i]
    mat[j,4]<-datasig$P_BOLT_LMM_INF[i]
    j<-j+1
  }
  
  print(i)
}

mat2<-na.omit(mat)
top<-dat2[which(dat2$SNP%in%mat2$SNP),]

dim(top)
#[1] 310  14

write.table(top, 'TopHits_PerBF_Females_GWAS.txt', quote=F, row.names=F, sep='\t')

###################################################################################################
###################################################################################################


    ##################################################################
    ### MAKE INPUT FILE FOR LDHUB ####################################
    ########################################################################
    dat2$SNP<-as.character(dat2[,1])
    dat2$A1<-as.character(dat2[,5])
    dat2$N<-457019
    dat2$A2<-as.character(dat2[,6])
    #dat2$OR<-exp(dat2$BETA)
    
    dat3<-subset(dat2,select=c("SNP","A1","A2","BETA","N","P_BOLT_LMM_INF"))
    dat3<-na.omit(dat3)
    names(dat3)[6]<-'P-value'
    names(dat3)[1]<-'snpid'

    #snpid   A1      A2      BETA      N       P-value
    
    write.table(dat3,file='LDHUB_PerBF_Females_GWAS_BOLTLMM_500K.txt', row.names=F, quote=F, sep='\t')
    
    #to zip: 
    system("zip LDHUB_PerBF_Females_GWAS_BOLTLMM_500K LDHUB_PerBF_Females_GWAS_BOLTLMM_500K.txt")
    ################################################################################################################################
    
    
    library(qqman)
    
    load('/xdisk/yann/mig2020/rsgrps/yann/Victoria/GWAS_BodyFatPer_SexCombo/Females/FullGWASresults_PerBF_Females_GWAS.rda')
    
    SNP<-as.character(dat2$SNP)
    CHR<-as.numeric(as.character(dat2$CHR))
    BP<-as.numeric(as.character(dat2$BP))
    P<-as.numeric(as.character(dat2$P_BOLT_LMM_INF))
    
    myGWAS<-data.frame(SNP, CHR, BP, P)
    myGWAS2<-na.omit(myGWAS)
    myGWAS3=myGWAS2[which(myGWAS2$P<0.01),]
    rm(myGWAS)
    
    setwd('/xdisk/yann/mig2020/rsgrps/yann/Victoria/GWAS_BodyFatPer_SexCombo/Females')
    
    png(file='manhattan_PerBF_Females_9.png', res = 300, width = 10, height = 6,units = 'in')
    manhattan(myGWAS3, col = c("tomato4", "tomato"), ylim=c(2,9))
    dev.off()
    #####
    ####
    png(file='manhattan_PerBF_Females_20.png', res = 300, width = 10, height = 6,units = 'in')
    manhattan(myGWAS3, col = c("tomato4", "tomato"), ylim=c(5,20))
    dev.off()
    
    png(file='manhattan_PerBF_Females_50.png', res = 300, width = 10, height = 6,units = 'in')
    manhattan(myGWAS3, col = c("tomato4", "tomato"), ylim=c(5,50))
    dev.off()
    
    png(file='QQ_PerBF_Females.png',res = 300, width = 10, height = 6,units = 'in')
    qq(myGWAS2$P, main='Q-Q plot of GWAS')
    dev.off()
    ####
    ####

    
###################################################################################################
###################################################################################################

    library(dplyr)
    
    #Find SNPs significant for favorable adiposity variants
    faSNPs <- c("rs11118306", "rs13389219", "rs2943653", "rs1801282", "rs2276936",
                 "rs40271", "rs632057", "rs998584", "rs972283", "rs2980888",
                 "rs7133378", "rs11045172", "rs7258937", "rs2267373")
    
    #use dplyr filter function to filter out favorable adiposity variants
    dat.fa <- filter(dat2, (SNP %in% faSNPs))

    write.table(dat.fa,file='FAVariants_PerBF_Females_GWAS.txt', row.names=F, quote=F)