####################################################
##     buildMatrix_plotCRISPR - October 2019      ##
## created by Berkley Gryder, gryderart@gmail.com ##
####################################################

### script for CRISPR data generated in house and processed through MAGIK

### 1. Set up folders and samples
  setwd("K:/projects/CRISPR")
  list.seq = "TKO_HGLib" 
  project = "MTC"
  project.folder = paste("projects/",project,"/",sep="") #define output location
  sample.all = list.dirs(path = "DATA/", full.names = F, recursive = F)
  sample.list = sample.all[grep(list.seq,sample.all)]
  sample.list = sample.list[grep(project,sample.list)]
  
### 2. Loop through the sample list and get sample counts
  
  #initiate dataframe before loop
  CRISPR.frame = read.table(paste("DATA/",sample.list[1],"/",sample.list[1],".count_normalized.txt",sep=""), sep="\t", header=T, quote="")
    CRISPR.frame = CRISPR.frame[order(CRISPR.frame$sgRNA),]
    CRISPR.frame = CRISPR.frame[,1:2]
    
  #loop through sample list to make a matrix  
  lapply(sample.list, function(x) {
    #read in txt file
    CRISPR.sample = read.table(paste("DATA/",x,"/",x,".count_normalized.txt",sep=""), sep="\t", header=T, quote="")
    #sort by sgRNA
    CRISPR.sample = CRISPR.sample[order(CRISPR.sample$sgRNA),]
    CRISPR.sample = as.data.frame(CRISPR.sample$read1)
    colnames(CRISPR.sample) = x
    CRISPR.frame <<- cbind(CRISPR.frame,CRISPR.sample)  #combine with others, writing outsite the loop
  })
  
  write.table(CRISPR.frame, file = paste(project.folder,list.seq,".count.norm.matrix.txt",sep=""), col.names = T, row.names = F, quote = F, sep = "\t")

  #subset to remove zeros
  CRISPR.frame.nozeros = subset(CRISPR.frame, !(CRISPR.frame$Sample_TKO_HGLib_PCR_MTC_VAN_IC70_D71_CRISPR_HFKMHBGXB == 0))
  CRISPR.frame.picks = subset(CRISPR.frame.nozeros, CRISPR.frame.nozeros$Gene %in% c("MYC","CTCF","YY1","GAPDH"))
  
### 3. Make rnk list
  CRISPR.frame.nozeros$DMSO_D25_v_D3 = log(x = CRISPR.frame.nozeros$Sample_TKO_HGLib_PCR_MTC_DMSO_D25_CRISPR_HLWGCBGX9+1, base = 2)-log(x = CRISPR.frame.nozeros$Sample_TKO_HGLib_PCR_MTC_DMSO_D3_CRISPR_HKLNMBGX9+1, 2)
  CRISPR.frame.nozeros$DMSO_D57_v_D3 = log(x = CRISPR.frame.nozeros$Sample_TKO_HGLib_PCR_MTC_DMSO_D57_CRISPR_HFKMHBGXB+1, base = 2)-log(x = CRISPR.frame.nozeros$Sample_TKO_HGLib_PCR_MTC_DMSO_D3_CRISPR_HKLNMBGX9+1, 2)
  
  dir.create(file.path(project.folder, "GSEA_ranklist"))
  
  CRISPR.frame.nozeros[,(ncol(CRISPR.frame)+1):ncol(CRISPR.frame.nozeros)]=round(CRISPR.frame.nozeros[,(ncol(CRISPR.frame)+1):ncol(CRISPR.frame.nozeros)], digits = 5)
  
  for (i in (ncol(CRISPR.frame)+1):ncol(CRISPR.frame.nozeros)) {
    Ranklist <- data.frame(CRISPR.frame.nozeros[, 2])
    Ranklist$L2FC_sgRNA <- CRISPR.frame.nozeros[, i]
    Ranklist = Ranklist[rev(order(Ranklist$L2FC_sgRNA)),]
    SampleName = colnames(CRISPR.frame.nozeros)[i]
    mytime <- format(Sys.time(), "%b_%d_%Y")
    myfile <- file.path(project.folder, "GSEA_ranklist", paste0(SampleName,"_",mytime,".rnk"))
    write.table(Ranklist, file = myfile, sep = "\t", row.names = FALSE, col.names = FALSE,
                quote = FALSE, append = FALSE)
  }
  
  
  
### 4. Plot the data
  
  #make long format
  library(tidyr)
  CRISPR.long  <- gather(CRISPR.frame, CRISPR_sample, count_norm, 3:(ncol(CRISPR.frame)))
  require(reshape)
  CRISPR.long = separate(data = CRISPR.long, col = CRISPR_sample, into = c("LIST_SEQ", "Sample"), sep = "_PCR_")
  CRISPR.long = CRISPR.long[!grepl("UISO",CRISPR.long$Sample),]
  
  #load in CR TF gene list
  genelist = read.table(file = paste(project.folder,"CR_TFs_VPMCC.genelist.txt",sep=""), header = F)
  CRISPR.long$CRTFs = CRISPR.long$Gene %in% genelist$V1
 
  library(ggplot2)
  ggplot(CRISPR.long, aes(x=Sample, y=log2(count_norm+1),fill=Sample))+theme_bw()+geom_violin(alpha=0.8)+geom_boxplot(outlier.shape = NA, alpha=0.2)+
    scale_fill_brewer(palette="RdPu",direction = 1)+facet_wrap(~CRTFs)
  
  #plots from custom hand picked gene set of interest
  CRISPR.long$CRTFs = CRISPR.long$Gene %in% c("MYC","ATOH1","ISL1","LHX3","POU4F3","SOX2","YY1","TGIF1","HES1","JUNB","INSM1","SMAD3")
  
  ggplot(data=subset(CRISPR.long,CRTFs=="TRUE"), aes(x=Sample, y=log2(count_norm+1), fill=Sample))+scale_fill_brewer(palette="RdPu",direction = 1)+
    theme_bw()+geom_violin()+geom_boxplot(outlier.shape = NA, alpha=0.2)+facet_wrap(~Gene)
  
  ggplot(data=subset(CRISPR.long,Gene=="INSM1"), aes(x=Sample, y=log2(count_norm+1), group=sgRNA, color=sgRNA))+scale_color_brewer(palette="Spectral",direction = 1)+
    theme_bw()+geom_line()

  