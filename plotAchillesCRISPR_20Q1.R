########################
# Achilles CRISPR data #
# Berkley Gryder 2020  #
# gryderart@gmail.com  #
########################

#########################
## 1. Setup environment # 
#########################

    setwd = "T:/Berkley/Experiments/"
    project.folder = "CRISPR/Achilles/Avana20Q1/"
    project.name = "RMS_TFs"
    
    #CRISPR data, gene effect level
          gene.effect = read.table(file = paste(project.folder,"Achilles_gene_effect.csv",sep=""),sep=",",header = T)  ###RUN THIS TO LOAD THE FIRST TIME ONLY, takes a while
          names(gene.effect)[1]<-"DepMap_ID"
          colnames(gene.effect) = gsub("\\..*","",colnames(gene.effect))
          duplicated.Achilles.genes = readLines(paste(project.folder,"duplicated.Achilles.genes.txt",sep=""))  ###gene list available in github
          gene.effect = gene.effect[ , -which(names(gene.effect) %in% duplicated.Achilles.genes)]
          
    #sample metadata (this was downloaded from https://depmap.org/portal/download/all/ as a text file, and edited in excel to subdivide Sarcomas)
          require(XLConnect)
          library(tidyr)
          wb = loadWorkbook(paste(project.folder,"sample_info_20Q1.xlsx",sep=""))
          sample.metadata = readWorksheet(wb, sheet = "sample_info", header = TRUE)
          #only keep metadata for cell lines present in Achilles data
          sample.metadata = subset(sample.metadata, sample.metadata$DepMap_ID %in% gene.effect$DepMap_ID)
          #calculate frequency per tumor type
          n.disease<-as.data.frame(table(sample.metadata$disease)); colnames(n.disease)=c("disease","Freq")
          library(plyr); sample.metadata = join(sample.metadata,n.disease, by = "disease")
          
          #make a VLOOKUP table for cell-line name to get $type
          metadata.lookup = data.frame(sample.metadata$DepMap_ID,sample.metadata$stripped_cell_line_name,sample.metadata$disease,sample.metadata$Freq)
          colnames(metadata.lookup) = c("DepMap_ID","stripped_cell_line_name","tumor_type","type.count")

          #visualize metadata and sample sizes
          library(ggplot2)
          ggplot(metadata.lookup, aes(x=tumor_type, color=tumor_type))+geom_bar(stat="count")+theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="none")+ 
            stat_count(binwidth=1, geom="text", aes(label=..count..), vjust=-0.5)+ggtitle("DepMap Cell lines with CRISPR data")
          
########################
## 2. DotPlot a subset #
########################
      ####SUBSET BASED ON SELECTED TUMOR TYPES
      samples.subset = subset(sample.metadata, sample.metadata$disease %in% c("FP Rhabdomyosarcoma","FN Rhabdomyosarcoma"))
      gene.effect.subset = subset(gene.effect, gene.effect$DepMap_ID %in% samples.subset$DepMap_ID)
      #replace gene.effect columns by reference
          gene.effect.subset = join(gene.effect.subset, metadata.lookup, by = "DepMap_ID")
          gene.effect.subset$DepMap_ID = gene.effect.subset$stripped_cell_line_name
          
      #make long format
        library(tidyr)
        gene.subset.long  <- gather(gene.effect.subset, gene, score, 2:(ncol(gene.effect)))
      
      #Subset by selection of genes
        gene.long.selectgenes = subset(gene.subset.long, gene.subset.long$gene %in% c("FOXO1","MYOD1","MYOG","PAX3","PAX7"))
        gene.long.fororder = subset(gene.subset.long, gene.subset.long$gene %in% c("MYOG"))
        #order samples by one gene score
        gene.long.selectgenes$DepMap_ID <- factor(gene.long.selectgenes$DepMap_ID, levels = gene.long.selectgenes$DepMap_ID[order(gene.long.fororder$score, decreasing = TRUE)])
        #order genes by score
        gene.long.selectgenes$gene = as.factor(gene.long.selectgenes$gene)
        gene.order = with(gene.long.selectgenes, reorder(gene, score, median))
        gene.long.selectgenes$gene <- factor(gene.long.selectgenes$gene, levels = gene.long.selectgenes$gene[order(gene.order, decreasing = TRUE)])
        
        #Ranked dot plot 
        library(ggplot2)
        ggplot(gene.long.selectgenes, aes(x=DepMap_ID, y=score, color = tumor_type))+geom_point(aes(size=0.5)) +facet_wrap(~gene,nrow=1)+theme(axis.text.x = element_text(angle = 90, hjust = 1))
        
#####################
## 3. RankBox plots #
#####################
        
      #Make ALL data into Long format
      library(tidyr)
      gene.long  <- gather(gene.effect, gene, score, 2:(ncol(gene.effect)))
      gene.long$score[is.na(gene.long$score)] <- 0
      
      #Subset through a list
      gene.picks = c("MYOD1","MYOG","PAX3","PAX7","SOX8","SNAI2")
      lapply(gene.picks, function(x) {
        #####
        gene.pick.long = subset(gene.long, gene.long$gene %in% x)
        colnames(gene.pick.long) = c("DepMap_ID","gene","score" )
        library(plyr)
        gene.pick.long.wmeta <- join(gene.pick.long, metadata.lookup, by = "DepMap_ID")
        gene.pick.long.wmeta$tumor_type = as.factor(gene.pick.long.wmeta$tumor_type)
        gene.pick.long.wmeta = subset(gene.pick.long.wmeta, gene.pick.long.wmeta$type.count > 2)
        
        #plot score, boxplot, y=score x=type
        library(ggplot2)
        output = paste(project.folder, x,".rankbox.pdf",sep = '')
          rank.box.plot = ggplot(gene.pick.long.wmeta, aes(x=reorder(tumor_type, -score, FUN=median), y=score))+geom_boxplot(outlier.shape = NA) +theme(axis.text.x = element_text(angle = 90, hjust = 1))+ggtitle(paste(x," cancer cell line dependencies",sep=""))+scale_y_continuous(limits = c(-2.2,1))
          ggsave(output,plot = rank.box.plot,width = 7, height = 5)
      })
        
      #next, make a stack of plots!
      
      plotdata <- function(x) {
        gene.pick.long = subset(gene.long, gene.long$gene %in% x)
        
        library(plyr)
        gene.pick.long.wmeta <- join(gene.pick.long, metadata.lookup, by = "DepMap_ID")
        gene.pick.long.wmeta$tumor_type = as.factor(gene.pick.long.wmeta$tumor_type)
        gene.pick.long.wmeta = subset(gene.pick.long.wmeta, gene.pick.long.wmeta$type.count > 1)
        
        #plot score, boxplot, y=score x=type
        library(ggplot2)
          ggplot(gene.pick.long.wmeta, aes(x=reorder(tumor_type, -score, FUN=median), y=score))+geom_boxplot(outlier.shape = NA) + scale_y_continuous(limits = c(-1.5,1))+
            theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.title.x=element_blank())+ggtitle(paste(x,"",sep=""))
      }
      
      library(gridExtra)
      output = paste(project.folder, project.name, ".grid.rankbox.pdf",sep = '')
      grid.rank.box.plot=do.call(grid.arrange,lapply(gene.picks, plotdata))
      ggsave(output,plot = grid.rank.box.plot,width = 12, height = 16)
      
      #plot summary across all cancers for gene picks
      gene.pick.long = subset(gene.long, gene.long$gene %in% gene.picks)
      library(plyr)
      gene.pick.long.wmeta <- join(gene.pick.long, metadata.lookup, by = "DepMap_ID")  #######PLEASE FIX: this mixed up tumors and their individual cell line in SNAI2 callout!
      gene.pick.long.wmeta$tumor_type = as.factor(gene.pick.long.wmeta$tumor_type)
      gene.pick.long.wmeta = subset(gene.pick.long.wmeta, gene.pick.long.wmeta$type.count > 1)
      ggplot(gene.pick.long.wmeta, aes(x=gene, y=score))+ geom_violin(scale = "width", width = 0.7)+geom_boxplot(outlier.shape = NA,width=0.3) + scale_y_continuous(limits = c(-1.8,1)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.title.x=element_blank())
      
############################  
## 4. HEAT MAP CRISPR data #
############################
  #for small gene sets, manually select columns
  library(pheatmap)
  library(RColorBrewer)
  gene.heated = data.frame(gene.effect$DepMap_ID, gene.effect$SOX8, gene.effect$PAX3, gene.effect$MYOD1, gene.effect$MYOG, gene.effect$MYCN)  ###this is completely manual and quite clunky
    gene.heated[is.na(gene.heated)] <- 0
    gene.heated.wmeta <- join(gene.heated, metadata.lookup, by = "DepMap_ID")
    gene.heated.wmeta <- subset(gene.heated.wmeta, gene.heated.wmeta$type.count > 2)
    gene.heated.matrix = as.matrix(gene.heated.wmeta[,c(2:(ncol(gene.heated.wmeta)-3))])
    row.names(gene.heated.matrix) = gene.heated.wmeta$tumor_type
    pheatmap(gene.heated.matrix,scale='none',cluster_rows = T,  color = colorRampPalette(c("red", "white", "lightblue"))(50),main=paste(project.name," heatmap",sep=""))
  
  #heatmap by type median
    gene.heated.type<-aggregate(gene.heated.wmeta[,c(2:(ncol(gene.heated.wmeta)-3))],by=list(gene.heated.wmeta$tumor_type),FUN=median)
      colnames(gene.heated.type)[1] <- "tumor_type"
      gene.heated.type.matrix = as.matrix(gene.heated.type[,c(2:(ncol(gene.heated.type)))])
      row.names(gene.heated.type.matrix) = gene.heated.type$tumor_type
      breaksList = seq(-5, 5, by = 0.2)
      pheatmap(gene.heated.type.matrix,scale='column',
               cluster_rows = T,  color = colorRampPalette(c("red", "white", "lightblue"))(50),
               main=paste(project.name," heatmap",sep=""))#breaks=breaksList
      

############################
## 5. Find selective genes #
############################
      library(plyr)
      gene.long.wmeta <- join(gene.long, metadata.lookup, by = "DepMap_ID")  #gene.long is made in step 3.RankBox plots
      gene.long.wmeta$tumor_type = as.factor(gene.long.wmeta$tumor_type)
      gene.long.wmeta = subset(gene.long.wmeta, gene.long.wmeta$type.count > 1)
      gene.long.wmeta$sample_gene = as.factor(do.call(paste, c(gene.long.wmeta[c("DepMap_ID","gene")], sep = "_"))) 
      gene.long.wmeta = gene.long.wmeta[!duplicated(gene.long.wmeta[,c('sample_gene')]),]  #remove if SAMPLE_ID has been duplicated
        #gene.spread.wmeta <- spread(gene.long.wmeta[,c(1,2,3,4)], DepMap_ID, score) #not enough memory!
      gene.long.wmeta$type_gene = as.factor(do.call(paste, c(gene.long.wmeta[c("tumor_type","gene")], sep = "_"))) 
        gene.aggregated<-aggregate(gene.long.wmeta[,c("score")],by=list(gene.long.wmeta$type_gene),FUN=mean)
        library(stringr)
        gene.aggregated[,3:4]<-str_split_fixed(gene.aggregated$Group.1,"_",2)
        gene.aggregated = gene.aggregated[,c(4,3,2)]
        colnames(gene.aggregated)= c("gene","tumor_type","mean_score")
        library(tidyr)
        gene.aggr.spread <- spread(gene.aggregated, tumor_type, mean_score)
        gene.aggr.spread$gene.mean<-rowMeans(gene.aggr.spread[,2:ncol(gene.aggr.spread)])
        gene.aggr.spread$FPRMS.v.mean = gene.aggr.spread$`FP Rhabdomyosarcoma` - gene.aggr.spread$gene.mean
        gene.aggr.spread.top = subset(gene.aggr.spread, gene.aggr.spread$FPRMS.v.mean < (-0.5))
        
        #heat the top!
        library(pheatmap)
        library(RColorBrewer)
        gene.aggr.spread.top.matrix = as.matrix(gene.aggr.spread.top[,c(2:(ncol(gene.aggr.spread.top)-2))])
        row.names(gene.aggr.spread.top.matrix) = gene.aggr.spread.top$gene
        breaksList = seq(-5, 5, by = 0.2)
        pheatmap(gene.aggr.spread.top.matrix,scale='row',cluster_rows = F,  color = colorRampPalette(c("red", "white", "lightblue"))(50),main=paste(project.name," heatmap",sep="")) #,breaks=breaksList

        #scatter
        ggplot(gene.aggr.spread, aes(x=`FP Rhabdomyosarcoma`, y=gene.mean))+geom_point(alpha=0.05)+
                geom_point(data=gene.aggr.spread.top, color="red")+scale_x_continuous(lim=c(-2.2,1))+scale_y_continuous(lim=c(-2.2,1))+
                geom_text(data=gene.aggr.spread.top,label=gene.aggr.spread.top$gene, hjust=1.1, vjust=0.5, size=4) + theme_bw()
        
################################
## 6. Rank order genes by type #
################################
        
        geneset.folder = "K:/projects/ChIP_seq/RNA_DATA/RNA_projects/Genesets/Qlucore_format/"
        geneset.name = "GRYDER_RH4_CR_TFs_CRISPRTop"
        geneset = read.table(file = paste(geneset.folder,geneset.name,".genelist.txt",sep=""),sep="\t",header = F)
        
        #give rank
        gene.aggr.spread = gene.aggr.spread[order(gene.aggr.spread$`FP Rhabdomyosarcoma`), ]
        gene.aggr.spread$rank = as.numeric(factor(gene.aggr.spread$`FP Rhabdomyosarcoma`,levels=unique(gene.aggr.spread$`FP Rhabdomyosarcoma`))) 
            
        #make subset
        gene.aggr.spread.geneset = subset(gene.aggr.spread, gene.aggr.spread$gene %in% geneset$V1)
        gene.aggr.spread$geneset = gene.aggr.spread$gene %in% geneset$V1
        
        #plots
        library(ggplot2)
        ggplot(gene.aggr.spread, aes(x=rank, y=`FP Rhabdomyosarcoma`)) + geom_bar(stat="identity",fill="grey90") +
          geom_line(color="black") + geom_bar(data=gene.aggr.spread.geneset,stat="identity",fill="red") + theme_bw()+ #coord_flip()+
          scale_y_continuous(position = "top") + #scale_x_continuous(lim=c(0,100))+
          labs(x = "Genes Ranked by CERES score") +labs(y = paste("CERES score",sep="")) +  geom_point(data=gene.aggr.spread.geneset)+
          geom_text(data=gene.aggr.spread.geneset,label=gene.aggr.spread.geneset$gene, hjust=-0.2, size=4) 
        
        ggplot(gene.aggr.spread.geneset,aes(x=rank)) + geom_histogram(color="black",fill="red",binwidth = 500)+theme_bw() + scale_x_continuous(lim=c(0,nrow(gene.aggr.spread)))
        
        #stats and boxplot
        geneset.p = signif(t.test(`FP Rhabdomyosarcoma` ~ geneset, data=gene.aggr.spread, var.equal=TRUE, conf.level=0.95)$p.value, digits = 2)
        ggplot(gene.aggr.spread,aes(y=`FP Rhabdomyosarcoma`,x=factor(geneset)))+geom_boxplot(width=0.5,outlier.shape=NA,aes(fill=factor(geneset))) +theme_bw() +
          annotate("text", x = 1.5, y = 1, label = paste("p = ",geneset.p))+ggtitle(label = "Subtype score by geneset", subtitle = geneset.name)
        
        
