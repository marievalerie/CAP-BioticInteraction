#Lepidostoma basale - Indoor functional enrichment of GO IDs (biological process terms)


##enrichment analysis
library(tidyr)
library(GO.db)
library(topGO)
##plotting
library(ggplot2)
library(scales)

setwd('C:/Users/mbras/OneDrive/Desktop/Doktorarbeit/Indoor_Genexpression/functional_annotation_eggnog/lepidostoma_trinity')

##gene2go must be prepared
GO_ids = read.csv('Trinity_Gene2GO.csv', sep=';', header = F)

##transform to longformat
long_GO <- gather(GO_ids, Gen, IDs, V2:V673)

##take out genes without GO terms
long_GO <- long_GO[which(long_GO$`IDs` != ""),] 

##remove variable column
long_GO <- long_GO[, c(1, 3)]

##sort by transcript/gene
gene.go <- long_GO[order(long_GO$V1), ]

##Create list with element for each gene, containing vectors with all terms for each gene
gene2GO <- tapply(gene.go$`IDs`, gene.go$V1, function(x)x)

head(gene2GO) #go IDs as strings

##load the results from the differential gene expression analysis
setwd('C:/Users/mbras/OneDrive/Desktop/Doktorarbeit/Indoor_Genexpression/DESeq/lep_numerical_pred')
pest <- as.data.frame(read.csv('DEG_lep_indoor_numerical_predictor_SVA_corrected_concentration.csv', header = T, row.names = 1))
pest_under_bi <- as.data.frame(read.csv('DEG_lep_indoor_numerical_predictor_SVA_corrected_pest_under_bi.csv', header = T, row.names = 1))
setwd('C:/Users/mbras/OneDrive/Desktop/Doktorarbeit/Indoor_Genexpression/functional_annotation_eggnog/lepidostoma_trinity/enrichment')

###insecticide effect without biotic interaction###

##Define vector that is 1 if gene is significantly DE (adj. p.val < chosen pcutoff) and 0 otherwise
DE <- pest
DE$X <- row.names(DE)
DE <- DE[, c(7,2,6)]
DE$padj[is.na(DE$padj)] <- 1

head(DE)
pcutoff = 0.05 #you can adjust this

##differentiate between up and downregulated 
DE_up <- DE
DE_up$up <-ifelse(DE$log2FoldChange < 0, 0, 1)
DE_up$padj <- ifelse(DE_up$padj < pcutoff, 1, 0)
tmp <- ifelse(DE_up$padj == 1 & DE_up$up == 1, 1, 0)
geneList_up <- tmp


DE_down <- DE
DE_down$down <-ifelse(DE_down$log2FoldChange > 0, 0, 1)
DE_down$padj <- ifelse(DE_down$padj < pcutoff, 1, 0)
tmp <- ifelse(DE_down$padj == 1 & DE_down$down == 1, 1, 0)
geneList_down <- tmp


##geneList need the same names as in the match for GO terms (gene2GO)
names(geneList_up) <- unlist(lapply(DE_up$X, function(x)x[1]))
names(geneList_down) <- unlist(lapply(DE_down$X, function(x)x[1]))

##Create topGOdata object: downregulated genes
GOdata_down <- new("topGOdata",
                   ontology = "BP",  #ontology criteria = biological process
                   allGenes = geneList_down, #gene universe = all genes which were tested for DE (i.e. used in DESeq)
                   geneSelectionFun = function(x)(x == 1), 
                   annot = annFUN.gene2GO, gene2GO = gene2GO) #gene ID-to-GO terms

##run enrichment test: here Fishers Exact Test
resultFisher.weight.down <- runTest(GOdata_down, algorithm = "weight01", statistic = "fisher")

##summarize in table
down <- GenTable(GOdata_down, weight01 = resultFisher.weight.down, topNodes = length(resultFisher.weight.down@score),  orderBy = 'weight01', numChar = 500)

##export results
write.csv(down, file = 'BP_enrichment_weight01_fisher_insecticide_downregulated.csv')

##upregulated genes
GOdata_up <- new("topGOdata",
                 ontology = "BP",  #ontology criteria = biological process
                 allGenes = geneList_up, #gene universe = all genes which were tested for DE (i.e. used in DESeq)
                 geneSelectionFun = function(x)(x == 1), 
                 annot = annFUN.gene2GO, gene2GO = gene2GO) #gene ID-to-GO terms

resultFisher.weight.up <- runTest(GOdata_up, algorithm = "weight01", statistic = "fisher")
up <- GenTable(GOdata_up, weight01 = resultFisher.weight.up,  topNodes = length(resultFisher.weight.up@score),  orderBy = 'weight01', numChar = 500)

##export results
write.csv(up, file = 'BP_enrichment_weight01_fisher_insecticide_upregulated.csv')


##plot enrichment results for the most significant terms
ntop <- 15 #15 most significant terms

##downregulated genes
ggdata <- down[1:ntop,]
ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term)) # fixes order
ggdata$weight01 <- as.numeric(ggdata$weight01)

pdf(file= 'GOenrichment_pesticide_downregulated_BP_top15.pdf', height = 9, width = 14)
gg1 <- ggplot(ggdata, aes(x = Term, y = Significant, size = Significant/Expected, fill = -log10(weight01))) +
  expand_limits(y = 7) +
  geom_point(shape = 21) +  #shape = 21 allows to fill with other color
  scale_size(range = c(1,10)) +#, breaks = c(1, 4, 10)) +
  scale_fill_continuous(low = 'yellow', high = 'darkblue') +
  guides(size=guide_legend("Enrichment ratio"))+ 
  
  xlab('') + ylab('Annotated genes') + 
  labs(
    title = 'Insecticide effect L. basale -\ndownregulated genes',
    subtitle = 'Top 15 Biological Process terms')+
  
  theme_bw(base_size = 28) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 19, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 15, face = 'bold', vjust = 1),
    
    axis.text.x = element_text(angle = 0, size = 15, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 15, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 15, face = 'bold'),
    axis.title.x = element_text(size = 15, face = 'bold'),
    axis.title.y = element_text(size = 15, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), #removes the border
    legend.key.size = unit(1, "cm"), #Sets overall area/size of the legend
    legend.text = element_text(size = 19, face = "bold"), # Text size
    title = element_text(size = 19, face = "bold")) +
  
  #change legend title
  labs(fill = '-log10(p-value)') +
  
  coord_flip()
gg1
dev.off()


#upregulated genes
ggdata <- up[1:ntop,]
ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term)) # fixes order
ggdata$weight01 <- as.numeric(ggdata$weight01)

pdf(file= 'GOenrichment_pesticide_upregulated_BP_top15.pdf',  height = 9, width = 14)
gg1 <- ggplot(ggdata,
              aes(x = Term, y = Significant, size = Significant/Expected, fill = -log10(weight01))) +
  expand_limits(y = 7) +
  geom_point(shape = 21) +  #shape = 21 allows to fill with other color
  scale_size(range = c(1,10)) +#, breaks = c(1, 4, 10)) +
  scale_fill_continuous(low = 'yellow', high = 'darkred') +
  guides(size=guide_legend("Enrichment ratio"))+ 
  
  xlab('') + ylab('Annotated genes') + 
  labs(
    title = 'Insecticide effect L. basale -\nupregulated genes',
    subtitle = 'Top 15 Biological Process terms')+
  
  theme_bw(base_size = 28) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 19, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 15, face = 'bold', vjust = 1),
    
    axis.text.x = element_text(angle = 0, size = 15, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 15, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 15, face = 'bold'),
    axis.title.x = element_text(size = 15, face = 'bold'),
    axis.title.y = element_text(size = 15, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 19, face = "bold"), # Text size
    title = element_text(size = 19, face = "bold")) +
  
  #change legend title
  labs(fill = '-log10(p-value)') +
  
  coord_flip()
gg1
dev.off()


###insecticide effect under biotic interaction###

##Define vector that is 1 if gene is significantly DE (adj. p.val < chosen pcutoff) and 0 otherwise
DE <- pest_under_bi
DE$X <- row.names(DE)
DE <- DE[, c(7,2,6)]
DE$padj[is.na(DE$padj)] <- 1

##differentiate between up and downregulated 
DE_up <- DE
DE_up$up <-ifelse(DE$log2FoldChange < 0, 0, 1)
DE_up$padj <- ifelse(DE_up$padj < pcutoff, 1, 0)
tmp <- ifelse(DE_up$padj == 1 & DE_up$up == 1, 1, 0)
geneList_up <- tmp

DE_down <- DE
DE_down$down <-ifelse(DE_down$log2FoldChange > 0, 0, 1)
DE_down$padj <- ifelse(DE_down$padj < pcutoff, 1, 0)
tmp <- ifelse(DE_down$padj == 1 & DE_down$down == 1, 1, 0)
geneList_down <- tmp


##geneList need the same names as in the match for GO terms (gene2GO)
names(geneList_up) <- unlist(lapply(DE_up$X, function(x)x[1]))
names(geneList_down) <- unlist(lapply(DE_down$X, function(x)x[1]))

##Create topGOdata object: downregulated genes
GOdata_down <- new("topGOdata",
                   ontology = "BP",  #ontology criteria = biological process
                   allGenes = geneList_down, #gene universe = all genes which were tested for DE (i.e. used in DESeq)
                   geneSelectionFun = function(x)(x == 1), 
                   annot = annFUN.gene2GO, gene2GO = gene2GO) #gene ID-to-GO terms

##run enrichment test: here Fishers Exact Test
resultFisher.weight.down <- runTest(GOdata_down, algorithm = "weight01", statistic = "fisher")

##summarize in table
down <- GenTable(GOdata_down, weight01 = resultFisher.weight.down, topNodes = length(resultFisher.weight.down@score),  orderBy = 'weight01', numChar = 500)

##export results
write.csv(down, file = 'BP_enrichment_weight01_fisher_insecticide_under_bi_downregulated.csv')

##upregulated genes
GOdata_up <- new("topGOdata",
                 ontology = "BP",  #ontology criteria = biological process
                 allGenes = geneList_up, #gene universe because all genes are present in the list; here only the ones with read abundance threshold are included because DESeq only worked with them
                 geneSelectionFun = function(x)(x == 1), ##function allows to use the genes which are de for treatment
                 annot = annFUN.gene2GO, gene2GO = gene2GO) #gene ID-to-GO terms

resultFisher.weight.up <- runTest(GOdata_up, algorithm = "weight01", statistic = "fisher")
up <- GenTable(GOdata_up, weight01 = resultFisher.weight.up, topNodes = length(resultFisher.weight.up@score), orderBy = 'weight01', numChar = 500)

##export results
write.csv(up, file = 'BP_enrichment_weight01_fisher_insecticide_under_bi_upregulated.csv')

##plot enrichment results for the most significant terms
ntop <- 15 #15 most significant terms

##downregulated genes
ggdata <- down[1:ntop,]
ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term)) # fixes order
ggdata$weight01 <- as.numeric(ggdata$weight01)

pdf(file= 'GOenrichment_pesticide_under_bi_downregulated_BP_top15.pdf', height = 9, width = 14)
gg1 <- ggplot(ggdata, aes(x = Term, y = Significant, size = Significant/Expected, fill = -log10(weight01))) +
  expand_limits(y = 7) +
  geom_point(shape = 21) +  #shape = 21 allows to fill with other color
  scale_size(range = c(1,10)) +#, breaks = c(1, 4, 10)) +
  scale_fill_continuous(low = 'yellow', high = 'darkblue') +
  guides(size=guide_legend("Enrichment ratio"))+ 
  
  xlab('') + ylab('Annotated genes') + 
  labs(
    title = 'Insecticide effect L. basale -\ndownregulated genes under biotic interaction',
    subtitle = 'Top 15 Biological Process terms')+
  
  theme_bw(base_size = 28) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 19, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 15, face = 'bold', vjust = 1),

    axis.text.x = element_text(angle = 0, size = 15, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 15, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 15, face = 'bold'),
    axis.title.x = element_text(size = 15, face = 'bold'),
    axis.title.y = element_text(size = 15, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), #removes the border
    legend.key.size = unit(1, "cm"), #Sets overall area/size of the legend
    legend.text = element_text(size = 19, face = "bold"), # Text size
    title = element_text(size = 19, face = "bold")) +
  
  #change legend title
  labs(fill = '-log10(p-value)') +
  
  coord_flip()
gg1
dev.off()


#upregulated genes
ggdata <- up[1:ntop,]
ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term)) # fixes order
ggdata$weight01 <- as.numeric(ggdata$weight01)

pdf(file= 'GOenrichment_pesticide_under_bi_upregulated_BP_top15.pdf',  height = 9, width = 14)
gg1 <- ggplot(ggdata,
              aes(x = Term, y = Significant, size = Significant/Expected, fill = -log10(weight01))) +
  expand_limits(y = 7) +
  geom_point(shape = 21) +  #shape = 21 allows to fill with other color
  scale_size(range = c(1,10)) +#, breaks = c(1, 4, 10)) +
  scale_fill_continuous(low = 'yellow', high = 'darkred') +
  guides(size=guide_legend("Enrichment ratio"))+ 

  xlab('') + ylab('Annotated genes') + 
  labs(
    title = 'Insecticide effect L. basale -\nupregulated genes under biotic interaction',
    subtitle = 'Top 15 Biological Process terms')+
  
  theme_bw(base_size = 28) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 19, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 15, face = 'bold', vjust = 1),

    axis.text.x = element_text(angle = 0, size = 15, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 15, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 15, face = 'bold'),
    axis.title.x = element_text(size = 15, face = 'bold'),
    axis.title.y = element_text(size = 15, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 19, face = "bold"), # Text size
    title = element_text(size = 19, face = "bold")) +
  
  #change legend title
  labs(fill = '-log10(p-value)') +
  
  coord_flip()
gg1
dev.off()
