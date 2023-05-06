#Epehemera danica - Indoor differential gene expression analysis with DESeq2


#analysis
library(tximport)
library(readr)
library(DESeq2)
library(sva)
library(PCAtools)
library(apeglm)
library(ashr)
library(stringr)
#plotting
library(devtools)
library(viridis)
library(reshape2)
library(RColorBrewer)
library(circlize)
library(ggplot2)
library(dplyr)


##set alpha for differential expression analysis
sig_lvl <- .05

##load file which points to individual salmon quantification files
setwd('/share/pool/mbrasseur/indoor/ephemera/salmon_quant/')
files = c(read.table('salmon.quant_files.txt'))
files <- c(files$V1)

##get file order to map library name to experimental replicate ID
order <- matrix(unlist(str_split(files, '/')), ncol = 3, byrow = TRUE)[,2]

##for ephemera seq run1: library name and experimental ID
channel <- read.table('/home/mbrasseur/ephemera_indoor/DESeq/genome_guided/indoor_seq_run1_ephemera_samples.txt', col.names = c('seq_batch', 'lib_name', 'channel'))
row.names(channel) <- channel$lib_name
channel <- channel[order,]

##apply experiment ID as file names
sample_names <- gsub('Eph-', '', channel$channel)
sample_names <- gsub('-Ind', '', sample_names)
names(files) <- sample_names

##check if files exist
all(file.exists(files))

##create tx2gene
tx2gene <- as.data.frame(read.table("/share/pool/mbrasseur/indoor/ephemera/stranded_trinity_ephemera_indoor.Trinity.fasta.gene_trans_map"))

##this must be in format 'tx ID, then gene ID' -> trinity-table must be reversed
tx2gene <- tx2gene[, c(2,1)]

##drop file of sample 24 (~35 % reads mapped to genome)
files <- files[-5]

##import data
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)


##change to DESeq wd
setwd('/share/pool/mbrasseur/indoor/Deseq/ephemera')

sampleData <- read.csv('../sample_treatments_indoor.csv', header = T, sep = ';', row.names = 1)
colnames(sampleData) <- c('Pesticide', 'Biotic_interaction', 'Treatment')

sampleData$Pesticide <- factor(sampleData$Pesticide, levels = c('control', 'low', 'medium', 'high')) #factor levels for plotting
sampleData$Biotic_interaction <- factor(sampleData$Biotic_interaction, levels = c('L', 'G', 'GL')) #this will not be included as covariate; only for exploratory analysis

##model expression with numerical predictors 
##beta for significant genes as the log fold change for each unit change in concentration

##conc. from day of sampling
sampleData$Concentration[sampleData$Pesticide == 'control'] <- 0
sampleData$Concentration[sampleData$Pesticide == 'low'] <- 0.3
sampleData$Concentration[sampleData$Pesticide == 'medium'] <- 4.93
sampleData$Concentration[sampleData$Pesticide == 'high'] <- 19.5

##order the data
sampleData <- sampleData[sample_names,]

##drop sample 24 from sampleData
sampleData <-  sampleData[-5,]

##be sure that samples in colData correspond to the right sample order
all(colnames(txi.salmon$counts) == rownames(sampleData))

dds1 <- DESeqDataSetFromTximport(txi.salmon, sampleData, design= ~ Concentration)

##no. of 'genes'
length(row.names(counts(dds1))) #700198

##filtering
dds1 <- estimateSizeFactors(dds1)
nc <- counts(dds1, normalized=TRUE)

filter2 <- rowSums(nc >= 12) >= 9 #Ephemera danica data set -> 9 replicates to estimate insecticide effect
dds1 <- dds1[filter2,]

##no. of 'genes' after filtering
length(row.names(counts(dds1))) #48610 genes

##run DESeq2
dds1 <- DESeq(dds1) 

##check if dispersion model fit is adequate
pdf(plot, file="gene_dispersion_eph_ind_numerical_predictor.pdf")
plot <- plotDispEsts(dds1)
dev.off()


##PCA based on normalized, vst() transformed count data
cols2 = c('grey', "#35B779FF", "#31688EFF","#440154FF")
shps = c(15, 16,17)

vsd <- vst(dds1)
p <- pca(assay(vsd), metadata = colData(dds1))

pdf('PCA_eph_ind_uncorrected.pdf', height = 6, width = 5)
biplot(p, showLoadings = F, sizeLoadingsNames = 2, boxedLoadingsNames = F, legendPosition = 'none', colby = 'Pesticide', colkey = cols2, shape = 'Biotic_interaction', shapekey = shps, title = 'E. danica', lab = NULL) 
dev.off()


##SVA analysis for batch correction based on normalized data
dat  <- counts(dds1, normalized = TRUE)

##specify full model
mod  <- model.matrix(~ Concentration, colData(dds1))

##specify reduced model
mod0 <- model.matrix(~   1, colData(dds1))
svseq <- svaseq(dat, mod, mod0) #estimates 7 surrogates

##perform frozen SVA to regress out latent factor (adjusted only for visualization)
newV = NULL #neccessary due to bug in the sva pacakge
dat_vst  <- vst(dds1) 
fsvaobj <- fsva(assay(dat_vst), mod, svseq, method = 'exact') 

##get adjusted data
data_adj <- fsvaobj$db

##PCA with adjusted data
p <- pca(data_adj, metadata = colData(dds1))

pdf('PCA_eph_ind_SVA_corrected.pdf', height = 6, width = 5)
biplot(p, showLoadings = F, sizeLoadingsNames = 2, boxedLoadingsNames = F, legendPosition = 'none', colby = 'Pesticide', colkey = cols2, shape = 'Biotic_interaction', shapekey = shps, lab = NULL, title = 'E. danica') 
dev.off()


##include latent factors as covariates in model
ddssva <- dds1
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
ddssva$SV3 <- svseq$sv[,3]
ddssva$SV4 <- svseq$sv[,4]   
ddssva$SV5 <- svseq$sv[,5]
ddssva$SV6 <- svseq$sv[,6]
ddssva$SV7 <- svseq$sv[,7]
design(ddssva) <- ~ SV1+SV2+SV3+SV4+SV5+SV6+SV7+Concentration

##reestimate parameters
ddssva <- DESeq(ddssva) #436 genes did not converge

##remove non-convergingn genes
ddsClean <- ddssva[which(mcols(ddssva)$betaConv),]

##check again dispersion
pdf(plot, file="gene_dispersion_SVA_corrected.pdf")
plot <- plotDispEsts(ddsClean)
dev.off()

##test for genes changing the expression as a function of insecticide concentration
conc <- results(ddsClean, name = "Concentration", lfcThreshold = 0, alpha = sig_lvl)

##effect size shrinkage
conc_ashr <- lfcShrink(ddsClean, type = 'ashr', res = conc)

##export results
write.csv(as.data.frame(conc), file= "DEG_eph_indoor_numerical_predictor_SVA_corrected_conc.csv")
write.csv(as.data.frame(conc_ashr), file= "DEG_eph_indoor_numerical_predictor_SVA_corrected_conc_ashr.csv")


#############Volcano plot#############

##read the data
conc <- as.data.frame(read.csv('DEG_eph_indoor_numerical_predictor_SVA_corrected_conc_ashr.csv', row.names = 1))

##define chr vector for color mapping
keyvals <- ifelse(
  conc$log2FoldChange < 0, 'royalblue','darkred')
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'darkred'] <- 'upregulated'
names(keyvals)[keyvals == 'black'] <- 'n.s.'
names(keyvals)[keyvals == 'royalblue'] <- 'downregulated'

##add expression info to df
conc <- conc %>% 
  mutate(
    Expression = case_when(log2FoldChange > 0 & padj < sig_lvl ~ "Upregulated\n(1422)",
                           log2FoldChange < 0 & padj < sig_lvl ~ "Downregulated\n(3279)",
                           TRUE ~ "n.s.\n(43473)"))

conc$Expression <- factor(conc$Expression, levels = c("Downregulated\n(3279)","Upregulated\n(1422)",'n.s.\n(43473)'))

##check if everything is right
conc %>% 
  count(Expression)

pdf('ephemera_volcano.pdf', width=4, height = 5)
p1 <- ggplot(conc, aes(log2FoldChange, -log(padj,10))) + # -log10 conversion  
  geom_point(aes(color = Expression), size = 1) + xlim(-0.8, 0.8)  +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"padj"))+
  scale_color_manual(values = c("dodgerblue3","firebrick3", "gray50" )) + theme_bw() +
  theme(legend.title = element_blank(), legend.text=element_text(size=12), legend.position = "bottom", axis.text=element_text(size=12), axis.title=element_text(size=14),
        panel.border = element_rect(linewidth = 1.5), panel.grid.minor = element_line(linewidth = 0.75), panel.grid.major = element_line(linewidth = 0.75))+
  guides(colour = guide_legend(override.aes = list(size=4)))
p1
dev.off()
