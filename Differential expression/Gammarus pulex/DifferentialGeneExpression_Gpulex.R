#Gammarus pulex - Indoor differential gene expression analysis with DESeq2


#analysis
library(tximport)
library(readr)
library(DESeq2)
library(sva)
library(PCAtools)
library(apeglm)
library(ashr)
library(DEGreport)
#plotting
library(devtools)
library(viridis)
library(reshape2)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(dendsort)
library(ggplot2)
library('UpSetR')
library('tidyverse')
library('grid')
library("cowplot")

##set alpha for differential expression analysis
sig_lvl <- .05

##load file which points to individual salmon quantification files
setwd('/share/pool/mbrasseur/indoor/gammarus/salmon_quant/')
files = c(read.table('salmon.quant_files.txt'))
files <- c(files$V1)

##apply experiment ID as file names
sample_names <- gsub('DE.*', '', files)
sample_names <- gsub('./Gam-', '', sample_names)
sample_names <- gsub('-Ind', '', sample_names)
sample_names <- gsub ('./Gam', '', sample_names)
names(files) <- sample_names

##check if files exist
all(file.exists(files))

##create tx2gene
tx2gene <- as.data.frame(read.table("/share/pool/mbrasseur/indoor/gammarus/stranded_spades_gammarus_indoor.fasta.clstr-98.tx2gene.map"))

##this must be in format 'tx ID, then gene ID'
tx2gene <- tx2gene[, c(2,1)]

##import data
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)

##change to DESeq wd
setwd('/share/pool/mbrasseur/indoor/Deseq/gammarus')

sampleData <- read.csv('../sample_treatments_indoor.csv', header = T, sep = ';', row.names = 1)
colnames(sampleData) <- c('Pesticide', 'Biotic_interaction', 'Treatment')

sampleData$Pesticide <- factor(sampleData$Pesticide, levels = c('control', 'low', 'medium', 'high'))
#levels(sampleData$Pesticide)
sampleData$Biotic_interaction <- factor(sampleData$Biotic_interaction, levels = c('G', 'L', 'GL'))


##model expression with numerical predictors 
##beta for significant genes as the log fold change for each unit change in concentration

##conc. from day of sampling
sampleData$Concentration[sampleData$Pesticide == 'control'] <- 0
sampleData$Concentration[sampleData$Pesticide == 'low'] <- 0.3
sampleData$Concentration[sampleData$Pesticide == 'medium'] <- 4.93
sampleData$Concentration[sampleData$Pesticide == 'high'] <- 19.5

##order the data
sampleData <- sampleData[sample_names,]

##be sure that samples in colData correspond to the right sample order
all(colnames(txi.salmon$counts) == rownames(sampleData))

##fit initial glms with DESeq 
dds1 <- DESeqDataSetFromTximport(txi.salmon, sampleData, design= ~ Biotic_interaction*Concentration)

##no. of 'genes'
length(row.names(counts(dds1))) # 760065 genes

##filtering
dds1 <- estimateSizeFactors(dds1)
nc <- counts(dds1, normalized=TRUE)

filter2 <- rowSums(nc >= 12) >= 6 #Gammarus pulex data set -> 6 replicates to estimate insecticide effect
dds1 <- dds1[filter2,]

##no. of 'genes' after filtering
length(row.names(counts(dds1))) #48640 genes

##run DESeq2
dds1 <- DESeq(dds1) 

##check if dispersion model fit is adequate
pdf(plot, file="gene_dispersion_gam_ind_numerical_predictor.pdf")
plot <- plotDispEsts(dds1)
dev.off()


##PCA based on normalized, vst() transformed count data
cols2 = c('grey', "#35B779FF", "#31688EFF","#440154FF")
shps = c(16, 17)

vsd <- vst(dds1)
p <- pca(assay(vsd), metadata = sampleData)

pdf('PCA_gam_ind_uncorrected.pdf', height = 6, width = 5)
biplot(p, showLoadings = F, sizeLoadingsNames = 2, boxedLoadingsNames = F, legendPosition = 'none', colby = 'Pesticide', colkey = cols2, shape = 'Biotic_interaction', shapekey = shps, lab = NULL, title = 'G. pulex') 
dev.off()


##SVA analysis for batch correction based on normalized data
dat  <- counts(dds1, normalized = TRUE)

##specify full model
mod  <- model.matrix(~ Biotic_interaction*Concentration, colData(dds1)) 

##specify null model
mod0 <- model.matrix(~   1, colData(dds1))
svseq <- svaseq(dat, mod, mod0) #estimates 5 surrogates

##perform frozen SVA to regress out latent factor (adjusted only for visualization)
newV = NULL #neccessary due to bug in the sva pacakge
dat_vst  <- vst(dds1) 
fsvaobj <- fsva(assay(dat_vst), mod, svseq, method = 'exact') 

##get adjusted data
data_adj <- fsvaobj$db


##PCA with adjusted data
p <- pca(data_adj, metadata = sampleData)

pdf('PCA_gam_ind_SVA_corrected.pdf', height = 6, width = 5)
biplot(p, showLoadings = F, sizeLoadingsNames = 2, boxedLoadingsNames = F, legendPosition = 'none', colby = 'Pesticide', colkey = cols2, shape = 'Biotic_interaction', shapekey = shps, lab = NULL, title = 'G. pulex') 
dev.off()

##include latent factors as covariates in model
ddssva <- dds1
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
ddssva$SV3 <- svseq$sv[,3]
ddssva$SV4 <- svseq$sv[,4]   
ddssva$SV5 <- svseq$sv[,5] 
design(ddssva) <- ~ SV1+SV2+SV3+SV4+SV5+Biotic_interaction*Concentration

##reestimate parameters
ddssva <- DESeq(ddssva) #1694 did not converge

##increase number of iterations to improve model convergence
ddssva <- nbinomWaldTest(ddssva, maxit=10000) #1559 rows did not converge

##remove non-convergingn genes
ddsClean <- ddssva[which(mcols(ddssva)$betaConv),]

##check again dispersion
pdf(plot, file="gene_dispersion_SVA_corrected.pdf")
plot <- plotDispEsts(ddsClean)
dev.off()

##test for genes changing the expression as a function of insecticide concentration when biotic interaction is at the refernce level (i.e., only G. pulex)
conc <- results(ddsClean, name ="Concentration", alpha = sig_lvl)

##effect size shrinkage
conc_ashr <- lfcShrink(ddsClean, type="ashr", res = conc)

###export results
write.csv(as.data.frame(conc), file= "DEG_gam_indoor_numerical_predictor_SVA_corrected_concentration.csv")
write.csv(as.data.frame(conc_ashr), file= "DEG_gam_indoor_numerical_predictor_SVA_corrected_concentration_ashr.csv")


##high pest effect under biotic interaction = GL
conc_under_bi <- results(ddsClean, contrast = list(c("Concentration", "Biotic_interactionGL.Concentration")), lfcThreshold = 0, alpha = sig_lvl)

##effect size shrinkage
conc_under_bi_ashr <- lfcShrink(ddsClean, type="ashr", res = conc_under_bi)

##export results
write.csv(as.data.frame(conc_under_bi), file= "DEG_gam_indoor_numerical_predictor_SVA_corrected_concentration_under_BI.csv")
write.csv(as.data.frame(conc_under_bi_ashr), file= "DEG_gam_indoor_numerical_predictor_SVA_corrected_concentration_under_BI_ashr.csv")



#############heatmap#############
##read the data with shrunken effect sizes
setwd('C:/Users/mbras/OneDrive/Desktop/Doktorarbeit/Indoor_Genexpression/DESeq/gam_numerical_pred')
conc <- as.data.frame(read.csv('DEG_gam_indoor_numerical_predictor_SVA_corrected_concentration_ashr.csv', row.names = 1))
conc_under_bi <- as.data.frame(read.csv('DEG_gam_indoor_numerical_predictor_SVA_corrected_concentration_under_BI_ashr.csv', row.names = 1))

##combine data
column <- c('conc', 'conc_under_bi')
colnms <- c('rmv')

res <- matrix(NA, ncol = 1, nrow = length(row.names(conc)))
colnames(res) <- 'rmv'
rownames(res) <- row.names(conc)
res <- as.data.frame(res)

for (i in column) {
  res1 <- eval(parse(text = i))
  res1$padj[is.na(res1$padj)] <- 1
  res1 <- res1[res1$padj < sig_lvl, ]
  res1 <- as.data.frame(res1)
  res <- merge(res, res1['log2FoldChange'],by ="row.names", all.x = T)
  rownames(res) <- res$Row.names
  res <- res[,-1]
  colnms <- c(colnms, i)
  colnames(res) <- colnms
  print(i)
  res <- as.data.frame(res)
}

res$rmv <- NULL ##remove dummy variable

res <- res[rowSums(is.na(res)) != ncol(res), ] ##remove non-significant transcripts
res_copy_sub <- res

##highest and lowest LFC
res[is.na(res)] <- 0
max(res)
min(res)

##use the ES from the treatment for clustering
res_hclust_sub <- res

res_hclust_sub$conc <- conc[row.names(res_hclust_sub), 'log2FoldChange']
res_hclust_sub$conc_under_bi <- conc_under_bi[row.names(res_hclust_sub), 'log2FoldChange'] 
row_dend_all = dendsort(hclust(dist(res_hclust_sub)))

##color vector
f1 = colorRamp2(c(-1,-0.5,-0.3, 0,0.3,0.5, 1), c("#003366", "#0066CC", '#99CCFF', "#FFFFFFFF", '#FF0000',"#CC0000", "#990000" ))

##plot heatmap
pdf('gam_only_pesticide_effects_ashr_heatmap.pdf', width = 5, height = 4.5)
Heatmap(as.matrix(res_copy_sub), na_col = 'lightgrey', col = f1, column_names_gp = gpar(fontsize = 11),
        column_labels = c(paste('\nInsecticide\n (', length(which(res_copy_sub$conc < 0)), '/', length(which(res_copy_sub$conc > 0)), ')', sep = ''),
                          paste('\nInsecticide under \nbiotic interaction (', length(which(res_copy_sub$conc_under_bi < 0)), '/', length(which(res_copy_sub$conc_under_bi > 0)), ')', sep = '')),
        column_names_rot = 0, 
        heatmap_legend_param = list(
          title = expression(log[2]~FC), at =c(-1,-0.5, 0,0.5, 1), title_gp = gpar(fontsize = 14),
          labels = c(-1,-0.5, 0,0.5, 1), labels_gp = gpar(fontsize(9)),legend_height = unit(2.5, "cm"),
          col_fun=f1), column_title = gt_render("***G. pulex***"), column_title_side = "bottom", show_heatmap_legend = TRUE, show_row_names = FALSE,
        cluster_columns = F, cluster_rows = row_dend_all)
dev.off()
