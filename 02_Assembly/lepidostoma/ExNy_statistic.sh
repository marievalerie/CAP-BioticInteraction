#!/bin/bash
#
#$ -S /bin/bash
#$ -N ExNy
#$ -cwd
#$ -j n
#$ -m e
#$ -M m.brasseur@leibniz-zfmk.de

#module load trinity/2.9.0
#module load R/3.6.2

module load trinity/2.13.2
module load R/4.1.2


#if salmon was used to generate abundances, use 
#'find . -maxdepth 2 -name "quant.sf" | tee salmon.quant_files.txt' 
#to prepare sample lists

#if clustered assembly; new gene trans map with:
$TRINITY_HOME/util/support_scripts/get_Trinity_gene_to_trans_map.pl ../stranded_trinity_lepidostoma_indoor.Trinity.fasta.clstr-98 > ../stranded_trinity_lepidostoma_indoor.Trinity.fasta.clstr-98.gene_trans_map


###lepidostoma abundance matrix
$TRINITY_HOME/util/abundance_estimates_to_matrix.pl --est_method salmon --gene_trans_map ../stranded_trinity_lepidostoma_indoor.Trinity.fasta.clstr-98.gene_trans_map --quant_files salmon.quant_files.txt --name_sample_by_basedir


####transcripts
##transcript counts as a functin of min TPM value
$TRINITY_HOME/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl salmon.isoform.TPM.not_cross_norm | tee trans_matrix.TPM.not_cross_norm.counts_by_min_TPM


##estimate ExN50 stats - transcripts
$TRINITY_HOME/util/misc/contig_ExN50_statistic.pl salmon.isoform.TMM.EXPR.matrix ../stranded_trinity_lepidostoma_indoor.Trinity.fasta.clstr-98 transcript | tee ExN50_trinity_assembly.stats


##plot ExN50 stats; tidyverse not available from module env; run this from head node
Rscript $TRINITY_HOME/util/misc/plot_ExN50_statistic.Rscript ExN50_trinity_assembly.stats


#genes
$TRINITY_HOME/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl salmon.gene.TPM.not_cross_norm | tee genes_matrix.TPM.not_cross_norm.counts_by_min_TPM
$TRINITY_HOME/util/misc/contig_ExN50_statistic.pl salmon.isoform.TMM.EXPR.matrix ../stranded_trinity_lepidostoma_indoor.Trinity.fasta.clstr-98 gene | tee ExN50.gene.stats


###new release of the script
/home/mbrasseur/bin/trinityrnaseq-v2.15.1/util/misc/contig_ExN50_statistic.pl salmon.isoform.TMM.EXPR.matrix ../stranded_trinity_lepidostoma_indoor.Trinity.fasta.clstr-98 transcript | tee ExN50.transcript.stats
{TRINITY_HOME}/util/misc/plot_ExN50_statistic.Rscript  ExN50.gene.stats 
