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

##gammarus abundance matrix
$TRINITY_HOME/util/abundance_estimates_to_matrix.pl --est_method salmon --gene_trans_map none --quant_files salmon.quant_files.txt --name_sample_by_basedir


##estimate ExN50 stats - transcripts
$TRINITY_HOME/util/misc/contig_ExN50_statistic.pl salmon.isoform.TMM.EXPR.matrix ../stranded_spades_gammarus_indoor.fasta.clstr-98 transcript | tee ExN50_spades_assembly.stats


##plot ExN50 stats; tidyverse not available from module env; run this from head node
Rscript $TRINITY_HOME/util/misc/plot_ExN50_statistic.Rscript ExN50_spades_assembly.stats
