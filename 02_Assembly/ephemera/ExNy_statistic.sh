#!/bin/bash
#
#$ -S /bin/bash
#$ -N ExNy
#$ -cwd
#$ -j n
#$ -m e
#$ -M m.brasseur@leibniz-zfmk.de

module load trinity/2.11.0
module load R/3.6.2

#if salmon was used to generate abundances, use 
#'find . -maxdepth 2 -name "quant.sf" | tee salmon.quant_files.txt' 
#to prepare sample lists

###ephemera abundance matrix
$TRINITY_HOME/util/abundance_estimates_to_matrix.pl --est_method salmon --gene_trans_map none --quant_files salmon.quant_files.txt --name_sample_by_basedir


##transcript counts as a functin of min TPM value
$TRINITY_HOME/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl salmon.isoform.TPM.not_cross_norm | tee trans_matrix.TPM.not_cross_norm.counts_by_min_TPM


##estimate ExN50 stats
$TRINITY_HOME/util/misc/contig_ExN50_statistic.pl salmon.isoform.TMM.EXPR.matrix ../stranded_trinity_ephemera_indoor.Trinity.fasta.clstr-98 | tee ExN50_trinity_assembly.stats


##plot ExN50 stats; module env is not working bc tidyverse is not installed -> use home R libs
Rscript $TRINITY_HOME/util/misc/plot_ExN50_statistic.Rscript ExN50_trinity_assembly.stats


##how many genes correspond to the Ex90 peak:
cat RSEM.isoform.TMM.EXPR.matrix.E-inputs | egrep -v ^\# | awk '$1 <= 90' | wc -l

####try with new release of the script for differences in transcript and gene level
/home/mbrasseur/bin/trinityrnaseq-v2.15.1/util/misc/contig_ExN50_statistic.pl salmon.isoform.TMM.EXPR.matrix ../stranded_trinity_ephemera_indoor.Trinity.fasta transcript | tee ExN50.transcript.stats

/home/mbrasseur/bin/trinityrnaseq-v2.15.1/util/misc/contig_ExN50_statistic.pl salmon.isoform.TMM.EXPR.matrix ../stranded_trinity_ephemera_indoor.Trinity.fasta gene | tee ExN50.gene.stats
