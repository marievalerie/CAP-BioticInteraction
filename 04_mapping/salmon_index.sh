#!/bin/bash
#
#$ -S /bin/bash
#$ -N salmon_idx
#$ -cwd
#$ -j n
#$ -m e
#$ -M m.brasseur@leibniz-zfmk.de

module load anaconda3/2022.05
conda activate salmon

salmon index -t ../gammarus/stranded_spades_gammarus_indoor.fasta.clstr-98 -i ../gammarus/spades_gam_indoor_clstr-98_salmon_idx -p 15

salmon index -t ../lepidostoma/stranded_trinity_lepidostoma_indoor.Trinity.fasta.clstr-98 -i ../lepidostoma/trinity_lep_indoor_clstr-98_salmon_idx -p 15

salmon index -t ../ephemera/stranded_trinity_ephemera_indoor.Trinity.fasta.clstr-98 -i ../ephemera/trinity_eph_indoor_clstr-98_salmon_idx -p 15
