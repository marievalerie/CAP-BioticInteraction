#!/bin/bash
#
#$ -S /bin/bash
#$ -N BUSCO
#$ -cwd
#$ -j n
#$ -m e
#$ -M m.brasseur@leibniz-zfmk.de

module load busco/4.0.6

##the output directory must not exist beforehand or use -f (force); no '/' allowed in output name

##lep
busco -i lepidostoma/stranded_trinity_lepidostoma_indoor.Trinity.fasta.clstr-98 -l BUSCO/insecta_odb10 -o lepidostoma_BUSCO_insecta_trin_ind -m transcriptome --offline -f --cpu 15

##eph
busco -i ephemera/stranded_trinity_ephemera_indoor.Trinity.fasta.clstr-98 -l BUSCO/insecta_odb10 -o ephemera_BUSCO_insecta_trin_ind -m transcriptome --offline -f --cpu 15

##gam
busco -i gammarus/stranded_spades_gammarus_indoor.fasta.clstr-98 -l BUSCO/arthropoda_odb10 -o gammarus_BUSCO_arthropoda_spad_ind -m transcriptome --offline -f --cpu 15
