#!/bin/bash
#
# -S /bin/bash
#$ -N CD-HIT
#$ -cwd
#$ -j n
#$ -m e
#$ -M m.brasseur@leibniz-zfmk.de

###gammarus
/home/mbrasseur/bin/cd-hit-v4.8.1-2019-0228/cd-hit-est -i ../gammarus/stranded_spades_gammarus_indoor/transcripts.fasta -o ../gammarus/stranded_spades_gammarus_indoor.fasta.clstr-98 -c 0.98 -M 20000 -T 10 -n 11 -d 0 -g 1 -sc 1 -sf 1 -r 0

###ephemera
/home/mbrasseur/bin/cd-hit-v4.8.1-2019-0228/cd-hit-est -i ../ephemera/stranded_trinity_ephemera_indoor.Trinity.fasta -o ../ephemera/stranded_trinity_ephemera_indoor.Trinity.fasta.clstr-98 -c 0.98 -M 20000 -T 10 -n 11 -d 0 -g 1 -sc 1 -sf 1 -r 0

###lepidostoma
/home/mbrasseur/bin/cd-hit-v4.8.1-2019-0228/cd-hit-est -i ../lepidostoma/stranded_trinity_lepidostoma_indoor.Trinity.fasta -o ../lepidostoma/stranded_trinity_lepidostoma_indoor.Trinity.fasta.clstr-98 -c 0.98 -M 20000 -T 10 -n 11 -d 0 -g 1 -sc 1 -sf 1 -r 0
