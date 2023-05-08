#!/bin/bash
#
#$ -S /bin/bash
#$ -N spades_gam_indoor
#$ -cwd
#$ -j n
#$ -m e
#$ -M m.brasseur@leibniz-zfmk.de

module load rnaspades/3.15.0

#rnaspades.py --dataset rnaSPAdes.yaml -m 800 -t $NSLOTS --ss-rf -o ../stranded_spades_gammarus_indoor #initial command, but more memory was required

rnaspades.py -m 1000 -t $NSLOTS -o ../stranded_spades_gammarus_indoor --restart-from last
