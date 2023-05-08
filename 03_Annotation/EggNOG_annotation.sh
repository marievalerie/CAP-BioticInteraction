#!/bin/bash
#
#$ -S /bin/bash
#$ -N eggnog
#$ -cwd 
#$ -j n
#$ -m e
#$ -M m.brasseur@leibniz-zfmk.de

module load anaconda3/2022.05
conda activate eggnog

emapper.py -i gammarus/stranded_spades_gammarus_indoor.fasta.clstr-98.transdecoder.pep -o gammarus/eggnog_mapper_results \
--cpu $NSLOTS --data_dir ~/bin/eggnog_db --pident 0.5 --evalue 0.00001 --score 50 \
--go_evidence all --excel --target_taxa 6656

emapper.py -i lepidostoma/stranded_trinity_lepidostoma_indoor.Trinity.fasta.clstr-98.transdecoder.pep -o lepidostoma/eggnog_mapper_results \
--cpu $NSLOTS --data_dir ~/bin/eggnog_db --pident 0.5 --evalue 0.00001 --score 50 \
--go_evidence all --excel --target_taxa 6656

emapper.py -i ephemera/stranded_trinity_ephemera_indoor.Trinity.clstr-98.fasta.transdecoder.pep -o ephemera/eggnog_mapper_results \
--cpu $NSLOTS --data_dir ~/bin/eggnog_db --pident 0.5 --evalue 0.00001 --score 50 \
--go_evidence all --excel --target_taxa 6656
