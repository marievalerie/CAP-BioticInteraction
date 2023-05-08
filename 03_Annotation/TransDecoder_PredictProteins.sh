#!/bin/bash
#
#$ -S /bin/bash
#$ -N transdecoder
#$ -cwd
#$ -j n
#$ -m e
#$ -M m.brasseur@leibniz-zfmk.de


dirname="/share/pool/mbrasseur/indoor" 

module load hmmer/3.3
module load blast+/2.9.0


##extract long open reading frames (min 100 AA); #strand specific?
~/bin/TransDecoder-TransDecoder-v5.5.0/TransDecoder.LongOrfs -t ${dirname}/lepidostoma/stranded_trinity_lepidostoma_indoor.Trinity.fasta.clstr-98
~/bin/TransDecoder-TransDecoder-v5.5.0/TransDecoder.LongOrfs -t ${dirname}/ephemera/stranded_trinity_ephemera_indoor.Trinity.clstr-98.fasta
~/bin/TransDecoder-TransDecoder-v5.5.0/TransDecoder.LongOrfs -t ${dirname}/gammarus/stranded_spades_gammarus_indoor.fasta.clstr-98


##include homology searches as ORF retention criteria,Swissprot updated
blastp -query ${dirname}/gammarus/stranded_spades_gammarus_indoor.fasta.clstr-98.transdecoder_dir/longest_orfs.pep -db ~/database/2023-uniprot_sprot.fasta -num_threads ${NSLOTS} -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > ${dirname}/gammarus/spades_transcriptome.blastp_for_homology.outfmt6 
blastp -query ${dirname}/ephemera/stranded_trinity_ephemera_indoor.Trinity.clstr-98.fasta.transdecoder_dir/longest_orfs.pep -db ~/database/2023-uniprot_sprot.fasta -num_threads ${NSLOTS} -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > ${dirname}/ephemera/trinity_transcriptome.blastp_for_homology.outfmt6 
blastp -query ${dirname}/lepidostoma/stranded_trinity_lepidostoma_indoor.Trinity.fasta.clstr-98.transdecoder_dir/longest_orfs.pep -db ~/database/2023-uniprot_sprot.fasta -num_threads ${NSLOTS} -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > ${dirname}/lepidostoma/trinity_transcriptome.blastp_for_homology.outfmt6 


###pfam database, updated
hmmsearch --cpu ${NSLOTS} --domtblout ${dirname}/ephemera/trinity_transcriptome_pfam_homology.domtblout ~/database/2023-Pfam-A.hmm ${dirname}/ephemera/stranded_trinity_ephemera_indoor.Trinity.clstr-98.fasta.transdecoder_dir/longest_orfs.pep
hmmsearch --cpu ${NSLOTS} --domtblout ${dirname}/lepidostoma/trinity_transcriptome_pfam_homology.domtblout ~/database/2023-Pfam-A.hmm ${dirname}/lepidostoma/stranded_trinity_lepidostoma_indoor.Trinity.fasta.clstr-98.transdecoder_dir/longest_orfs.pep
hmmsearch --cpu ${NSLOTS} --domtblout ${dirname}/gammarus/spades_transcriptome_pfam_homology.domtblout ~/database/2023-Pfam-A.hmm ${dirname}/gammarus/stranded_spades_gammarus_indoor.fasta.clstr-98.transdecoder_dir/longest_orfs.pep

##predict proteins
~/bin/TransDecoder-TransDecoder-v5.5.0/TransDecoder.Predict -t ${dirname}/lepidostoma/stranded_trinity_lepidostoma_indoor.Trinity.fasta.clstr-98 --retain_pfam_hits ${dirname}/lepidostoma/trinity_transcriptome_pfam_homology.domtblout --retain_blastp_hits ${dirname}/lepidostoma/trinity_transcriptome.blastp_for_homology.outfmt6
~/bin/TransDecoder-TransDecoder-v5.5.0/TransDecoder.Predict -t ${dirname}/ephemera/stranded_trinity_ephemera_indoor.Trinity.clstr-98.fasta --retain_pfam_hits ${dirname}/ephemera/trinity_transcriptome_pfam_homology.domtblout --retain_blastp_hits ${dirname}/ephemera/trinity_transcriptome.blastp_for_homology.outfmt6
~/bin/TransDecoder-TransDecoder-v5.5.0/TransDecoder.Predict -t ${dirname}/gammarus/stranded_spades_gammarus_indoor.fasta.clstr-98 --retain_pfam_hits ${dirname}/gammarus/spades_transcriptome_pfam_homology.domtblout --retain_blastp_hits ${dirname}/gammarus/spades_transcriptome.blastp_for_homology.outfmt6


##finally perform blastp with predicted proteins 
blastp -query ${dirname}/lepidostoma/stranded_trinity_lepidostoma_indoor.Trinity.fasta.clstr-98.transdecoder.pep -db ~/database/2023-uniprot_sprot.fasta -num_threads ${NSLOTS} -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > ${dirname}/lepidostoma/trinity_transcriptome_final_annotation.outfmt6
blastp -query ${dirname}/ephemera/stranded_trinity_ephemera_indoor.Trinity.clstr-98.fasta.transdecoder.pep -db ~/database/2023-uniprot_sprot.fasta -num_threads ${NSLOTS} -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > ${dirname}/ephemera/trinity_transcriptome_final_annotation.outfmt6
blastp -query ${dirname}/gammarus/stranded_spades_gammarus_indoor.fasta.clstr-98.transdecoder.pep -db ~/database/2023-uniprot_sprot.fasta -num_threads ${NSLOTS} -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > ${dirname}/gammarus/spades_transcriptome_final_annotation.outfmt6

