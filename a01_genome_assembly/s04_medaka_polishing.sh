#!/bin/sh
### Account information
#PBS -W group_list=cu_10081 -A cu_10081
### Number of nodes#PBS -l nodes=1:ppn=4
### Memory
#PBS -l mem=180gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds>
#PBS -l walltime=00:04:00:00

module load tools anaconda3/4.4.0 medaka/1.2.0

cd /home/people/knunie/cu_10081/projects/m_acridum_pg2021/generateddata/ont_genome_asm

for id in 01 02 03 04 05 06; do

BASECALLS=${id}_Reads/${id}_rawReads_no_mt_v7.m4kb95p80x.fq.gz
DRAFT=${id}_Flye_95p_10kbOvl/assembly.fasta

## First Medaka round (uncomment if needed)
medaka_consensus -i ${BASECALLS} -d ${DRAFT} -o ${id}_FlyeMedaka_r1 -t 4 -m r941_min_high_g360

## Second Medaka round
medaka_consensus \
    -i ${BASECALLS} \
    -d ${id}_FlyeMedaka_r1/*fasta \
    -o ${id}_FlyeMedaka_r2 \
    -t 4 \
    -m r941_min_high_g360

rm ${id}_FlyeMedaka_r1/calls_to_draft.bam*
rm ${id}_FlyeMedaka_r1/consensus.fasta.mmi
rm ${id}_FlyeMedaka_r2/calls_to_draft.bam*

qzip ${id}_FlyeMedaka_r2/consensus.fasta

done
