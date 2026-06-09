#!/bin/sh
### Account information
#PBS -W group_list=cu_10081 -A cu_10081
### Number of nodes#PBS -l nodes=1:ppn=40
### Memory
#PBS -l mem=180gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds>
#PBS -l walltime=00:48:00:00

module load tools samtools/1.12 pigz/2.3.4 seqkit/0.13.2 seqtk/1.3

cd /home/people/knunie/cu_10081/projects/m_acridum_pg2021/generateddata/ont_genome_asm
for id in 01 02 03 04 05 06; do

cd ${id}_Reads

### List mitochondrial reads
#samtools view -F 0x904 -q 40 ${id}_sort.bam CP058939.1 ${id} | cut -f 1 | sort | uniq | sed 's/^/^@/g' > ${id}_mt_Rname.lst

### Raw reads
reads=~/cu_10081/projects/m_acridum_pg2021/data/barcode${id}/barcode${id}.merged.fq.gz

### Exclude mitochondrial reads from the original raw read file. NOTE: grep does not work with lines > 8 kb
#cat ${reads} | paste - - - - | sed -f <(sed 's|.*|/&/d|' ${id}_mt_Rname.lst) | tr "\t" "\n" | gzip - > ${id}_rawReads_no_mt_v7.fq.gz

cat ${id}_mt_Rname.lst | sed 's|.*|/&/d|' > {id}_mt_Rname.sed_delete
zcat ${reads} | paste - - - - | sed -f {id}_mt_Rname.sed_delete | tr "\t" "\n" | gzip - > ${id}_rawReads_no_mt_v7.fq.gz

zcat ${id}_rawReads_no_mt_v7.fq.gz | seqkit seq -m 8000 | gzip - > ${id}_rawReads_NoMt_m8kb.fq.gz
zcat ${id}_rawReads_no_mt_v7.fq.gz | seqkit seq -m 10000 | gzip - > ${id}_rawReads_NoMt_m10kb.fq.gz
zcat ${id}_rawReads_no_mt_v7.fq.gz | seqkit seq -m 12000 | gzip - > ${id}_rawReads_NoMt_m12kb.fq.gz

cd ..

done
