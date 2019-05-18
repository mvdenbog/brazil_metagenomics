#!/bin/bash

#PBS -o /data3/projects/cisbi-0136/stdout
#PBS -e /data3/projects/cisbi-0136/stderr

## INFILE=`ls ./Dakar-2510_GGACTCCT-TATCCTCT_L001_R1_001.fastq`
## INFILE=$1
echo INFILE is ${INFILE} >> /data3/projects/cisbi-0136/stdout
## Q=$2
## Q="26"

cwd=${PBS_O_WORKDIR}
## echo cwd is ${cwd} >> /data3/projects/cisbi-0136/stdout
INBASE=`basename ${INFILE}`
INFILE=${cwd}/${INBASE}
INDIR=`dirname ${INFILE}`
## echo INFILE is ${INFILE} >> /data3/projects/cisbi-0136/stdout

BASEDIR=/data3/projects/cisbi-0136

## COMMENT: to start, make sure that the input file (fastq) does not contain spaces:
# use Streaming Editor (sed) to replace spaces with underscores (_):
sed -i "s/ /_/g" ${INFILE}

## COMMENT: 1. execute fastqc on my input file (fastq):
perl ${BASEDIR}/Tools/FastQC/fastqc ${INFILE}

## COMMENT: 2. execute prinseq below on input file (fastq):

perl ${BASEDIR}/Tools/prinseq-lite-0.20.4/prinseq-lite.pl -derep 1 -out_format 3 -out_good ${INDIR}/${INBASE}_prinseq_good -out_bad ${INDIR}/${INBASE}_prinseq_bad  -min_qual_score ${Q} -fastq ${INDIR}/${INBASE}

## COMMENT: 3. execute DeconSeq on the output of prinseq:

perl ${BASEDIR}/Tools/DeconSeq/deconseq-standalone-0.4.3/deconseq.pl  -f ${INDIR}/${INBASE}_prinseq_good.fastq  -dbs human -dbs_retain viral -c 95 -i 94 -group 1 -out_dir ${INDIR}

## rename deconseq result files:
DECONSEQ_CONT=`ls ${INDIR}/*_cont.fq`
DECONSEQ_CLEAN=`ls ${INDIR}/*_clean.fq`

mv -i ${DECONSEQ_CONT} ${INDIR}/${INBASE}_cont.fq
mv -i ${DECONSEQ_CLEAN} ${INDIR}/${INBASE}_clean.fq

### COMMENT: 4. execute Kraken2 classification on resulting deconseq output:

${BASEDIR}/Tools/Kraken2/kraken2/kraken2 --threads 10 -db ${BASEDIR}/DB/Kraken2/ --unclassified-out ${INDIR}/${INBASE}_unclassified-out --classified-out  ${INDIR}/${INBASE}_classified-out --output  ${INDIR}/${INBASE}_output --report-zero-counts --report  ${INDIR}/${INBASE}_report ${INDIR}/${INBASE}_clean.fq

### Associate complete taxonomy to each taxid:

export INFILE_TAXID=${INDIR}/${INBASE}_output
echo infile_taxid is ${INFILE_TAXID} >> ${INDIR}/log
/data3/projects/cisbi-0136/Tools/Kraken2Krona/taxid2taxo2krona.pbs

### Next step

# Identify all viral sequences from Kraken's classification

echo grepping virus in ${INFILE_TAXID}.taxo >> ${INDIR}/log
grep Virus ${INFILE_TAXID}.taxo | cut -f 2 > ${INFILE_TAXID}.taxo.virus_list

perl /data3/projects/cisbi-0136/Tools/get_fastq_from_name_seq.pl ${INFILE_TAXID}.taxo.virus_list ${INFILE}_clean.fq ${INFILE_TAXID}.taxo.virus_list.fastq

# Assemble all viral sequences using SPAdes

python /data3/projects/cisbi-0136/Tools/Spades/SPAdes-3.13.0-Linux/bin/spades.py -t 4 -k 21,31,51,71 -s ${INFILE_TAXID}.taxo.virus_list.fastq --only-assembler -o ${INDIR}/spades_out

# BlastN the resulting contigs:

/data1/biotools/ncbi-blast-2.4.0+/bin/blastn -query ${INDIR}/spades_out/scaffolds.fasta -db /data1/public_data/nt/nt -outfmt "6 qseqid sseqid length pident mismatch gapopen qstart qend sstart send evalue bitscore" -evalue 1e-3 -max_target_seqs 1 -task blastn -out ${INDIR}/spades_out/scaffolds.blastn

# BlastN all unclassified sequences from Kraken, and associate with taxonomie (in order not to loose any possible viral sequences in the unclassifies sequences)

perl /data3/projects/cisbi-0136/fastq2fasta.pl ${INDIR}/${INBASE}_unclassified-out ${INDIR}/${INBASE}_unclassified-out.fasta 

# before blasting, reduce the number of sequences, using homology with the cd-hit tool:

/data1/biotools/cdhit/cd-hit-est -c 0.99 -d 0 -T 4 -M 0 -i ${INDIR}/${INBASE}_unclassified-out.fasta  -o ${INDIR}/${INBASE}_unclassified-out.cdhit.fasta

/data1/biotools/ncbi-blast-2.4.0+/bin/blastn -query ${INDIR}/${INBASE}_unclassified-out.cdhit.fasta -db /data1/public_data/nt/nt -outfmt "6 qseqid sseqid length pident mismatch gapopen qstart qend sstart send evalue bitscore" -evalue 1e-3 -max_target_seqs 1 -task blastn -out ${INDIR}/${INBASE}_unclassified-out.cdhit.blastn


