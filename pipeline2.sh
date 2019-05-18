#!/bin/bash

INFILE=$1
Q=$2
echo INFILE is ${INFILE} >> /data3/projects/cisbi-0136/stdout

INBASE=$(basename ${INFILE})
INFILE=$(readlink -f ${INFILE})
INDIR=$(dirname ${INFILE})

SCRIPTDIR=/data3/projects/cisbi-0136/Pipeline

BASEDIR=/data3/projects/cisbi-0136

## COMMENT: to start, make sure that the input file (fastq) does not contain spaces:
# use Streaming Editor (sed) to replace spaces with underscores (_):

job_sed=$(qsub -v INFILE=${INFILE} ${SCRIPTDIR}/sed.pbs); 

## COMMENT: 1. wait for sed to finish, and execute fastqc on my input file (fastq):
job_fastqc=$(qsub -W depend=afterok:$job_sed -v INFILE=${INFILE},BASEDIR=${BASEDIR} ${SCRIPTDIR}/fastqc.pbs); 

## COMMENT: 2. wait for fastqc to finish, and execute prinseq below on input file (fastq):

job_prinseq=$(qsub -W depend=afterok:$job_fastqc -v INFILE=${INFILE},INDIR=${INDIR},INBASE=${INBASE},BASEDIR=${BASEDIR},Q=${Q} ${SCRIPTDIR}/prinseq.pbs); 

## COMMENT: 3. wait for prinseq to finish, and execute DeconSeq on the output of prinseq:

job_deconseq=$(qsub -W depend=afterok:$job_prinseq -v INFILE=${INFILE},INDIR=${INDIR},INBASE=${INBASE},BASEDIR=${BASEDIR} ${SCRIPTDIR}/deconseq.pbs); 



### COMMENT: 4. wait for deconseq to finish, and execute Kraken2 classification on resulting deconseq output:

job_kraken=$(qsub -W depend=afterok:$job_deconseq -v INDIR=${INDIR},INBASE=${INBASE}_clean.fq,BASEDIR=${BASEDIR} -l nodes=1:ppn=1,mem=35gb ${SCRIPTDIR}/kraken.pbs); 


### COMMENT: 5. wait for kraken to finish, and do 2 things:
   #### 1. Associate complete taxonomy to each taxid:
   #### 2. output krona html file

job_assoc_taxo_taxid=$(qsub -W depend=afterok:$job_kraken -v INFILE=${INFILE}_clean.fq_output /data3/projects/cisbi-0136/Tools/Kraken2Krona/taxid2taxo2krona.pbs); 


# COMMENT 6. Wait for taxo task to finish, and identify all viral sequences from Kraken's classification

job_id_viral_from_kraken=$(qsub -W depend=afterok:$job_assoc_taxo_taxid -v INFILE=${INFILE}_clean.fq_output.taxo ${SCRIPTDIR}/id_viral_from_kraken.pbs);

# COMMENT: 7. wait for the previous job to finish, and Assemble all viral sequences using SPAdes

job_spades=$(qsub -W depend=afterok:$job_id_viral_from_kraken -l nodes=1:ppn=4,mem=10gb -v INFILE=${INFILE}_clean.fq_output.taxo.virus_list.fastq,INDIR=${INDIR} ${SCRIPTDIR}/spades.pbs);


# COMMENT 8. wait for assembly to finish, and BlastN the resulting contigs/scaffolds:

job_blastn_contigs=$(qsub -W depend=afterok:$job_spades -v INDIR=${INDIR} ${SCRIPTDIR}/blastn_contigs.pbs);

# COMMENT 9: wait for blastn to finish, and BlastN all unclassified sequences from Kraken, and associate with taxonomy (in order not to loose any possible viral sequences in the unclassifies sequences). We have to convert all the fastq sequences to fasta first:

## 9.1 convert fastq to fasta:

job_fastq2fasta=$(qsub -W depend=afterok:$job_blastn_contigs -v INFILE=${INFILE}_clean.fq_unclassified-out ${SCRIPTDIR}/fastq2fasta.pbs);

## 9.2 before blasting, reduce the number of sequences, using homology with the cd-hit tool:

job_cdhit=$(qsub -W depend=afterok:$job_fastq2fasta -v INFILE=${INFILE}_clean.fq_unclassified-out.fasta ${SCRIPTDIR}/cdhit.pbs); 

## 9.3 perform blastN on the resulting cdhit collapsed sequences:
job_blastn_unclass=$(qsub -W depend=afterok:$job_cdhit -v INFILE=${INFILE}_clean.fq_unclassified-out.cdhit.fasta,INDIR=${INDIR} ${SCRIPTDIR}/blastn_unclass.pbs);

## 9.4 associate species-ID and full taxonomy to each Accession number we got from blastn, and also make a krona representation for the results:

job_assoc_blastn_taxo=$(qsub -W depend=afterok:$job_blastn_unclass -v INFILE=${INFILE}_clean.fq_unclassified-out.cdhit.blastn /data3/projects/cisbi-0136/Pipeline/blastn2taxo2krona.pbs);


