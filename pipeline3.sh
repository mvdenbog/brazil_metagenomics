#!/bin/bash

INFILE=$1
Q=$2
## echo INFILE is ${INFILE} >> /data3/projects/cisbi-0136/stdout

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

job_prinseq=$(qsub -W depend=afterok:$job_fastqc -v INFILE=${INFILE},INDIR=${INDIR},INBASE=${INBASE},BASEDIR=${BASEDIR},Q=${Q} -l nodes=1:ppn=1,mem=10gb,walltime=10:0:0 ${SCRIPTDIR}/prinseq.pbs); 

## COMMENT: 3. wait for prinseq to finish, and execute DeconSeq on the output of prinseq:

### NEW PARALLELIZED VERSION:

job_deconseq_prepare=$(qsub -W depend=afterok:$job_prinseq -v INFILE=${INFILE},Q=${Q},NB_FILE=200 ${SCRIPTDIR}/pgpDeconseq_prepare.pbs); 

job_deconseq_submission=$(qsub -W depend=afterok:$job_deconseq_prepare -v INFILE=${INFILE},INDIR=${INDIR},INBASE=${INBASE},BASEDIR=${BASEDIR},Q=${Q},NB_FILE=200,DB_CONT="human",DB_RETAIN="viral",COVERAGE=95,IDENTITY=94,GROUP=1 ${SCRIPTDIR}/pgpDeconseq_submit.pbs); 

job_deconseq_cleanup=$(qsub -W depend=afterokarray:$job_deconseq_submission -v INFILE=${INFILE},INDIR=${INDIR},INBASE=${INBASE},BASEDIR=${BASEDIR},Q=${Q} ${SCRIPTDIR}/pgpDeconseq_cleanup.pbs); 



#### sh ${BASEDIR}/Tools/DeconSeq/pgpDeconseq_sge.sh -f ${INFILE}_prinseq_good.fastq  -d human -r viral -c 95 -i 94 -g 1 -n 200


##### OLD SEQUENTIAL VERSION: job_deconseq=$(qsub -W depend=afterok:$job_prinseq -v INFILE=${INFILE},INDIR=${INDIR},INBASE=${INBASE},BASEDIR=${BASEDIR} ${SCRIPTDIR}/deconseq.pbs); 



### COMMENT: 4. wait for deconseq to finish, and execute Kraken2 classification on resulting deconseq output:

KRAKEN_IN=`basename ${INBASE} .fastq`_clean.fastq
KRAKEN_OUT=${INDIR}/`basename ${INBASE} .fastq`_clean.fastq_output
job_kraken=$(qsub -W depend=afterok:$job_deconseq_cleanup -v INDIR=${INDIR},INBASE=${KRAKEN_IN},BASEDIR=${BASEDIR} -l nodes=1:ppn=1,mem=35gb,walltime=10:0:0 ${SCRIPTDIR}/kraken.pbs); 


### COMMENT: 5. wait for kraken to finish, and do 2 things:
   #### 1. Associate complete taxonomy to each taxid:
   #### 2. output krona html file


job_assoc_taxo_taxid=$(qsub -W depend=afterok:$job_kraken -v INFILE=${KRAKEN_OUT} -l walltime=10:0:0 /data3/projects/cisbi-0136/Tools/Kraken2Krona/taxid2taxo2krona.pbs); 


# COMMENT 6. Wait for taxo task to finish, and identify all viral sequences from Kraken's classification

job_id_viral_from_kraken=$(qsub -W depend=afterok:$job_assoc_taxo_taxid -v INFILE=${KRAKEN_OUT}.taxo ${SCRIPTDIR}/id_viral_from_kraken.pbs);

# COMMENT: 7. wait for the previous job to finish, and Assemble all viral sequences using SPAdes

job_spades=$(qsub -W depend=afterok:$job_id_viral_from_kraken -l nodes=1:ppn=8,mem=30gb,walltime=15:0:0  -v INFILE=${KRAKEN_OUT}.taxo.virus_list.fastq,INDIR=${INDIR} ${SCRIPTDIR}/spades.pbs);


# COMMENT 8. wait for assembly to finish, and BlastN the resulting contigs/scaffolds:

job_blastn_contigs=$(qsub -N blastn_contigs -W depend=afterok:$job_spades -v INDIR=${INDIR} -l nodes=1:ppn=1,mem=30gb,walltime=23:59:0 ${SCRIPTDIR}/blastn_contigs.pbs);

# COMMENT 9: wait for blastn to finish, and BlastN all unclassified sequences from Kraken, and associate with taxonomy (in order not to loose any possible viral sequences in the unclassifies sequences). We have to convert all the fastq sequences to fasta first:

## 9.1 convert fastq to fasta:

IN=${INDIR}/`basename ${INBASE} .fastq`_clean.fastq_unclassified-out
job_fastq2fasta=$(qsub -W depend=afterok:$job_blastn_contigs -v INFILE=${IN} ${SCRIPTDIR}/fastq2fasta.pbs);

## 9.2 before blasting, reduce the number of sequences, using homology with the cd-hit tool:

job_cdhit=$(qsub -W depend=afterok:$job_fastq2fasta -v INFILE=${IN}.fasta ${SCRIPTDIR}/cdhit.pbs); 

## 9.3 perform blastN on the resulting cdhit collapsed sequences:
### job_blastn_unclass_prepare=$(qsub -W depend=afterok:$job_cdhit -v INFILE=${IN}.cdhit.fasta,INDIR=${INDIR} ${SCRIPTDIR}/blastn_unclass_prepare.pbs);

job_blastn_unclass_prepare=$(qsub -W depend=afterok:$job_cdhit -v INFILE=${IN}.cdhit.fasta,INDIR=${INDIR} ${SCRIPTDIR}/blastn_unclass_prepare.pbs);

### job_blastn_unclass_submit=$(echo /data3/projects/cisbi-0136/Tools/Pblastall/PgpBlastall2/sge_blastall ${IN}.cdhit.blastn pgp_blastall_prepare.list -db /data3/projects/cisbi-0136/DB/Blast/Gbk_vrl/genbank_release_vrl  -evalue 1e-3 -max_target_seqs 1 -task blastn   | qsub -W depend=afterok:$job_blastn_unclass_prepare -l nodes=1:ppn=1,mem=30gb -t 1-1000 )

job_blastn_unclass_submit=$(echo /data3/projects/cisbi-0136/Tools/Pblastall/PgpBlastall2/sge_blastall ${IN}.cdhit.blastn ${INDIR}/pgp_blastall_prepare.list -db /data1/public_data/nt/nt -evalue 1e-3 -max_target_seqs 1 -task blastn   | qsub -N blast_unclassd -W depend=afterok:$job_blastn_unclass_prepare -l nodes=1:ppn=1,mem=30gb,walltime=23:59:0 -t 1-1000);

job_assoc_blastn_taxo=$(qsub -W depend=afterokarray:$job_blastn_unclass_submit -v INFILE=${IN}.cdhit.blastn /data3/projects/cisbi-0136/Pipeline/blastn2taxo2krona.pbs); 

