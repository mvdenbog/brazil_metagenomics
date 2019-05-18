#!/bin/bash

INFILE=$1
Q=$2
## echo INFILE is ${INFILE} >> /data3/projects/cisbi-0136/stdout

INBASE=$(basename ${INFILE})
INFILE=$(readlink -f ${INFILE})
INDIR=$(dirname ${INFILE})

SCRIPTDIR=/data3/projects/cisbi-0136/Pipeline

BASEDIR=/data3/projects/cisbi-0136

# BLAST_DB=/data3/projects/cisbi-0136/DB/Blast/Nt/nt
BLAST_DB=/data3/projects/cisbi-0136/DB/Blast/Gbk_vrl/genbank_release_vrl

jobid=$(echo "echo START" | qsub)

if [[ `head -n 1 ${INFILE}` =~ " " ]]; then
echo "submitting sed"
jobid=$(qsub -W depend=afterok:$jobid -v INFILE=${INFILE} ${SCRIPTDIR}/sed.pbs); 
fi;

# IN is  ${INFILE}
fastqc_out=${INDIR}/`basename ${INBASE} .fastq`_fastqc.zip
if [[ ! -e ${fastqc_out} ]]; then
echo "submitting fastqc"
jobid=$(qsub -W depend=afterok:$jobid -v INFILE=${INFILE},BASEDIR=${BASEDIR} ${SCRIPTDIR}/fastqc.pbs); 
fi

# IN is ${INDIR}/${INBASE}
prinseq_out=${INDIR}/${INBASE}_prinseq_good.fastq
if [[ ! -e ${prinseq_out} ]]; then
echo "submitting prinseq"
jobid=$(qsub -W depend=afterok:$jobid -v INFILE=${INFILE},INDIR=${INDIR},INBASE=${INBASE},BASEDIR=${BASEDIR},Q=${Q} -l nodes=1:ppn=1,mem=10gb,walltime=10:0:0 ${SCRIPTDIR}/prinseq.pbs); 
fi;


deconseq_out=${INDIR}/`basename ${INBASE} .fastq`_clean.fastq
if [[ ! -e ${deconseq_out} ]]; then
echo "submitting deconseq"
jobid=$(qsub -W depend=afterok:$jobid -v INFILE=${INFILE},Q=${Q},NB_FILE=200 ${SCRIPTDIR}/pgpDeconseq_prepare.pbs); 
jobid=$(qsub -W depend=afterok:$jobid -v INFILE=${INFILE},INDIR=${INDIR},INBASE=${INBASE},BASEDIR=${BASEDIR},Q=${Q},NB_FILE=200,DB_CONT="human",DB_RETAIN="viral",COVERAGE=95,IDENTITY=94,GROUP=1 ${SCRIPTDIR}/pgpDeconseq_submit.pbs); 
jobid=$(qsub -W depend=afterokarray:$jobid -v INFILE=${INFILE},INDIR=${INDIR},INBASE=${INBASE},BASEDIR=${BASEDIR},Q=${Q} ${SCRIPTDIR}/pgpDeconseq_cleanup.pbs); 
fi

KRAKEN_IN=`basename ${INBASE} .fastq`_clean.fastq
KRAKEN_OUT=${INDIR}/`basename ${INBASE} .fastq`_clean.fastq_output
if [[ ! -e ${KRAKEN_OUT} ]]; then
echo "submitting kraken"
jobid=$(qsub -W depend=afterok:$jobid -v INDIR=${INDIR},INBASE=${KRAKEN_IN},BASEDIR=${BASEDIR} -l nodes=1:ppn=1,mem=35gb,walltime=10:0:0 ${SCRIPTDIR}/kraken.pbs); 
fi;

taxo_out=${KRAKEN_OUT}.taxo
if [[ ! -e ${taxo_out} ]]; then
echo "submitting assoc "
jobid=$(qsub -W depend=afterok:$jobid -v INFILE=${KRAKEN_OUT} -l walltime=10:0:0 /data3/projects/cisbi-0136/Tools/Kraken2Krona/taxid2taxo2krona.pbs); 
fi;


kraken_viral_out=${taxo_out}.virus_list.fastq
if [[ ! -e ${kraken_viral_out} ]]; then
echo "submitting viral id from kraken"
jobid=$(qsub -W depend=afterok:$jobid -v INFILE=${KRAKEN_OUT}.taxo ${SCRIPTDIR}/id_viral_from_kraken.pbs);
fi

if [[ ! -d spades_out ]]; then
echo "submitting spades"
jobid=$(qsub -W depend=afterok:$jobid -l nodes=1:ppn=8,mem=30gb,walltime=15:0:0  -v INFILE=${KRAKEN_OUT}.taxo.virus_list.fastq,INDIR=${INDIR} ${SCRIPTDIR}/spades.pbs);
fi;


#OUT is spades_out/scaffolds.blastn
if [[ ! -e spades_out/scaffolds.blastn ]]; then
echo "submitting blastn contigs"
jobid=$(qsub -N blastn_contigs -W depend=afterok:$jobid -v INDIR=${INDIR} -l nodes=1:ppn=1,mem=30gb,walltime=23:59:0  ${SCRIPTDIR}/blastn_contigs.pbs);
fi;

IN=${INDIR}/`basename ${INBASE} .fastq`_clean.fastq_unclassified-out
fasta_out=${IN}.fasta
if [[ ! -e ${IN}.fasta ]]; then
echo "submitting fastq2fasta "
jobid=$(qsub -W depend=afterok:$jobid -v INFILE=${IN} ${SCRIPTDIR}/fastq2fasta.pbs);
fi;

#OUT is ${IN}.cdhit.fasta
if [[ ! -e ${IN}.cdhit.fasta ]]; then
echo "submitting cdhit"
jobid=$(qsub -W depend=afterok:$jobid -v INFILE=${IN}.fasta ${SCRIPTDIR}/cdhit.pbs); 
fi;

#OUT is  ${IN}.cdhit.blastn
if [[ ! -e ${IN}.cdhit.blastn ]]; then
echo "submitting blastn unclass"
jobid=$(qsub -W depend=afterok:$jobid -v INFILE=${IN}.cdhit.fasta,INDIR=${INDIR} ${SCRIPTDIR}/blastn_unclass_prepare.pbs);
jobid=$(echo /data3/projects/cisbi-0136/Tools/Pblastall/PgpBlastall2/sge_blastall ${IN}.cdhit.blastn ${INDIR}/pgp_blastall_prepare.list -db ${BLAST_DB} -evalue 1e-3 -max_target_seqs 1 -task blastn   | qsub -N blast_unclassd -W depend=afterok:$jobid -l nodes=1:ppn=10,mem=30gb,walltime=23:59:0 -t 1-300);
else
echo "submitting blastn unclass bogus array"
jobid=$(echo sleep 1 | qsub -W depend=afterok:$jobid -t 1)
fi

#OUT is ${IN}.cdhit.blastn.taxo
if [[ ! -e ${IN}.cdhit.blastn.taxo ]]; then
echo "submitting blastn taxo"
jobid=$(qsub -W depend=afterokarray:${jobid} -v INFILE=${IN}.cdhit.blastn -l walltime=23:59:0 /data3/projects/cisbi-0136/Pipeline/blastn2taxo2krona.pbs); 
fi;
