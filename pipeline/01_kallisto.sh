#!/bin/bash
#SBATCH --nodes 1 --ntasks 8 --mem 4G --time 2:00:00 --out logs/kallisto.%a.log

module load kallisto
CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi

N=${SLURM_ARRAY_TASK_ID}
INDIR=deve_reads
OUTDIR=results/kallisto
SAMPLEFILE=samples.csv
TRANSCRIPTS=cds/Phycomycesb_Phyb12.kallisto.idx 

mkdir -p $OUTDIR

if [ ! $N ]; then
 N=$1
fi

if [ ! $N ]; then
 echo "cannot run without a number provided either cmdline or --array in sbatch"
 exit
fi
IFS=,
tail -n +2 $SAMPLEFILE | sed -n ${N}p | while read SAMPLE GENOTYPE CONDITION REP

do
 OUTNAME=$GENOTYPE.$CONDITION.r${REP}
 if [ ! -f $OUTDIR/$OUTNAME/abundance.h5 ]; then
  kallisto quant -i $TRANSCRIPTS -o $OUTDIR/$OUTNAME -t $CPU --bias $INDIR/${SAMPLE}_1.fq.gz $INDIR/${SAMPLE}_2.fq.gz 
 fi
done
