#! /bin/bash

#$ -cwd
#$ -e ErrFiles/star.$TASK_ID.err
#$ -o ErrFiles/star.$TASK_ID.log
#$ -q broad
#$ -P regevlab
#$ -l h_vmem=60g
#$ -l h_rt=60:00:00
#$ -l os="RedHat7"

#$ -t 1-4

source /broad/software/scripts/useuse

source ~/ForPyth.sh
conda activate MASPip2

#mkdir m64297e_220708_013849
#cd m64297e_220708_013849

use .java-jdk-1.8.0_181-x86-64
use .bcftools-1.15.1

SEEDFILE=bams.txt
bam=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $1}')
#vcf_col=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $2}')
vcf_col=BN0009_P86_A1_BN0009_P86_A1
vcf=/stanley/levin_asap_storage/612-eqtl/GenotypeData/Clean_vcf/Combined/comb_new.no.chr.vcf.gz ##wrong vcf
outdir=samp_$SGE_TASK_ID
mkdir $outdir
cd $outdir

#bam=/stanley/levin_asap_storage/MASseq/m64297e_220708_013849/m64297e_220708_013849.bam

nextflow=/stanley/levin_asap/ssimmons/eQTL/Code/AlleleExpressionPipeline/version2.0/ASE_pipeline/nextflow
pipeline=/stanley/levin_asap/ssimmons/eQTL/Code/AlleleExpressionPipeline/MASSeq/FromBam/ASEFromBam.nf
$nextflow $pipeline --bam $bam --input_vcf $vcf --vcf_col $vcf_col

