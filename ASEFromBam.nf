params.bam //input bam
params.gtf="$projectDir/ref/genes.gtf" //input gtf
params.outdir="output" //output directory
params.pyScript="$projectDir/scripts/GeneAnnotate.py" //python script for annotating 
params.matpy="$projectDir/scripts/MakeMat.py" //makes umi_tools output into MM matrix
params.AlleleCountJar="$projectDir/scripts/AlleleMASSeq.jar"
params.input_vcf //vcf used
params.vcf_col // sample name in the vcf to use
params.umi_tag="Jr" //tag in the bam used for UMI
params.cbc_tag="Jp" //tag in the bam used for CBC
params.method="Long"
params.strand=2

//prepare VCF for downstream
process PrepVCF
{
    publishDir "${params.outdir}/new_vcf", mode: 'rellink'

    input:
    env vcf_col from params.vcf_col
    path input_vcf, stageAs:"input.vcf.gz" from params.input_vcf

    output:
    path "new.vcf" into new_vcf_ch




    '''
    tabix -p vcf input.vcf.gz
    bcftools view -H -O v -s $vcf_col input.vcf.gz | grep -v "0|0" | grep -v "1|1" > new.vcf
    '''


}


//Uses feature counts to annotate bam in CORE format (since long read)
process FeatureCount_intron
{
    publishDir "${params.outdir}/FC_intron", mode: 'rellink'

    input:
    path gtf, stageAs:"genes.gtf" from params.gtf
    path bam, stageAs:"Aligned.sortedByCoord.out.bam" from params.bam
    env strand from params.strand

    output:
    path "gene_assigned" into gene_intron_assign
    path "gene_assigned.summary" into gene_inton_assign_summary
    path "Aligned.sortedByCoord.out.bam.featureCounts" into read_gene_ann

    '''
    featureCounts -a genes.gtf -L -t gene -g gene_name -o gene_assigned -s $strand -R CORE Aligned.sortedByCoord.out.bam
    '''
}


process FeatureCount_exon
{
    publishDir "${params.outdir}/FC_exon", mode: 'rellink'

    input:
    path gtf, stageAs:"genes.gtf" from params.gtf
    path bam, stageAs:"Aligned.sortedByCoord.out.bam" from params.bam
    env strand from params.strand

    output:
    path "gene_assigned" into gene_exon_assign
    path "gene_assigned.summary" into gene_exon_assign_summary

    '''
    featureCounts -a genes.gtf -L -t exon -g gene_name -o gene_assigned -s $strand Aligned.sortedByCoord.out.bam
    '''
}


//Combines CORE format annotation and bam file, adding gene info as tag
process AnnotBam
{
    publishDir "${params.outdir}/Bam", mode: 'rellink'
    
    input:
    path pyScript, stageAs:"GeneAnnotate.py" from params.pyScript
    path bam, stageAs:"Aligned.sortedByCoord.out.bam" from params.bam
    path out, stageAs:"Aligned.sortedByCoord.out.bam.featureCounts" from read_gene_ann
    env method from params.method

    output:
    path "ann.bam" into ann_bam

    '''
    python GeneAnnotate.py Aligned.sortedByCoord.out.bam Aligned.sortedByCoord.out.bam.featureCounts ann.bam $method
    '''
}


//Uses UMITools to count the number of UMI per gene per cell
process UMITools
{
    publishDir "${params.outdir}/UMITools", mode: 'rellink'
    input:
    path ann_bam, stageAs:"ann.bam" from ann_bam
    env umi_tag from params.umi_tag
    env cbc_tag from params.cbc_tag

    output:
    path "counts.txt" into counts_umi

    '''
    samtools index ann.bam
    umi_tools count --umi-tag=$umi_tag --cell-tag=$cbc_tag --per-gene --gene-tag=XT --skip-tags-regex=NA --per-cell -I ann.bam -S counts.txt --extract-umi-method=tag
    '''
}

//makes the output of UMI tools a MM format matrix 
process MakeMM
{
    publishDir "${params.outdir}/Matrix", mode: 'rellink'

    input:
    path counts, stageAs:"counts.txt" from counts_umi
    path matpy, stageAs:"MakeMat.py" from params.matpy

    output:
    path "matrix.mtx.gz" into mat
    path "features.tsv.gz" into feat
    path "barcodes.tsv.gz" into cell

    '''
    python MakeMat.py counts.txt
    '''

}

//This process runs a java based method to count number of alleles, giving results as gene level
process CountAlleles
{
    publishDir "${params.outdir}/AlleleCounts", mode: 'rellink'

    input:
    path ann_bam, stageAs:"ann.bam" from ann_bam
    path vcf, stageAs:"new.vcf" from new_vcf_ch 
    path AlleleCountJar, stageAs:"AlleleMASSeq.jar" from params.AlleleCountJar
    env umi_tag from params.umi_tag
    env cbc_tag from params.cbc_tag

    output:
    path "counts.txt" into gene_counts_ch

    '''
    echo "Count Alleles UMIs"
    java -jar AlleleMASSeq.jar ann.bam new.vcf counts.txt $cbc_tag $umi_tag
    '''

}


