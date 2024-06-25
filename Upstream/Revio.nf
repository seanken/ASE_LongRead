#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.Primers="$projectDir/seqs/10x_3kit_primers.fasta"
params.CBC="$projectDir/seqs/3M-february-2018-REVERSE-COMPLEMENTED.txt.gz"
params.refFA="$projectDir/seqs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
//params.refBed
params.inputSkeraBam
params.outdir="output"
params.uniqBampy="$projectDir/RemoveDuplicates.py"

//Use Revio conda package 


workflow{
    lima=RunLima(params.inputSkeraBam,params.Primers)
    tag=RunTag(lima)
    refine=RunRefine(tag,params.Primers)
    correct=RunCorrect(refine,params.CBC)
    fa=MakeFasta(correct)
    bamfile=RunMinimap(fa,params.refFA) //,params.refBed)
    //sortBam=sortBamName(bamfile)
    //uniqueBam=uniqBam(bamfile,params.uniqBampy)
}

process RunLima
{
    input:
    path "movie.skera.bam"
    path "10x_3kit_primers.fasta"

    output:
    path "movie.lima.5p--3p.bam"

    '''
    lima --isoseq movie.skera.bam 10x_3kit_primers.fasta movie.lima.bam
    '''

}


process RunTag
{
    input:
    path "movie.lima.bam"

    output:
    path "movie.tagged.bam"

    '''
    isoseq tag --design T-12U-16B movie.lima.bam movie.tagged.bam
    '''
}

process RunRefine
{
    input:
    path "movie.tagged.bam"
    path "10x_3kit_primers.fasta"

    output:
    path "movie.refined.bam"

    '''
    isoseq refine --require-polya movie.tagged.bam 10x_3kit_primers.fasta movie.refined.bam
    '''
}

process RunCorrect
{
    input:
    path "movie.refined.bam"
    path "3M-february-2018-REVERSE-COMPLEMENTED.txt.gz"

    output:
    path "movie.corrected.bam"

    '''
    isoseq correct --barcodes 3M-february-2018-REVERSE-COMPLEMENTED.txt.gz movie.refined.bam movie.corrected.bam
    '''
}

//Should we add dedup?
//process RunDeDup
//{
//    input:
//    output:
//
//    '''
//    '''
//}

process MakeFasta
{
    input:
    path "movie.corrected.bam"

    output:
    path "movie.fasta"

    '''
    samtools fasta -T CB,XM movie.corrected.bam > movie.fasta
    '''
}


process RunMinimap
{
    publishDir "${params.outdir}/MappedBam", mode: 'move'

    input:
    path "movie.fasta"
    path "GRCm38.primary_assembly.genome.fa.gz"
    //path "gencode.vM10.annotation.bed"

    output:
    path "mm2.mapped.bam"

    //got rid of secondary=no and the splice bed

    '''
    minimap2 --secondary=no -ayY --MD --eqx -x splice:hq --cs=long GRCm38.primary_assembly.genome.fa.gz movie.fasta > mm2.mapped.sam
    samtools view -bS mm2.mapped.sam > mm2.mapped.bam
    '''
}

process sortBamName
{
    input:
    path "input.bam"

    output:
    path "input.sort.bam"

    '''
    samtools sort -n -o input.sort.bam input.bam
    '''
}

process uniqBam{
    publishDir "${params.outdir}/UniqueBam", mode: 'move'

    input:
    path "input.bam"
    path "RemoveDuplicates.py"

    output:
    path "out.bam"

    '''
    python RemoveDuplicates.py input.bam out.bam
    '''
}