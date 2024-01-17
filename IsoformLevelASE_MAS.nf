params.ref
params.fa
params.gtf
params.bam

params.outdir="output" //output directory
params.pyScript="$projectDir/scripts/GeneAnnotate_isoform.py" //python script for annotating by isoform
params.matpy="$projectDir/scripts/MakeMat.py" //makes umi_tools output into MM matrix
params.AlleleCountJar="$projectDir/scripts/AlleleMASSeq.jar"
params.vcf //vcf used
params.vcf_col // sample name in the vcf to use
params.umi_tag="Jr" //tag in the bam used for UMI
params.cbc_tag="Jp" //tag in the bam used for CBC
params.method="MAS"
params.picardjar="/seq/software/picard-public/current/picard.jar"



if(!params.ref && !(params.fa && params.gtf))
{
    println("Need either a pregenerated reference or a gtf and fasta")
    System.exit(0)
}



workflow{
    //If it doesn't exist make minimap2 reference
    if(!params.ref)
    {
        GenerateTransFA(params.fa,params.gtf)
        GenerateMinimap2Ref(GenerateTransFA.out)
    }

    minimap2Ref=params.ref ? params.ref : GenerateMinimap2Ref.out

    if(!params.bam)
    {
        println("No bam, only generate reference")
        System.exit(0)
    }
    if(!params.vcf)
    {
        println("No VCF, will not run ASE step")
    }

    if(!params.vcf)
    {
        println("No VCF column, will not run ASE step")
    }

    //make fastq and map to transcriptome
    BamToFastq(params.bam,params.picardjar)
    MapBamToTrans(BamToFastq.out,minimap2Ref)

    //Take results and annotated genome mapped bam file with it
    AnnotBam(params.pyScript,params.bam,MapBamToTrans.out,params.method)
    
    //Get isoform level expression for each cell
    UMITools(AnnotBam.out,params.umi_tag,params.cbc_tag)
    MakeMM(UMITools.out,params.matpy)
    
    //If you have a vcf cleans up vcf and runs ASE script
    if(params.vcf && params.vcf_col)
    {
        CleanVCF(params.vcf,params.vcf_col)
        CountAlleles(AnnotBam.out,CleanVCF.out,params.AlleleCountJar,params.umi_tag,params.cbc_tag)
    }
    


}


////////////
//Section for ASE related processes
//////////
//prepare VCF for downstream
process CleanVCF
{
    publishDir "${params.outdir}/new_vcf", mode: 'rellink'

    input:
    path "input.vcf.gz"
    env vcf_col

    output:
    path "new.vcf"


    '''
    tabix -p vcf input.vcf.gz
    bcftools view -H -O v -s $vcf_col input.vcf.gz | grep -v "0|0" | grep -v "1|1" > new.vcf
    '''


}



//This process runs a java based method to count number of alleles, giving results as gene level
process CountAlleles
{
    publishDir "${params.outdir}/AlleleCounts", mode: 'rellink'

    input:
    path "ann.bam" 
    path "new.vcf"
    path "AlleleMASSeq.jar"
    env umi_tag 
    env cbc_tag 

    output:
    path "counts.txt" 

    '''
    echo "Count Alleles UMIs"
    java -jar AlleleMASSeq.jar ann.bam new.vcf counts.txt $cbc_tag $umi_tag
    '''

}



/////
//Section for counting isoform level expression
////
//
process UMITools
{
    publishDir "${params.outdir}/UMITools", mode: 'rellink'
    input:
    path "ann.bam"
    env umi_tag 
    env cbc_tag 

    output:
    path "counts.txt"

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
    path "counts.txt"
    path "MakeMat.py"

    output:
    path "matrix.mtx.gz"
    path "features.tsv.gz"
    path "barcodes.tsv.gz"

    '''
    python MakeMat.py counts.txt
    '''

}


//////////////
//Steps for combining trans mapped bam and genome mapped bam
//////////////


//Combines transcriptome annotation and bam file, adding isform info as tag
process AnnotBam
{
    publishDir "${params.outdir}/Bam", mode: 'rellink'
    
    input:
    path "GeneAnnotate.py"
    path "Aligned.sortedByCoord.out.bam"
    path "trans.bam"
    env method 

    output:
    path "ann.bam"

    '''
    echo hello
    python GeneAnnotate.py Aligned.sortedByCoord.out.bam trans.bam ann.bam $method
    '''
}



////////////
//Prepare input bam for remapping
/////////
//Converts mapped bam to fastq
process BamToFastq
{
    input:
    path "input.bam"
    path "PICARD.jar"

    output:
    path "unmapped.fq"

    '''
    java -Xmx2g -jar PICARD.jar SamToFastq I=input.bam FASTQ=unmapped.fq
    '''
    //bedtools bamtofastq -i input.bam -fq unmapped.fq

}

//Takes a fastq file and maps to the transriptome reference using minimap2
//Might need to improve later
process MapBamToTrans 
{
    input:
    path "reads.fq"
    path "ref.mmi"

    output:
    path "alignment.bam"


    '''
    minimap2 -ax map-pb ref.mmi reads.fq > alignment.sam
    samtools view -Sb alignment.sam > alignment.bam
    '''

}


//////////////////
////Steps to make reference
/////////////
//Generate a transcriptome fasta (including premrna) from the genome fasta and a gtf
process GenerateTransFA
{
    input:
    path "genome.fa"
    path "genes.gtf"

    output:
    path "trans.fa"

    '''
    awk \'{if(\$3=="gene"){print \$0}}\' genes.gtf > gene.only.gtf
    gtfToGenePred genes.gtf test.genePhred 
    genePredToBed test.genePhred trans.bed 
    rm test.genePhred
    bedtools getfasta -s -fi genome.fa -bed gene.only.gtf > trans1.fa
    bedtools getfasta -split -s -fi genome.fa -bed trans.bed -name > trans2.fa
    cat trans1.fa trans2.fa > trans.fa
    rm trans1.fa
    rm trans2.fa
    '''

}
        
//Turns transcriptome fasta into a minimap2 reference
process GenerateMinimap2Ref
{
    publishDir "${params.outdir}/Ref", mode: 'rellink'

    input:
    path "trans.fa"

    output:
    path "ref.mmi"

    '''
    minimap2 -d ref.mmi trans.fa
    '''

}

