# Pipeline for processing long reads for ASE

This pipeline is built to take in long read single cell or single nucleus 10X data (ONT, PacBio, or MAS-Seq) and producing ASE specific counts at the gene or isoform level.

## Upstream processing

The pipeline assumes certain processing upstream. In particular, it assumes a bam file mapped to the reference genome (preferably with minimap2 though should work with others), with tags corresponding to UMI and CBC.

For upstream processing there are numerous choices. For MAS-Seq data we have an example upstream processing pipeline that takes the output of the skera tool, see in the Upstream directory. One needs to pass it a list of the MAS-Seq primers (Primers argument), a list of CBC (CBC argument), a fasta of the reference genome (refFA), and the input bam from Skera (inputSkeraBam). Can also specify the output directory (outdir). For example:

```
nextflow /path/to/Revio.nf --Primers $primers --CBC $cbc --refFA $fa --outdir output --inputSkeraBam $bam
```

which would output a processed bam in output.

## ASE Pipeline

For extracting gene level ASE from long read data we use ASEFromBam.nf (ASEFromBam.Revio.nf is the same, except has strand hard coded to 1). The options are:

`--bam`: The input bam. Must be mapped to the reference and sorted by genomic position, and have CBC/UMI information as tags.

`--gtf`: A gtf with information about the transcriptome.

`--outdir`: The output directory.

`--input_vcf`: The vcf to use.

`--vcf_col`: The column in the vcf with information about the current sample.

`--umi_tag`: The tag in the bam with umi information.

`--cbc_tag`: The tag in the bam with cbc information.

`--method`: The method used, either MAS or Long. Same, except with MAS uses the rq tag to filter out low qualty reads. 

`--strand`: The strand to use (default 2, with the MAS-Seq pipeline we shared should use 1 instead).

Will produce many ouputs in the outdirectory, including an annotated bam (annotated with gene), gene level count matrices (both sparse (in Matrix) and not sparse (in UMITools)) and gene level allele counts. Can be used with out scAlleleExpression package.

## Isoform level ASE pipeline

If instead of gene level results on wants isoform level resutls (including premrna as a seperate isoform for each gene) one can use ASEFromBamTrans.nf. Same arguments, except no strand argument (ignores strand).

## Required packages to run

In order to run the above ASE pipelines, need the following on the path:
1) bcftools
2) tabix
3) python with pysam for gene level analysis
4) python with pysam, intervaltree, pandas, numpy, and scipy installed for isoform level analysis

For the Revio upstream pipeline need lima, isoseq, samtools, and minimap2.