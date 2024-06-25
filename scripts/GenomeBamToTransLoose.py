import pysam
from intervaltree import Interval,IntervalTree
import sys

##
##The main function
##
def main(gtf,bamfile,outfile):
    print("Load GTF")
    gtfRet=loadGTF(gtf)
    gtfExons=gtfRet[0]
    gtfTrans=gtfRet[1]


    print("Process Reads")

    print("First read through")
    dct_score={} ##Used to remove alternative alignments so only one per read
    genefil=pysam.AlignmentFile(bamfile, "rb")
    for seq in genefil:
        if seq.is_supplementary:
            continue;
        readname=seq.query_name
        if not seq.has_tag("AS"):
            continue;
        curScore=seq.get_tag("AS")
        oldScore=dct_score.get(readname,0)
        if oldScore>curScore+.25:
            continue; 
        if curScore>oldScore+.25:
            dct_score[readname]=curScore
            continue;
        dct_score[readname]=curScore+.5 ##Only non integer if multiple equally good alignments
    genefil.close()

    print("Second")
    numReads=0;
    numAssigned=0;
    numOverlap=0;
    inSam=pysam.AlignmentFile(bamfile, "rb")
    outSam=pysam.AlignmentFile(outfile,"wb",template=inSam)	
    for read in inSam:
        if read.is_supplementary:
            continue;
        numReads=numReads+1
        if numReads % 100000==0:
            print(numReads)
            print(numAssigned)
            print(numOverlap)
            print(numOverlap/numReads)
            print(numAssigned/numReads)
            print(" ")

        if not read.has_tag("AS"):
            continue;
        curScore=read.get_tag("AS")
        readname=read.query_name
        oldScore=dct_score.get(readname,0)
        if curScore+.25 < oldScore:
            continue; 
        
        chrom=read.reference_name
        start=read.reference_start
        end=read.reference_end
        exons=ReadToExons(read)
        #print(read.cigarstring)
        if not chrom in gtfTrans:
            continue
        tree=gtfTrans[chrom]
        overlapTrans=[interval.data for interval in tree.overlap(start,end)] ##get transcripts overlapping region covered by read
        if len(overlapTrans)==0:
            continue;
        overlap=[CompareExons(exons,gtfExons[trans],start,end,trans) for trans in overlapTrans]
        overlap=[stat for stat in overlap if stat[2]/stat[1]>.95 and stat[2]/stat[0]>.95]
        if len(overlap)==0:
            continue;
        numOverlap=numOverlap+1
        maxOverlap=max([stat[2] for stat in overlap])
        overlap=[stat[3] for stat in overlap if stat[2]>maxOverlap-10]
        #print(overlap)
        

        ##Check if all same gene, if not continue
        genes=[trans.split("_")[1] for trans in overlap]
        if len([g for g in genes if g==genes[0]])!=len(genes):
            continue;
        gene=genes[0]
        transNames=[trans.split("_")[0] for trans in overlap]
        trans="spliced_"+gene
        if "premrna" in transNames:
            trans="ambig_"+gene
        if len(overlap)==1:
            trans=overlap[0]
        read.set_tag("XT",trans,"Z")
        outSam.write(read)
        numAssigned=numAssigned+1
    inSam.close()
    outSam.close()




def loadGTF(gtf,addChr=True):
    fil=open(gtf,"r")
    gtfExons={}
    gtfTrans={}
    for line in fil:
        if line[0]=="#":
            continue;
        s=line.strip().split()
        chrom=s[0]
        if addChr:
            chrom="chr"+chrom
        regionType=s[2]
        start=int(s[3])
        end=int(s[4])
        transID=""
        transNames=[s[i+1].strip(";").strip("\"") for i in range(0,len(s)-1) if s[i]=="transcript_id"]
        geneNames=[s[i+1].strip(";").strip("\"") for i in range(0,len(s)-1) if s[i]=="gene_id"]
        if len(transNames)==1 and len(geneNames)==1:
            transID=transNames[0]+"_"+geneNames[0]
        if regionType=="gene":
            if len(geneNames)==1:
                transID="premrna_"+geneNames[0]

        if len(transID)==0:
            continue;
        #transID=transID.strip(";").strip("\"")
        if regionType=="transcript" or regionType=="gene":
            regions=gtfTrans.get(chrom,[])
            regions.append([start,end,transID])
            gtfTrans[chrom]=regions
        if regionType=="exon" or regionType=="gene":
            exons=gtfExons.get(transID,[])
            exons.append([start,end])
            gtfExons[transID]=exons
    fil.close()
    ##Change to interval tree
    print("make tree")
    for chrom in gtfTrans.keys():
        listTrans=gtfTrans[chrom]
        tree=IntervalTree([Interval(lst[0],lst[1],lst[2]) for lst in listTrans])
        gtfTrans[chrom]=tree
    ##Sort exons
    print("sort")
    for trans in gtfExons.keys():
        gtfExons[trans]=sorted(gtfExons[trans])
    ret=[gtfExons,gtfTrans]
    return(ret)


##
##Given the exon regions from a long read and a gtf entry compares them
##
def CompareExons(readExons,gtfExons,startRead,endRead,trans):
    #print(gtfExons)
    #print(readExons)
    gtfExons=[[max(exon[0],startRead),min(exon[1],endRead)] for exon in gtfExons if exon[0]<endRead and exon[1]>startRead]
    if startRead>endRead-10 or len(gtfExons)<1:
        return([1,1,0])
    lenRead=sum([(exon[1]-exon[0]+1) for exon in readExons])
    lenGTF=sum([(exon[1]-exon[0]+1) for exon in gtfExons])
    interExons=intersectExons(gtfExons,readExons)
    if len(interExons)<1:
        return([1,1,0])
    lenInter=sum([(exon[1]-exon[0]+1) for exon in interExons])
    return([lenRead,lenGTF,lenInter,trans])


##
##Overlaps 2 sets of exons, assumes both are ordered
##
def intersectExons(exons1,exons2):
    i1=0
    i2=0
    interExons=[]
    while i1<len(exons1) and i2<len(exons2):
        exon1=exons1[i1]
        exon2=exons2[i2]
        start=max(exon1[0],exon2[0])
        end=min(exon1[1],exon2[1])
        if start<end:
            interExons.append([start,end])
        if exon1[1]>exon2[1]:
            i2=i2+1
        else:
            i1=i1+1
    return(interExons)


##
##Goes from a long read to exon coverage
##
def ReadToExons(read):
    blocks=read.get_blocks()
    cleanBlocks=[]
    if len(blocks)<1:
        return([]);
    start=-1
    end=-1
    for block in blocks:
        startNew=block[0]
        endNew=block[1]
        if start<0:
            start=startNew
            end=endNew
            continue;
        if startNew<end+5:
            end=endNew
            continue;
        cleanBlocks.append([start,end])
        start=startNew
        end=endNew
    cleanBlocks.append([start,end])
    return(cleanBlocks)


if __name__=="__main__":
    gtf=sys.argv[1]
    bamfile=sys.argv[2]
    outfile=sys.argv[3]
    main(gtf,bamfile,outfile)

