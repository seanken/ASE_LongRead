import sys
import pysam

def AnnotBam(bam,gene_ann,outfile,verbose=False):
	samfile=pysam.AlignmentFile(bam, "rb")
	outbam=pysam.AlignmentFile(outfile,"wb",template=samfile)	
	genefil=open(gene_ann,"r")
	geneInfo=""
	iterat=0
	for seq in samfile:
		iterat=iterat+1;
		if iterat % 100000==0:
			if verbose:
				print("Iteration: "+ str(iterat));
		geneInfo=genefil.readline().strip().split();
		rq=seq.get_tag("rq") ##read quality check
		if rq<.99:
			continue;
		gene=geneInfo[3]
		numAssign=int(geneInfo[2])
		assign=geneInfo[1]
		nam=geneInfo[0]
		readnam=seq.query_name
		if nam!=readnam:
			print("Misaligned reads!")
			print(nam)
			print(readnam)
			break;
		seq.set_tag("XT",gene,"Z")
		seq.set_tag("XS",assign,"Z")
		seq.set_tag("XN",numAssign,"i")
		outbam.write(seq)
	outbam.close()
	samfile.close()


if __name__=="__main__":
	args=sys.argv
	bam=args[1]
	gene_ann=args[2]
	outbam=args[3]

	AnnotBam(bam,gene_ann,outbam,verbose=True)

