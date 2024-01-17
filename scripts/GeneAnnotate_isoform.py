import sys
import pysam


def AnnotBam(bam,gene_ann,outfile,useRQ=False,verbose=False):
	samfile=pysam.AlignmentFile(bam, "rb")
	outbam=pysam.AlignmentFile(outfile,"wb",template=samfile)	
	genefil=pysam.AlignmentFile(gene_ann, "rb")
	dct_read_iso={}
	dct_score={}
	for seq in genefil:
		if seq.is_supplementary:
			continue;
		trans=seq.reference_name
		readname=seq.query_name
		curTrans=dct_read_iso.get(readname,trans)
		if not seq.has_tag("AS"):
			continue;
		curScore=seq.get_tag("AS")
		oldScore=dct_score.get(readname,0)
		if curTrans==trans:
			dct_read_iso[readname]=trans
			dct_score[readname]=max(oldScore,curScore)
		else:
			if oldScore>curScore:
				continue;
			if oldScore==curScore:
				dct_read_iso[readname]="MULTI"
				continue;
			if curScore>oldScore:
				dct_score[readname]=curScore
				dct_read_iso[readname]=trans	
	genefil.close()
	
	iterat=0
	for seq in samfile:
		iterat=iterat+1;
		if iterat % 100000==0:
			if verbose:
				print("Iteration: "+ str(iterat));

		if useRQ:
			rq=seq.get_tag("rq"); ##read quality check
			if rq<.99:
				continue;
		readnam=seq.query_name
		gene=dct_read_iso.get(readnam,"MULTI")
		if not gene=="MULTI":
			seq.set_tag("XT",gene,"Z")
			outbam.write(seq)
	outbam.close()
	samfile.close()


if __name__=="__main__":
	args=sys.argv
	bam=args[1]
	gene_ann=args[2]
	outbam=args[3]
	method=args[4]
	useRQ=False
	if method=="MAS":
		useRQ=True
	AnnotBam(bam,gene_ann,outbam,useRQ,verbose=True)

