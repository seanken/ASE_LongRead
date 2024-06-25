import pysam

##inBam is input bam, should be sorted by read name
def removeDup(inBam,outBam):
    samfile = pysam.AlignmentFile(bamfile, "rb")
    outfile = pysam.AlignmentFile(bamfile, "wb",template=samfile) ##Add header part
    seqLast=None
    matched=False
    for seq in samfile:
        if seqLast is None:
            seqLast=seq
            continue;
        namLast=seqLast.query_name
        nam=seq.query_name
        if nam!=namLast and not matched:
            outfile.write(seqLast)
        matched=True
        if nam!=namLast:
            matched=False
        seqLast=seq
    if not matched:
        outfile.write(seqLast)
    samfile.close()
    outfile.close()
