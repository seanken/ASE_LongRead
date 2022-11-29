import numpy as np
import pandas
from scipy.sparse import csr_matrix
import sys
import os
from scipy.io import mmwrite

def MakeSparseMatrix(dat):
	rows_dct={}
	cols_dct={}
	rows=[]
	cols=[]
	vals=[]
	numRow=0;
	numCol=0;
	print("Prep for matrix")
	rowName={}
	colName={}
	print(dat.shape)
	for i in range(0,dat.shape[0]):
		gene=dat.iat[i,0]
		cell=dat.iat[i,1]
		count=dat.iat[i,2]
		vals.append(count)
		geneNum=rows_dct.get(gene,numRow)
		if geneNum==numRow:
			rows_dct[gene]=geneNum
			numRow=numRow+1
		rowName[geneNum]=gene
		cellNum=cols_dct.get(cell,numCol)
		if cellNum==numCol:
			cols_dct[cell]=cellNum
			numCol=numCol+1
		colName[cellNum]=cell
		rows.append(geneNum)
		cols.append(cellNum)
	print("Make sparse matrix!")
	mat=csr_matrix((np.array(vals), (np.array(rows), np.array(cols))),shape=(numRow,numCol))
	print("Save!")
	mmwrite("matrix.mtx",mat)
	fil=open("features.tsv","w");
	for i in range(0,numRow):
		fil.write(rowName[i]+"\t"+rowName[i]+"\tGene Expression\n")
	fil.close()
	fil=open("barcodes.tsv","w");
	for i in range(0,numCol):
		fil.write(colName[i]+"\n")
	fil.close()
	print("compress")
	os.system("gzip matrix.mtx")
	os.system("gzip features.tsv")
	os.system("gzip barcodes.tsv")

if __name__=="__main__":
	args=sys.argv
	infil=args[1]
	dat=pandas.read_csv(infil,sep="\t")
	MakeSparseMatrix(dat)
