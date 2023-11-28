import pandas as pd
import scanpy as sc
from gtfparse import read_gtf
import polars as pl

genecodepath="../data/genecode/v44/all/ungz/gencode.v44.chr_patch_hapl_scaff.annotation.gtf"
lncRNApath="../data/lncipedia/lncRNAs.txt"
genecode=read_gtf(genecodepath)
lncRNAs=pd.read_table(lncRNApath,header=None)
celldata=sc.read_h5ad("../data/raw/gse212890/h5ad/all.h5ad")
data=genecode.filter((pl.col("gene_type")=="miRNA" )| (pl.col("gene_type")=="lncRNA") | (pl.col("gene_type")=="protein_coding"))

genes=celldata.var["gene"]
gene_type=[]
genelist=[]
for g in genes:
    if g in data["gene_name"].to_list():
        gene_type.append(data.filter(pl.col("gene_name")==g)["gene_type"][0])
        genelist.append(g)
    else:
        if g in lncRNAs:
            gene_type.append("lncRNA")
            genelist.append(g)

df=pd.DataFrame(columns=["gene_type"],data=gene_type)
df.to_csv("../data/raw/gse212890/h5ad/genetype.csv",index=None)
newdata=celldata[:,genelist]
newdata.write("../data/raw/gse212890/h5ad/all.h5ad".replace("all","newAll"))
newdata.var["gene_type"]=gene_type
newdata.write("../data/raw/gse212890/h5ad/all.h5ad".replace("all","newAll"))

