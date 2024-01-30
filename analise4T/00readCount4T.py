import os

import numpy as np
import pandas as pd
import polars as pl
import scanpy as sc
from scipy.sparse import csr_matrix
from gtfparse import read_gtf

def read(countpath,metapath):
    data=sc.read(countpath)
    metainfo=pd.read_table(metapath)
    cellinfo=metainfo[metainfo["cellID"].isin(data.var.index)]
    adata=sc.AnnData(data.X.T,var=data.obs,obs=cellinfo)
    adata.X = csr_matrix(adata.X)
    adata.write(countpath.replace("/ungz/","/h5ad/").replace(".counts.txt",".h5ad"))

def addGenetype(h5adpath):
    adata=sc.read_h5ad(h5adpath)


    adata.var["gene"]=adata.var.index
    genes=adata.var["gene"]
    genelist=[]
    genetypes=[]
    for gene in genes:
        if gene in data["gene_name"].to_list():
            genetypes.append(data[data["gene_name"]==gene]["gene_type"].to_list()[0])
            genelist.append(gene)

        elif gene in lncRNAs:
            genetypes.append("lncRNA")
            genelist.append(gene)
    newadata = adata[:, genelist]
    newadata.var["gene_type"]=genetypes
    newadata.write(h5adpath.replace("/h5ad/","/newh5ad/"))

# countpath="../data/raw/gse156728/CD8/ungz/GSE156728_BC_10X.CD8.counts.txt"
# metapath="../data/raw/gse156728/metadata/ungz/GSE156728_metadata.txt"
# path="../data/raw/gse156728/CD8/ungz/"
# names=os.listdir(path)
# for n in names:
#     read(path+n,metapath)
# h5adpath="../data/raw/gse156728/CD8/h5ad/GSE156728_BC_10X.CD8.h5ad"
# addGenetype(h5adpath)
path="../data/raw/gse156728/CD8/h5ad/"
names=os.listdir(path)
genecodepath = "../data/genecode/v44/all/ungz/gencode.v44.chr_patch_hapl_scaff.annotation.gtf"
lncRNApath = "../data/lncipedia/lncRNA268.txt"
genecode = read_gtf(genecodepath)
lncRNAs = pd.read_table(lncRNApath, header=None)[0].to_list()
data = genecode[genecode["gene_type"].isin(["lncRNA","protein_coding"])]
i=0
for n in names:
    print(path+n)
    addGenetype(path+n)
    print(">>>>>>>>>>>>>>>>>" + str(i))
    i += 1