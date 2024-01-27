import pandas as pd
import scanpy as sc
def read(cellinfo_path,geneinfo_path,counts_path):
    cellinfo = pd.read_csv(cellinfo_path, index_col=0)
    geneinfo = pd.read_csv(geneinfo_path, index_col=0)
    adata = sc.read(counts_path)
    adata = sc.AnnData(adata.X, obs=cellinfo, var=geneinfo)
    adata.obs_names.name = 'cell'
    adata.var.index = adata.var['gene']
    adata.var_names.name = 'gene'
    adata.write(cellinfo_path.split("ungz/")[0]+"h5ad/all.h5ad")
    cd56hcd16l=adata[adata.obs["Majortype"]=="CD56highCD16low"]
    cd56lcd16h=adata[adata.obs["Majortype"]=="CD56lowCD16high"]
    cd56lcd16h.write(cellinfo_path.split("ungz/")[0]+"h5ad/cd56lcd16h.h5ad")
    cd56hcd16l.write(cellinfo_path.split("ungz/")[0]+"h5ad/cd56hcd16l.h5ad")
cellinfo_path="../data/raw/gse212890/ungz/GSE212890_NK_metadata.csv"
geneinfo_path="../data/raw/gse212890/ungz/GSE212890_NK.genes.csv"
counts_path="../data/raw/gse212890/ungz/GSE212890_NK_counts.mtx"

read(cellinfo_path,geneinfo_path,counts_path)