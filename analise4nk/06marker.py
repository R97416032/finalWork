import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn
from gtfparse import read_gtf
import polars as pl
cd56hcd16l="../data/raw/gse212890/h5ad/tau/cd56hcd16l_tau_var.csv"
cd56lcd16h="../data/raw/gse212890/h5ad/tau/cd56lcd16h_tau_var.csv"

def see_distribution(path,mean):

    adata=pd.read_csv(path)

    if mean:
        seaborn.distplot(adata[adata["gene_type"]=="protein_coding"]['tau_usemean'],bins=50)
        seaborn.distplot(adata[adata["gene_type"]=="lncRNA"]['tau_usemean'],bins=50)
        plt.show()
    else:
        seaborn.distplot(adata[adata["gene_type"] == "protein_coding"]['tau_usemedian'], bins=50)
        seaborn.distplot(adata[adata["gene_type"] == "lncRNA"]['tau_usemedian'], bins=50)
        plt.show()
# see_distribution(cd56lcd16h,True)
# see_distribution(cd56hcd16l,True)

def get_ranks(path,n):
    adata = sc.read_h5ad(path)
    sc.tl.rank_genes_groups(adata,groupby="celltype",n_genes=n)
    sc.pl.rank_genes_groups(adata)
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    df = pd.DataFrame(
        {group + '_' + key: result[key][group]
         for group in groups for key in ['names', 'logfoldchanges', 'pvals', 'pvals_adj']})
    savepath=path.replace(".h5ad",".csv").replace("/h5ad/tau/","/ranks_gene/")
    df.to_csv(savepath)

def get_geneset4type(tau_path,ranks_path,threshold):
    taudata=pd.read_csv(tau_path)
    ranksdata=pd.read_csv(ranks_path)
    tau_genes=taudata[taudata["tau_usemean"]>=threshold][["gene","tau_usemean"]]
    types=[x for x in ranksdata.columns if x.endswith("names")]

    dict={}
    for type in types:
        rankgene=ranksdata[type]
        genes=list(set(tau_genes["gene"])&set(rankgene))
        dict[type.replace("_names","")]=genes
    df = pd.DataFrame.from_dict(dict, orient='index').transpose()
    print(df)
    df.to_csv(ranks_path.replace("/ranks_gene/","/type_genes/").replace("_tau","type_genes"),index=False)

def see_typegenes(typegenes_path,h5adpath):
    adata=sc.read_h5ad(h5adpath)
    print(adata)
    typegenes=pd.read_csv(typegenes_path).fillna("AA")

    dict =typegenes.to_dict(orient='list')
    marker_genes_dict={}
    print(marker_genes_dict)
    for key,value in dict.items():
        marker_genes_dict[key]=[x for x in value[0:20] if x!="AA"]

    print(marker_genes_dict)
    sc.pl.dotplot(adata, marker_genes_dict ,groupby="celltype")


cd56hcd16lh5ad="../data/raw/gse212890/h5ad/tau/cd56hcd16l_tau.h5ad"
cd56lcd16hh5ad="../data/raw/gse212890/h5ad/tau/cd56lcd16h_tau.h5ad"
#获取ranks
# get_ranks(cd56hcd16lh5ad,150)
get_ranks(cd56lcd16hh5ad,300)


#获取marker
# ranks_path="../data/raw/gse212890/ranks_gene/cd56hcd16l_tau.csv"
# get_geneset4type(cd56hcd16l,ranks_path,0.70)
# ranks_path="../data/raw/gse212890/ranks_gene/cd56lcd16h_tau.csv"
# get_geneset4type(cd56lcd16h,ranks_path,0.70)

# typespath_hl="../data/raw/gse212890/type_genes/cd56hcd16ltype_genes.csv"
# see_typegenes(typespath_hl,cd56hcd16lh5ad)

typespath_lh="../data/raw/gse212890/type_genes/cd56lcd16htype_genes.csv"
see_typegenes(typespath_lh,cd56lcd16hh5ad)





