import os

import numpy as np
import pandas as pd
import scanpy as sc
import seaborn
from matplotlib import pyplot as plt


def tau(path):
    adata=sc.read_h5ad(path)
    print(adata)
    # sc.pl.highest_expr_genes(adata, n_top=20, )
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata.raw = adata
    print(adata)
    tau_mean=[]
    tau_median=[]
    tau_cluster=[]
    genes=adata.var["gene"]
    clusters=list(set(adata.obs["meta.cluster"]))
    clusters.sort()
    print(clusters)
    vars = []
    for g in genes:
        mean=[]
        median=[]

        vars.append(np.var(adata[:,adata.var["gene"]==g].X.toarray()))
        for c in clusters:
            mean.append(np.mean(adata[adata.obs["meta.cluster"]==c,adata.var["gene"]==g].X.toarray()))
            median.append(np.median(adata[adata.obs["meta.cluster"] == c, adata.var["gene"] == g].X.toarray()))
        a=1/(len(clusters)-1)
        maxmean=np.max(mean)
        # maxmedian=np.max(mean)
        tau_cluster.append(clusters[np.argmax(mean)])
        if maxmean!=0:
            tau_mean.append(a*(len(clusters)-np.sum(mean)/maxmean))
        else:
            tau_mean.append(0)
        # if maxmedian!=0:
        #     tau_median.append(a*(len(clusters)-np.sum(median)/maxmedian))
        #
        # else:
        #     tau_median.append(0)
    adata.var["tau_usemean"]=tau_mean
    adata.var["tau_usemedian"]=tau_median
    adata.var["tau_cluster"] = tau_cluster
    adata.var["tau_usemean"].replace(-np.inf,0,inplace=True)
    adata.var["tau_usemedian"].replace(-np.inf,0,inplace=True)
    # print(vars)
    print(adata)
    print(len(vars))
    adata.var["vars"]=vars
    print(adata)
    adata.var.to_csv(path.replace('.h5ad', '_tau_var.csv').replace("/qch5ad/","/qch5ad_tau/"),index=None)
    adata.write(path.replace(".h5ad", "_tau.h5ad").replace("/qch5ad/","/qch5ad_tau/"))


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

def find(tau,h5ad,threshold):
    tau=pd.read_csv(tau)
    adata=sc.read_h5ad(h5ad)
    tau_genes = tau[tau["tau_usemean"] >= threshold][["gene", "tau_usemean"]]
    celltype=list(set(adata.obs["meta.cluster"]))
    celltype.sort()
    dict={}
    for c in celltype:
        dict.update({c:[]})
    for g in tau_genes["gene"]:
        temp=""
        maxexp=0
        for type in celltype:
            if np.mean(adata[adata.obs["meta.cluster"] == type][:,g].X.toarray().reshape(-1))>maxexp :
                if sum(adata[adata.obs["meta.cluster"] == type][:,g].X.toarray().reshape(-1)>0)/len(adata[adata.obs["meta.cluster"] == type][:,g].X.toarray().reshape(-1))>=0.05:
                    temp=type
                    maxexp=np.mean(adata[adata.obs["meta.cluster"] == type][:,g].X.toarray().reshape(-1))
        if temp in celltype:
            mark=dict.get(temp)
            mark.append(g)
            dict.update({temp:mark})
    df = pd.DataFrame.from_dict(dict, orient='index').transpose()
    df.to_csv(h5ad.replace("/h5ad/","/marker/").replace("_tau.h5ad","_marker.csv"))





