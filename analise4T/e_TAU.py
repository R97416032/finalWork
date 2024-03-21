import os

import numpy as np
import pandas as pd
import scanpy as sc
import seaborn
from matplotlib import pyplot as plt


def tau(path):
    adata=sc.read_h5ad(path)
    # sc.pl.highest_expr_genes(adata, n_top=20, )
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata.raw = adata
    tau_mean=[]
    tau_mean_nonzero=[]
    tau_median=[]
    tau_cluster4mean=[]
    tau_cluster4mean_nonzero=[]
    tau_cluster4median=[]
    maxzero_per=[]
    medianzero_per=[]
    genes=adata.var["gene"]
    clusters=list(set(adata.obs["meta.cluster"]))
    clusters.sort()
    vars = []
    for g in genes:
        mean=[]
        median=[]
        mean_nonzero=[]
        zeros=[]
        vars.append(np.var(adata[:,adata.var["gene"]==g].X.toarray()))
        for c in clusters:
            temp=adata[adata.obs["meta.cluster"]==c,adata.var["gene"]==g].X.toarray()
            mean.append(np.mean(temp))
            if len(temp[np.nonzero(temp)])>0:
                mean_nonzero.append(np.mean(temp[np.nonzero(temp)]))
            else:
                mean_nonzero.append(0)
            median.append(np.median(temp))
            zeros.append(len(np.nonzero(temp)[0])/len(temp))


        a=1/(len(clusters)-1)
        maxmean=np.max(mean)
        maxmedian=np.max(median)



        maxmean_nonzero=np.max(mean_nonzero)
        tau_cluster4mean.append(clusters[np.argmax(mean)])
        tau_cluster4mean_nonzero.append(clusters[np.argmax(mean_nonzero)])
        tau_cluster4median.append(clusters[np.argmax(median)])
        maxzero_per.append(np.max(zeros))
        medianzero_per.append(np.median(zeros))


        if maxmean!=0:
            tau_mean.append(a*(len(clusters)-np.sum(mean)/maxmean))
        else:
            tau_mean.append(0)
        if maxmean_nonzero!=0:
            tau_mean_nonzero.append(a*(len(clusters)-np.sum(mean_nonzero)/maxmean_nonzero))
        else:
            tau_mean_nonzero.append(0)
        if maxmedian!=0:
            median=np.array(median)
            minmedian = np.min(median[np.nonzero(median)[0]])
            if minmedian!=maxmedian:
                tau_median.append(a * (len(clusters) - (np.sum(median) - len(np.nonzero(median)[0]) * minmedian) / (
                            maxmedian - minmedian)))
            else:
                tau_median.append(a * (len(clusters) - np.sum(median) /
                            maxmedian ))

        else:
            tau_median.append(0)
    adata.var["tau_usemean"]=tau_mean
    adata.var["tau_usemean_nonzero"]=tau_mean_nonzero
    adata.var["tau_usemedian"]=tau_median
    adata.var["maxzero_per"]=np.around(maxzero_per, 4)
    adata.var["medianzero_per"]=np.around(medianzero_per, 4)

    adata.var["tau_cluster4mean"] = tau_cluster4mean
    adata.var["tau_cluster4mean_nonzero"] = tau_cluster4mean_nonzero
    adata.var["tau_cluster4median"] = tau_cluster4median

    adata.var["tau_usemean"].replace(-np.inf,0,inplace=True)
    adata.var["tau_usemean_nonzero"].replace(-np.inf,0,inplace=True)
    adata.var["tau_usemedian"].replace(-np.inf,0,inplace=True)
    # print(vars)
    adata.var["vars"]=vars
    adata.var.to_csv(path.replace('.h5ad', '_tau_var.csv').replace("/qch5ad/","/qch5ad_tau/csv/"),index=None)
    adata.write(path.replace(".h5ad", "_tau.h5ad").replace("/qch5ad/","/qch5ad_tau/h5ad/"))


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





