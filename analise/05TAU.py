
import scanpy as sc
import numpy as np
from matplotlib import pyplot as plt

cd56lcd16h_path="../data/raw/gse212890/h5ad/cd56lcd16h.h5ad"
cd56hcd16l_path="../data/raw/gse212890/h5ad/cd56hcd16l.h5ad"
def tau(path):
    adata=sc.read_h5ad(path)
    tau_mean=[]
    tau_median=[]
    genes=adata.var["gene"]
    clusters=list(set(adata.obs["celltype"]))
    clusters.sort()
    print(clusters)
    vars = []
    for g in genes:
        mean=[]
        median=[]
        vars.append(np.var(adata[:,adata.var["gene"]==g].X.toarray()))
        for c in clusters:
            mean.append(np.mean(adata[adata.obs["celltype"]==c,adata.var["gene"]==g].X.toarray()))
            median.append(np.median(adata[adata.obs["celltype"] == c, adata.var["gene"] == g].X.toarray()))
        a=1/(len(clusters)-1)
        maxmean=np.max(mean)
        maxmedian=np.max(mean)
        if maxmean!=0:
            tau_mean.append(a*(len(clusters)-np.sum(mean)/maxmean))
        else:
            tau_mean.append(0)
        if maxmedian!=0:
            tau_median.append(a*(len(clusters)-np.sum(median)/maxmedian))

        else:
            tau_median.append(0)
    adata.var["tau_usemean"]=tau_mean
    adata.var["tau_usemedian"]=tau_median
    adata.var["tau_usemean"].replace(-np.inf,0,inplace=True)
    adata.var["tau_usemedian"].replace(-np.inf,0,inplace=True)
    print(vars)
    print(adata)
    print(len(vars))
    adata.var["vars"]=vars


    adata.var.to_csv(path.replace('.h5ad', '_tau_var.csv').replace("/h5ad/","/h5ad/tau/"))
    adata.write(path.replace(".h5ad", "_tau.h5ad").replace("/h5ad/","/h5ad/tau/"))
tau(cd56lcd16h_path)
tau(cd56hcd16l_path)