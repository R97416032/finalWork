import numpy as np
import scanpy as sc

def tau(path):
    adata=sc.read_h5ad(path)
    print(adata)
    sc.pl.highest_expr_genes(adata, n_top=20, )
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata.raw = adata
    print(adata)
    tau_mean=[]
    tau_median=[]
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


    # adata.var.to_csv(path.replace('.h5ad', '_tau_var.csv').replace("/h5ad/","/h5ad/tau/"))
    # adata.write(path.replace(".h5ad", "_tau.h5ad").replace("/h5ad/","/h5ad/tau/"))


path="../data/raw/gse156728/CD8/qch5ad/GSE156728_BC_10X.CD8_qc.h5ad"
#注意路径需要qc后的新路径
tau(path)