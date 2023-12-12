
import scanpy as sc
import numpy as np
from matplotlib import pyplot as plt

cd56lcd16h_path="../data/raw/gse212890/h5ad/cd56lcd16h_gauss_knn_neigh100_pc10.h5ad"
cd56hcd16l_path="../data/raw/gse212890/h5ad/cd56hcd16l_gauss_noknn_neigh100_pc10.h5ad"
def tau(path):
    adata=sc.read_h5ad(path)
    tau=[]
    genes=adata.var["gene"]
    clusters=list(set(adata.obs["leiden"]))
    clusters.sort()
    print()
    for g in genes:
        mean=[]
        for c in clusters:
            mean.append(np.mean(adata[adata.obs["leiden"]==c,adata.var["gene"]==g].X))
        a=1/(len(clusters)-1)
        maxi=np.max(mean)
        tau.append(a*(len(clusters)-np.sum(mean)/maxi))

    adata.var["tau"]=tau
    adata.var["tau"].replace(-np.inf,0,inplace=True)

    adata.var.to_csv(path.replace('.h5ad','_var.csv'))
    adata.write(path.replace(".h5ad","_tau.h5ad"))
tau(cd56lcd16h_path)