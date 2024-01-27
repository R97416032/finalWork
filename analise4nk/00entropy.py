import math

import numpy as np
import scanpy as sc


cd56hcd16l_min0="../data/raw/gse212890/h5ad/cd56hcd16l_color_min0.h5ad"
cd56hcd16l_min_not0="../data/raw/gse212890/h5ad/cd56hcd16l_color_min_not0.h5ad"
cd56lcd16_min0="../data/raw/gse212890/h5ad/cd56lcd16h_color_min0.h5ad"
cd56lcd16_min_not0="../data/raw/gse212890/h5ad/cd56lcd16h_color_min_not0.h5ad"
def cul_entropy(path):
    adata=sc.read_h5ad(path)
    entropy=[]
    n_obs=len(adata.obs_names)
    for g in adata.var["gene"]:
        colors=adata[:,g].X.toarray().reshape(-1)
        clist=list(set(colors))
        clist.sort()
        ans=0
        for c in clist:
            ans=ans-(sum(colors==c)/n_obs)*math.log2(sum(colors==c)/n_obs)
        entropy.append(ans)
    adata.var["entropy"]=entropy
    adata.write(path.replace(".h5ad","_entropy.h5ad").replace("/h5ad/","/h5ad/entropy/"))

cul_entropy(cd56hcd16l_min0)
cul_entropy(cd56hcd16l_min_not0)
cul_entropy(cd56lcd16_min0)
cul_entropy(cd56lcd16_min_not0)
