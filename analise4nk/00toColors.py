import numpy as np
import scanpy as sc


cd56lcd16h="../data/raw/gse212890/h5ad/cd56lcd16h.h5ad"
cd56hcd16l="../data/raw/gse212890/h5ad/cd56hcd16l.h5ad"
def color(v,res):
    if v >= 380 and v < 420:
        res.append(7)
    elif v >= 420 and v < 450:
        res.append(6)
    elif v >= 450 and v < 490:
        res.append(5)
    elif v >= 490 and v < 560:
        res.append(4)
    elif v >= 560 and v < 590:
        res.append(3)
    elif v >= 590 and v < 620:
        res.append(2)
    elif v >= 620 and v <= 780:
        res.append(1)
    else:
        res.append(0)
def tocolors(path):
    adata = sc.read_h5ad(path)

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # 将线粒体基因标记为 mt
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)
    n_vars = len(adata.var_names)
    n_obs =len(adata.obs_names)
    X = np.zeros(n_obs * n_vars).reshape(n_obs, n_vars)
    XX = np.zeros(n_obs * n_vars).reshape(n_obs, n_vars)
    print(adata)
    cells=adata.obs["cellID"]
    print(cells)
    for c,i in zip(cells,range(n_obs)):
        counts=adata[c,:].X.toarray()[0]
        color0=[]
        color1=[]
        maxcount=max(counts)
        mincount=min(counts)
        min_notzero=min(filter(lambda x: x > 0, counts))
        for count in counts:
            wavelen0=(count-mincount)*(780-380)/(maxcount-mincount)+380
            wavelen1=(count-min_notzero)*(780-380)/(maxcount-min_notzero)+380
            color(wavelen0,color0)
            color(wavelen1,color1)
        X[i,:]=color0
        XX[i,:]=color1

    adata.X = X
    adata.write(path.replace(".h5ad", "_color_min0.h5ad"))
    adata.X=XX
    adata.write(path.replace(".h5ad","_color_min_not0.h5ad"))
# tocolors(cd56hcd16l)
tocolors(cd56lcd16h)