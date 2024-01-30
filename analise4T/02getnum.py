import os

import numpy as np
import pandas as pd
import scanpy as sc

def getALlcluster(path):

    names = os.listdir(path)
    types = []
    for n in names:
        allclusters(path + n, types)

    types.sort()
    with open(path.split("/")[-2]+"_clusters.txt", "w") as file:
        # 将list中的每个元素写入到文件中
        for item in types:
            file.write("%s\n" % item)
    return types
def allclusters(path,types):
    adata = sc.read_h5ad(path)
    for i in list(set(adata.obs["meta.cluster"].values.tolist())):
        if i not in types:
            types.append(i)

def getnum417cluster(path):
    types=getALlcluster(path)
    datanames=os.listdir(path)
    datanames.sort()
    df=pd.DataFrame()

    df["names"]=[x.replace(".h5ad","").replace(x.split("_")[0]+"_","").replace("_qc","") for x in datanames]

    for t in types:
        num = []
        for n in datanames:
            adata=sc.read_h5ad(path+n)
            num.append(sum(adata.obs["meta.cluster"]==t))
        df[t]=num
    df.to_csv(path.replace("/qch5ad/","/nums/")+"clusters.csv",index=None)





path = "../data/raw/gse156728/CD8/qch5ad/"
getnum417cluster(path)