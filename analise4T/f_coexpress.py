import os

import pandas as pd
import scanpy as sc
import numpy as np
import scanpy as sc
import seaborn
from matplotlib import pyplot as plt
from scipy.spatial import distance




def getNum(path,threshold):
    adata=sc.read_h5ad(path)
    protein_nums=adata[:,adata.var["gene_type"]=="protein_coding"].var["tau_usemean"].values
    protein_nums.sort()
    lncRNA_nums = adata[:, adata.var["gene_type"] == "lncRNA"].var["tau_usemean"].values
    lncRNA_nums.sort()
    return protein_nums[-int(len(protein_nums)*threshold)],lncRNA_nums[-int(len(lncRNA_nums)*threshold)]

def get_genelist(path,pt,lt):
    adata=sc.read_h5ad(path)
    proteins=adata[:,(adata.var["gene_type"]=="protein_coding") & (adata.var["tau_usemean"]>=pt)&(adata.var["vars"] <= 0.5)]
    lncRNAs=adata[:,(adata.var["gene_type"]=="lncRNA") & (adata.var["tau_usemean"]>=lt)&(adata.var["vars"] <= 0.5)]
    proteins.write(path.replace("/h5ad/","/enrichment_analysis_genes/h5ad/").replace(".h5ad","_proteins.h5ad"))
    lncRNAs.write(path.replace("/h5ad/", "/enrichment_analysis_genes/h5ad/").replace(".h5ad", "lncRNAs.h5ad"))
    p_clusters=list(set(proteins.var["tau_cluster"].values.tolist()))
    p_clusters.sort()
    l_clusters = list(set(lncRNAs.var["tau_cluster"].values.tolist()))
    l_clusters.sort()
    pdata={}
    for pc in p_clusters:
        pdata[pc]=proteins.var[proteins.var["tau_cluster"]==pc]["gene"].values.tolist()
    df=fillDataFrame(pdata)
    df.to_csv(path.replace("/h5ad/", "/enrichment_analysis_genes/csv/").replace(".h5ad", "_protein.csv"), index=None)
    ldata = {}
    for lc in l_clusters:
        ldata[lc] = lncRNAs.var[lncRNAs.var["tau_cluster"] == lc]["gene"].values.tolist()
    df = fillDataFrame(ldata)
    df.to_csv(path.replace("/h5ad/", "/enrichment_analysis_genes/csv/").replace(".h5ad", "_lncRNAs.csv"), index=None)

def fillDataFrame(data):
    max_length = max(len(v) for v in data.values())
    for key in data:
        data[key] += [''] * (max_length - len(data[key]))

    return pd.DataFrame(data)
def getCoexpress(path,genelist):
    n=len(genelist)
    adata=sc.read_h5ad(path)
    res=[]
    for i in range(n):
        for j in range(i+1,n):
            g1=adata[:,genelist[i]].X.todense()
            g2=adata[:,genelist[j]].X.todense()
            bool1=np.where(g1>0)
            bool2=np.where(g2>0)
            nonzero=list(set(list(bool1[0])+list(bool2[0])))
            nonzero.sort()
            sumlength = len(g1)
            dis=distance.euclidean(np.array(g1[nonzero]).reshape(-1),np.array(g2[nonzero]).reshape(-1))
            dis=(len(nonzero)/sumlength)*dis-(sumlength-len(nonzero))/sumlength
            ans=genelist[i]+"\t"+genelist[j]+"\t"+str(dis)+"\n"
            res.append(ans)
    with open(path.replace("/qch5ad_tau/h5ad/","/coexpress/").replace("_qc_tau.h5ad",".txt"), 'w') as f:
        f.writelines("x\ty\tdis\n")
        f.writelines(res)


# copath = "../data/raw/gse156728/CD8/coexpress/GSE156728_BC_10X.CD8.txt"
# data=pd.read_table(copath)
# newdata=data[data["dis"]>1]
# print()
# print(data["dis"][np.where(data["dis"]>1)])
# seaborn.distplot(newdata["dis"],bins=20)
# plt.show()