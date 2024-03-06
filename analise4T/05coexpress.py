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
    proteins=adata[:,adata.var["gene_type"]=="protein_coding"].var["tau_usemean"]>=pt
    proteins_var = adata[:, adata.var["gene_type"] == "protein_coding"].var["vars"] <= 0.5
    lncRNAs=adata[:, adata.var["gene_type"] == "lncRNA"].var["tau_usemean"]>=lt
    lncRNAs_var = adata[:, adata.var["gene_type"] == "lncRNA"].var["vars"] <=0.5
    ps=proteins[np.where(proteins)[0]].index
    psv=proteins_var[np.where(proteins_var)[0]].index
    ls=lncRNAs[np.where(lncRNAs)[0]].index
    lsv=lncRNAs_var[np.where(lncRNAs_var)[0]].index
    resp=list(set(ps.values.tolist())&set(psv.values.tolist()))
    resl=list(set(ls.values.tolist())&set(lsv.values.tolist()))

    return resp+resl
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
path = "../data/raw/gse156728/CD8/qch5ad_tau/h5ad/"
names=os.listdir(path)
for n in names:
    print(path+n)
    print(">>>>>>>>>>>>>>")
    a, b = getNum(path+n, 0.10)
    listgene = get_genelist(path+n, a, b)
    print(len(listgene))
    # getCoexpress(path+n, listgene)
    # exit()

# copath = "../data/raw/gse156728/CD8/coexpress/GSE156728_BC_10X.CD8.txt"
# data=pd.read_table(copath)
# newdata=data[data["dis"]>1]
# print()
# print(data["dis"][np.where(data["dis"]>1)])
# seaborn.distplot(newdata["dis"],bins=20)
# plt.show()