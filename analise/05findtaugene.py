import numpy as np
import pandas as pd
import scanpy as sc
#tau路径
tauhl="../data/raw/gse212890/h5ad/tau/cd56hcd16l_tau_var.csv"
taulh="../data/raw/gse212890/h5ad/tau/cd56lcd16h_tau_var.csv"

#数据路径
h5adhl="../data/raw/gse212890/h5ad/tau/cd56hcd16l_tau.h5ad"
h5adlh="../data/raw/gse212890/h5ad/tau/cd56lcd16h_tau.h5ad"

#marker路径
markerhl="../data/raw/gse212890/marker/cd56hcd16l_marker.csv"
markerlh="../data/raw/gse212890/marker/cd56lcd16h_marker.csv"
def find(tau,h5ad,threshold):
    tau=pd.read_csv(tau)
    adata=sc.read_h5ad(h5ad)
    tau_genes = tau[tau["tau_usemean"] >= threshold][["gene", "tau_usemean"]]
    celltype=list(set(adata.obs["celltype"]))
    celltype.sort()
    dict={}
    for c in celltype:
        dict.update({c:[]})
    for g in tau_genes["gene"]:
        temp=""
        maxexp=0
        for type in celltype:
            if np.mean(adata[adata.obs["celltype"] == type][:,g].X.toarray().reshape(-1))>maxexp :
                if sum(adata[adata.obs["celltype"] == type][:,g].X.toarray().reshape(-1)>0)/len(adata[adata.obs["celltype"] == type][:,g].X.toarray().reshape(-1))>=0.05:
                    temp=type
                    maxexp=np.mean(adata[adata.obs["celltype"] == type][:,g].X.toarray().reshape(-1))
        if temp in celltype:
            mark=dict.get(temp)
            mark.append(g)
            dict.update({temp:mark})
    df = pd.DataFrame.from_dict(dict, orient='index').transpose()
    df.to_csv(h5ad.replace("/h5ad/tau/","/marker/").replace("_tau.h5ad","_marker.csv"))

def cul_percent(h5ad,marker,name):
    adata=sc.read_h5ad(h5ad)
    markers=pd.read_csv(marker,index_col=0).fillna("AA")
    print(markers)
    celltype=markers.columns
    for ctype in celltype:
        mean=[]
        per=[]
        genes=[]
        for g in markers[ctype]:
            if g!="AA":
                genes.append(g)
        for g in genes   :
            per.append(sum(adata[adata.obs["celltype"] == ctype][:, g].X.toarray().reshape(-1) > 0) / len(
                adata[adata.obs["celltype"] == ctype][:, g].X.toarray().reshape(-1)))
            mean.append(np.mean(adata[adata.obs["celltype"] == ctype][:,g].X.toarray().reshape(-1)))
        dict={ctype:genes,"mean":mean,"per":per}
        df=pd.DataFrame.from_dict(dict,orient="index").transpose()
        df.to_csv(marker.replace("/marker/","/marker/"+name+"/").replace("marker.csv",ctype+"_marker.csv"))




find(tauhl,h5adhl,0.75)
find(taulh,h5adlh,0.75)
hlname="cd56hcd16l"
lhname="cd56lcd16h"
cul_percent(h5adlh,markerlh,lhname)
cul_percent(h5adhl,markerhl,hlname)
