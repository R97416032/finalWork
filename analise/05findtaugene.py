import pandas as pd
import scanpy as sc
#tau路径
tauhl="../data/raw/gse212890/h5ad/tau/cd56hcd16l_tau_var.csv"
taulh="../data/raw/gse212890/h5ad/tau/cd56lcd16h_tau_var.csv"

#数据路径
h5adhl="../data/raw/gse212890/h5ad/tau/cd56hcd16l_tau.h5ad"
h5adlh="../data/raw/gse212890/h5ad/tau/cd56lcd16h_tau.h5ad"

def find(tau,h5ad,threshold):
    tau=pd.read_csv(tau)
    adata=sc.read_h5ad(h5ad)
    tau_genes = tau[tau["tau_usemean"] >= threshold][["gene", "tau_usemean"]]
    for g in tau_genes["gene"]:
        ax=sc.pl.dotplot(adata, g ,groupby="celltype",show=False,save=g+".png")
        print(ax)
        print(ax['mainplot_ax'])
        exit()

find(tauhl,h5adhl,0.75)