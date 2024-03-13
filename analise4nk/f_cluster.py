import re

import numpy as np
import scanpy as sc
from matplotlib import pyplot as plt
cd56lcd16h_path="../data/raw/gse212890/h5ad/cd56lcd16h.h5ad"
cd56hcd16l_path="../data/raw/gse212890/h5ad/cd56hcd16l.h5ad"


def cluster(path):
    adata=sc.read_h5ad(path)
    print(adata)
    # sc.pl.highest_expr_genes(adata, n_top=20,show=False)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # 将线粒体基因标记为 mt
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    # sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
    #              jitter=0.4, multi_panel=True,show=False)

    # sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt',show=False)
    # sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',show=False)

    sc.pp.normalize_total(adata,target_sum=1e6)
    sc.pp.log1p(adata)



    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    # sc.pl.highly_variable_genes(adata,show=False)
    print(adata)

    adata.raw = adata
    # 获取只有特异性基因的数据集
    # adata = adata[:, adata.var.highly_variable]
    # 回归每个细胞的总计数和表达的线粒体基因的百分比的影响。
    # sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    # 将每个基因缩放到单位方差。阈值超过标准偏差 10。
    # sc.pp.scale(adata,zero_center=False, max_value=10)
    print(adata)
    # 绘制 PCA 图
    sc.tl.pca(adata, svd_solver='arpack')
    # sc.pl.pca(adata, color=["FGFBP2", "CCL3", "CX3CR1", "GNLY", "NKG7", "GZMH", "PRF1","CREM","RGS1"],show=False)
    sc.pl.pca_variance_ratio(adata, log=False,show=False)

    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10,method='gauss',knn=True)
    # exit()

    # 如果设置了 adata 的 .raw 属性时，下图显示了“raw”（标准化、对数化但未校正）基因表达矩阵。
    # sc.pl.umap(adata, color=['CCL3',"CREM","RGS1","IL7R","GZMH"],show=False)

    # sc.pl.umap(adata, color=['CST3', 'CCL3',"RGS1","IL7R","GZMH"], use_raw=False,show=False)
    sc.tl.umap(adata)
    # 计算
    sc.tl.leiden(adata)
    print(adata)
    # 绘制
    sc.pl.umap(adata, color=['leiden'],legend_fontsize=12, legend_loc='on data',show=False)

    sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
    sc.pl.rank_genes_groups(adata, n_genes=40, sharey=False)
    print(adata)

    print(adata.obs["leiden"])
    print(adata.uns["rank_genes_groups"])
    adata.write(path.replace(".h5ad","_gauss_knn_neigh100_pc10.h5ad"))


cluster(cd56lcd16h_path)