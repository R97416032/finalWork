import os

import numpy as np
import scanpy as sc


def qc(path):
    adata=sc.read_h5ad(path)
    print(adata)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # 将线粒体基因标记为 mt
    # mitochondrial genes

    adata.var['mt'] = adata.var_names.str.startswith('MT-')

    # hemoglobin genes. 血红蛋白基因

    adata.var['hb'] = adata.var_names.str.contains('^HB[^P]')

    # ribosomal genes

    adata.var['ribo'] = adata.var_names.str.startswith('RPS', 'RPL')

    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'hb', 'ribo'], percent_top=None, log1p=False, inplace=True)

    # n_genes_by_counts：每个细胞中，有表达的基因的个数；
    #
    # total_counts：每个细胞的基因总计数（总表达量, umi数）；
    #
    # pct_counts_mt：每个细胞中，线粒体基因表达量占该细胞所有基因表达量的百分比
    #
    # pct_counts_hb: 每个细胞中，血红蛋白基因表达量占该细胞所有基因表达量的百分比
    #
    # pct_counts_ribo: 每个细胞中，核糖体RNA基因表达量占该细胞所有基因表达量的百分比
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo', 'pct_counts_hb'],
                 jitter=0.4, multi_panel=True, show=False,save="_"+path.split("/")[-1].replace(path.split("/")[-1].split("_")[0]+"_","").replace(".h5ad","_")+"QC.pdf")

    mito_genes = adata.var_names.str.startswith('MT-')

    malat1 = adata.var_names.str.startswith('MALAT1')

    hb_genes = adata.var_names.str.contains('^HB[^(P)]')

    remove = np.add(mito_genes, malat1)

    remove = np.add(remove, hb_genes)

    keep = np.invert(remove)

    adata = adata[:, keep]
    sc.pp.normalize_total(adata, target_sum=1e4)  ##标准化

    print("--------------------------------")
    print(adata)
    print(sum(adata.var["gene_type"]=="lncRNA"))
    print(sum(adata.var["gene_type"]=="protein_coding"))

    adata.write(path.replace("/newh5ad/","/qch5ad/").replace(".h5ad","_qc.h5ad"))
    print("--------------------------------")

path="../data/raw/gse156728/CD8/newh5ad/"
names=os.listdir(path)
for n in names:
    qc(path+n)

