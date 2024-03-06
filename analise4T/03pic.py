import os

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


def getcellnumpic(path):
    df=pd.read_csv(path)
    df.index=df["names"]

    names=df["names"]
    types=df.columns[1:]

    col_n=[]
    col_t=[]
    col_nums=[]
    for n in names:
        nums=sum(df.loc[n].tolist()[1:])
        print(nums)
        for t in types:
            col_n.append(n)
            col_t.append(t)
            col_nums.append(df[df["names"]==n][t].values[0]/nums)
    data=pd.DataFrame()
    data["names"]=col_n
    data["types"]=col_t
    data["nums"]=col_nums
    data.to_csv(path.split("/")[-1].replace(".csv","_new.csv"),index=None)


# path="../data/raw/gse156728/CD8/nums/"
# names=os.listdir(path)
# names.remove('newclusters.csv')
# print(names)
# for n in names:
#     getcellnumpic(path+n)


# path="../data/raw/gse156728/CD8/qch5ad_tau/csv2/"
# names=os.listdir(path)
# ns=[]
# taus=[]
# types=[]
# for n in names:
#     print(n)
#     df=pd.read_csv(path+n)
#     ns.extend([n.split("_")[1]]*df.shape[0])
#     taus.extend(df["tau_usemean"])
#     types.extend(df["gene_type"])
# df=pd.DataFrame()
#
# df["names"]=ns
# df["type"]=types
# df["tau"]=taus
# df.to_csv("tau.csv",index=None)
#
df=pd.read_csv("./tau.csv")
# g = sns.FacetGrid(df, hue="type",col="names",col_wrap=4)
# g.map(sns.kdeplot, "tau")
# g.add_legend()
# plt.show()
sns.displot(df, x="tau",hue="type",col="names",col_wrap=4,kind='kde')
plt.show()

# path="../data/raw/gse156728/CD8/coexpress/"
# names=os.listdir(path)
# ns=[]
# nums=[]
# nums1=[]
# for n in names:
#     print(n)
#     ns.append(n.split("_")[1])
#     df=pd.read_table(path+n)
#     # nums.append(len(df[df["dis"]>1]))
#     # nums1.append(len(df["dis"]))
#
#     ax=sns.displot(df[df["dis"]>0],x="dis",kind='kde')
#     plt.show()
#     ax=sns.displot(df[df["dis"] > 1], x="dis",kind='kde')
#
#     plt.show()
# df=pd.DataFrame()
# df["names"]=ns
# df["nums"]=nums
# df["nums1"]=nums1
# df["per"]=np.array(nums)/np.array(nums1)
# df.to_csv("conums.csv",index=None)