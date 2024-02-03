import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

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


path="../data/raw/gse156728/CD8/nums/"
names=os.listdir(path)
names.remove('newclusters.csv')
print(names)
for n in names:
    getcellnumpic(path+n)