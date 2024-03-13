import os

from e_TAU import tau, find

path="../data/raw/gse156728/CD8/qch5ad/"
taupath="../data/raw/gse156728/CD8/qch5ad_tau/"



#注意路径需要qc后的新路径,计算各个tau以及表达最大的类别
names=os.listdir(path)
taucsvs=os.listdir(taupath)
for n in names:
    print(n)
    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
    tau(path+n)
    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")

#看分布
# for n in taucsvs:
#     if(".h5ad" in n ):
#         continue
#     else:
#         print(n)
#         see_distribution(taupath+n,True)

#获取每个类的标记基因，范围太小 目前不用 2024-03
# csv_path="../data/raw/gse156728/CD8/qch5ad_tau/csv/"
# h5ad_path="../data/raw/gse156728/CD8/qch5ad_tau/h5ad/"
# csv_list=os.listdir(csv_path)
# h5ad_list=os.listdir(h5ad_path)
# csv_list.sort()
# h5ad_list.sort()
#
# for c,h in zip(csv_list,h5ad_list):
#     print(c)
#     print(h)
#     print(">>>>>>>>>>>>>>>>>>>>>>")
#     find(csv_path+c,h5ad_path+h,0.8)

