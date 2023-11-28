import scanpy as sc
datapath="../data/raw/gse212890/h5ad/newAll.h5ad"
allData=sc.read_h5ad(datapath)

cd56hcd16l=allData[allData.obs["Majortype"]=="CD56highCD16low"]
cd56lcd16h=allData[allData.obs["Majortype"]=="CD56lowCD16high"]
cd56lcd16h.write(datapath.replace("newAll","cd56lcd16h"))
cd56hcd16l.write(datapath.replace("newAll","cd56hcd16l"))