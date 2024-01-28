import os

from utils.decompress import un_gz

# un_gz("raw/gse212890/GSE212890_NK.genes.csv.gz")
# un_gz("raw/gse212890/GSE212890_NK.barcodes.csv.gz")
# un_gz("raw/gse212890/GSE212890_NK_counts.mtx.gz")
# un_gz("raw/gse212890/GSE212890_NK_metadata.csv.gz")

path="raw/gse156728/metadata/"

def ugz(path):
    paths=os.listdir(path)
    for p in paths:

        un_gz(path+"/"+p)

ugz(path)