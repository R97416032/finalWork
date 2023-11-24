import gzip
import os
import tarfile


def un_tar(file_name):
    tar = tarfile.open(file_name)  
    names = tar.getnames()
    file_name=file_name.replace(".tar","")
    if os.path.isdir(file_name):
        pass
    else:  
        os.mkdir(file_name)
    #由于解压后是许多文件，预先建立同名文件夹  
    for name in names:  
        tar.extract(name, file_name + "/")
    tar.close()


def un_gz(file_name):
    """ungz zip file"""
    f_name = file_name.replace(".gz", "").split("/")[-1]
    # 获取文件的名称，去掉
    g_file = gzip.GzipFile(file_name)
    print(f_name)
    if os.path.isdir(file_name.replace(".gz", "").replace(f_name,"ungz")):
        pass
    else:
        os.mkdir(file_name.replace(".gz", "").replace(f_name,"ungz"))
    # 创建gzip对象
    open(file_name.replace(".gz", "").replace(f_name,"ungz/")+f_name, "wb+").write(g_file.read())
    # gzip对象用read()打开后，写入open()建立的文件中。
    g_file.close()
    # 关闭gzip对象
