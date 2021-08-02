# bedtools
mv bedtools.static.binary bedtools
chmod a+x bedtools

# macs2
wget https://files.pythonhosted.org/packages/e2/61/85d30ecdd34525113e28cb0c5a9f393f93578165f8d848be5925c0208dfb/MACS2-2.2.7.1.tar.gz
tar -xzvf
python setup.py install

# deeptools
wget https://files.pythonhosted.org/packages/8c/9b/92fbfc413548a5d690061ef533c7521af71ebce733ce50e4b6eea3d69d4a/deepTools-3.5.1-py3-none-any.whl

# go ATACdemultiplex
go get -v -u gitlab.com/Grouumf/ATACdemultiplex/...