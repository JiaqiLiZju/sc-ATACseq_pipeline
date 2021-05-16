# Python (both python2 and python3)
# pysam
# h5py
# numpy
# pybedtools
# snaptools
# deeptools
# macs2

# install using pip
pip install snaptools

# Install snaptools from source code
git clone https://github.com/r3fang/snaptools.git
cd snaptools
pip install -e .
./bin/snaptools

# install deeptools
pip install deeptools

# Install deeptools from source code
git clone https://github.com/deeptools/deepTools.git
wget https://github.com/deeptools/deepTools/archive/1.5.12.tar.gz
tar -xzvf
python setup.py install --prefix /User/Tools/deepTools2.0

# install macs2
pip install MACS2

# Install macs2 from source code
python setup.py install --prefix=$HOME
