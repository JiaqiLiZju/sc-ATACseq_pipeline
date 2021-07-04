import h5py
import pandas as pd
from sys import argv

h5file = argv[1]
# <KeysViewHDF5 ['AM', 'BD', 'FM', 'HD']>
# <KeysViewHDF5 ['CM', 'PE', 'PL', 'PP', 'SA', 'SE', 'TN', 'UM', 'UQ', 'US', 'name']>

fh = h5py.File(h5file, 'r')
df = pd.DataFrame({  
            "CM":fh["BD/CM"][:], 
            "PE":fh["BD/PE"][:],
            "PL":fh["BD/PL"][:],
            "PP":fh["BD/PP"][:],
            "SA":fh["BD/SA"][:],
            "SE":fh["BD/SE"][:],
            "TN":fh["BD/TN"][:],
            "UM":fh["BD/UM"][:],
            "UQ":fh["BD/UQ"][:],
            "US":fh["BD/US"][:],
            },
            index=map(lambda x: x.decode(), fh["BD/name"][:]))

df.to_csv(h5file.replace("snap", 'metadata.csv'))

fh.close()
