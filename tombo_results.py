import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from nuc_tool import CalcNucPositions
from pathlib import Path
from scipy.ndimage import median_filter

def read_fasta(filename, name=None):
    with open(filename, 'r') as f:
        content = f.readlines()
        for line in range(len(content)):
            if content[line][0] == '>':
                seq = content[line + 1]
                return (seq)


def read_stats(filename, seq, length=0, add_nuc=True):
    blocks = [0, 1]
    result = []
    for block in blocks:
        df = pd.read_hdf(filename, key=f'Statistic_Blocks/Block_{block}/block_stats')
        for id in df['read_id'].unique():
            mod = np.zeros((len(seq)))
            df_id = df[df['read_id'] == id]
            mod[np.asarray(df_id['pos'])] = np.asarray(df_id['stat'])
            if (df_id['pos'].max() - df_id['pos'].min()) > length:
                result.append(mod)

    if add_nuc:
        nuc = CalcNucPositions(seq, mu=-8)
        nuc -= 0.5
        nuc *= 2
        for i in range(25):
            result.append(np.zeros((len(seq))))
        for i in range(100):
            result.append(nuc)

    return np.asarray(result)


filename = r'/media/noort/Data/users/kuijntjes/20220816_barcode10_selection/read_stats.5mC.tombo.per_read_stats'
fasta = r'/media/noort/Data/users/kuijntjes/20220816_barcode08_selection/LinearReference_CP130.fasta'

seq = read_fasta(fasta).upper()

result = read_stats(filename, seq, 3000)
print(len(result))

plt.suptitle(Path(filename))
colorscale = 0.5

plt.imshow(result, cmap='bwr', vmin=-colorscale, vmax=colorscale, origin='lower')
plt.xlim((0, len(seq)))
plt.show()
