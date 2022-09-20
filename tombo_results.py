import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
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
    blocks = [1, 0]
    result = []
    for block in blocks:
        df = pd.read_hdf(filename, key=f'Statistic_Blocks/Block_{block}/block_stats')
        for id in df['read_id'].unique():
            mod = np.zeros((len(seq)))
            df_id = df[df['read_id'] == id]
            mod[np.asarray(df_id['pos'])] = np.asarray(df_id['stat'])
            if (df_id['pos'].max() - df_id['pos'].min()) > length:
                result.append(mod)

    return np.asarray(result)

def find_motif(seq, motif):
    return seq.lower().replace(motif.lower(), motif.upper())

def add_seq(data, seq, ref_seq = None, scale=2):
    if not ref_seq:
        nuc = CalcNucPositions(seq, mu=-8)
    else:
        nuc = find_motif(seq, ref_seq)
        nuc = np.asarray([l == l.upper() for l in nuc]).astype(float)

    nuc *= scale

    data = np.concatenate((data, [0*nuc]*30))
    data = np.concatenate((data, [nuc]*80))
    data = np.concatenate((data, [0*nuc]*30))
    return data

def create_plots(filename, fasta, seq601):
    seq = read_fasta(fasta).upper()

    data = read_stats(filename, seq, 3000)
    colorscale = 1.0

    data = add_seq(data, seq, ref_seq=seq601, scale=colorscale)

    plt.figure(figsize=(8, 2))
    plt.text(50, len(data) + 50, filename.split(r'/')[-2])
    ax = plt.gca()
    im = ax.imshow(data, cmap='bwr', vmin=-colorscale, vmax=colorscale, origin='lower')
    plt.xlim((0, len(seq)))
    plt.xlabel('i')
    plt.ylabel(f'LLR {filename.split(".")[-3]}')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="1%", pad=0.1)
    plt.colorbar(im, cax=cax)
    plt.savefig(filename.replace('.tombo.per_read_stats', '.jpg'), dpi=1200)
    plt.show()

if __name__ == '__main__':
    filename = r'/media/noort/Data/users/noort/20220816_barcode10_selection/read_stats.MOD.tombo.per_read_stats'
    fasta = r'/media/noort/Data/users/noort/20220816_barcode08_selection/LinearReference_CP130.fasta'
    seq601 = 'CTGGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGCAAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCAGATATATACATCCTGT'

    for mod in ['6mA', '5mC']:
        create_plots(filename.replace('MOD', mod), fasta, seq601)

