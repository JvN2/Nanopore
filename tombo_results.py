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
    strand =[]
    for block in blocks:
        df = pd.read_hdf(filename, key=f'Statistic_Blocks/Block_{block}/block_stats')
        for id in df['read_id'].unique():
            mod = np.zeros((len(seq)))
            df_id = df[df['read_id'] == id]
            mod[np.asarray(df_id['pos'])] = np.asarray(df_id['stat'])
            if (df_id['pos'].max() - df_id['pos'].min()) > length:
                result.append(mod)
                strand.append(block)

    return np.asarray(result), np.asarray(strand)


def find_motif(seq, motif, context=None, format='string'):
    found_bps = seq.lower().replace(motif.lower(), motif.upper())
    if context is not None:
        context_bps = seq.lower().replace(context.lower(), context.upper())
        found_bps = ''.join([a if np.char.isupper(b) else b for a,b in zip(context_bps, found_bps)])
    if format == 'string':
        return found_bps
    if format == 'numeric':
        return np.asarray([1.0 if bp == bp.upper() else 0 for bp in found_bps])

def add_seq(data, seq, ref_seq=None, scale=2):
    if not ref_seq:
        nuc = np.roll(CalcNucPositions(seq, mu=-9), -54)
    else:
        nuc = find_motif(seq, ref_seq)
        nuc = np.asarray([l == l.upper() for l in nuc]).astype(float)

    nuc = 2 * (nuc - 0.5)
    nuc *= scale

    data = np.concatenate((data, [0 * nuc] * 30))
    data = np.concatenate((data, [nuc] * 80))
    data = np.concatenate((data, [0 * nuc] * 30))
    return data


def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'same') / w


def create_plots(filename, fasta, seq601):
    seq = read_fasta(fasta).upper()

    llr, strands = read_stats(filename, seq, 3000)

    CpG = find_motif(seq, 'C', 'CG', format='numeric')
    GpC = find_motif(seq, 'G', 'GC', format='numeric')
    for i, strand in enumerate(strands):
        if strand == 0:
            llr[i] *= CpG
        else:
            llr[i] *= GpC


    llr = np.asarray([moving_average(d, 10) for d in llr])

    colorrange = 1.0

    #llr = add_seq(llr, seq, ref_seq=None, scale=colorrange)
    llr = add_seq(llr, seq, ref_seq=seq601, scale=colorrange)

    plt.figure(figsize=(8, 2))
    title = filename.split(r'/')[-2] + f'\n\n{filename.split(".")[-3]}'
    #    plt.text(50, len(data) + 50, title)

    plt.text(0, 1.0, title, horizontalalignment='left', verticalalignment='bottom', transform=plt.gca().transAxes)

    ax = plt.gca()
    im = ax.imshow(llr, cmap='bwr', vmin=-colorrange, vmax=colorrange, origin='lower')
    plt.xlim((0, len(seq)))
    plt.xlabel('i')
    plt.ylabel(f'molecule#')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="1%", pad=0.1)
    cbar = plt.colorbar(im, cax=cax)
    cbar.set_label('LLR', rotation=270)
    plt.savefig(filename.replace('.tombo.per_read_stats', '.jpg'), dpi=1200)
    plt.show()


if __name__ == '__main__':
    filename = r'/media/noort/Data/users/noort/20220816_barcode10_selection/read_stats.MOD.tombo.per_read_stats'
    fasta = r'/media/noort/Data/users/noort/20220816_barcode08_selection/LinearReference_CP130.fasta'
    seq601 = 'CTGGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGCAAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCAGATATATACATCCTGT'
    mods = ['6mA', '5mC']


    for mod in mods:
        create_plots(filename.replace('MOD', mod), fasta, seq601)
