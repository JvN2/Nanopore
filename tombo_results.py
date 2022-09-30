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
    strand = []
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
        found_bps = ''.join([a if np.char.isupper(b) else b for a, b in zip(context_bps, found_bps)])
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


def create_plots(filename, fasta, seq601, label=''):
    seq = read_fasta(fasta).upper()

    llr, strands = read_stats(filename, seq, 3000)

    if '5mC' in filename:
        CpG = find_motif(seq, 'C', 'CG', format='numeric')
        GpC = find_motif(seq, 'G', 'GC', format='numeric')
        for i, strand in enumerate(strands):
            if strand == 0:
                llr[i] *= CpG
            else:
                llr[i] *= GpC

    mean_llr = np.mean(llr, axis=0)
    mean_llr[np.isclose(mean_llr,0)] = np.nan

    llr = np.asarray([moving_average(d, 10) for d in llr])

    colorrange = 1.0

    # llr = add_seq(llr, seq, ref_seq=None, scale=colorrange)
    llr = add_seq(llr, seq, ref_seq=seq601, scale=colorrange)

    plt.figure(figsize=(18, 2))
    title = filename.split(r'/')[-2] + f'\n\n{filename.split(".")[-3]} {label}'
    plt.text(0, 1.0, title, horizontalalignment='left', verticalalignment='bottom', transform=plt.gca().transAxes)

    ax = plt.gca()

    if True:
        plt.plot(-mean_llr, linewidth=0, marker='o', color='red', markerfacecolor='none')
        nucs = find_motif(seq, seq601, format='numeric')
        nucs[-1]=0
        plt.plot(-1.5*(nucs-0.5), color='blue')
#        plt.fill(nucs*3, facecolor='blue', alpha=0.5)
#        plt.fill((1-nucs)*-3, facecolor='red', alpha=0.5)
        plt.ylabel(f'LLR')
    else:
        im = ax.imshow(llr, cmap='bwr', vmin=-colorrange, vmax=colorrange, origin='lower')
        plt.ylabel(f'molecule#')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="1%", pad=0.1)
        cbar = plt.colorbar(im, cax=cax)
        cbar.set_label('LLR', rotation=270)

    plt.xlim((0, len(seq)))
    plt.xlabel('i')

    plt.savefig(filename.replace('.tombo.per_read_stats', '.jpg'), dpi=1200)
    plt.show()


if __name__ == '__main__':
    barcode = 10
    experiment = r'/media/noort/Data/users/noort/20220816_1950_MN30914_AJF795_9344cc69/sample_description.xlsx'
    filename = fr'/media/noort/Data/users/noort/20220816_barcode{barcode:02d}_selection/read_stats.MOD.tombo.per_read_stats'
    fasta = rf'/media/noort/Data/users/noort/20220816_barcode{barcode:02d}_selection/LinearReference_CP130.fasta'
    seq601 = 'CTGGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGCAAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCAGATATATACATCCTGT'
    mods = ['6mA', '5mC']

    df = pd.read_excel(experiment, index_col='Barcode')
    label=df.at[barcode, 'Enzymes']
    print(df)

    for mod in mods:
        create_plots(filename.replace('MOD', mod), fasta, seq601, label)
