import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from nuc_tool import CalcNucPositions
from pathlib import Path
from scipy.ndimage import median_filter
import h5py, re, os
from scipy import stats
import re
import tqdm
import glob
from tombo.tombo_stats import TomboModel


def conjugate_strand(watson, reverse=True):
    conjugate = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    crick = ''.join([conjugate[base] for base in watson.upper()])
    if not reverse:
        return crick
    else:
        return crick[::-1]


def read_fasta(filename):
    with open(filename, 'r') as f:
        content = f.readlines()
        for line in range(len(content)):
            if content[line][0] == '>':
                seq = content[line + 1]
    seq = ''.join([s for s in seq.upper() if s in ['A', 'C', 'G', 'T']])
    return seq


def read_squigle(fast5_filename, df=None, ref_fasta_filename=None):
    if df is None:
        if ref_fasta_filename is None:
            ref_fasta_filename = glob.glob(str(Path(fast5_filename).parent) + '/*.fasta')[0]
        df = pd.DataFrame(list(read_fasta(ref_fasta_filename)), columns=['base'])
    try:
        events = pd.read_hdf(fast5_filename, key=f'Analyses/RawGenomeCorrected_000/BaseCalled_template/Events')
        with h5py.File(fast5_filename, 'r') as hdf:
            read_id = list(hdf['Raw']['Reads'].keys())[0]
            attributes = dict(hdf['Analyses/RawGenomeCorrected_000/BaseCalled_template/Alignment'].attrs.items())
        events.index += attributes['mapped_start']
        if attributes["mapped_strand"] == '-':
            events.index = events.index[::-1]
        df[f'{attributes["mapped_strand"]}{read_id}'] = events['norm_mean']
        return df
    except KeyError:
        print(f'Squigle not found in: {fast5_filename}')
        return df


def compute_squigle(seq, strand='+', data=None):
    model_filename = os.getcwd() + r'/Data/r9.4_180mv_450bps_6mer.model'
    model = TomboModel(model_filename, is_text_model=True)
    seq = seq.upper()
    if strand == '-':
        seq = conjugate_strand(seq)
    mean, sd = model.get_exp_levels_from_seq(seq)

    if strand == '-':
        mean = mean[::-1]
        sd = sd[::-1]
        offset = 3
    else:
        offset = 2

    df = pd.DataFrame(mean, columns=['mean'])
    df['sd'] = sd
    df.index += offset
    df['normalized'] = (mean - np.mean(mean)) / np.std(mean)
    if data is not None:
        df['data'] = data
        fit = np.polyfit(df.dropna()['mean'], df.dropna()['data'], 1)
        df['fitted'] = np.polyval(fit, mean)
    return df


if __name__ == '__main__':
    data_dir = r'/mnt/c/tmp/depc'
    filenames = glob.glob(data_dir + r'/*.fast5')
    df = None
    for f in filenames:
        df = read_squigle(f, df=df)
    df.sort_index(axis=1, inplace=True)
    seq = ''.join(df['base'])

    data = []
    for read in df.columns:
        if read[0] in ['+', '-']:
            model = compute_squigle(seq, strand=read[0], data=df[read])
            data.append(df[read] - model['fitted'])
            data.append(df[read])
            plt.figure(figsize=(25, 2))
            plt.plot(df[read], 'r')
            plt.plot(model['fitted'][df[read].notna()], 'k')
            plt.ylim((-4, 4))
            plt.xlim((0, model.index[-1]))
            plt.title(read)
            plt.ylabel('normalized current')
            plt.tight_layout()
            plt.savefig(data_dir + f'/squigle{read}.jpg', dpi=1200)
            plt.close()


    def bins_to_centers(bins):
        return (bins[1:] + bins[:-1]) / 2

    # print(np.shape(data))
    # plt.imshow(data, cmap='bwr', aspect=100, vmin=-1, vmax=1)
    # plt.plot(np.asarray(data).T)
    # plt.ylim((-3,3))
    # plt.plot(model['fitted'], color = 'black')

    for d in data:
        his, bins = np.histogram(d[~np.isnan(d)], bins=np.linspace(-3, 3, 500))
        his = his / np.sum(his)
        plt.plot(bins_to_centers(bins), his)
        # plt.yscale("log")
    plt.show()
