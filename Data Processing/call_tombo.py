import h5py, mappy, glob
import numpy as np
from tombo import tombo_helper, tombo_stats, resquiggle
import pandas as pd
import matplotlib.pyplot as plt

def test_resquiggle(fast5_fn):
    fast5_data = h5py.File(fast5_fn, 'r')
    seq_samp_type = tombo_helper.get_seq_sample_type(fast5_data)

    # prep aligner, signal model and parameters
    aligner = mappy.Aligner(reference_fn, preset=str('map-ont'), best_n=1)
    std_ref = tombo_stats.TomboModel(seq_samp_type=seq_samp_type)
    rsqgl_params = tombo_stats.load_resquiggle_parameters(seq_samp_type)

    # extract data from FAST5
    map_results = resquiggle.map_read(fast5_data, aligner, std_ref)
    all_raw_signal = tombo_helper.get_raw_read_slot(fast5_data)['Signal'][:]
    if seq_samp_type.rev_sig:
        all_raw_signal = all_raw_signal[::-1]
    map_results = map_results._replace(raw_signal=all_raw_signal)
    print(map_results.genome_seq)


def read_tombo_mods(folder):
    filenames = glob.glob(folder + r'/*.per_read_stats')
    all_mods = pd.DataFrame()
    for filename in filenames:
        with h5py.File(filename, 'r') as hf:
            blocks = [str(block) for block in hf[f'Statistic_Blocks']]
        for block in blocks:
            mods = pd.read_hdf(filename, f'Statistic_Blocks/{block}/block_stats')
            with h5py.File(filename, 'r') as hf:
                read_pars = dict(hf[f'Statistic_Blocks/{block}'].attrs)
                index = [read.decode() for read in hf[f'Statistic_Blocks/{block}/read_ids']]
            mods['read_id'] = np.asarray(index)[mods['read_id']]
            mods['strand'] = read_pars['strand']
            mods['chrm'] = read_pars['chrm']
            all_mods = all_mods.append(mods, ignore_index=True)

    all_mods = pd.pivot_table(all_mods,
                              values=['chrm', 'strand', 'pos', 'stat'],
                              index='read_id',
                              aggfunc={'strand': 'first', 'chrm': 'first', 'pos': list, 'stat': list})
    length = [(np.max(row['pos']) - np.min(row['pos'])) for i, row in all_mods.iterrows()]
    all_mods['length'] = length

    return all_mods


def print_mods(df, read_nr=0):
    chrm = df.iloc[read_nr]['chrm']
    stat = df.iloc[read_nr]['stat']
    pos = np.asarray(df.iloc[read_nr]['pos']).astype(int)
    pos = pos[np.asarray(stat) > 0]
    seq = list(mappy.Aligner(reference_fn).seq(chrm).lower())
    seq = [s.upper() if i in pos else s for i, s in enumerate(seq)]
    seq = ''.join(seq)
    print(seq)


def plot_mods(df, reference_fn, chrms=None, min_length=0):
    aligner = mappy.Aligner(reference_fn)
    if chrms is None:
        chrms = aligner.seq_names

    for chrm in chrms:
        df1 = df[df['chrm'] == chrm]
        df1 = df1.sort_values('strand')
        df1 = df1[df1['length'] > min_length]

        if len(df1) > 0:
            im = np.zeros((len(df1), len(aligner.seq(chrm))))
            for i, _ in enumerate(im):
                im[i, df1.iloc[i]['pos']] = df1.iloc[i]['stat']
            plt.imshow(im, vmin=-4, vmax=4, origin='lower', cmap=plt.get_cmap('seismic'),
                       aspect='auto')  # , interpolation='none') #
            plt.title(chrm)
            plt.xlabel('i (bp)')
            plt.ylabel('read nr')
            plt.show()
        else:
            print('No reads found that match criteria')


# set file locations
fast5_fn = rf'/media/noort/Data/users/noort/test20/tombo/0/00a2b3c7-6add-44f2-b171-3bd51e1ec637.fast5'
reference_fn = '/media/noort/Data/users/noort/ref_files/combined.fa'
tombo_folder = '/media/noort/Data/users/noort/test20'

df = read_tombo_mods(tombo_folder)
chrms = df['chrm'].unique()
plot_mods(df, reference_fn, chrms=chrms, min_length = 4000)
# print_mods(df, 20)
