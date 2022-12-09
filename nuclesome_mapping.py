import EditChromatin
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from sklearn import decomposition


def test():
    seq = ''.join(np.random.choice(['A', 'C', 'G', 'T'], 1000))
    print(seq)
    w = 147
    mu = -19

    nucleosome_occupancy = EditChromatin.calc_nucleosome_positions(seq, w, mu, flank_size=0)

    fig = plt.figure(figsize=(13, 2))
    plt.tight_layout()
    fig.subplots_adjust(bottom=0.2, top=0.8, left=0.05, right=0.98)
    plt.ylim((0, 1))
    plt.xlim((0, len(seq) - 1))
    plt.ylabel(r'Occupancy')
    plt.plot(nucleosome_occupancy)
    plt.show()


def read_stats(filename, seq, length=0, range=None, motifs=None):
    strands = ['-', '+']
    result = []
    strand = []
    names=[]
    for block, s in enumerate(strands):
        df = pd.read_hdf(filename, key=f'Statistic_Blocks/Block_{block}/block_stats')
        for id in tqdm(df['read_id'].unique(), desc=f'Reading strand: {block}'):
            mod = np.zeros((len(seq)))
            df_id = df[df['read_id'] == id]
            mod[np.asarray(df_id['pos'])] = np.asarray(df_id['stat'])
            names.append(f'{s}{id}')
            if range is None:
                if (df_id['pos'].max() - df_id['pos'].min()) > length and df_id['pos'].min():
                    result.append(mod)
                    strand.append(block)
            else:
                if (df_id['pos'].max() >= range[1] and df_id['pos'].min()) <= range[0]:
                    result.append(mod)
                    strand.append(block)
    result = np.asarray(result)

    # if motifs is not None:
    #     filter = np.zeros(len(seq))
    #     for motif in motifs:
    #         index = find_motif(seq, motif, format='index')
    #         offset = [i for i, c in enumerate(motif) if c.isupper()]
    #         filter[index + offset[0]] = 1
    #     result = np.multiply(result, filter)
    df_out = pd.DataFrame(result.T, columns=names)
    df_out.index.name = 'i'
    df_out = df_out.reindex(sorted(names, key=int), axis=1)
    return df_out


def read_fasta(filename, name=None):
    with open(filename, 'r') as f:
        content = f.readlines()
        for line in range(len(content)):
            if content[line][0] == '>':
                seq = content[line + 1]
    seq = ''.join([s for s in seq.upper() if s in ['A', 'C', 'G', 'T']])
    return seq

def sort_traces(df, n_components = 0):
    traces = df.values.T
    if n_components:
        pca = decomposition.PCA(n_components=n_components)
        pca.fit(traces)
        data = pca.transform(traces)
        index = np.argsort(data[:, 0])
    else:
        index = np.argsort(np.sum(traces, axis=1))
    new_columns = [df.columns[i] for i in index]
    df = df.reindex(new_columns, axis=1)
    return df

def test2():
    filename = r'/media/noort/Data/users/noort/20220816_barcode10_selection/read_stats.5mC.tombo.per_read_stats'
    fast_filename = r'/media/noort/Data/users/noort/20220816_barcode10_selection/LinearReference_CP130.fasta'

    seq = read_fasta(fast_filename).upper()
    tmp = EditChromatin.calc_nucleosome_positions(seq, 80, -7)

    llr, strands = read_stats(filename, seq)
    colorrange = 2
    # plt.imshow(llr, cmap='bwr', vmin=-colorrange, vmax=colorrange, origin='lower')
    # plt.plot(llr[700])
    # plt.show()

    # df2 = pd.read_csv(tombo_filename)
    # print(df2)
    llr /= 2
    mu = 2

    if False:
        i = 700
        #for mu in np.linspace(2, 6, 4):
        for i in range(5):
            i+=705
            energy = EditChromatin.moving_average(llr[i], 5)
            dyad_probability = EditChromatin.vanderlick(energy, mu)
            nucleosome_occupancy = np.convolve(dyad_probability, np.ones(147), mode='same')
            plt.plot(nucleosome_occupancy, label=f'{i}')
        #    plt.plot(tmp, color = 'k', linewidth=2, label='model')


        plt.fill_between(range(0, len(tmp)), tmp, tmp * 0, alpha=0.2, color='k', linewidth=0.4, label='model')
        plt.ylim((0, 1.5))
        plt.xlim((0, len(tmp)))
        plt.legend(loc="upper left")
        plt.show()

    occupancy=[]
    for l in tqdm(llr, desc='mapping nucs '):
        energy = EditChromatin.moving_average(l, 5)
        dyad_probability = EditChromatin.vanderlick(energy, mu)
        nucleosome_occupancy = np.convolve(dyad_probability, np.ones(147), mode='same')
        occupancy.append(nucleosome_occupancy)

    occupancy = sort_traces(occupancy, strands, 'pca')
    ax = plt.gca()
    im = ax.imshow(occupancy, cmap='bwr', vmin=0.5, vmax=1, origin='lower')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="1%", pad=0.1)
    cbar = plt.colorbar(im, cax=cax)
    plt.show()

def plot_traces(llr, occupancy, model):
    light_curve = model

    fig = plt.figure(1, figsize=(10, 4))

    ax1 = plt.subplot(2, 1, 1)
    ax1.imshow(llr, cmap='bwr', origin='lower', vmin = -2, vmax = 2)
    ax11 = ax1.twinx()
    ax11.plot(light_curve, color = 'black')
    ax11.set_ylim((0,1))
    ax11.set_xlim(ax1.get_xlim())
    ax1.set_aspect( "auto")


    ax2 = plt.subplot(2, 1, 2)
    ax2.imshow(occupancy, cmap='bwr', origin='lower', vmin = 0.5, vmax = 0.75)
    ax21 = ax2.twinx()
    ax21.plot(light_curve, color = 'black')
    ax21.set_ylim((0,1))
    ax21.set_xlim(ax1.get_xlim())
    ax2.set_aspect( "auto")

    plt.show()

if __name__ == '__main__':
    filename = r'/media/noort/Data/users/noort/20220816_barcode10_selection/read_stats.5mC.tombo.per_read_stats'
    fast_filename = r'/media/noort/Data/users/noort/20220816_barcode10_selection/LinearReference_CP130.fasta'
    seq = read_fasta(fast_filename).upper()
    model = EditChromatin.calc_nucleosome_positions(seq, 80, -7)

    df = read_stats(filename, seq)
    df = df.filter(regex=r'\+')

    df = sort_traces(df, n_components=4)


    # ax = plt.gca()
    #
    #
    # im = ax.imshow(df.values.T, cmap='bwr', vmin=-2, vmax=2, origin='lower')
    # divider = make_axes_locatable(ax)
    # cax = divider.append_axes("right", size="1%", pad=0.1)
    # cbar = plt.colorbar(im, cax=cax)
    #
    # ax2 = ax.twinx()
    # ax2.set_aspect("equal")
    # ax2.plot(tmp, color = 'black')
    #
    # plt.show()





    mu = 2
    occupancy=[]
    flanking_bps = 500
    for llr in tqdm(df.values.T, desc='mapping nucs '):
        energy = EditChromatin.moving_average(llr/2, 5)
        energy = np.concatenate((np.zeros(flanking_bps), energy, np.zeros(flanking_bps)))
        dyad_probability = EditChromatin.vanderlick(energy, mu)
        nucleosome_occupancy = np.convolve(dyad_probability, np.ones(147), mode='same')
        occupancy.append(nucleosome_occupancy[flanking_bps:-flanking_bps])

    #occupancy = sort_traces(occupancy)
    plot_traces(df.values.T, occupancy, model)
    breakpoint()

    ax = plt.gca()
    im = ax.imshow(occupancy, cmap='bwr', vmin=0.5, vmax=0.75, origin='lower')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="1%", pad=0.1)
    cbar = plt.colorbar(im, cax=cax)
    plt.show()