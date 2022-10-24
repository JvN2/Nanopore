import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tqdm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from nuc_tool import CalcNucPositions
from pathlib import Path
from scipy.ndimage import median_filter
import h5py, re
from scipy import stats


def read_fasta(filename, name=None):
    with open(filename, 'r') as f:
        content = f.readlines()
        for line in range(len(content)):
            if content[line][0] == '>':
                seq = content[line + 1]
    seq = ''.join([s for s in seq.upper() if s in ['A', 'C', 'G', 'T']])
    return seq


def conjugate_strand(watson, reverse=True):
    conjugate = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    crick = ''.join([conjugate[base] for base in watson.upper()])
    if not reverse:
        return crick
    else:
        return crick[::-1]


def read_model(seq, strand='+'):
    import tombo.tombo_stats
    modelfile = '/home/noort/PycharmProjects/Nanopore/temp/r9.4_180mv_450bps_6mer.model'
    model = tombo.tombo_stats.TomboModel(modelfile, is_text_model=True)
    seq = seq.upper()
    if strand == '-':
        seq = conjugate_strand(seq)
    mean, sd = model.get_exp_levels_from_seq(seq)
    squigle = np.empty(len(seq))
    squigle[:] = np.nan
    if strand == '+':
        squigle[2:2 + len(mean)] = mean
    else:
        squigle[0: len(mean)-6] = mean[::-1][6:]
    return squigle

def cross_correlate(a,b):
    from scipy import signal
    a[~np.isfinite(a)] = 0
    b[~np.isfinite(b)] = 0
    n = len(a)
    cc = signal.correlate(b, a, mode='same') / np.sqrt(
        signal.correlate(a, a, mode='same')[int(n / 2)] * signal.correlate(b, b, mode='same')[int(n / 2)])
    delay_arr = np.linspace(-0.5*n, 0.5*n, n)
    delay = delay_arr[np.argmax(cc)]
    return cc, delay

def read_squiggles(filename, fasta):
    import glob
    from ont_fast5_api.fast5_file import Fast5File

    path = str(Path(filename).parent) + r'/tombo/0/*.fast5'
    files = glob.glob(path)

    seq = read_fasta(fasta)
    index601 = find_motif(seq, seq601, format='index')
    seq = find_motif(seq, seq601, format='string')

    fit = [0.10309092, -9.42492226]
    display_nr = -11

    plus_squigles = []
    min_squigles = []
    min_squigles_quality = []
    for f in tqdm.tqdm(files[:500]):
        squigle = np.empty(len(seq))
        squigle[:] = np.nan
        try:
            df = pd.read_hdf(f, key=f'Analyses/RawGenomeCorrected_000/BaseCalled_template/Events')
            with Fast5File(f, mode="r") as f5:
                read_strand = dict(f5.get_analysis_attributes(r'RawGenomeCorrected_000/BaseCalled_template/Alignment'))[
                    'mapped_strand']
                offset = dict(f5.get_analysis_attributes(r'RawGenomeCorrected_000/BaseCalled_template/Alignment'))[
                    'mapped_start']

            df.index += offset
            squigle[df.index] = df['norm_mean']
            if read_strand == '-':
                squigle = squigle[::-1]

            ref = read_model(seq.upper(), read_strand)
            selection = np.isfinite(ref) * np.isfinite(squigle)
            fit = np.polyfit(ref[selection], squigle[selection], 1)
            if fit[0] > 0.05:
                if read_strand == '+':
                    plus_squigles.append(squigle)
                else:
                    min_squigles.append(squigle)
            if len(plus_squigles) == display_nr:
                ref = read_model(seq.upper(), '+')
                plt.plot(np.polyval(fit, ref), color = 'k')
                plt.plot(squigle, color = 'r')
                plt.title('+ strand')
                plt.show()
            #cc, shift = cross_correlate(squigle, np.polyval(fit, ref))
            #min_squigles_quality.append([fit[0], fit[1]])

        except KeyError:
            pass

    #min_squigles_quality=np.asarray(min_squigles_quality).T
    #plt.scatter(min_squigles_quality[0], min_squigles_quality[1])
    #plt.show()

    xlabels = np.asarray([s + f'\n:\n{i+1}' if (i+1)%10 == 0 else s for i, s in enumerate(seq)])
    x = np.arange(len(seq))
    i = index601[3]
    linkerlength = 50
    selection = np.arange(147 + linkerlength) + i

    for strand, color, squigles in zip(['+', '-'], ['red', 'green'], [plus_squigles, min_squigles]):
        fig = plt.figure(figsize=(25, 3))
        plt.tight_layout()
        plt.xticks(x[selection] - x[selection][0] + 1, xlabels[selection])
        plt.ylim((-5.5, 5.5))
        plt.ylabel(r'nomalized current')
        fig.subplots_adjust(bottom=0.2, top=0.8, left=0.05, right=0.98)
        title = filename.split(r'/')[-2] + f'\n\n{label}'
        plt.text(0, 1.0, title, horizontalalignment='left', verticalalignment='bottom', transform=plt.gca().transAxes)
        squigles = np.asarray(squigles)

        cleaned_squigles = [s[~np.isnan(s)] if len(s[~np.isnan(s)]) > 0 else [-1e6] for s in
                            squigles[:, selection + 1].T]

        violin_parts = plt.violinplot(cleaned_squigles, showextrema=False, widths=0.5)
        for pc in violin_parts['bodies']:
            pc.set_facecolor(color)
            pc.set_edgecolor(color)
        ref = read_model(seq, strand)

        #    fit = np.polyfit( ref, median_squigle, 1)
        #    print(fit)
        fit = [0.10309092, -9.42492226]
        ref = np.polyval(fit, ref)

        plt.text(1, 0.97, f'\nmodel ', horizontalalignment='right', verticalalignment='top',
                 transform=plt.gca().transAxes)
        plt.text(1, 0.97, f'N={len(cleaned_squigles[0])}, {strand} strand ', horizontalalignment='right', verticalalignment='top',
                 transform=plt.gca().transAxes, color=color)

        width = 1.5
        median_cleaned_squigles = np.asarray([np.median(s) for s in cleaned_squigles])
        plt.scatter(x[selection] + 1 - x[selection][0], median_cleaned_squigles, marker="_", color=color,
                    linewidths=width)
        plt.scatter(x[selection] - x[selection][0]+1, ref[selection+1], marker="_", color='k', linewidths=width)

        # Difference plot
        plt.plot(x[selection] + 1 - x[selection][0], 0*median_cleaned_squigles - 5, color = 'k', linewidth=0.5)
        plt.scatter(x[selection] + 1 - x[selection][0], ref[selection+1] - median_cleaned_squigles - 5, marker="_",
                    color=color, linewidths=width)
        plt.text(-1, -5, 'difference',  horizontalalignment='right', verticalalignment='center')

        plt.xlim((0, x[selection][-1] - x[selection][0] + 1.5))
        plt.savefig(filename.replace('MOD.tombo.per_read_stats', f'squigle{strand}.jpg'), dpi=1200)
        # plt.close()

    plt.show()
    return


def read_stats(filename, seq, length=0, range=None, motifs=None):
    print(filename)
    strands = [1, 0]
    result = []
    strand = []
    for block in strands:
        df = pd.read_hdf(filename, key=f'Statistic_Blocks/Block_{block}/block_stats')
        for id in df['read_id'].unique():
            mod = np.zeros((len(seq)))
            df_id = df[df['read_id'] == id]
            mod[np.asarray(df_id['pos'])] = np.asarray(df_id['stat'])
            if range is None:
                if (df_id['pos'].max() - df_id['pos'].min()) > length and df_id['pos'].min():
                    result.append(mod)
                    strand.append(block)
            else:
                if (df_id['pos'].max() >= range[1] and df_id['pos'].min()) <= range[0]:
                    result.append(mod)
                    strand.append(block)
    result = np.asarray(result)

    if motifs is not None:
        filter = np.zeros(len(seq))
        for motif in motifs:
            index = find_motif(seq, motif, format='index')
            offset = [i for i, c in enumerate(motif) if c.isupper()]
            filter[index + offset[0]] = 1
        result = np.multiply(result, filter)
    return result, np.asarray(strand)


def find_motif(seq, motif, context=None, format='string'):
    seq = seq.lower()
    motif = motif.lower()
    if context is not None:
        context = context.lower()
    found_bps = seq.replace(motif, motif.upper())
    if context is not None:
        context_bps = seq.replace(context, context.upper())
        found_bps = ''.join([a if np.char.isupper(b) else b for a, b in zip(context_bps, found_bps)])
    if format == 'string':
        return found_bps
    if format == 'numeric':
        return np.asarray([1.0 if bp == bp.upper() else 0 for bp in found_bps])
    if format == 'index':
        return np.asarray([m.start() for m in re.finditer(motif, seq)])


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


def create_plots(filename, fasta, seq601, label='', type='mean'):
    seq = read_fasta(fasta).upper()

    if '5mC' in filename:
        motifs = ['Cg', 'Gc', 'cG', 'gC']
    elif '6mA' in filename:
        motifs = ['A', 'T']
    else:
        motifs = None

    llr, strands = read_stats(filename, seq, length=5800, motifs=motifs)
    llr = sort_traces(llr, strands)

    mean_llr = np.mean(llr, axis=0)
    mean_llr[np.isclose(mean_llr, 0)] = np.nan

    # llr = np.asarray([moving_average(d, 10) for d in llr])

    colorrange = 1.0

    # llr = add_seq(llr, seq, ref_seq=None, scale=colorrange)
    # llr = add_seq(llr, seq, ref_seq=seq601, scale=colorrange)

    plt.figure(figsize=(25, 2))
    title = filename.split(r'/')[-2] + f'\n\n{filename.split(".")[-3]} {label}'
    plt.text(0, 1.0, title, horizontalalignment='left', verticalalignment='bottom', transform=plt.gca().transAxes)

    ax = plt.gca()

    if type == 'mean':
        plt.plot(-mean_llr, linewidth=0, marker='o', color='red', markerfacecolor='none')
        nucs = find_motif(seq, seq601, format='numeric')
        nucs[-1] = 0
        plt.plot(-1.5 * (nucs - 0.5), color='blue')
        plt.ylabel(f'LLR')
        plt.xlim((0, len(seq)))
        plt.xlabel('i')
    elif type == 'motif':
        import matplotlib.patches as mpatches

        index601 = find_motif(seq, seq601, format='index')
        seq = find_motif(seq, seq601, format='string')
        xlabels = [s + f'\n:\n{i}' if i % 10 == 0 else s for i, s in enumerate(seq)]

        i = index601[3]
        linkerlength = 25
        plotrange = [-linkerlength, 147 + linkerlength]
        xlabels = xlabels[i + plotrange[0]: i + plotrange[1]]

        labels = []
        for s, color in enumerate(['red', 'blue']):
            subllr = llr[strands == s, i + plotrange[0]: i + plotrange[1]]
            violin_parts = plt.violinplot(-subllr, showextrema=False)
            for pc in violin_parts['bodies']:
                pc.set_facecolor(color)
                pc.set_edgecolor('black')

            labels.append((mpatches.Patch(color=color), f'strand {s}'))

        plt.legend(*zip(*labels), loc=2)
        plt.xticks(np.arange(len(xlabels)) + 1, xlabels)
        plt.ylim((-10, 10))
        plt.ylabel(f'-LLR')
        plt.tight_layout()
    else:
        im = ax.imshow(-llr, cmap='bwr', vmin=-colorrange, vmax=colorrange, origin='lower')
        plt.ylabel(f'molecule#')
        plt.xlim((0, len(seq)))
        plt.xlabel('i')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="1%", pad=0.1)
        cbar = plt.colorbar(im, cax=cax)
        cbar.set_label('LLR', rotation=270)

    plt.savefig(filename.replace('.tombo.per_read_stats', '.jpg'), dpi=1200)
    plt.show()


def show_squigle(filename, fasta):
    import hdf5plugin
    import tombo
    def printall(name, obj):
        print(name, dict(obj.attrs))

    with h5py.File(filename, 'r') as hf:
        read = list(hf['Raw']['Reads'].keys())[0]
        print(read)
        data = hf['Raw']['Reads'][read]['Signal']

        print(data[()])
        return

        # print(np.asarray(hf[reads[0]]['Analyses/Basecall_1D_000/BaseCalled_template/ModBaseProbs']))

        plt.plot((np.asarray(data)))
        plt.show()


def read_tombo(fast5_fn, reference_fn):
    import h5py, mappy
    from tombo import tombo_helper, tombo_stats, resquiggle
    import hdf5plugin
    # set read values

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

    # run full re-squiggle
    rsqgl_results = resquiggle.resquiggle_read(
        map_results, std_ref, rsqgl_params, all_raw_signal=all_raw_signal)


def sort_traces(traces, strands, method='sum'):
    from sklearn import decomposition
    new_traces = np.empty_like(traces)
    for strand in np.unique(strands):
        selection = traces[strands == strand]
        if method == 'pca':
            pca = decomposition.PCA(n_components=1)
            pca.fit(traces[strands == strand])
            data = pca.transform(selection)
            i = np.argsort(data[:, 0])
        else:
            i = np.argsort(np.sum(selection, axis=1))
        new_traces[strands == strand] = selection[i]
    return new_traces


if __name__ == '__main__':
    barcode = 8
    experiment = r'/media/noort/Data/users/noort/20220816_1950_MN30914_AJF795_9344cc69/sample_description.xlsx'
    tombo_filename = fr'/media/noort/Data/users/noort/20220816_barcode{barcode:02d}_selection/read_stats.MOD.tombo.per_read_stats'
    fasta = rf'/media/noort/Data/users/noort/20220816_barcode{barcode:02d}_selection/LinearReference_CP130.fasta'
    seq601 = 'CTGGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGCAAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCAGATATATACATCCTGT'
    mods = ['6mA', '5mC']
    fast5_filename = r'/media/noort/Data/users/noort/20220816_barcode08_selection/tombo/0/0bc22d60-5c34-48c6-9cfd-163aa82be630.fast5'

    df = pd.read_excel(experiment, index_col='Barcode')
    label = df.at[barcode, 'Enzymes']
    print(df)

    for mod in mods[:1]:
        #        create_plots(tombo_filename.replace('MOD', mod), fasta, seq601, label, type='motif')
        read_squiggles(tombo_filename, fasta)
