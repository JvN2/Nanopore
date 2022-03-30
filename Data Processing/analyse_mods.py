import glob

import h5py
import mappy
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

def read_guppy_mods(folder, reference_fn):
    aligner = mappy.Aligner(reference_fn)
    print(aligner.seq_names)
    filenames = glob.glob(folder + r'/*.fast5')
    all_mods = pd.DataFrame(columns=['chrm', 'strand', 'pos', 'stat'])
    for filename in filenames:
        with h5py.File(filename, 'r') as hf:
            for read_id in hf.keys():
                # print(read_id)
                read_mods = np.asarray(hf[read_id]['Analyses/Basecall_1D_000/BaseCalled_template/ModBaseProbs'])
                # print(read_mods)
                seq = np.asarray(hf[read_id]['Analyses/Basecall_1D_000/BaseCalled_template/Fastq']).astype(str).item()
                seq = seq.split('\n')[1]
                print(seq)
                # all_mods['read_id']=['chrm', '+/-', read_mods[0], read_mods[2]/read_mods[1], read_mods[4]/read_mods[3]]
                # print(all_mods)
                try:
                    # enumerate all alignments to avoid mappy memory leak
                    alignment = list(aligner.map(seq))[0]
                except IndexError:
                    raise print('Alignment not produced')

                chrm = alignment.ctg
                print(chrm)
                # subtract one to put into 0-based index
                ref_start = alignment.r_st
                ref_end = alignment.r_en
                strand = '+' if alignment.strand == 1 else '-'
                num_match = alignment.mlen
                num_ins, num_del, num_aligned = 0, 0, 0
                for op_len, op in alignment.cigar:
                    if op == 1:
                        num_ins += op_len
                    elif op in (2, 3):
                        num_del += op_len
                    elif op in (0, 7, 8):
                        num_aligned += op_len
                    elif op == 6:
                        pass
                    else:
                        # soft and hard clipping are not reported in the
                        # mappy cigar
                        raise th.TomboError('Invalid cigar operation')

                # store number of clipped bases relative to read sequence
                if strand == '+':
                    num_start_clipped_bases = alignment.q_st
                    num_end_clipped_bases = len(seq) - alignment.q_en
                else:
                    num_start_clipped_bases = len(seq) - alignment.q_en
                    num_end_clipped_bases = alignment.q_st


                # extract genome sequence from mappy aligner
                # expand sequence to get model levels for all sites (need to handle new
                # sequence coordinates downstream)
                dnstrm_bases = std_ref.kmer_width - std_ref.central_pos - 1
                if ((seq_samp_type.name == RNA_SAMP_TYPE and strand == '+') or
                        (seq_samp_type.name == DNA_SAMP_TYPE and strand == '-' and
                         USE_START_CLIP_BASES) or
                        (seq_samp_type.name == DNA_SAMP_TYPE and strand == '+' and
                         not USE_START_CLIP_BASES)):
                    if ref_start < std_ref.central_pos:
                        ref_start = std_ref.central_pos
                    ref_seq_start = ref_start - std_ref.central_pos
                    ref_seq_end = ref_end + dnstrm_bases
                else:
                    if ref_start < dnstrm_bases:
                        ref_start = dnstrm_bases
                    ref_seq_start = ref_start - dnstrm_bases
                    ref_seq_end = ref_end + std_ref.central_pos
                genome_seq = aligner.seq(chrm, ref_seq_start, ref_seq_end)
                if genome_seq is None or genome_seq == '':
                    raise th.TomboReads('Invalid mapping location')

                if sys.version_info[0] < 3:
                    genome_seq = genome_seq.decode()
                if strand == '-':
                    genome_seq = th.rev_comp(genome_seq)
                # discordant mapping to sequence extraction is due to reads mapping up to
                # the end of a seqeunce record (and don't need to carry around record lens),
                # so don't error on these discordant lengths here
                # if len(genome_seq) != ref_end - ref_start + std_ref.kmer_width - 1:
                #    raise th.TomboError('Discordant mapped position and sequence')
                genome_loc = th.genomeLocation(ref_start, strand, chrm)

                # store sequence at the end of the read without an adapter
                # for simpler read start identification (start of RNA genomic sequence
                # end of DNA genomic sequence)
                start_clip_bases = None
                if USE_START_CLIP_BASES:
                    start_clip_bases = seq_data.seq[alignment.q_en:][::-1]

                return th.resquiggleResults(
                    align_info=align_info, genome_loc=genome_loc, genome_seq=genome_seq,
                    mean_q_score=seq_data.mean_q_score, start_clip_bases=start_clip_bases)
                break


def read_megalodon_mods(filename):
    all_mods = pd.read_csv(filename, sep = '\t')
    all_mods['LLR'] = all_mods['mod_log_prob']-all_mods['can_log_prob']
    all_mods = pd.pivot_table(all_mods,
                                  values=['chrm', 'strand', 'pos', 'LLR'],
                                  index='read_id',
                                  aggfunc={'strand': 'first', 'chrm': 'first', 'pos': list, 'LLR': list})
    all_mods['length'] = [(np.max(row['pos']) - np.min(row['pos'])) for i, row in all_mods.iterrows()]
    return all_mods

def read_tombo_events(filename):
    events = pd.read_hdf(filename,'Analyses/RawGenomeCorrected_000/BaseCalled_template/Events')
    with h5py.File(filename, 'r') as hf:
        read_nr = list(hf['Raw/Reads'].keys())[0]
        squiggle = np.asarray(hf[f'Raw/Reads/{read_nr}/Signal'])
    return events, squiggle


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
    all_mods['length'] = [(np.max(row['pos']) - np.min(row['pos'])) for i, row in all_mods.iterrows()]
    return all_mods

def print_mods(df, reference_fn, read_nr=0):
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
                im[i, df1.iloc[i]['pos']] = df1.iloc[i]['LLR']
            plt.imshow(im, vmin=-4, vmax=4, origin='lower', cmap=plt.get_cmap('seismic'),
                       aspect='auto')  # , interpolation='none') #
            plt.title(chrm)
            plt.xlabel('i (bp)')
            plt.ylabel('read nr')
            plt.show()
        else:
            print('No reads found that match criteria')

megalodon_filename = '/home/noort/data/Analysed_2020-05-12/megalodon/per_read_modified_base_calls.txt'
tombo_folder = '/home/noort/data/Analysed_2022-01-18_12samplemethylationtest/tombo'
tombo_filename = '/media/noort/Data/users/noort/Analysed_2020-05-12/tombo/0/0a3fcfc7-628b-4397-aa6f-4a3e262407f0.fast5'
guppy_folder = '/home/noort/data/Analysed_2022-01-18_12samplemethylationtest/guppy/workspace/barcode02'
reference_filename = '/media/noort/Data/users/noort/ref_files/combined.fa'

# read_guppy_mods(guppy_folder, reference_filename)
all_mods = read_megalodon_mods(megalodon_filename)
plot_mods(all_mods, reference_filename, min_length=3000)

df = read_tombo_mods(tombo_folder,)
# events, squiggle = read_tombo_events(tombo_filename)
# print(events)
# squiggle_model = []
# norm_mean = events['norm_mean'].values
# length = events['length'].values
# squiggle_model= []
# for nm, l in zip(norm_mean, length):
#     squiggle_model.append(nm*np.ones(l))
# squiggle_model = np.asarray(squiggle_model).reshape((1,-1))
#
# print(squiggle_model)
# plt.plot(squiggle_model)
# plt.show()
