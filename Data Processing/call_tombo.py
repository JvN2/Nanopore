import h5py, mappy, glob, sys
import numpy as np
from tombo import tombo_stats, resquiggle, _plot_commands
import pandas as pd
import matplotlib.pyplot as plt

# import tombo modules/functions
from tombo import tombo_stats as ts
from tombo import tombo_helper as th

from tombo._default_parameters import (
    EXTRA_SIG_FACTOR, MASK_FILL_Z_SCORE,
    MASK_BASES, DEL_FIX_WINDOW, MAX_DEL_FIX_WINDOW,
    MIN_EVENT_TO_SEQ_RATIO, MAX_RAW_CPTS, SHIFT_CHANGE_THRESH,
    SCALE_CHANGE_THRESH, SIG_MATCH_THRESH, DNA_SAMP_TYPE, RNA_SAMP_TYPE,
    USE_RNA_EVENT_SCALE, RNA_SCALE_NUM_EVENTS, RNA_SCALE_MAX_FRAC_EVENTS,
    START_CLIP_PARAMS, STALL_PARAMS, COLLAPSE_RNA_STALLS, COLLAPSE_DNA_STALLS)
START_CLIP_PARAMS = th.startClipParams(*START_CLIP_PARAMS)
DEFAULT_STALL_PARAMS = th.stallParams(**STALL_PARAMS)
USE_START_CLIP_BASES = False

def map_read(
        fast5_data, aligner, std_ref,
        seq_samp_type=th.seqSampleType(DNA_SAMP_TYPE, False),
        bc_grp='Basecall_1D_000', bc_subgrp='BaseCalled_template',
        map_thr_buf=None, q_score_thresh=0, seq_len_rng=None, read_id = None):
    """Extract read sequence from the Fastq slot providing useful error messages

    Args:
        fast5_data (:class:`tombo.tombo_helper.readData`): read information
        aligner (mappy.Aligner): aligner object
        std_ref (:class:`tombo.tombo_stats.TomboModel`): canonical model (in order to extract extended genomic sequence)
        seq_samp_type (:class:`tombo.tombo_helper.seqSampleType`): sequencing sample type (default: DNA)
        bc_grp (str): group location containing read information (optional; default: 'Basecall_1D_000')
        bc_subgrp (str): sub-group location containing read information (optional; default: 'BaseCalled_template')
        map_thr_buf (mappy.ThreadBuffer): mappy thread buffer object (optional; default: None)
        q_score_thresh (float): basecalling mean q-score threshold (optional; default: 0/no filtering)
        seq_len_rng (tuple): allowed mapped sequence length range (optional; default: None/no filtering)

    Returns:
        :class:`tombo.tombo_helper.resquiggleResults` containing valid mapping values (signal to sequence assignment attributes will be ``None``)
    """
    seq_data = resquiggle.get_read_seq(
        fast5_data, bc_grp, bc_subgrp, seq_samp_type, q_score_thresh)

    try:
        # enumerate all alignments to avoid mappy memory leak
        alignment = list(aligner.map(str(seq_data.seq), buf=map_thr_buf))[0]
    except IndexError:
        raise th.TomboError('Alignment not produced')

    chrm = alignment.ctg
    # subtract one to put into 0-based index
    ref_start = alignment.r_st
    ref_end = alignment.r_en
    if not (seq_len_rng is None or
            seq_len_rng[0] < ref_end - ref_start < seq_len_rng[1]):
        raise th.TomboError(
            'Mapped location not within --sequence-length-range')
    strand = '+' if alignment.strand == 1 else '-'
    num_match = alignment.mlen
    num_ins, num_del, num_aligned = 0, 0, 0
    for op_len, op in alignment.cigar:
        if op == 1: num_ins += op_len
        elif op in (2,3): num_del += op_len
        elif op in (0,7,8): num_aligned += op_len
        elif op == 6: pass
        else:
            # soft and hard clipping are not reported in the
            # mappy cigar
            raise th.TomboError('Invalid cigar operation')

    # store number of clipped bases relative to read sequence
    if strand == '+':
        num_start_clipped_bases = alignment.q_st
        num_end_clipped_bases = len(seq_data.seq) - alignment.q_en
    else:
        num_start_clipped_bases = len(seq_data.seq) - alignment.q_en
        num_end_clipped_bases = alignment.q_st

    if not read_id:
        read_id = seq_data.id.decode()

    align_info = th.alignInfo(
        read_id, bc_subgrp, num_start_clipped_bases,
        num_end_clipped_bases, num_ins, num_del, num_match,
        num_aligned - num_match)

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
    #if len(genome_seq) != ref_end - ref_start + std_ref.kmer_width - 1:
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

def resquiggle_multifast5(fast5_fn):
    print('\n@JvN: multifast5')
    seq_samp_type = th.seqSampleType('DNA', False)
    # prep aligner, signal model and parameters
    aligner = mappy.Aligner(reference_fn, preset=str('map-ont'), best_n=1)
    std_ref = tombo_stats.TomboModel(seq_samp_type=seq_samp_type)
    rsqgl_params = tombo_stats.load_resquiggle_parameters(seq_samp_type)

    with h5py.File(fast5_fn, 'r') as hf:
        for i, fast5_data in enumerate(hf.keys()):
            fast5_data = hf[fast5_data]
            map_results = map_read(fast5_data, aligner, std_ref, read_id=str(fast5_data))
            # all_raw_signal = tombo_helper.get_raw_read_slot(fast5_data)['Signal'][:]
            print(map_results.align_info)
            print(map_results.segs)
            # map_results = map_results._replace(raw_signal=all_raw_signal)
            # print(map_results.align_info)

            if i > 3:
                break


def resquiggle_fast5(fast5_fn):
    print('\n@JvN: fast5')
    fast5_data = h5py.File(fast5_fn, 'r')

    # seq_samp_type = tombo_helper.get_seq_sample_type(fast5_data)
    seq_samp_type = th.seqSampleType('DNA', False)
    # prep aligner, signal model and parameters
    aligner = mappy.Aligner(reference_fn, preset=str('map-ont'), best_n=1)
    std_ref = tombo_stats.TomboModel(seq_samp_type=seq_samp_type)
    rsqgl_params = tombo_stats.load_resquiggle_parameters(seq_samp_type)

    # extract data from FAST5
    map_results = resquiggle.map_read(fast5_data, aligner, std_ref)
    all_raw_signal = th.get_raw_read_slot(fast5_data)['Signal'][:]
    # if seq_samp_type.rev_sig:
    #     all_raw_signal = all_raw_signal[::-1]
    map_results = map_results._replace(raw_signal=all_raw_signal)
    print(map_results.align_info)
    print(map_results.segs)

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
fast5_fn = rf'/home/noort/data/test1/tombo/0/00a2b3c7-6add-44f2-b171-3bd51e1ec637.fast5'
reference_fn = '/media/noort/Data/users/noort/ref_files/combined.fa'
fast5_fn = '/home/noort/data/Analysed_2020-05-12/tombo/2/0a4e18b3-e06e-4cf5-a338-06c3aba8bf33.fast5'

# tombo_folder = '/media/noort/Data/users/noort/test1/tombo_tmp
# df = read_tombo_mods(tombo_folder)
# chrms = df['chrm'].unique()
# plot_mods(df, reference_fn, chrms=chrms, min_length = 4000)
# print_mods(df, 20)

fast5_fn = '/home/noort/data/Analysed_2020-05-12/FAL22238_pass_1d8860a5_0.fast5'
resquiggle_multifast5(fast5_fn)
fast5_fn = '/home/noort/data/Analysed_2020-05-12/tombo/0/09120485-9b3b-4897-b64d-4452de65574b.fast5'
try:
    resquiggle_fast5(fast5_fn)
except tombo_helper.TomboError:
    print(f'{fast5_fn}: None')

# _plot_commands.plot_single_read(fast5_fn)
# tombo detect_modifications alternative_model --fast5-basedirs $path_to_fast5s --statistics-file-basename