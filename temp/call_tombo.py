import h5py, mappy
import numpy as np
from tombo import tombo_helper, tombo_stats, resquiggle
import show_squiggle
from Bio import pairwise2
import matplotlib.pyplot as plt
import glob

def call_tombo(fast5_fn, reference_fn, show = False):
    # set read values
    # fast5_fn, reference_fn = fast5_file, fasta_file
    fast5_data = h5py.File(fast5_fn, 'r')
    seq_samp_type = tombo_helper.get_seq_sample_type(fast5_data)

    # prep aligner, signal model and parameters
    aligner = mappy.Aligner(reference_fn, preset=str('map-ont'), best_n=1)
    std_ref = tombo_stats.TomboModel(seq_samp_type=seq_samp_type)
    rsqgl_params = tombo_stats.load_resquiggle_parameters(seq_samp_type)

    # extract data from FAST5
    map_results = resquiggle.map_read(fast5_data, aligner, std_ref)
    if map_results:
        all_raw_signal = tombo_helper.get_raw_read_slot(fast5_data)['Signal'][:]
        if seq_samp_type.rev_sig:
            all_raw_signal = all_raw_signal[::-1]
        map_results = map_results._replace(raw_signal=all_raw_signal)

        # run full re-squiggle
        rsqgl_results = resquiggle.resquiggle_read(
            map_results, std_ref, rsqgl_params, all_raw_signal=all_raw_signal)

        # or run individual steps
        num_events = tombo_stats.compute_num_events(
            all_raw_signal.shape[0], len(map_results.genome_seq),
            rsqgl_params.mean_obs_per_event)
        valid_cpts, norm_signal, scale_values = resquiggle.segment_signal(
            map_results, num_events, rsqgl_params)
        event_means = tombo_stats.compute_base_means(norm_signal, valid_cpts)
        dp_results = resquiggle.find_adaptive_base_assignment(
            valid_cpts, event_means, rsqgl_params, std_ref, map_results.genome_seq)
        norm_signal = norm_signal[
                      dp_results.read_start_rel_to_raw:
                      dp_results.read_start_rel_to_raw + dp_results.segs[-1]]
        segs = resquiggle.resolve_skipped_bases_with_raw(
            dp_results, norm_signal, rsqgl_params)


        fast5 = show_squiggle.read_fast5(fast5_file)
        print(f'Guppy:  {len(fast5["sequence"])} bp')
        print(f'Tombo: {len(rsqgl_results.genome_seq)} bp')

        seq_601 = 'CTGGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGCAAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCAGATATATACATCCTGT'

        alignments = pairwise2.align.globalxx(fast5['sequence'], rsqgl_results.genome_seq)
        # alignments = pairwise2.align.globalxx(fast5['sequence'], seq_601)
        print(alignments[0][0])
        print(alignments[0][1])
        match = (np.asarray(list(alignments[0][0])) == np.asarray(list(alignments[0][1])))
        print(f'Matched {np.sum(match)} bp\n')
        if show:
            plt.plot(np.cumsum(match))
            plt.show()


fasta_file = r'/media/noort/Data/users/noort/test2/S_CP130_pUC18_16x197.fasta'
# fasta_file = r'/media/noort/Data/users/noort/test2/HS_BAC.fasta'
fast5_file = r'/media/noort/Data/users/noort/test4/barcode02/called/0/0a873c73-a0cf-4654-ab8f-236667df61de.fast5'
fast5_files = glob.glob(r'/media/noort/Data/users/noort/test4/barcode02/called/0/*.fast5')
# fast5_files = glob.glob(r'/media/noort/Data/users/noort/test4/called04/single_fast5/0/*.fast5')

for fast5_file in fast5_files:
    print(fast5_file)
    call_tombo(fast5_file, fasta_file)
