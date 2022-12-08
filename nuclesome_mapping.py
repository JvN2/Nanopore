import EditChromatin
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt

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
    strands = [1, 0]
    result = []
    strand = []
    for block in strands:
        df = pd.read_hdf(filename, key=f'Statistic_Blocks/Block_{block}/block_stats')
        for id in tqdm(df['read_id'].unique(), desc=f'Reading strand: {block}'):
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

    # if motifs is not None:
    #     filter = np.zeros(len(seq))
    #     for motif in motifs:
    #         index = find_motif(seq, motif, format='index')
    #         offset = [i for i, c in enumerate(motif) if c.isupper()]
    #         filter[index + offset[0]] = 1
    #     result = np.multiply(result, filter)
    return result, np.asarray(strand)

def read_fasta(filename, name=None):
    with open(filename, 'r') as f:
        content = f.readlines()
        for line in range(len(content)):
            if content[line][0] == '>':
                seq = content[line + 1]
    seq = ''.join([s for s in seq.upper() if s in ['A', 'C', 'G', 'T']])
    return seq


if __name__ == '__main__':
    filename = r'/media/noort/Data/users/noort/20220816_barcode10_selection/read_stats.5mC.tombo.per_read_stats'
    fast_filename = r'/media/noort/Data/users/noort/20220816_barcode10_selection/LinearReference_CP130.fasta'

    seq = read_fasta(fast_filename).upper()
    tmp = EditChromatin.calc_nucleosome_positions(seq, 100, -4 )

    llr, strands = read_stats(filename, seq)
    colorrange = 2
    #plt.imshow(llr, cmap='bwr', vmin=-colorrange, vmax=colorrange, origin='lower')
    # plt.plot(llr[700])
    # plt.show()

    #df2 = pd.read_csv(tombo_filename)
    #print(df2)

    mu= 10
    dyad_probability = EditChromatin.vanderlick(llr[700], mu)
    nucleosome_occupancy = np.convolve(dyad_probability, np.ones(146), mode='same')
    plt.plot(nucleosome_occupancy)
    plt.plot(tmp)
    plt.show()