import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import nuc_tool as nt


def vanderlick(Energy, mu, convolve = True):
    E_out = Energy - mu
    footprint = 147

    forward = np.zeros(len(Energy))
    for i in range( len(Energy) ):
        tmp = sum(forward[ max( i - footprint , 0):i])
        forward[i] = np.exp(E_out[i] - tmp)

    backward = np.zeros(len(Energy))
    r_forward = forward[::-1]
    for i in range( len(Energy) ):
        backward[i] = 1 - sum(r_forward[max(i - footprint , 0):i]* backward[max(i - footprint,0):i])

    P = forward * backward[::-1]

    if convolve:
        P = np.convolve(P,np.ones(146), mode = 'same')
    return P

def megalodon_command(data_folder, config_file, ref_file, output_folder):
    print('<<<megalodon command>>>:')
    cmd = f'megalodon {data_folder}'
    cmd += f' --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server'
    cmd += f' --guppy-config {config_file}'
    cmd += f' --remora-modified-bases dna_r9.4.1_e8 fast 0.0.0 5hmc_5mc CG 0'
    # cmd += f' --remora-modified-bases res_dna_r941_min_modbases-all-context_v001.cfg'
    # cmd += f' --outputs basecalls mod_mappings mods'
    cmd += f' --outputs basecalls mod_mappings'
    cmd += f' --reference {ref_file}'
    # cmd += f' --devices 0 1 --processes 20'
    cmd += f' --output-directory {output_folder}'
    cmd += f' --overwrite'
    cmd += f' --write-mods-text'
    print(cmd, '\n')


def read_ref(filename):
    txt = open(filename, "rt").read()
    txt = txt.split('>')
    txt = [i.split('\n')[:2] for i in txt if i]
    df = pd.DataFrame(txt, columns=['chrm', 'seq'])
    df.set_index('chrm', inplace=True)
    return df

def plot_mods(ref_file, mod_base_file, chrm):
    refs = read_ref(ref_file)
    seq = refs.loc[chrm]['seq'].upper()

    all_mods = pd.read_csv(mod_base_file, sep='\t')
    mods = all_mods[all_mods['chrm'] == chrm]
    mods.insert(len(mods.columns), 'LLR', mods['mod_log_prob']- mods['can_log_prob'], True)

    df1 = mods.groupby('read_id').agg({'strand': 'first', 'pos' : np.ptp})
    df1.rename(columns = {'pos':'range'}, inplace = True)

    for item in ['pos','mod_base', 'LLR']:
        df1 = pd.merge(df1, mods.groupby('read_id')[item].apply(np.asarray).reset_index(name=item).set_index('read_id'),
                       left_index=True, right_index=True, how='outer')

    df1 = df1[df1['range'] > 3500]
    df1.sort_values(by=['strand'], inplace=True)

    #only keep most likely modification
    for read_id, row in df1.iterrows():
        n_mods = len(np.unique(row['mod_base']))
        mod_index = [row == np.max(row) for row in np.reshape(row['LLR'], (-1, n_mods))]
        mod_index = np.reshape(mod_index, (-1))
        for item in ['pos', 'mod_base', 'LLR']:
            df1.at[read_id,item] = row[item][mod_index]

    # fill image with LLR
    im = np.zeros( (len(df1), len(refs.loc[chrm]["seq"])))
    for i, _ in enumerate(im):
        im[i, df1.iloc[i]['pos']] = -df1.iloc[i]['LLR']

    # convert LLR to nucleosome occupancy
    # mu = -5
    # plt.plot(nt.CalcNucPositions(seq, mu=mu))
    # for i, e in enumerate(im[:2]):
    #     # im[i] = vanderlick(e, 10)
    #     print(i)
    #     im[i] = nt.CalcNucPositions(seq, mu=mu, penalty= im[i])
    #     plt.plot(im[i])
    # plt.show()
    # return
    if im.any():
        plt.imshow(im, vmin= -2, vmax= 2,  origin = 'lower', cmap=plt.get_cmap('seismic'), aspect = 'auto')#, interpolation='none') #
        plt.show()
    else:
        print('No reads found that match criteria')



data_folder = r'/media/noort/Data/users/noort/test4/'
# data_folder = r'/media/noort/Data/users/noort/test4/barcode04'
output_folder = r'megalodon_data4/'
config_file = r'dna_r9.4.1_450bps_modbases_5mc_hac.cfg'
ref_file = r'/media/noort/Data/users/noort/ref_files/ref_601s.fasta'
mod_base_file = data_folder + output_folder + 'per_read_modified_base_calls.txt'

# megalodon_command(data_folder, config_file, ref_file, output_folder)

mod_base_file = r'/media/noort/Data/users/noort/test8/megalodon_results//per_read_modified_base_calls.txt'
plot_mods(ref_file, mod_base_file, 'pCP130')