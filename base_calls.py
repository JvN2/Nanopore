import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def megalodon_command(data_folder, config_file, ref_file, output_folder):
    print('<<<megalodon command>>>:')
    cmd = f'megalodon {data_folder}'
    cmd += f' --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server'
    cmd += f' --guppy-config {config_file}'
    cmd += f' --remora-modified-bases dna_r9.4.1_e8 fast 0.0.0 5hmc_5mc CG 0'
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


data_folder = r'/media/noort/Data/users/noort/test3/'
output_folder = r'megalodon_data4/'
config_file = r'dna_r9.4.1_450bps_modbases_5mc_hac.cfg'
ref_file = r'/media/noort/Data/users/noort/ref_files/ref_601s.fasta'
mod_base_file = data_folder + output_folder + 'per_read_modified_base_calls.txt'

try:
    refs = read_ref(ref_file)
    # for chr, row in refs.iterrows():
    #     print(chr, len(row['seq']))

    all_mods = pd.read_csv(mod_base_file, sep='\t')
    # print(all_mods.head(3))

    for chrm in ['pCP130']:
        mods = all_mods[all_mods['chrm'] == chrm]
        seq = refs.loc[chrm]["seq"].upper()
        read_ids = list(set(mods['read_id']))
        print(f'{len(read_ids)} reads for {chrm} ({len(seq)} bp)')
        print(f'{seq.count("CG"):4}', seq)
        seq+='*'
        im =[]
        for read_id in read_ids[:50]:
            # mod_string = list(' ') * len(seq)
            mod_string = ['.' if seq[i:i+2] == 'CG' else ' ' for i,_ in enumerate(list(seq)) ]
            for _, row in mods[mods['read_id'] == read_id].iterrows():
                mod_string[row['pos']] = row['mod_base']
            print(f'{np.sum(mods["read_id"] == read_id):4}', ''.join(mod_string))
            im.append(np.asarray(mod_string) == 'm')
        im = np.asarray(im).astype(int)
        plt.imshow(im)
        plt.show()

except FileNotFoundError:
    print('File(s) not foond, run Megalodon first.')
    megalodon_command(data_folder, config_file, ref_file, output_folder)
