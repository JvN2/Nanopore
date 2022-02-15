from ont_fast5_api.fast5_interface import get_fast5_file
import matplotlib.pyplot as plt
import h5py
import numpy as np


def read_info(filename):
    def printall(name, obj):
        print(name, dict(obj.attrs))

    with h5py.File(filename, 'r') as hf:
        reads = list(hf.keys())
        hf[reads[0]].visititems(printall)
        # plt.plot((np.asarray(hf[reads[0]]['Raw']['Signal'])))
        # plt.show()
        print(np.asarray(hf[reads[0]]['Analyses/Basecall_1D_000/BaseCalled_template/ModBaseProbs']))


def read_modifications(filename):
    def to_ascii(text):
        ascii_values = [ord(character) for character in text]
        return ascii_values

    with get_fast5_file(filename, mode="r") as f5:
        for read_id in list(f5.get_read_ids())[:3]:
            read = f5.get_read(read_id)
            latest_basecall = read.get_latest_analysis('Basecall_1D')
            mod_base_table = read.get_analysis_dataset(latest_basecall, 'BaseCalled_template/ModBaseProbs')
            called_base_table = read.get_analysis_dataset(latest_basecall, 'BaseCalled_template/Fastq').split('\n')
            seq = called_base_table[1]
            # for line in called_base_table:
            #     print(line)
            # plt.plot(np.asarray(to_ascii(called_base_table[3])).astype(float)/255)
            # plt.show()

            mod_metadata = read.get_analysis_attributes(f'{latest_basecall}/BaseCalled_template/ModBaseProbs')
            called_metadata = read.get_analysis_attributes(f'{latest_basecall}/BaseCalled_template/Fastq')

            # for mods in mod_base_table:
            #     print(mods)
            alphabet = read.get_analysis_attributes(f'{latest_basecall}/BaseCalled_template/ModBaseProbs')[
                'output_alphabet']
            print(alphabet)
            print(seq)
            methylated = 1 -(mod_base_table[:, 1] / 255.)
            res = ''
            for i in methylated:
                if i > 0.95:
                    res += '1'
                else:
                    res += '0'
            print(res)
            plt.plot(methylated, '.')
            plt.show()


filename = r'/home/john/Documents/Nanopore/barcode04/called/workspace/ACR467_pass_barcode04_6749c10b_0.fast5'
# filename = r'/home/john/Documents/Nanopore/202006/called/workspace/FAL22238_pass_1d8860a5_10.fast5.temp'
read_info(filename)
read_modifications(filename)

