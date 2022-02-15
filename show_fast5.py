import linecache

import h5py
import matplotlib.pyplot as plt
import numpy as np
from ont_fast5_api.fast5_interface import get_fast5_file


def ascii_to_value(text):
    values = np.asarray([ord(character) for character in text]).astype(float)
    return values / 255


def read_info(filename, show=False):
    def printall(name, obj):
        print(name, dict(obj.attrs))

    with h5py.File(filename, 'r') as hf:
        reads = list(hf.keys())
        hf[reads[0]].visititems(printall)

        print(np.asarray(hf[reads[0]]['Analyses/Basecall_1D_000/BaseCalled_template/ModBaseProbs']))
        if show:
            plt.plot((np.asarray(hf[reads[0]]['Raw']['Signal'])))
            plt.show()


def read_fast5(filename, read_nr=0):
    with get_fast5_file(filename, mode="r") as f5:
        read_id = list(f5.get_read_ids())[read_nr]
        read = f5.get_read(read_id)
        latest_basecall = read.get_latest_analysis('Basecall_1D')
        called_base_table = read.get_analysis_dataset(latest_basecall, 'BaseCalled_template/Fastq').split('\n')

        seq = called_base_table[1]
        quality = ascii_to_value(called_base_table[3])

        try:
            mod_base_table = read.get_analysis_dataset(latest_basecall, 'BaseCalled_template/ModBaseProbs')
            mod = 1 - (mod_base_table[:, 1] / 255.)
            mod_alphabet = \
                read.get_analysis_attributes(f'{latest_basecall}/BaseCalled_template/ModBaseProbs')['output_alphabet']
        except:
            mod = None
            mod_alphabet = None


def read_fastq(filename, read_nr=0):
    seq = linecache.getline(filename, (read_nr * 4) + 2)
    quality = ascii_to_value(linecache.getline(filename, (read_nr * 4) + 4))

    print(seq)
    print(quality)
    print(len(seq), len(quality))
    plt.plot(quality)
    plt.show()


filename = r'/home/john/Documents/Nanopore/barcode04/called/workspace/ACR467_pass_barcode04_6749c10b_0.fast5'
# filename = r'/home/john/Documents/Nanopore/202006/called/workspace/FAL22238_pass_1d8860a5_10.fast5.temp'
# read_info(filename)
# read_fast5(filename)

filename = r'C:\Users\noort\Downloads\fastq_runid_1d8860a5a693832b45cbe93cd13b8886c3f0b5f4_0_0.fastq'
read_fastq(filename, 19)
