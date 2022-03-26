import h5py
import matplotlib.pyplot as plt
from ont_fast5_api.fast5_interface import get_fast5_file
import numpy as np
from tombo import resquiggle

def printall(name, obj):
    print(name, dict(obj.attrs))

def read_fast5(filename, print_all = False):
    with h5py.File(filename, 'r') as hdf:
        if print_all:
            hdf.visititems(printall)

        read_nr = list(hdf['Raw/Reads'])[0]
        dict = {'signal': [], 'readID': '', 'digitisation': 0.0, 'offset': 0.0, 'range': 0.0, 'sampling_rate': 0.0}
        readID = hdf['Raw/Reads'][read_nr].attrs['read_id'].decode()
        dict['readID'] = readID
        dict['digitisation'] = hdf['UniqueGlobalKey/channel_id'].attrs['digitisation']
        dict['offset'] = hdf['UniqueGlobalKey/channel_id'].attrs['offset']
        dict['range'] = float("{0:.2f}".format(hdf['UniqueGlobalKey/channel_id'].attrs['range']))
        dict['sampling_rate'] = hdf['UniqueGlobalKey/channel_id'].attrs['sampling_rate']
        dict['sequence'] = resquiggle.get_read_seq(hdf).seq

    with get_fast5_file(filename, mode="r") as f5:
        for read in f5.get_reads():
                dict['signal'] = read.get_raw_data()
    return dict

if __name__ == '__main__':
    filename = '/home/noort/data/data_to_docker/ACR467_pass_barcode04_6749c10b_0.fast5'
    filename = '/home/noort/data/Analysed_2022-01-18_12samplemethylationtest/guppy/workspace/barcode04/ACR467_pass_barcode04_6749c10b_0.fast5'
    # filename = '/media/noort/Data/users/noort/test2/called0/workspace/201bc07d-38fa-4566-b458-7cac4fd26964.fast5'
    # filename = '/media/noort/Data/users/noort/test2/called0/workspace/1c8d0eaf-9dad-4b2b-bfe2-a12aa8952561.fast5'
    # filename =r'/media/noort/Data/users/noort/test4/called04/single_fast5/0/97331af7-b49b-4676-97fc-6480624f6697.fast5'
    dict = read_fast5(filename, print_all=True)



    for key in dict:
        print(key,'>', dict[key])

    print('length >',len(dict['sequence']))
    t = np.arange(len(dict['signal']))/dict['sampling_rate']
    s = dict['signal']*dict['range']/dict['digitisation'] + dict['offset']

    plt.plot(t, s)
    plt.xlabel('t (s)')
    plt.ylabel('I (pA)')
    plt.show()
