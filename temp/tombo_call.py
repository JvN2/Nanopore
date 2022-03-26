import os
from multiprocessing import Pool

from matplotlib import pyplot as plt
from ont_fast5_api import multi_fast5
from ont_fast5_api.conversion_tools.conversion_utils import get_fast5_file_list, get_progress_bar
from ont_fast5_api.conversion_tools.multi_to_single_fast5 import convert_multi_to_single
from ont_fast5_api.fast5_interface import get_fast5_file

def print_all_raw_data(fast5_filepath):
    with get_fast5_file(fast5_filepath, mode="r") as f5:
        for read in f5.get_reads():
            raw_data = read.get_raw_data()
            print(read.read_id, raw_data)

    return raw_data



def batch_convert_multi_files_to_single(input_path, output_folder, threads, recursive, follow_symlinks):
    pool = Pool(threads)
    file_list = get_fast5_file_list(input_path, recursive, follow_symlinks=follow_symlinks)
    pbar = get_progress_bar(len(file_list))

    def update(result):
        input_file = result[0]
        with open(os.path.join(output_folder, "filename_mapping.txt"), 'a') as output_table:
            for filename in result[1]:
                output_table.write("{}\t{}\n".format(input_file, filename))
        pbar.update(pbar.currval + 1)

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    results_array = []
    for batch_num, filename in enumerate(file_list):
        results_array.append(pool.apply_async(convert_multi_to_single,
                                              args=(filename, output_folder,
                                                    str(batch_num)),
                                              callback=update))

    pool.close()
    pool.join()
    pbar.finish()

filename = "/home/noort/data/test21/tombo/0/00a2b3c7-6add-44f2-b171-3bd51e1ec637.fast5"
# This can be a single- or multi-read file
sq = print_all_raw_data(filename)

from tombo import resquiggle

resquiggle()