from ont_fast5_api import multi_fast5
from ont_fast5_api.fast5_interface import get_fast5_file

def print_all_raw_data(fast5_filepath):
    with get_fast5_file(fast5_filepath, mode="r") as f5:
        for read in f5.get_reads():
            raw_data = read.get_raw_data()
            print(read.read_id, raw_data)


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


./.local/lib/python3.8/site-packages (4.0.0)/multi_to_single_fast5 --input_path /media/noort/Data/users/noort/test/pass/fastq_runid_6749c10bc43747f33508adb64ddbcc4a8b108508_0_0.fastq --save_path /media/noort/Data/users/noort/test/multi_to_single --recursive --t 20

filename = "/media/noort/Data/users/noort/test/pass/fastq_runid_6749c10bc43747f33508adb64ddbcc4a8b108508_0_0.fastq"
# This can be a single- or multi-read file
print_all_raw_data(filename)
multi_fast5.Fast5Read.