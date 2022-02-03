def guppy_command(data_folder, output_folder, config_file):
    print('<<<guppy command>>>:')
    cmd = f'guppy_basecaller --input_path {data_folder} --save_path {output_folder} --config {config_file} --fast5_out'
    print(cmd, '\n')


def tombo_command(fast5_folder, fasta_file):
    print('<<<tombo command>>>:')
    cmd = f'tombo resquiggle {fast5_folder} {fasta_file} --processes 4 --num-most-common-errors 5'
    print(cmd, '\n')
    # cmd = f'tombo preprocess annotate_raw_with_fastqs --fast5 -basedir <fast5s-base-directory> --fastq-filenames <reads.fastq>'
    # print(cmd, '\n')


def multi_to_single_command(data_folder, output_folder):
    print('<<<multi_to_single command>>>:')
    cmd = f'multi_to_single_fast5 -i {data_folder} -s {output_folder} --recursive'
    print(cmd, '\n')

if __name__ == '__main__':
    data_folder = r'/media/noort/Data/users/noort/test4/barcode02'
    output_folder = r'/media/noort/Data/users/noort/test4/called2'
    config_file = r'/media/noort/Data/users/noort/test2/dna_r9.4.1_450bps_hac.cfg'


    config_file = r'/opt/ont/guppy/data/dna_r9.4.1_450bps_modbases_5mc_hac.cfg'
    data_folder = r'/media/noort/Data/users/noort/test4/barcode04'
    output_folder = r'/media/noort/Data/users/noort/test4/called04'

    guppy_command(data_folder, output_folder, config_file)

    # data_folder = r'/media/noort/Data/users/noort/test4/barcode02'
    # data_folder = r'/media/noort/Data/users/noort/test3/called1/workspace'
    data_folder = output_folder
    output_folder = output_folder + r'/single_fast5'
    multi_to_single_command(data_folder, output_folder)


    fasta_file = r'/media/noort/Data/users/noort/test2/S_CP130_pUC18_16x197.fasta'
    ont_fast5_api=r'/media/noort/Data/users/noort/test2/HS_BAC.fasta'
    fast5_folder = r'/media/noort/Data/users/noort/test3/called1'
    tombo_command(fast5_folder, fasta_file)