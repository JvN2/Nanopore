def guppy_command(data_folder, output_folder, config_file):
    print('<<<guppy command>>>:')
    cmd = f'guppy_basecaller --input_path {data_folder} --save_path {output_folder} --config {config_file} --fast5_out'
    print(cmd)


def tombo_command(fast5_folder, fasta_file):
    print('<<<tombo command>>>:')
    cmd = f'tombo resquiggle {fast5_folder} {fasta_file} --processes 4 --num-most-common-errors 5'
    # cmd = f'tombo preprocess annotate_raw_with_fastqs --fast5-basedir <fast5s-base-directory> --fastq-filenames <reads.fastq>'
    print(cmd)



if __name__ == '__main__':
    data_folder = r'/media/noort/Data/users/noort/test2/reads0'
    output_folder = r'/media/noort/Data/users/noort/test2/called0'
    config_file = r'/media/noort/Data/users/noort/test2/dna_r9.4.1_450bps_hac.cfg'
    guppy_command(data_folder, output_folder, config_file)

    fasta_file = r'/media/noort/Data/users/noort/test2/S_CP130_pUC18_16x197.fasta'
    fast5_folder = r'/media/noort/Data/users/noort/test2/called0/workspace'
    tombo_command(fast5_folder, fasta_file)