def guppy_command(data_folder, output_folder, config_file):
    print('<<<guppy command>>>:')
    cmd = f'guppy_basecaller' \
           f' --input_path {data_folder}' \
           f' --save_path {output_folder}' \
           f' --config {config_file}' \
           f' --fast5_out'
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


def megalodon_command(data_folder, config_file, ref_file):
    print('<<<megalodon_command command>>>:')
    cmd = f'megalodon {data_folder}'
    cmd += f' --guppy-config {config_file}'
    cmd += f' --remora-modified-bases dna_r9.4.1_e8 fast 0.0.0 5hmc_5mc CG 0'
    cmd += f' --outputs basecalls mappings mod_mappings mods'
    cmd += f' --reference {ref_file}'
    # cmd += f' --devices 0 1 --processes 20'
    cmd += f' --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server'
    cmd += f' --overwrite'

    print(cmd, '\n')


if __name__ == '__main__':
    config_file = r'/opt/ont/guppy/data/dna_r9.4.1_450bps_hac.cfg'
    # config_file = r'/opt/ont/guppy/data/dna_r9.4.1_450bps_modbases_5mc_hac.cfg'
    data_folder = r'/media/noort/Data/users/noort/test/barcode04'
    output_folder = data_folder + r'/called_nomethylation'
    guppy_command(data_folder, output_folder, config_file)

    config_file = r'dna_r9.4.1_450bps_modbases_5mc_hac.cfg'
    ref_file = r'//media/noort/Data/users/noort/test/S_CP130_pUC18_16x197.fasta'
    megalodon_command(output_folder, config_file, ref_file)

    if False:
        data_folder = output_folder
        output_folder = output_folder + r'/single_fast5'
        multi_to_single_command(data_folder, output_folder)

    if False:
        fasta_file = r'/media/noort/Data/users/noort/test2/S_CP130_pUC18_16x197.fasta'
        ont_fast5_api=r'/media/noort/Data/users/noort/test2/HS_BAC.fasta'
        fast5_folder = r'/media/noort/Data/users/noort/test3/called1'
        tombo_command(fast5_folder, fasta_file)

