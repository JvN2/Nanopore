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
    print('<<<megalodon command>>>:')
    cmd = f'megalodon {data_folder}'
    cmd += f' --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server'
    cmd += f' --guppy-config {config_file}'
    cmd += f' --remora-modified-bases dna_r9.4.1_e8 fast 0.0.0 5hmc_5mc CG 0'
    # cmd += f' --outputs basecalls mod_mappings mods'
    cmd += f' --outputs basecalls mod_mappings'
    cmd += f' --reference {ref_file}'
    # cmd += f' --devices 0 1 --processes 20'
    cmd += f' --output-directory {data_folder+ r"megalodon"}'
    cmd += f' --overwrite'
    cmd += f' --write-mods-text'

    print(cmd, '\n')

def megalodon_alphabet():
    print('<<<megalodon_aphabet command>>>:')
    cmd = r'megalodon_extras modified_bases describe_alphabet'
    cmd += f' --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server'
    print(cmd, '\n')

if __name__ == '__main__':
    if True:
        config_file = r'/opt/ont/guppy/data/dna_r9.4.1_450bps_hac.cfg'
        config_file = r'/opt/ont/guppy/data/dna_r9.4.1_450bps_modbases_5mc_hac.cfg'
        data_folder = r'/media/noort/Data/users/noort/test'

        config_file = r'/opt/ont/guppy/data/dna_r9.4.1_450bps_modbases_5mc_hac.cfg'
        data_folder = r'/media/noort/Data/users/noort/test2/20220118_1825_MN30914_ACR467_7759a37d/fast5_pass/All_data_combined_singlefast5/0'

        output_folder = data_folder + r'/called_with_methylation'
        output_folder = r'/media/noort/Data/users/noort/test2'
        guppy_command(data_folder, output_folder, config_file)

    if True:
        data_folder = r'/media/noort/Data/users/noort/test3/'
        config_file = r'dna_r9.4.1_450bps_modbases_5mc_hac.cfg'
        ref_file = r'/media/noort/Data/users/noort/ref_files/combined.fasta'
        ref_file = r'/media/noort/Data/users/noort/ref_files/GAL_locus+14xPP7+2xRS_circle_V2.fasta'
        ref_file = r'/media/noort/Data/users/noort/ref_files/GCF_000146045.2_R64_genomic.fna'
        ref_file = r'/media/noort/Data/users/noort/ref_files/GCF_000146045.2_R64_genomic+601_arrays.fna'
        ref_file = r'/media/noort/Data/users/noort/ref_files/ref_601s.fasta'
        megalodon_command(data_folder, config_file, ref_file)
        megalodon_alphabet()
    if False:
        data_folder = output_folder
        output_folder = output_folder + r'/single_fast5'
        multi_to_single_command(data_folder, output_folder)

    if False:
        fasta_file = r'/media/noort/Data/users/noort/test2/S_CP130_pUC18_16x197.fasta'
        ont_fast5_api=r'/media/noort/Data/users/noort/test2/HS_BAC.fasta'
        fast5_folder = r'/media/noort/Data/users/noort/test3/called1'
        tombo_command(fast5_folder, fasta_file)


    import subprocess

    cmd = ['echo', 'More output']
    cmd = ['megalodon', '-h']
    # cmd = ['guppy_basecaller', '-h']
    process = subprocess.Popen(cmd ,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    # print(stdout.decode(), stderr.decode())
