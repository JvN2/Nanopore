from cli_guppy import cmd_guppy, cmd_docker


def cmd_multi_to_single(fast5_dir, single=True):
    cmd = [rf'multi_to_single_fast5']
    cmd.append(rf'-i {fast5_dir}/guppy/workspace/')
    cmd.append(rf'-s  {fast5_dir}/tombo/')
    if single: cmd = ' \\\n'.join(cmd)
    return cmd


def cmd_resquigle(fast5_dir, reference_genome):
    cmd = [rf'tombo resquiggle {fast5_dir} {reference_genome}']
    cmd.append(rf'--num-most-common-errors 5')
    cmd.append(rf'--overwrite')
    cmd.append(rf'--processes 18')
    cmd.append(rf'--skip-index')
    return ' \\\n'.join(cmd)


def cmd_detect_mods(fast5_dir, model='alternative', filename_base=''):
    if model == 'alternative':
        cmd = [rf'tombo detect_modifications alternative_model']
        cmd.append(rf'--alternate-bases 5mC 6mA')
    else:
        cmd = [rf'tombo detect_modifications de_novo']
    cmd.append(rf'--statistics-file-basename {filename_base}.stats')
    cmd.append(rf'--per-read-statistics-basename {filename_base}.read_stats')
    cmd.append(rf'--processes 18')
    cmd.append(rf'--fast5-basedirs {fast5_dir}/tombo')
    return ' \\\n'.join(cmd)

if __name__ == '__main__':
    config = 'docker'
    local_base_folder = rf'/home/kuijntjes/Desktop/2022-10-12_WholeCellExtractGalLocusCindy/no_sample/20221012_1615_MN30914_AJA380_b85344d5/fast5_pass/barcode01'
    if config == 'local':
        fast5_dir = local_base_folder
        rerio_models = r'/home/noort/Downloads/rerio-master/basecall_models/'
        reference_genome = local_base_folder + r'/LinearReference115_130.fasta'
        guppy_config = r'res_dna_r941_min_modbases-all-context_v001.cfg'
        # guppy_bin = r'/opt/ont/guppy/bin/guppy_basecall_server'
    else:
        fast5_dir = r'/home/data'
        reference_genome = r'/home/data/LinearGalInplasmid.fasta'
        rerio_models = r'/home/rerio/basecall_models/'
        guppy_config = r'res_dna_r941_min_modbases-all-context_v001.cfg'
        # guppy_bin = r'/usr/bin/guppy_basecall_server'

    data_dir = r'/home/kuijntjes/Desktop/2022-10-12_WholeCellExtractGalLocusCindy/no_sample/20221012_1615_MN30914_AJA380_b85344d5/fast5_pass'
    data_dir = r'c:/tmp/DEPC'

    fast5_dir = 'barcode03'

    fast5_dir = r'/home/data'
    #fast5_dir = '/'
    print(cmd_docker(data_dir), '\n')
    print(cmd_guppy(fast5_dir, guppy_config, rerio_models), '\n')
    print(cmd_multi_to_single(fast5_dir), '\n')
    print(cmd_resquigle(fast5_dir, reference_genome), ' \n')
    print(cmd_detect_mods(fast5_dir, filename_base=fast5_dir), '\n')
