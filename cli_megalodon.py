from cli_guppy import cmd_guppy
from megalodon.backends import ModelInfo

def cmd_megalodon_extras(guppy_config, guppy_bin):
    cmd = [rf'megalodon_extras modified_bases describe_alphabet']
    cmd.append(rf'--guppy-server-path {guppy_bin}')
    cmd.append(rf'--guppy-config {guppy_config}')
    return ' \\\n'.join(cmd)


def cmd_megalodon(fast5_dir, reference_genome, guppy_bin, guppy_config, rerio_models):
    cmd = [rf'megalodon {fast5_dir}']
    cmd.append(rf'--guppy-server-path {guppy_bin}')
    # cmd.append(rf'--guppy-config {guppy_config}')
    # cmd.append(rf'--guppy-params "-d {rerio_models}"')
    cmd.append(rf'--reference {reference_genome}')
    cmd.append(rf'--remora-modified-bases dna_r9.4.1_e8 fast 0.0.0 5hmc_5mc CG 0')
    cmd.append(rf'--outputs per_read_mods')
    cmd.append(rf'--output-directory {fast5_dir}megalodon')
    cmd.append(rf'--write-mods-text')
    cmd.append(rf'--processes 18')
    cmd.append(rf'--overwrite ')
    return ' \\\n'.join(cmd)


def cmd_start_docker():
    image = 'badcad167147'
    cmd = [rf'docker run --name box -it {image} /bin/bash']
    return ' \\\n'.join(cmd)

def cmd_copy_to_docker():
    cmd = [rf'docker cp /media/noort/Data/users/noort/docker/to_container/. box:/home/data']
    return ' \\\n'.join(cmd)

if __name__ == '__main__':

    config = 'docker'
    if config == 'local':
        fast5_dir = r'/home/noort/data/test21'
        rerio_models = r'/home/noort/Downloads/rerio-master/basecall_models/'
        reference_genome = r'/media/noort/Data/users/noort/ref_files/combined.fa'
        guppy_config = r'res_dna_r941_min_modbases-all-context_v001.cfg'
        guppy_bin = r'/opt/ont/guppy/bin/guppy_basecall_server'
    else:
        fast5_dir = r'/home/data/'
        reference_genome = r'/home/data/combined.fa'
        rerio_models = r'/home/rerio/basecall_models/'
        guppy_config = r'res_dna_r941_min_modbases-all-context_v001.cfg'
        guppy_bin = r'/usr/bin/guppy_basecall_server'


    print(cmd_start_docker())
    print(cmd_copy_to_docker(), '\n')
    print(cmd_guppy(fast5_dir, guppy_config, rerio_models), '\n')
    print(cmd_megalodon(fast5_dir, reference_genome, guppy_bin, guppy_config, rerio_models), '\n')