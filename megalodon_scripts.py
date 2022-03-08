config = 'docker'

if config == 'local':
    FAST5_FILES = r'/media/noort/Data/users/noort/test6/'
    RERIO_MODELS = r'/home/noort/Downloads/rerio-master/basecall_models/'
    REFRENCE_GENOME = r'/media/noort/Data/users/noort/ref_files/combined.fasta'
    GUPPY_BIN = r'/opt/ont/guppy/bin/guppy_basecall_server'
    GUPPY_CONFIG = r'res_dna_r941_min_modbases-all-context_v001.cfg'
else:
    FAST5_FILES = r'/home/data/'
    REFRENCE_GENOME = r'/home/data/combined.fasta'
    RERIO_MODELS = r'/home/rerio/basecall_models/'
    GUPPY_BIN = r'/usr/bin/guppy_basecall_server'
    GUPPY_CONFIG = r'res_dna_r941_min_modbases-all-context_v001.cfg'
    # GUPPY_CONFIG = r'dna_r9.4.1_450bps_fast.cfg'

def cmd_megalodon_extras():
    cmd = [rf'megalodon_extras modified_bases describe_alphabet']
    cmd.append(rf'--guppy-server-path {GUPPY_BIN}')
    return ' \\\n'.join(cmd)

def cmd_megalodon():
    cmd = [rf'megalodon {FAST5_FILES}']
    cmd.append(rf'--guppy-server-path {GUPPY_BIN}')
    cmd.append(rf'--guppy-config {GUPPY_CONFIG}')
    cmd.append(rf'--guppy-params "-d {RERIO_MODELS}"')
    cmd.append(rf'--reference {REFRENCE_GENOME}')
    # cmd.append(rf'--remora-modified-bases dna_r9.4.1_e8 fast 0.0.0 5hmc_5mc CG 0')
    cmd.append(rf'--outputs mod_basecalls')
    cmd.append(rf'--output-directory {FAST5_FILES}results_megalodon')
    cmd.append(rf'--write-mods-text')
    cmd.append(rf'--processes 18')
    cmd.append(rf'--overwrite ')
    return ' \\\n'.join(cmd)

def cmd_guppy():
    cmd = [rf'guppy_basecaller']
    cmd.append(rf'-i {FAST5_FILES}')
    cmd.append(rf'-s {FAST5_FILES}results_guppy')
    cmd.append(rf'-d {RERIO_MODELS}')
    cmd.append(rf'-c {GUPPY_CONFIG}')
    cmd.append(rf'--cpu_threads_per_caller 18')
    return ' \\\n'.join(cmd)

print(cmd_megalodon_extras())
print()
print(cmd_megalodon())
print()
print(cmd_guppy())