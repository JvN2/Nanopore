config = 'docker'

if config == 'local':
    LOCATION_OF_FAST5_FILES = r'/media/noort/Data/users/noort/test6/'
    PATH_TO_RERIO_MODELS = r'/home/noort/Downloads/rerio-master/basecall_models/'
    REFRENCE_GENOME = r'/media/noort/Data/users/noort/ref_files/combined.fasta'
    LOCATION_OF_GUPPY_BIN = r'/opt/ont/guppy/bin/guppy_basecall_server'
    GUPPY_CONFIG = r'res_dna_r941_min_modbases-all-context_v001.cfg'
else:
    LOCATION_OF_FAST5_FILES = r'/home/data/'
    PATH_TO_RERIO_MODELS = r'/home/rerio/basecall_models/'
    REFRENCE_GENOME = r'/home/data/combined.fasta'
    LOCATION_OF_GUPPY_BIN = r'/usr/bin/guppy_basecall_server'
    GUPPY_CONFIG = r'dna_r9.4.1_450bps_hac.cfg'

def cmd_megalodon():
    cmd = [rf'megalodon {LOCATION_OF_FAST5_FILES}']
    cmd.append(rf'--guppy-server-path {LOCATION_OF_GUPPY_BIN}')
    cmd.append(rf'--guppy-config res_dna_r941_min_modbases-all-context_v001.cfg ')
    cmd.append(rf'--guppy-params "-d {PATH_TO_RERIO_MODELS}"')
    cmd.append(rf'--reference {REFRENCE_GENOME}')
    cmd.append(rf'--outputs basecalls mod_basecalls')
    cmd.append(rf'--output-directory {LOCATION_OF_FAST5_FILES}megalodon_results')
    cmd.append(rf'--write-mods-text')
    cmd.append(rf'--overwrite ')
    return ' \\\n'.join(cmd)

def cmd_guppy():
    cmd = [rf'guppy_basecaller']
    cmd.append(rf'-i {LOCATION_OF_FAST5_FILES}')
    cmd.append(rf'-s {LOCATION_OF_FAST5_FILES}megalodon_results')
    cmd.append(rf'-d {PATH_TO_RERIO_MODELS}')
    cmd.append(rf'-c {GUPPY_CONFIG}')
    return ' \\\n'.join(cmd)

print(cmd_megalodon())
print()
print(cmd_guppy())