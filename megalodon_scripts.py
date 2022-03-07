
config = 'docker'

if config == 'local':
    LOCATION_OF_FAST5_FILES = r'/media/noort/Data/users/noort/test6/'
    PATH_TO_RERIO_MODELS = r'/home/noort/Downloads/rerio-master/basecall_models/'
    REFRENCE_GENOME = r'/media/noort/Data/users/noort/ref_files/combined.fasta'
    OUTPUT_DIR = LOCATION_OF_FAST5_FILES + r'megalodon_results/'
    LOCATION_OF_GUPPY_BIN = r'/opt/ont/guppy/bin/guppy_basecall_server'
    GUPPY_CONFIG = r'res_dna_r941_min_modbases-all-context_v001.cfg'
else:
    LOCATION_OF_FAST5_FILES = r'/home/data/'
    PATH_TO_RERIO_MODELS = r'/home/rerio/basecall_models/'
    REFRENCE_GENOME = r'/home/data/combined.fasta'
    OUTPUT_DIR = LOCATION_OF_FAST5_FILES + r'megalodon_results/'
    LOCATION_OF_GUPPY_BIN = r'/usr/bin/guppy_basecall_server'
    GUPPY_CONFIG = r'dna_r9.4.1_450bps_hac.cfg'

def cmd_megalodon():
    cmd = rf'megalodon {LOCATION_OF_FAST5_FILES} '
    cmd += rf'--guppy-server-path {LOCATION_OF_GUPPY_BIN} '
    cmd += rf'--guppy-config res_dna_r941_min_modbases-all-context_v001.cfg '
    cmd += rf'--guppy-params "-d {PATH_TO_RERIO_MODELS}" '
    cmd += rf'--reference {REFRENCE_GENOME} '
    cmd += rf'--outputs basecalls mod_basecalls '
    cmd += rf'--output-directory {OUTPUT_DIR} '
    cmd += rf'--write-mods-text '
    cmd += rf'--overwrite '
    return cmd

def cmd_guppy():
    cmd = rf'guppy_basecaller'
    cmd += rf' -i {LOCATION_OF_FAST5_FILES}'
    cmd += rf' -s guppy_results'
    cmd += rf' -d {PATH_TO_RERIO_MODELS}'
    cmd += rf' -c {GUPPY_CONFIG}'
    return cmd

print(cmd_megalodon().replace('--','\\\n--'))
print()
print(cmd_guppy().replace(' -',' \\\n-'))