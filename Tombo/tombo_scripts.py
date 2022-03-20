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


# multi-to single
def cmd_to_single(single=True):
    cmd = [rf'multi_to_single_fast5']
    cmd.append(rf'-i {FAST5_FILES}')
    cmd.append(rf'-s  {FAST5_FILES}tombo/')
    if single: cmd = ' \\\n'.join(cmd)
    return cmd


def cmd_preprocess(single=True):
    cmd = [rf'tombo preprocess annotate_raw_with_fastqs']
    cmd.append(rf'--fast5-basedir {FAST5_FILES}')
    cmd.append(
        rf'--fastq-filenames {FAST5_FILES}guppy/pass/*.fastq')
    cmd.append(rf'--overwrite')
    if single: cmd = ' \\\n'.join(cmd)
    return cmd


# re-squiggle raw reads
def cmd_resquigle():
    cmd = [rf'tombo resquiggle {FAST5_FILES}tombo {REFRENCE_GENOME}']
    # cmd.append(rf'--processes = 4')
    cmd.append(rf'--num-most-common-errors 5')
    return ' \\\n'.join(cmd)


# run modified base detection
def cmd_detect_mods():
    cmd = [rf'tombo detect_modifications alternative_model']
    cmd.append(rf'--fast5-basedirs {FAST5_FILES}')
    cmd.append(rf'--statistics-file-basename cpg_testing')
    cmd.append(rf'--alternate-bases CpG')
    cmd.append(rf'--processes 4')
    return ' \\\n'.join(cmd)


def cmd_guppy():
    cmd = [rf'guppy_basecaller']
    cmd.append(rf'-i {FAST5_FILES}')
    cmd.append(rf'-s {FAST5_FILES}guppy')
    cmd.append(rf'-d {RERIO_MODELS}')
    cmd.append(rf'-c {GUPPY_CONFIG}')
    cmd.append(rf'--cpu_threads_per_caller 18')
    cmd.append(rf'--fast5_out')
    return ' \\\n'.join(cmd)


print(cmd_guppy())
print()
print(cmd_preprocess())
print()
print(cmd_resquigle())
print()
print(cmd_to_single())
