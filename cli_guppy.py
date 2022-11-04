def cmd_docker(data_dir, target_dir=r'/home/data'):
    #image = rf'jvn2/nanopore:0.1'
    image = rf'sidizhao/tombo'
    #data_dir = rf'/media/noort/Data/users/noort/20220816_1950_MN30914_AJF795_9344cc69'
    # cmd = [rf'docker run --name box -it {image} /bin/bash']
    cmd = [rf'docker run -it --mount src={data_dir},target={target_dir},type=bind {image} sh']
    return ' \\\n'.join(cmd)


def cmd_guppy(fast5_dir, guppy_config, rerio_models):
    cmd = [rf'guppy_basecaller']
    cmd.append(rf'-i {fast5_dir}')
    cmd.append(rf'-s {fast5_dir}/guppy')
    cmd.append(rf'-d {rerio_models}')
    cmd.append(rf'-c {guppy_config}')
    cmd.append(rf'--cpu_threads_per_caller 18')
    cmd.append(rf'--fast5_out')
    cmd.append(rf'-r')
    return ' \\\n'.join(cmd)


if __name__ == '__main__':
    config = 'docker'
    if config == 'local':
        fast5_dir = r'/home/noort/data/test21'
        rerio_models = r'/home/noort/Downloads/rerio-master/basecall_models/'
        guppy_config = r'res_dna_r941_min_modbases-all-context_v001.cfg'
        # guppy_bin = r'/opt/ont/guppy/bin/guppy_basecall_server'
    else:
        fast5_dir = r'/home/kuijntjes/Desktop/2022-10-12_WholeCellExtractGalLocusCindy/no_sample/20221012_1615_MN30914_AJA380_b85344d5/fast5_pass/barcode01'
        rerio_models = r'/home/rerio/basecall_models/'
        guppy_config = r'res_dna_r941_min_modbases-all-context_v001.cfg'
        # guppy_bin = r'/usr/bin/guppy_basecall_server'

    print(cmd_docker(fast5_dir), '\n')
    print(cmd_guppy('/home/data', guppy_config, rerio_models), '\n')
