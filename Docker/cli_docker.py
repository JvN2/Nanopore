def cmd_docker(data_dir, docker_ref):
    cmd = [f'docker run -it']
    cmd.append(rf'--mount src={data_dir},target=/home/data,type=bind')
    cmd.append(rf'{docker_ref} sh')
    return ' \\\n'.join(cmd)


data_dir = r'/media/noort/Data/users/kuijntjes/2022-08-16_RecsInAmsterdam_BandsOnGel/CP115_CP130_reconstitutions_methylations/20220816_1950_MN30914_AJF795_9344cc69/fast5_pass/barcode09'
docker_ref = 'jvn2/nanopore:0.1'
docker_ref = 'sidizhao/tombo'

print(cmd_docker(data_dir, docker_ref))
