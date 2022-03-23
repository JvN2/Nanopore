import subprocess

exec_file = r'/opt/ont/guppy/bin/guppy_basecaller'
inputfile = r'D/media/noort/Data/users/kuijntjes/2022-01-18_12samplemethylationtest/no_sample/20220118_1825_MN30914_ACR467_7759a37d/fast5_pass/barcode04/ACR467_pass_barcode04_6749c10b_0.fast5'
savefile = r'/media/noort/Data/users/kuijntjes/2022-01-18_12samplemethylationtest/no_sample/20220118_1825_MN30914_ACR467_7759a37d/fast5_pass/barcode04/ACR467_pass_barcode04_6749c10b_0'
flowcell = r'FLO-FLG001'
kit = r'SQK-LSK109'

# cmd = f'{exec_file} --input_path {inputfile} --save_path  {savefile} --flowcell {flowcell} --kit {kit}'
# cmd = ['which', 'guppy_basecall_server']
# p = subprocess.Popen(cmd, stdout=subprocess.PIPE,cwd=r'/usr/bin/', shell=True)
# print(p.communicate()[0].decode())

#Follow this protocol:
# https://community.nanoporetech.com/protocols/Guppy-protocol/v/gpb_2003_v1_revab_14dec2018/run-guppy-on-linux
# Call this line from command prompt:
# guppy_basecaller --input_path /media/noort/Data/users/noort/2022-01-18_12samplemethylationtest/no_sample/20220118_1825_MN30914_ACR467_7759a37d/fast5_pass/barcode04 --save_path /media/noort/Data/users/noort/2022-01-18_12samplemethylationtest/no_sample/20220118_1825_MN30914_ACR467_7759a37d/fast5_pass/barcode04/test --flowcell FLO-FLG001 --kit SQK-LSK109