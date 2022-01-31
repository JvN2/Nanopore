from pyguppyclient import GuppyBasecallerClient, yield_reads
from pyguppy_client_lib import helper_functions

# server_args = ["--log_path", "/home/noort/guppy_server_logs",
#                "--config", "dna_r9.4.1_450bps_fast.cfg",
#                "--port", '-auto']
# helper_functions.run_server(server_args, "/opt/ont/guppy/bin")


config = "dna_r9.4.1_450bps_fast"
read_file = "/media/noort/Data/users/kuijntjes/2022-01-18_12samplemethylationtest/no_sample/20220118_1825_MN30914_ACR467_7759a37d/fast5_pass/barcode04/ACR467_pass_barcode04_6749c10b_0.fast5"

with GuppyBasecallerClient(config_name=config) as client:
    for read in yield_reads(read_file):
        called = client.basecall(read)
        print(read.read_id, called.seq[:50], called.move)