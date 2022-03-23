from pyguppy_client_lib import helper_functions

# A basecall server requires:
#  * A location to put log files (on your PC)
#  * An initial config file to load
#  * A port to run on
server_args = ["--log_path", "/media/noort/Data/users/noort/guppy_server_logs",
               "--config", "dna_r9.4.1_450bps_fast.cfg",
               "--port", "5556"]

# The second argument is the directory where the
# guppy_basecall_server executable is found. Update this as
# appropriate.
helper_functions.run_server(server_args, r'/opt/ont/guppy/bin/')


server_args = [ '-p', 'auto', '-l',
 '/media/noort/Data/users/noort/test10//megalodon_data2/guppy_log', '-c',
 '/home/noort/Downloads/rerio-master/basecall_models/res_dna_r941_min_modbases-all-context_v001.cfg', '--quiet']

