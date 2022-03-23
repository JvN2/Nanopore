from pyguppy_client_lib import helper_functions
from pyguppy_client_lib.pyclient import PyGuppyClient

# A basecall server requires:
#  * A location to put log files (on your PC)
#  * An initial config file to load
#  * A port to run on
server_args = ["--log_path", "/home/myuser/guppy_server_logs",
               "--config", "dna_r9.4.1_450bps_fast.cfg",
               "--port", '-auto']
# The second argument is the directory where the
# guppy_basecall_server executable is found. Update this as
# appropriate.

server_args = ['--version']
helper_functions.run_server(server_args, "/opt/ont/guppy/bin")

# helper_functions.get_server_stats(r'/opt/ont/guppy/bin/guppy_basecall_server', 1)


# client = PyGuppyClient(
#     "127.0.0.1:5555",
#     "dna_r9.4.1_450bps_fast"
# )
# client.connect()