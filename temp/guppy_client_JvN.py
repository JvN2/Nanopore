from pyguppy_client_lib import helper_functions

# A basecall server requires:
#  * A location to put log files (on your PC)
#  * An initial config file to load
#  * A port to run on
server_args = ["--log_path", "/home/myuser/guppy_server_logs",
               "--config", "dna_r9.4.1_450bps_fast.cfg",
               "--port", "5556"]
# The second argument is the directory where the
# guppy_basecall_server executable is found. Update this as
# appropriate.
helper_functions.run_server(server_args, "/opt/ont/guppy/bin")