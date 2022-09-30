import docker
import os

os.system('echo hello world form os')


os.system('whereis systemctl')
os.system('systemctl start docker')

# client = docker.from_env()
# client.containers.run("hello-world")




# import subprocess
#
# logfile = "/tmp/output.log"
# with open(logfile, "a") as output:
#     subprocess.call("docker -v", shell=True, stdout=output, stderr=output)
#
# with open(logfile) as lines:
#     for l in lines:
#         print(l)

