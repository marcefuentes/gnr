#! /usr/bin/env python

import os
import subprocess
import logging

# Purpose: resubmit unfinished jobs
# Usage: python resubmit.py

hours = 23
queues = ["clk", "epyc"]
executable = "/home/ulc/ba/mfu/code/gnr/bin/gnr"
mail_user = "marcelinofuentes@gmail.com"
input_file_extension = ".glo"
output_file_extensions = [".csv", ".frq", ".gl2"]

log_file = "/home/ulc/ba/mfu/code/gnr/results/resubmit.log"
logging.basicConfig(filename=log_file,
                    level=logging.DEBUG,
                    format="%(asctime)s %(levelname)s: %(message)s")
blue = "\033[94m"
cyan = "\033[96m"
green = "\033[32m"
red = "\033[91m"
yellow = "\033[33m"
bold = "\033[1m"
reset_format = "\033[0m"


def get_free_slots(queue):
    command = ["sacctmgr", "-p", "show", "qos", "format=name,maxwall"]
    output = subprocess.check_output(command).decode().strip()
    qos_name = queue + "_short"
    for line in output.split("\n"):
        if line.startswith(qos_name):
            fields = line.strip().split("|")
            maxwall = int(fields[1].split(":")[0])
            break
    if maxwall is None:
        print(f"{red}QOS {qos_name} not found{reset_format}")
        logging.error(f"QOS {qos_name} not found")
        exit()
    if hours >= maxwall:
        qos_name = queue + "_medium"

    command = ["sacctmgr", "-p", "show", "qos", "format=name,maxsubmit"]
    output = subprocess.check_output(command).decode().strip()
    for line in output.split("\n"):
        if line.startswith(qos_name):
            fields = line.strip().split("|")
            total_slots = int(fields[1])
            break

    output = subprocess.check_output(f"squeue -t RUNNING,PENDING -r -o '%j' | grep -E '^{queue}' | wc -l",
                                     shell=True)
    filled_slots = int(output.decode().strip())
    free_slots = total_slots - filled_slots
    print(f"\n{filled_slots} queued jobs, {cyan}{free_slots}{reset_format} free slots in {queue}")
    return free_slots

names = [name[:-4] for name in os.listdir() if name.endswith(input_file_extension)]
job_array = []
for name in names:
    if not os.path.isfile(name + output_file_extensions[0]):
        job_array.append(int(name))
    else:
        with open(name + ".csv") as f:
            if sum(1 for line in f) < 10:
                job_array.append(int(name))
                for extension in [".csv", ".frq", ".gl2"]:
                    file = name + extension
                    os.remove(file)

if len(job_array) == 0:
    print(f"{yellow}{bold}All jobs completed{reset_format}")
    exit()

path = os.getcwd()
path_folders = path.split("/")
new_path = "/".join(path_folders[8:])
logging.info(f"Submitting failed jobs in {new_path}")

for queue in queues:
    free_slots = get_free_slots(queue) 
    while free_slots > 0 and len(job_array) > 0:
        num_jobs_to_submit = min(free_slots, len(job_array))
        first_job = job_array[0]
        last_job = job_array[num_jobs_to_submit - 1]
        job_name = f"{queue}-{last_job}"
        job_time = f"{hours}:59:00"
        cmd = ["sbatch",
               "--job-name", job_name,
               "--output", f"{job_name}.%j.out",
               "--constraint", queue,
               "--nodes=1",
               "--tasks=1",
               "--time", job_time,
               "--mem=4MB",
               "--mail-type=begin,end",
               "--mail-user", mail_user,
               "--array", ",".join(map(str, job_array[:num_jobs_to_submit])),
               "--wrap", f"srun {executable} ${{SLURM_ARRAY_TASK_ID}}"]
        output = subprocess.run(cmd, stdout=subprocess.PIPE, text=True).stdout.strip()
        print(output)
        logging.info(output)
        print(job_array[:num_jobs_to_submit])
        logging.info(f"{given_print}/{job_array[0]}-{job_array[num_jobs_to_submit - 1]} to {queue}")
        del job_array[:num_jobs_to_submit]
        free_slots -= num_jobs_to_submit 

if len(job_array) > 0:
    print(f"{red}{len(job_array)} jobs not submitted{reset_format}")
else:
    print(f"{bold}{green}All jobs submitted{reset_format}")
    print(f"{cyan}{free_slots}{reset_format} free slots\n")

print("")

