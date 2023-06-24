#! /usr/bin/env python

import os
import subprocess
import logging

# Purpose: resubmit jobs that failed
# Usage: python resubmit.py

hours = 15
queues = ["clk", "epyc"]
executable = "/home/ulc/ba/mfu/code/gnr/bin/gnr"
mail_user = "marcelinofuentes@gmail.com"

log_file = "/home/ulc/ba/mfu/submit.log"
logging.basicConfig(filename=log_file,
                    level=logging.DEBUG,
                    format="%(asctime)s %(levelname)s: %(message)s")
blue = "\033[94m"
cyan = "\033[96m"
red = "\033[91m"
yellow = "\033[33m"
bold = "\033[1m"
reset_format = "\033[0m"

def get_qos_max_submit(queue):
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
        logging.error(f"QOS "{qos_name}" not found")
        exit()
    if hours >= maxwall:
        qos_name = queue + "_medium"
    command = ["sacctmgr", "-p", "show", "qos", "format=name,maxsubmit"]
    output = subprocess.check_output(command).decode().strip()
    for line in output.split("\n"):
        if line.startswith(qos_name):
            fields = line.strip().split("|")
            return int(fields[1])
    print(f"{red}QOS {qos_name} not found{reset_format}")

job_array = []

for file in os.listdir("."):
    if file.endswith(".csv"):
        if sum(1 for line in open(file)) < 10:
            base = os.path.splitext(file)[0]
            job_array.append(base)

if len(job_array) == 0:
    print(f"{yellow}All jobs completed{reset_format}")
    exit()

path = os.getcwd()
path_folders = path.split("/")
new_path = "/".join(path_folders[8:])
logging.info(f"Submitting failed jobs in {new_path}")

for queue in queues:
    print(f"{bold}{cyan}\n{queue}:{reset_format}")
    output = subprocess.check_output(f"squeue -t RUNNING,PENDING -r -o "%j" | grep -E "^{queue}" | wc -l", shell=True)
    num_jobs_in_queue = int(output.decode().strip())
    print(f"{blue}{num_jobs_in_queue} jobs in queue{reset_format}")
    maxsubmit = get_qos_max_submit(queue)
    available_slots = maxsubmit - num_jobs_in_queue 
    print(f"{blue}{available_slots} slots available{reset_format}")

    while available_slots > 0 and len(job_array) > 0:

        num_jobs_to_submit = min(available_slots, len(job_array))
        first_job = job_array[0]
        last_job = job_array[num_jobs_to_submit - 1]
        job_name = f"{queue}-{os.getcwd().split("/")[-1]}"
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
               "--array", ",".join(job_array[:num_jobs_to_submit]),
               "--wrap", f"srun {executable} ${{SLURM_ARRAY_TASK_ID}}"]
        print(f"{blue}Submitting jobs {first_job} to {last_job}{reset_format}")
        logging.info(f"Submitting jobs {first_job} to {last_job}")
        result = subprocess.run(cmd, stdout=subprocess.PIPE)
        print(result.stdout.decode().strip())
        logging.info(result.stdout.decode().strip())
        for base in job_array[:num_jobs_to_submit]:
            for extension in [".csv", ".frq", ".gl2"]:
                file = base + extension
                os.remove(file)
        del job_array[:num_jobs_to_submit]
        available_slots -= num_jobs_to_submit 

if len(job_array) > 0:
    print(f"{red}{len(job_array)} jobs not submitted{reset_format}")
else:
    print(f"{yellow}All jobs submitted{reset_format}")
    print(f"{blue}{available_slots} slots available{reset_format}")

print("")

