#! /usr/bin/env python

import os
import subprocess
import logging

hours = 23
mechanisms = ["none", "p", "pr", "r", "p8", "pr8", "r8"]
givens = ["given000", "given050", "given095", "given100"]
queues = ["clk", "epyc"]
executable = "/home/ulc/ba/mfu/code/gnr/bin/gnr"
mail_user = "marcelinofuentes@gmail.com"
input_file_extension = ".glo"

last_job_file = "/home/ulc/ba/mfu/code/gnr/results/last_submitted_job.tmp"
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
    print(f"{blue}{filled_slots} queued jobs{reset_format}")
    free_slots = total_slots - filled_slots
    print(f"{blue}{free_slots} free slots{reset_format}")
    return free_slots

def get_job_min(path):
    job_min = 9999
    for file in os.listdir(path):
        if file.endswith(input_file_extension):
            basename = int(file.split(".")[0])
            if basename < job_min:
                job_min = basename
    return job_min

def get_job_max(path):
    job_max = 0
    for file in os.listdir(path):
        if file.endswith(input_file_extension):
            basename = int(file.split(".")[0])
            if basename > job_max:
                job_max = basename
    return job_max

def submit_jobs(path, last_job, free_slots):
    os.chdir(path)
    job_max = get_job_max(path)
    if last_job == job_max:
        print(f"{yellow}No jobs to submit in this folder{reset_format}")
        logging.info(f"No jobs to submit in this folder")
        return 0, free_slots
    if last_job == 0:
        job_min = get_job_min(path)
    else:
        job_min = last_job + 1
    if os.path.isfile(os.path.join(path, str(job_min) + ".csv")):
        print(f"{red}Found {str(job_min)}.csv. Skipping this folder{reset_format}")
        return 0, free_slots
    num_jobs_to_submit = min(free_slots, job_max - job_min + 1)
    last_job = job_min + num_jobs_to_submit - 1
    print(f"{blue}Submitting jobs {job_min} to {last_job}{reset_format}")
    logging.info(f"Submitting jobs {job_min} to {last_job}")
    job_name = f"{queue}-{os.getcwd().split('/')[-1]}"
    job_array = f"{job_min}-{last_job}"
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
           "--array", job_array,
           "--wrap", f"srun {executable} ${{SLURM_ARRAY_TASK_ID}}"]
    result = subprocess.run(cmd, stdout=subprocess.PIPE)
    print(result.stdout.decode().strip())
    logging.info(result.stdout.decode().strip())
    free_slots -= num_jobs_to_submit
    if free_slots and last_job == job_max:
        last_job = 0
    return last_job, free_slots

if os.path.isfile(last_job_file):
    with open(last_job_file, "r") as f:
        path, last_job = f.read().strip().split(",")
    last_job = int(last_job)
else:
    user_input = input("Submit jobs in current folder? (y/n): ")
    if user_input.lower() == "n":
        exit()
    path = os.getcwd()
    last_job = 0
path_folders = path.split("/")
mechanism, given = path_folders[-2], path_folders[-1]
mechanism_index = mechanisms.index(mechanism)
given_index = givens.index(given)

if last_job:
    path_print = "/".join(path_folders[-3:])
    print(f"{blue}\nLast job submitted: {path_print}/{last_job}{reset_format}")

for queue in queues:
    print(f"{bold}{cyan}\n{queue}:{reset_format}")
    free_slots = get_free_slots(queue) 

    while free_slots and mechanism_index < len(mechanisms):
        path_folders[-2] = mechanisms[mechanism_index]
        while free_slots and given_index < len(givens):
            path_folders[-1] = givens[given_index]
            path = '/'.join(path_folders)
            path_print = '/'.join(path_folders[-3:])
            if os.path.isdir(path):
            print(f"{blue}Working in {path_print}{reset_format}")
            logging.info(f"Working in {path_print}")
                last_job, free_slots = submit_jobs(path, last_job, free_slots)
            else:
                print(f"{red}Skipping {path_print}{reset_format}")
                logging.info(f"Skipping {path_print}")
            if free_slots:
                given_index += 1
        if given_index == len(givens):
            given_index = 0
        if free_slots:
            mechanism_index += 1

    if free_slots and mechanism_index == len(mechanisms):
        print(f"{bold}{yellow}All jobs submitted{reset_format}")
        logging.info("All jobs submitted")
        print(f"{blue}{free_slots} slots free{reset_format}\n")
        os.remove(last_job_file)
        exit()
    else:
        with open(last_job_file, "w") as f:
            f.write(f"{path},{last_job}")

print("")

