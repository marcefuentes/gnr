#! /usr/bin/env python

import os
import subprocess
import logging

hours = 23
queues = ["clk", "epyc"]
executable = "/home/ulc/ba/mfu/code/gnr/bin/gnr"
mail_user = "marcelinofuentes@gmail.com"
input_file_extension = ".glo"

last_job_file = "/home/ulc/ba/mfu/code/gnr/results/last_submitted_job.tmp"
log_file = "/home/ulc/ba/mfu/code/gnr/results/submit.log"
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
yesno = f"[{bold}{green}Yes{reset_format}/{bold}{red}No{reset_format}]"

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

def folder_list(path):
    folders = []
    for item in os.listdir(path):
        item_path = os.path.join(path, item)
        if os.path.isdir(item_path):
            folders.append(item_path)
    folders.sort(key=lambda x: os.path.getctime(x))
    return folders

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

def submit_jobs(free_slots, given, last_job):
    os.chdir(given)
    job_max = get_job_max(given)
    if last_job == 0:
        job_min = get_job_min(given)
    else:
        job_min = last_job + 1
    given_folders = given.split("/")
    given_print = "/".join(given_folders[-3:])
    if os.path.isfile(os.path.join(given, str(job_min) + ".csv")):
        print(f"{red}{given_print}/{str(job_min)}.csv already exists{reset_format}")
        exit()
    num_jobs_to_submit = min(free_slots, job_max - job_min + 1)
    last_job = job_min + num_jobs_to_submit - 1
    job_name = f"{queue}-{last_job}"
    job_array = f"{job_min}-{last_job}"
    cmd = ["sbatch",
           "--job-name", job_name,
           "--output", f"{job_name}.%j.out",
           "--constraint", queue,
           "--nodes=1",
           "--tasks=1",
           "--time", f"{hours}:59:00",
           "--mem=4MB",
           "--mail-type=begin,end",
           "--mail-user", mail_user,
           "--array", job_array,
           "--wrap", f"srun {executable} ${{SLURM_ARRAY_TASK_ID}}"]
    output = subprocess.run(cmd, stdout=subprocess.PIPE, text=True).stdout.strip()
    print(output)
    logging.info(output)
    print(f"{given_print}/{job_array}")
    logging.info(f"{given_print}/{job_array} to {queue}")
    free_slots -= num_jobs_to_submit
    if last_job == job_max:
        last_job = 0
        mechanism = os.path.dirname(given)
        givens = folder_list(mechanism)
        given_index = givens.index(given) + 1
        if given_index < len(givens):
            given = givens[given_index]
        else:
            mechanisms = folder_list(os.path.dirname(mechanism))
            mechanism_index = mechanisms.index(mechanism) + 1
            if mechanism_index < len(mechanisms):
                mechanism = mechanisms[mechanism_index]
                givens = folder_list(mechanism)
                given = givens[0]
            else:
                print(f"{bold}{green}All jobs submitted{reset_format}")
                print(f"{cyan}{free_slots}{reset_format} free slots\n")
                exit()

    return free_slots, given, last_job

for queue in queues:
    free_slots = get_free_slots(queue) 
    while free_slots:
        if os.path.isfile(last_job_file):
            with open(last_job_file, "r") as f:
                given, last_job = f.read().strip().split(",")
            os.remove(last_job_file)
            last_job = int(last_job)
        else:
            mechanisms = folder_list(os.getcwd())
            givens = folder_list(mechanisms[0])
            given = givens[0]
            given_folders = given.split("/")
            given_print = "/".join(given_folders[-3:])
            print(f"\n{bold}Submit jobs in {given_print}?{reset_format} {yesno} ", end="")
            user_input = input()
            if user_input.lower() == "n":
                exit()
            last_job = 0
        free_slots, given, last_job = submit_jobs(free_slots, given, last_job)
        with open(last_job_file, "w") as f:
            f.write(f"{given},{last_job}")

print()

