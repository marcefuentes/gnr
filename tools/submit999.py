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
        logging.error(f"QOS '{qos_name}' not found")
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

def submit_jobs_in_folder(path, last_job, available_slots):
    os.chdir(path)
    job_max = get_job_max(path)
    if last_job == job_max:
        print(f"{yellow}No jobs to submit in this folder{reset_format}")
        logging.info(f"No jobs to submit in this folder")
        return 0, available_slots
    if last_job == 0:
        job_min = get_job_min(path)
    else:
        job_min = last_job + 1
    if os.path.isfile(os.path.join(path, str(job_min) + ".csv")):
        print(f"{red}Found {str(job_min)}.csv. Skipping this folder{reset_format}")
        return 0, available_slots
    num_jobs_to_submit = min(available_slots, job_max - job_min + 1)
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
    available_slots -= num_jobs_to_submit
    if available_slots and last_job == job_max:
        last_job = 0
    return last_job, available_slots

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

path_list = "/".join(path_folders[:-2])
folders = os.listdir(path_list)
folders.sort()
folder = path_folders[-2]
folder_index = folders.index(folder)

path_list = "/".join(path_folders[:-1])
subfolders = os.listdir(path_list)
subfolders.sort()
subfolder = path_folders[-1]
subfolder_index = subfolders.index(subfolder)

if last_job:
    path_print = "/".join(path_folders[-3:])
    print(f"{blue}\nLast job submitted: {path_print}/{last_job}{reset_format}")

for queue in queues:
    print(f"{bold}{cyan}\n{queue}:{reset_format}")
    output = subprocess.check_output(f"squeue -t RUNNING,PENDING -r -o '%j' | grep -E '^{queue}' | wc -l",
                                     shell=True)
    num_jobs_in_queue = int(output.decode().strip())
    print(f"{blue}{num_jobs_in_queue} jobs in queue{reset_format}")
    maxsubmit = get_qos_max_submit(queue)
    available_slots = maxsubmit - num_jobs_in_queue 
    print(f"{blue}{available_slots} slots available{reset_format}")

    while available_slots and folder_index < len(folders):
        path_folders[-2] = folders[folder_index]
        path_list = '/'.join(path_folders[:-1])
        subfolders = os.listdir(path_list)
        while available_slots and subfolder_index < len(subfolders):
            path_folders[-1] = subfolders[subfolder_index]
            path_print = '/'.join(path_folders[-3:])
            print(f"{blue}Working in {path_print}{reset_format}")
            logging.info(f"Working in {path_print}")
            path = '/'.join(path_folders)
            last_job, available_slots = submit_jobs_in_folder(path, last_job, available_slots)
            if available_slots:
                subfolder_index += 1
        if subfolder_index == len(subfolders):
            subfolder_index = 0
        folder_index += 1

    if available_slots and folder_index == len(folders):
        print(f"{bold}{yellow}All jobs submitted{reset_format}")
        logging.info("All jobs submitted")
        print(f"{blue}{available_slots} slots available{reset_format}\n")
        os.remove(last_job_file)
        exit()
    else:
        with open(last_job_file, "w") as f:
            f.write(f"{path},{last_job}")

print("")

