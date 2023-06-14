#! /usr/bin/env python

import os
import subprocess
import logging

hours = 23
queues = ["clk", "epyc"]
executable = "/home/ulc/ba/mfu/code/gnr/bin/gnr"
mail_user = "marcelinofuentes@gmail.com"
input_file_extension = ".glo"

next_job_file = "/home/ulc/ba/mfu/code/gnr/results/last_submitted_job.tmp"
log_file = "/home/ulc/ba/mfu/submit.log"
logging.basicConfig(filename=log_file,
                    level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s: %(message)s')
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

def submit_jobs_in_folder(path, job_min):
    os.chdir(path)
    job_max = get_job_max(path)
    if job_min >= job_max:
        print(f"{yellow}No jobs to submit in this folder{reset_format}")
        logging.info(f"No jobs to submit in this folder")
        job_min = 0
        return job_min, available_slots
    if job_min == 0:
        job_min = get_job_min(path)
    if os.path.isfile(os.path.join(path, str(job_min) + '.csv')):
        print(f"{red}Found unexpected {str(job_min)}.csv{reset_format}")
        exit()
    num_jobs_to_submit = min(available_slots, job_max - job_min + 1)
    job_last = job_min + num_jobs_to_submit - 1
    job_name = f"{queue}-{os.getcwd().split('/')[-1]}"
    job_array = f"{job_min}-{job_last}"
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
    print(f"{blue}Submitting jobs {job_min} to {job_last}{reset_format}")
    logging.info(f"Submitting jobs {job_min} to {job_last}")
    result = subprocess.run(cmd, stdout=subprocess.PIPE)
    print(result.stdout.decode().strip())
    logging.info(result.stdout.decode().strip())
    available_slots -= num_jobs_to_submit
    job_min = job_last + 1
    return job_min, available_slots

for queue in queues:
    print(f"{bold}{cyan}\n{queue}:{reset_format}")
    output = subprocess.check_output(f"squeue -t RUNNING,PENDING -r -o '%j' | grep -E '^{queue}' | wc -l",
                                     shell=True)
    num_jobs_in_queue = int(output.decode().strip())
    print(f"{blue}{num_jobs_in_queue} jobs in queue{reset_format}")
    maxsubmit = get_qos_max_submit(queue)
    available_slots = maxsubmit - num_jobs_in_queue 
    print(f"{blue}{available_slots} slots available{reset_format}")

    if available_slots:
        if os.path.isfile(next_job_file):
            with open(next_job_file, "r") as f:
                path, job_min = f.read().strip().split(",")
        else:
            user_input = input("Submit jobs in current folder? (y/n): ")
            if user_input.lower() == 'n':
                exit()
            path = os.getcwd()
            job_min = get_job_min(path)

        path_folders = path.split('/')
        path_list = '/'.join(path_folders[:-2])
        folders = os.listdir(path_list)
        folder, subfolder = path_folders[-2], path_folders[-1]
        folder_index = 0
        while available_slots and folder_index < len(folders):
            path_list = '/'.join(path_folders[:-1])
            subfolders = os.listdir(path_list)
            subfolder_index = 0
            while available_slots and subfolder_index < len(subfolders):
                path_folders[-2:] = folders[folder_index],
                                    subfolders[subfolder_index]
                path = '/'.join(path_folders)
                path_print = '/'.join(path_folders[-3:])
                print(f"{blue}Working in {path_print}{reset_format}")
                logging.info(f"Working in {path_print}")
                job_min, available_slots = submit_jobs_in_folder(path, job_min)
                subfolder_index += 1
                if available_slots:
                    job_min = 0
            folder_index += 1
            subfolder_index = 0

        if folder_index == len(folders):
            print(f"{bold}{yellow}All jobs submitted{reset_format}")
            logging.info("All jobs submitted")
            print(f"{blue}{available_slots} slots available{reset_format}\n")
            os.remove(next_job_file)
            exit()
        else:
            with open(next_job_file, "w") as f:
                f.write(f"{path},{job_min}")

print("")

