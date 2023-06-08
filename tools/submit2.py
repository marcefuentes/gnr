#! /usr/bin/env python

import os
import subprocess
import logging

hours = 23
folders = ['none', 'p', 'p8', 'p8r', 'pr', 'r', 'r8']
subfolders = ['given000', 'given050', 'given095', 'given100']
queues = ['clk', 'epyc']
executable = "/home/ulc/ba/mfu/code/gnr/bin/gnr"
mail_user = "marcelinofuentes@gmail.com"

job_min = 100
job_max = 541
job_file = "last_submitted_job.tmp"
folder_file = "/home/ulc/ba/mfu/code/gnr/results/active_folder.tmp"
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

if os.path.isfile(folder_file):
    with open(folder_file, "r") as f:
        path = f.read().strip()
    path_folders = path.split('/')
    path_print = '/'.join(path_folders[11:])
    os.chdir(path)
    if os.path.isfile(job_file):
        with open(job_file, "r") as f:
            last_job = int(f.read().strip())
            print(f"\n{blue}Last submitted job is {last_job} in {path_print}{reset_format}")
            logging.info(f"Last submitted job is {last_job} in {path_print}")
    else:
        print(f"\n{red}{job_file} does not exist in {path_print}{reset_format}")
        logging.error(f"{job_file} does not exist in {path_print}")
        exit()
else:
    user_input = input("Submit jobs in current folder? (y/n): ")
    if user_input.lower() == 'n':
        exit()
    last_job = job_min
    path = os.getcwd()
    path_folders = path.split('/')
    path_print = '/'.join(path_folders[11:])
    print(f"{blue}Submitting jobs in {path_print}{reset_format}")
    logging.info(f"Submitting jobs in {path_print}")
    with open(folder_file, "w") as f:
        f.write(path)

folder, subfolder = path_folders[-2], path_folders[-1]
folder_index = folders.index(folder)
subfolder_index = subfolders.index(subfolder)

for queue in queues:
    print(f"{bold}{cyan}\n{queue}:{reset_format}")
    output = subprocess.check_output(f"squeue -t RUNNING,PENDING -r -o '%j' | grep -E '^{queue}' | wc -l", shell=True)
    num_jobs_in_queue = int(output.decode().strip())
    print(f"{blue}{num_jobs_in_queue} jobs in queue{reset_format}")
    maxsubmit = get_qos_max_submit(queue)
    available_slots = maxsubmit - num_jobs_in_queue 
    print(f"{blue}{available_slots} slots available{reset_format}")

    if available_slots:
        should_break = False
        for folder in folders[folder_index:]:
            for subfolder in subfolders[subfolder_index:]:
                path_folders[-2:] = folder, subfolder
                path = '/'.join(path_folders)
                path_print = '/'.join(path_folders[11:])
                if not os.path.isdir(path):
                    print(f"{red}Skipping {path_print}{reset_format}")
                    logging.info(f"Skipping {path_print}")
                else:
                    os.chdir(path)
                    print(f"{blue}Working in {path_print}{reset_format}")
                    logging.info(f"Working in {path_print}")
                    if os.path.isfile(os.path.join(path, str(last_job + 1) + '.csv')):
                        print(f"{red}Found unexpected {path_print}/{str(last_job + 1)}.csv{reset_format}")
                        exit()
                    else:
                        if last_job == job_max and os.path.isfile(job_file):
                            os.remove(job_file)
                            last_job = job_min
                        else:
                            num_jobs_to_submit = min(available_slots, job_max - last_job)
                            job_name = f"{queue}-{os.getcwd().split('/')[-1]}"
                            first_job = last_job + 1
                            last_job = last_job + num_jobs_to_submit
                            job_array = f"{first_job}-{last_job}"
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
                            print(f"{blue}Submitting jobs {first_job} to {last_job}{reset_format}")
                            logging.info(f"Submitting jobs {first_job} to {last_job}")
                            available_slots -= num_jobs_to_submit 
                            if available_slots == 0:
                                with open(job_file, "w") as f:
                                    f.write(str(last_job))
                                with open(folder_file, "w") as f:
                                    f.write(path)
                                should_break = True
                                subfolder_index = subfolders.index(subfolder)
                                folder_index = folders.index(folder)
                                break
                            else:
                                if last_job == job_max:
                                    last_job = job_min
            if should_break:
                break
            if subfolder_index == len(subfolders) - 1:
                subfolder_index = 0
        if not should_break:
            print(f"{bold}{yellow}All jobs submitted{reset_format}")
            logging.info("All jobs submitted")
            print(f"{blue}{available_slots} slots available{reset_format}\n")
            os.remove(folder_file)
            exit()

print("")

