#! /usr/bin/env python

import os
import subprocess
import logging

hours = 18
folders = ['none', 'p', 'p8', 'p8r', 'pr', 'r']
queues = ['clk', 'epyc']
executable = "/home/ulc/ba/mfu/code/gnr/bin/gnr"

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
    new_path = '/'.join(path_folders[11:])
    os.chdir(path)
    if os.path.isfile(job_file):
        with open(job_file, "r") as f:
            last_job = int(f.read().strip())
            print(f"\n{blue}Last submitted job is {last_job} in {new_path}{reset_format}")
            logging.info(f"Last submitted job is {last_job} in {new_path}")
    else:
        print(f"\n{red}{job_file} does not exist{reset_format}")
        logging.error(f"{job_file} does not exist")
        exit()
else:
    user_input = input("Submit jobs in current folder? (y/n): ")
    if user_input.lower() == 'n':
        exit()
    last_job = job_min
    path = os.getcwd()
    path_folders = path.split('/')
    new_path = '/'.join(path_folders[8:])
    logging.info(f"Submitting jobs in {new_path}")
    with open(folder_file, "w") as f:
        f.write(path)

for queue in queues:
    print(f"{bold}{cyan}\n{queue}:{reset_format}")
    output = subprocess.check_output(f"squeue -t RUNNING,PENDING -r -o '%j' | grep -E '^{queue}' | wc -l", shell=True)
    num_jobs_in_queue = int(output.decode().strip())
    print(f"{blue}{num_jobs_in_queue} jobs in queue{reset_format}")
    maxsubmit = get_qos_max_submit(queue)
    available_slots = maxsubmit - num_jobs_in_queue 
    print(f"{blue}{available_slots} slots available{reset_format}")

    while available_slots > 0:

        if last_job == job_max:
            os.remove(job_file)
            folder_name = os.path.basename(path)
            folder_index = folders.index(folder_name)
            changed_dir = False
            for folder in folders[folder_index + 1:]:
                next_folder = os.path.join('../', folder)
                if os.path.isdir(next_folder):
                    os.chdir(next_folder)
                    changed_dir = True
                    last_job = job_min
                    path = os.getcwd()
                    path_folders = path.split('/')
                    new_path = '/'.join(path_folders[11:])
                    print(f"{blue}Moving to folder {new_path}{reset_format}")
                    logging.info(f"Moving to folder {new_path}")
                    with open(job_file, 'w') as f:
                        f.write(str(job_min))
                    with open(folder_file, 'w') as f:
                        f.write(path)
                    break
            if not changed_dir:
                print(f"{bold}{yellow}All jobs completed{reset_format}")
                logging.info("All jobs completed")
                print(f"{blue}{available_slots} slots available{reset_format}\n")
                os.remove(folder_file)
                exit()

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
               "--mail-user=marcelinofuentes@gmail.com",
               "--array", job_array,
               "--wrap", f"srun {executable} ${{SLURM_ARRAY_TASK_ID}}"]
        print(f"{blue}Submitting jobs {first_job} to {last_job}{reset_format}")
        logging.info(f"Submitting jobs {first_job} to {last_job}")
        result = subprocess.run(cmd, stdout=subprocess.PIPE)
        print(result.stdout.decode().strip())
        logging.info(result.stdout.decode().strip())
        with open(job_file, "w") as f:
            f.write(str(last_job))
        available_slots -= num_jobs_to_submit 

print("")

