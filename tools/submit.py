#! /usr/bin/env python

import os
import subprocess
import logging


folders = ['none', 'p', 'p8', 'p8r', 'pr', 'r']
queues = ['epyc', 'clk']
maxsubmits = [200, 250]
job_min = 100
job_max = 541
job_file = "last_submitted_job.tmp"
folder_file = "/home/ulc/ba/mfu/code/gnr/results/active_folder.tmp"
slurm_file = "/home/ulc/ba/mfu/code/gnr/results/job.sh"
log_file = "/home/ulc/ba/mfu/submit.log"
logging.basicConfig(filename=log_file,
                    level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s: %(message)s')

if os.path.isfile(folder_file):
    with open(folder_file, "r") as f:
        path = f.read().strip()
    path_folders = path.split('/')
    new_path = '/'.join(path_folders[11:])
    print(f"\nActive folder is {new_path}")
    logging.info(f"Active folder is {new_path}")
    os.chdir(path)
    if os.path.isfile(job_file):
        with open(job_file, "r") as f:
            last_job = int(f.read().strip())
            print(f"Last submitted job is {last_job}")
            logging.info(f"Last submitted job is {last_job}")
    else:
        print(f"{job_file} does not exist")
        logging.error(f"{job_file} does not exist")
        exit()
else:
    user_input = input("There is no active folder. Continue? (y/n): ")
    if user_input.lower() == 'n':
        exit()
    last_job = job_min
    path = os.getcwd()
    path_folders = path.split('/')
    new_path = '/'.join(path_folders[8:])
    print(f"\nActive folder is {new_path}")
    logging.info(f"Active folder is {new_path}")
    with open(folder_file, "w") as f:
        f.write(path)

for queue, maxsubmit in zip(queues, maxsubmits):

    print(f"\n\033[96m{queue}:\033[0m")
    logging.info(f"{queue}:")
    output = subprocess.check_output(f"squeue -t RUNNING,PENDING -r -o '%j' | grep -E '^{queue}' | wc -l", shell=True)
    num_jobs_in_queue = int(output.decode().strip())
    print(f"{num_jobs_in_queue} jobs in queue")
    logging.info(f"{num_jobs_in_queue} jobs in queue")
    available_slots = maxsubmit - num_jobs_in_queue 
    print(f"{available_slots} slots available (estimate)")
    logging.info(f"{available_slots} slots available (estimate)")

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
                    print(f"Moving to folder {new_path}")
                    logging.info(f"Moving to folder {new_path}")
                    with open(job_file, 'w') as f:
                        f.write(str(job_min))
                    with open(folder_file, 'w') as f:
                        f.write(path)
                    break
            if not changed_dir:
                print("All jobs completed")
                logging.info("All jobs completed")
                os.remove(folder_file)
                exit()

        num_jobs_to_submit = min(available_slots, job_max - last_job)
        job_name = f"{queue}-{os.getcwd().split('/')[-1]}"
        first_job = last_job + 1
        last_job = last_job + num_jobs_to_submit
        array = f"{first_job}-{last_job}"
        subprocess.run(["sbatch",
                        "--job-name", job_name,
                        "--output", f"{job_name}.%j.out",
                        "-C", queue,
                        "--array", array,
                        slurm_file])
        with open(job_file, "w") as f:
            f.write(str(last_job))
        print(f"with jobs {first_job} to {last_job}")
        logging.info(f"with jobs {first_job} to {last_job}")

        output = subprocess.check_output(f"squeue -t RUNNING,PENDING -r -o '%j' | grep -E '^{queue}' | wc -l", shell=True)
        num_jobs_in_queue = int(output.decode().strip())
        print(f"{num_jobs_in_queue} jobs in queue")
        logging.info(f"{num_jobs_in_queue} jobs in queue")
        available_slots = maxsubmit - num_jobs_in_queue
        print(f"{available_slots} slots available (estimate)")
        logging.info(f"{available_slots} slots available (estimate)")

print("")

