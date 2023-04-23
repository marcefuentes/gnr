#! /usr/bin/env python

import os
import subprocess

folders = ['none', 'p', 'p8', 'p8r', 'pr', 'r']
queues = ['epyc', 'clk']
maxsubmits = [200, 400]
job_min = 100
job_max = 541
job_file = "last_job_submitted.tmp"
folder_file = "/home/ulc/ba/mfu/code/gnr/results/active_folder.tmp"

if os.path.isfile(folder_file):
    with open(folder_file, "r") as f:
        path = f.read().strip()
    print(f"\nActive folder is {path}")
    os.chdir(path)
    if os.path.isfile(job_file):
        with open(job_file, "r") as f:
            last_job = int(f.read().strip())
            print(f"Last job submitted is {last_job}")
    else:
        print("There is no file with last job submitted")
        exit()
else:
    user_input = input("There is no active folder. Continue? (y/n): ")
    if user_input.lower() == 'n':
        exit()
    last_job = job_min
    path = os.getcwd()
    print(f"Active folder is now {path}")
    with open(folder_file, "w") as f:
        f.write(path)

for queue, maxsubmit in zip(queues, maxsubmits):

    print(f"\n\033[96m{queue}:\033[0m")
    output = subprocess.check_output(f"squeue -t RUNNING,PENDING -r -o '%j' | grep -E '^{queue}' | wc -l", shell=True)
    num_jobs_in_queue = int(output.decode().strip())
    print(f"{num_jobs_in_queue} jobs in queue")
    available_slots = maxsubmit - num_jobs_in_queue

    while available_slots > 0:

        print(f"{available_slots} slots available")
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
                    print(f"Active folder is now {path}")
                    with open(folder_file, 'w') as f:
                        f.write(path)
                    break
            if not changed_dir:
                print("All jobs completed")
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
                        "job.sh"])
        print(f"with jobs {first_job} to {last_job}")
        available_slots = available_slots - num_jobs_to_submit
        num_jobs_in_queue = num_jobs_in_queue + num_jobs_to_submit
        print(f"{num_jobs_in_queue} jobs in queue")

    print(f"{available_slots} slots available")

print("")

with open(job_file, "w") as f:
    f.write(str(last_job))
