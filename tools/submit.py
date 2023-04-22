#! /usr/bin/env python

import os
import subprocess

queues = ['epyc', 'clk']
maxsubmits = [200, 400]
job_min = 100
job_max = 541
job_file = "last_job_submitted.tmp"
folder_file = "/home/ulc/ba/mfu/code/gnr/results/active_folder.tmp"

if os.path.isfile(folder_file):
    with open(folder_file, "r") as f:
        path = f.read().strip()
    os.chdir(path)
    with open(job_file, "r") as f:
        last_job = int(f.read().strip())
else:
    last_job = job_min
    cwd = os.getcwd()
    with open(folder_file, "w") as f:
        f.write(cwd)

for queue, maxsubmit in zip(queues, maxsubmits):

    print(f"\n\033[96m{queue}:\033[0m")
    output = subprocess.check_output(f"squeue -t RUNNING,PENDING -r -o '%j' | grep -E '^{queue}' | wc -l", shell=True)
    num_jobs_in_queue = int(output.decode().strip())
    print(f"{num_jobs_in_queue} jobs in queue")
    available_slots = maxsubmit - num_jobs_in_queue
    print(f"{available_slots} slots available")
    num_jobs_to_submit = min(available_slots, job_max - last_job)

    if num_jobs_to_submit > 0:
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
    else:
        print("No jobs to submit")

    available_slots = available_slots - num_jobs_to_submit
    print(f"{available_slots} slots still available")

print("")

if last_job == job_max:
    os.remove(job_file)
    os.remove(folder_file)
else:
    with open(job_file, "w") as f:
        f.write(str(last_job))
