#! /usr/bin/env python

import os
import subprocess

queues = ['epyc', 'clk']
batchlims = [200, 400]
job_min = 100
job_max = 541
filename = "last_job_submitted.tmp"

last_job = job_min
if os.path.isfile(filename):
    with open(filename, "r") as f:
        last_job = int(f.read().strip())

for queue, batchlim in zip(queues, batchlims):

    output = subprocess.check_output(f"squeue -t RUNNING,PENDING -r -o '%j' | grep -E '^{queue}' | wc -l", shell=True)
    num_jobs_in_queue = int(output.decode().strip())
    print(f"{num_jobs_in_queue} jobs in {queue}")
    available_slots = batchlim - num_jobs_in_queue
    print(f"{available_slots} slots available in {queue}")

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
    print(f"{available_slots} slots still available in {queue}")

if last_job == job_max:
    os.remove(filename)
else:
    with open(filename, "w") as f:
        f.write(str(last_job))
