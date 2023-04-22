#! /usr/bin/env python

import os
import subprocess
import sys

# Get the queue name from the command line argument
queue = sys.argv[1]

# Define the batch limit based on the queue
batchlim = 400 if queue == "clk" else 200

# Find the number of running and queued jobs in the queue
output = subprocess.check_output(f"squeue -t RUNNING,PENDING -r -o '%j' | grep -E '^{queue}' | wc -l", shell=True)
num_jobs_in_queue = int(output.decode().strip())
available_slots = batchlim - num_jobs_in_queue
print(f"{num_jobs_in_queue} jobs in {queue}")
print(f"{available_slots} slots available in {queue}")

# Find last job submitted
last_job_submitted = 100
if os.path.isfile("last_job_submitted.txt"):
    with open("last_job_submitted.txt", "r") as f:
        last_job_submitted = int(f.read().strip())

# Submit jobs
num_jobs_to_submit = min(available_slots, 541 - last_job_submitted)
if num_jobs_to_submit > 0:
    job_name = f"{queue}-{os.getcwd().split('/')[-1]}"
    first_job = last_job_submitted + 1
    last_job = last_job_submitted + num_jobs_to_submit
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

# Update last_job_submitted.txt
if last_job == 541:
    os.remove("last_job_submitted.txt")
else:
    with open("last_job_submitted.txt", "w") as f:
        f.write(str(last_job))

# Print the number of available slots in the batch queue
available_slots = available_slots - num_jobs_to_submit
print(f"{available_slots} slots still available in {queue}")

