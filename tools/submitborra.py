#! /usr/bin/env python

import os
import subprocess
import logging

hours = 23
givens = ["given000", "given050", "given095", "given100"]
queues = ["clk", "epyc"]
executable = "/home/ulc/ba/mfu/code/gnr/bin/gnr"
mail_user = "marcelinofuentes@gmail.com"
input_file_extension = ".glo"

last_job_file = "/home/marcelino/code/gnr/results/last_submitted_job.tmp"
log_file = "/home/marcelino/code/gnr/results/submit.log"
logging.basicConfig(filename=log_file,
                    level=logging.DEBUG,
                    format="%(asctime)s %(levelname)s: %(message)s")
blue = "\033[94m"
cyan = "\033[96m"
red = "\033[91m"
yellow = "\033[33m"
bold = "\033[1m"
reset_format = "\033[0m"

def get_free_slots(queue):
    return 250

def get_free_slots2(queue):
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
    print(f"{blue}{filled_slots} queued jobs{reset_format}")
    free_slots = total_slots - filled_slots
    print(f"{blue}{free_slots} free slots{reset_format}")
    return free_slots

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

def submit_jobs(finished, given, last_job, free_slots):
    os.chdir(given)
    job_max = get_job_max(given)
    if last_job == 0:
        job_min = get_job_min(given)
    else:
        job_min = last_job + 1
    if os.path.isfile(os.path.join(given, str(job_min) + ".csv")):
        print(f"{red}Found {str(job_min)}.csv. Exiting{reset_format}")
        exit()
    num_jobs_to_submit = min(free_slots, job_max - job_min + 1)
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
    #result = subprocess.run(cmd, stdout=subprocess.PIPE)
    #print(result.stdout.decode().strip())
    #logging.info(result.stdout.decode().strip())
    free_slots -= num_jobs_to_submit
    if last_job == job_max:
        last_job = 0
        mechanism = os.path.dirname(given)
        givens = []
        for item in os.listdir(mechanism):
            item_path = os.path.join(mechanism, item)
            if os.path.isdir(item_path):
                givens.append(item_path)
        givens.sort(key=lambda x: os.path.getctime(x))
        given_index = givens.index(given)
        if given_index + 1 < len(givens):
            given = givens[given_index + 1]
        else:
            mechanisms = []
            parent = os.path.dirname(mechanism)
            for item in os.listdir(parent):
                item_path = os.path.join(parent, item)
                if os.path.isdir(item_path):
                    mechanisms.append(item_path)
            mechanisms.sort(key=lambda x: os.path.getctime(x))
            mechanism_index = mechanisms.index(mechanism)
            if mechanism_index + 1 < len(mechanisms):
                mechanism = mechanisms[mechanism_index + 1]
                givens = []
                for item in os.listdir(mechanism):
                    item_path = os.path.join(mechanism, item)
                    if os.path.isdir(item_path):
                        givens.append(item_path)
                givens.sort(key=lambda x: os.path.getctime(x))
                given = givens[0]
            else:
                finished = 1

    return finished, given, last_job, free_slots

for queue in queues:
    print(f"{bold}{cyan}\n{queue}:{reset_format}")
    free_slots = get_free_slots(queue) 
    finished = 0

    while free_slots and not finished:
        if os.path.isfile(last_job_file):
            with open(last_job_file, "r") as f:
                given, last_job = f.read().strip().split(",")
            os.remove(last_job_file)
            last_job = int(last_job)
        else:
            user_input = input("Submit jobs in current folder? (y/n): ")
            if user_input.lower() == "n":
                exit()
            path = os.getcwd()
            mechanisms = []
            for item in os.listdir(path):
                item_path = os.path.join(path, item)
                if os.path.isdir(item_path):
                    mechanisms.append(item_path)
            mechanisms.sort(key=lambda x: os.path.getctime(x))
            givens = []
            for item in os.listdir(mechanisms[0]):
                item_path = os.path.join(mechanisms[0], item)
                if os.path.isdir(item_path):
                    givens.append(item_path)
            givens.sort(key=lambda x: os.path.getctime(x))
            given = givens[0]
            last_job = 0
        given_folders = given.split("/")
        given_print = "/".join(given_folders[-3:])
        print(f"{blue}Working in {given_print}{reset_format}")
        logging.info(f"Working in {given_print}")
        finished, given, last_job, free_slots = submit_jobs(finished, given, last_job, free_slots)

        if finished:
            print(f"{bold}{yellow}All jobs submitted{reset_format}")
            logging.info("All jobs submitted")
            print(f"{blue}{free_slots} slots free{reset_format}\n")
            exit()
        else:
            with open(last_job_file, "w") as f:
                f.write(f"{given},{last_job}")

print("")

