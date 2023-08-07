#! /usr/bin/env python

import subprocess

queues = ["clk", "epyc"]

blue = "\033[94m"
cyan = "\033[96m"
green = "\033[32m"
red = "\033[91m"
yellow = "\033[33m"
bold = "\033[1m"
reset_format = "\033[0m"

for queue in queues:
    command = ["sacctmgr", "-p", "show", "qos", "format=name,maxwall"]
    output = subprocess.check_output(command).decode().strip().split("\n")
    qos_name = queue + "_short"
    for line in output:
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
    output = subprocess.check_output(command).decode().strip().split("\n")
    for line in output:
        if line.startswith(qos_name):
            fields = line.strip().split("|")
            total_slots = int(fields[1])
            break

    command = ["squeue", "-t", "RUNNING,PENDING", "-r", "-o", "%j"]
    grep_command = ["grep", "-E", f"^{queue}"]
    wc_command = ["wc", "-l"]
    output = subprocess.check_output(command + grep_command + wc_command).decode().strip()
    filled_slots = int(output)
    free_slots = total_slots - filled_slots
    print(f"\n{filled_slots} queued jobs, {cyan}{free_slots}{reset_format} free slots in {queue}")
