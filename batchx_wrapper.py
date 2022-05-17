#!/usr/bin/python3

"""
Parse batchx user json for launching armadillo
"""
import json
import os
import subprocess
from sys import exit

print("Starting checks...", flush=True)
## Set environment
defaults = {
    "identity": 90,
    "lendiff": 15,
    "mlen": 100,
    "outputName": "armadillo_data"}
        
## Parse json
with open("/batchx/input/input.json", "r") as inputFile:
    inputJson = inputFile.read()
parsedJson = json.loads(inputJson)
parsedJson = {**defaults, **parsedJson}

os.symlink(parsedJson["reference"], "/tmp/"+os.path.basename(parsedJson["reference"]))
parsedJson["reference"] = "/tmp/"+os.path.basename(parsedJson["reference"])

if "referenceIdx" not in parsedJson:
    print("Indexing reference file", flush = True)
    p1 = subprocess.Popen(f'samtools faidx {parsedJson["reference"]}', shell = True)
else:
    os.symlink(parsedJson["referenceIdx"], "/tmp/"+os.path.basename(parsedJson["referenceIdx"]))
    parsedJson["referenceIdx"] = "/tmp/"+os.path.basename(parsedJson["referenceIdx"])

run_args = ["identity", "lendiff", "mlen", "outputName"]
cmd_run = ''

for attribute, value in parsedJson.items():
    cmd_run += f'--{attribute} {value} ' if attribute in run_args else ''

cmd_final = f'python3 /armadillo-dataprep/data_prep.py -i {parsedJson["blast8"]} -g {parsedJson["reference"]} {cmd_run}'

## Prepare output
outputJson = {}
try:
    p1.wait()
except NameError:
    pass

try:
    print("Everything OK, starting reference building!", flush=True)
    r = subprocess.run(cmd_final, shell = True, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    if r.returncode == 0:
        print("The run was successful! Closing the environment", flush = True)
    else:
        print("Error during the run. Exiting now" , flush = True)
        print(r.stdout)
        print(r.stderr.decode())
        exit(1)

    outputJson["armadillo_data"] = f'/batchx/output/{parsedJson["outputName"]}.tar.gz'
    
except subprocess.CalledProcessError as e:
        print(e)
        exit(e.returncode)

with open('/batchx/output/output.json', 'w+') as json_file:
    json.dump(outputJson, json_file)
