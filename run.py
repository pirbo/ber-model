#!/usr/local/bin/python

import os
import shutil
from string import Template
import argparse as ap
import simplejson as js
import random
from copy import deepcopy
import string
import synth
import Queue
import threading
import subprocess

random.seed(124357)

def merge(d1, d2):
    d3 = dict()
    for (k, v) in d1.items(): d3[k] = v
    for (k, v) in d2.items(): d3[k] = v
    return d3

def orc(f): #open read close
    f = open(f, "r")
    ret = f.read()
    f.close()
    return ret

def package(l, n):
    for i in range(0, len(l), n): yield(tuple(l[i:i+n]))

from pprint import *
def prepare_fields_for_sim(data):
    cstr = data["constraints"]
    defaults = cstr["defaults"]
    vol = defaults["volume"]
    av =  6.022e23
    kD_constraints = cstr["kD"]
    for k in data.keys():
        if data[k] == None: data[k] = defaults["k_on"]
        if type(data[k]) in [str, unicode] and data[k] in defaults: data[k] = defaults[data[k]]
    for k_on, k_off, v in kD_constraints:
        if k_on in data: data[k_off] = (v  * av * vol) *  data[k_on]
        else:
            data[k_on] = defaults["k_on"]
            data[k_off] = data[k_on] / (v  * 6.022e23 * vol)
    pprint(data)

def fix_xtra(xtra):
    return dict(package([string.lstrip(s, "-") for s in xtra], 2))

def gen_experiment(data):
    constants = data["constants"]
    fields = data["fields"]
    expr = dict()
    todo = []
    for (k, v) in fields.items():
        typ, val = v["type"], v["value"]
        if typ == "constant": expr[k] = val
        if typ == "random_uniform":
            if isinstance(val, str) and val in constants: val = constants[val]
            if len(val) < 3: val.insert(0, "float")
            t, lower, upper = val
            t = {"float": float, "int": int}[t]
            expr[k] = t(random.uniform(lower, upper))
        else: todo.append((k, v))
    for (k, v) in todo:
        typ, val = v["type"], v["value"]
        if typ == "Kd_constraint":
            field, ratio = val
            expr[k] = expr[field]/(ratio * 6.022e23 * constants["volume"])
    return expr

def synth_dna(args, xtra):
    settings = js.loads(orc(os.path.join("src", "default_synth_dna.json")))
    meta = settings["meta"]
    data = settings["data"]
    if args.descriptor:
        loaded = js.loads(orc(args.descriptor))
        meta = merge(meta, loaded["meta"])
        data = merge(data, loaded["data"])
    meta = merge(meta, xtra)
    data = merge(meta, data)
    
    if not os.path.exists(meta["temp_dir"]): os.makedirs(meta["temp_dir"])
    if not os.path.exists(meta["out_dir"]): os.makedirs(meta["out_dir"])
    
    output = os.path.join(meta["out_dir"], meta["out"])
    data["out"] = output

    raw = "\n\n".join([synth.gen_strand(data["bp_per_strand"], data["weights"]) for i in range(data["strands"]) ])
    with open(output, "w", 0) as f: f.write(raw)
    

def damage_dna(opts, xtra): pass #TODO
def simulate_dna(opts, xtra):
    settings = js.loads(orc(os.path.join("src", "default_simulate_dna.json")))
    meta = settings["meta"]
    data = settings["data"]
    if args.descriptor:
        loaded = js.loads(orc(args.descriptor))
        meta = merge(meta, loaded["meta"])
        data = merge(data, loaded["data"])
    meta = merge(meta, xtra)
    if args.sanity: meta["sanity"] = True
    data = merge(meta, data)
    if not os.path.exists(meta["temp_dir"]): os.makedirs(meta["temp_dir"])
    if not os.path.exists(meta["out_dir"]): os.makedirs(meta["out_dir"])
    
    output = meta["out"] + ".out"
    parsed = os.path.join(meta["temp_dir"], meta["out"] + ".ka")

    prepare_fields_for_sim(data)
    
    src = orc(os.path.join("src", "initial.ka")) + "\n" + \
          "\n".join([orc(os.path.join("src", s + ".ka")) for s in meta["model"]]) + \
          "\n" + orc(os.path.join("src", "observables.ka"))
    if "sanity" in meta and meta["sanity"]:
        src = src + "\n" + orc(os.path.join("src", "sanity.ka"))
        
    src = Template(src)
    with open(parsed, "w", 0) as f:
        f.write(src.safe_substitute(data))

    print ("Running simulation with")
    if "time" in meta:
        com = "KaSim -i src/sig/dna.ka -i {} -i {} -t {} -p {} -o {} -d {}".format(
            parsed, meta["dna_file"], meta["time"], meta["plot_points"], output, meta["out_dir"])
    else:
        com = "KaSim -i src/sig/dna.ka -i {} -i {} -e {} -p {} -o {} -d {}".format(
            parsed, meta["dna_file"], meta["events"], meta["plot_points"], output, meta["out_dir"])
    if "seed" in meta: com = com + " -seed {}".format(meta["seed"])
    print (com)
    os.system(com)

def worker(q):
    while True:
        item = q.get()
        com = "python run.py simulate -descriptor {}".format(item)
        print "Running \n{}".format(com)
        print threading.active_count()
        print threading.current_thread()
        subprocess.call(com, shell=True)
        q.task_done()

def experiment_dna(args, xtra):
    settings = js.loads(orc(args.descriptor))
    meta = settings["meta"]
    data = settings["data"]
    meta = merge(meta, xtra)
    data = merge(meta, data)
    direct = os.path.join("experiments", meta["experiment_name"])
    if not os.path.exists(direct): os.makedirs(direct)

    q = Queue.Queue()
    for i in range(meta["runs"]):
        d = gen_experiment(data)
        m = deepcopy(meta)
        m["out_dir"] = direct
        m["out"] = "{}_{}".format(meta["experiment_name"], i)
        s = {"meta": m, "data": d}
        descriptor = os.path.join(direct, m["out"]+".json")
        q.put(descriptor)
        f = open(descriptor, "w", 0)
        js.dump(s, f, sort_keys=True, indent=4)
        f.close()
    for i in range(meta["parallel"]):
        t = threading.Thread(target=worker, args=[q])
        t.daemon = True
        t.start()
    q.join()
        

parser = ap.ArgumentParser(description = "Task runner for DNA BER pathway simulation")
subparsers = parser.add_subparsers(dest = "type")

parser_synth = subparsers.add_parser("synth", help = "synthesis of DNA strands")
parser_damage = subparsers.add_parser("damage", help = "simulate damage of DNA data")
parser_simulate = subparsers.add_parser("simulate", help = "simulate BER pathway")
parser_simulate.add_argument("-sanity", action='store_true', default=False)
parser_experiment = subparsers.add_parser("experiment", help = "generates/runs series of experiments")

for p in [parser_synth, parser_damage, parser_simulate]:
    p.add_argument("-descriptor", help = "mode descriptor", default=None)

parser_experiment.add_argument("descriptor", help = "experiment description")

args, xtra = parser.parse_known_args()
xtra = fix_xtra(xtra)
{
    "synth": synth_dna,
    "damage": damage_dna,
    "simulate": simulate_dna,
    "experiment": experiment_dna
}[args.type](args, xtra)
