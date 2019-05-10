#!/usr/bin/env python3

import os
import re
import sys
import gzip
import atexit
import shutil
import sqlite3
import tempfile
import subprocess

EXPECTED_HEADER = "seqnames\tstart\tend\tstrand\tREF\tALL\tNS\tAN\tAC_A\tAC_C\tAC_G\tAC_T\thg\tpanTro\tgorGor\tponAbe\trheMac\troot\tponAbe-else\thg-pan-gor\thg-panTro\n"

SCRIPT_DIR = os.path.dirname(__file__)
RES_DIR = os.path.join(SCRIPT_DIR, "../results/snptable")

DB_PATH = os.path.join(SCRIPT_DIR, "../results/snpdb.sqlite3")

def get_results():
    res = os.listdir(RES_DIR)
    res = [ x for x in res if re.search("\.log$", x) ]
    res = [ re.sub("\.log$", "", x) for x in res ]

    ## Ensure all tasks have finished
    for id in res:
        with open(os.path.join(RES_DIR, id + ".log")) as logfile:
            log = logfile.readlines()
            if not any([ re.search("Finished", x) for x in log ]):
                sys.exit("Task {} did not finish".format(id))

    ## Check uniqueness
    res = [[id, *id.split("_", 1)] for id in res]
    res = sorted(res, key = lambda x : x[0])
    check_uniq = dict()
    for id, rg, date in res:
        if rg in check_uniq:
            check_uniq[rg] = check_uniq[rg] + 1
            sys.exit("Found duplicated records")
        else:
            check_uniq[rg] = 1

    return res

def db_pipe():
    ## Check sqlite3 exists
    if not shutil.which("sqlite3"):
        sys.exit("sqlite3 program not found")
    ## Delete original DB file
    if os.path.isfile(DB_PATH):
        os.unlink(DB_PATH)

    db_program = [
        "sqlite3",
        DB_PATH,
        "-init",
        os.path.join(SCRIPT_DIR, "build_db.sql")
    ]

    r, w = os.pipe()
    r = open(r, "r")
    w = open(w, "w")
    process = subprocess.Popen(db_program, stdin = r,
                               stdout = sys.stdout,
                               stderr = sys.stderr)
    return (w, process)

def write_pipe(res, pipe):
    ## Write header
    ## FIXME: If we later use CREATE statement, we should not write header
    #pipe.write(EXPECTED_HEADER)

    for id, rg, date in res:
        gzfile = os.path.join(RES_DIR, id+".tsv.gz")
        tbl_input = gzip.open(gzfile, "rt")
        header = tbl_input.readline()                ## Read the header
        if header == "":                             ## In case that file is empty
            continue
        if header != EXPECTED_HEADER:
            print("Header of {} does not match with expected header".format(id))
        while True:
            line = tbl_input.readline()
            if line == "": break
            pipe.write(line)
    return True

def main():
    res = get_results()
    pipe, process = db_pipe()
    #res = res[0:3]
    write_pipe(res, pipe)
    pipe.close()
    process.wait()

if __name__ == "__main__":
    main()

# Local Variables:
# compile-command: "python3 build_db.py"
# End:
