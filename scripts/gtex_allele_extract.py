import os
import re
import sys
import itertools

try:
    SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
except NameError:
    if os.path.exists("scripts/gtex_allele_extract.py"):
        SCRIPT_DIR = os.path.join(".", "scripts")
    else:
        raise Exception("Can not locate script dir")
PROJECT_ROOT = os.path.join(SCRIPT_DIR, "..")
DB_PATH = os.path.join(PROJECT_ROOT, "results/snpdb.sqlite3")

sys.path.insert(0, os.path.join(PROJECT_ROOT, "lib"))

from snpdb import SNPDB

db = SNPDB(DB_PATH)

GTEx_FILE = os.path.join(PROJECT_ROOT, "raw_data/GTEx_v7_Binaryfile_4919.txt")

def read_gtex_file():
    with open(GTEx_FILE, "r", newline="") as gtex_file:
        headers = next(gtex_file).split()
        for line in gtex_file:
            rsid = line.split()[0]
            rsid_num = re.sub("rs", "", rsid)
            yield rsid_num

RES_FILE = os.join.path(PROJECT_ROOT, "results/GTEx_v7_Binaryfile_4919_ancestral_allele.txt")

def main():
    res_file = open(RES_FILE, "w")
    firstline = True
    for row in db.get_rsidnum_batch(read_gtex_file()):
        if firstline:
            res_file.write("\t".join(row.keys()))
            res_file.write("\n")
        line = (str(v) if v is not None else "NA" for v in tuple(row))
        line = "\t".join(line)
        res_file.write(line)
        res_file.write("\n")

if __name__ == "__main__":
    main()

# for each in db.get_rsidnum_batch(itertools.islice(read_gtex_file(), 10)):
#     print(each)

