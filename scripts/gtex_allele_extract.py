import os
import sys

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

#def read_gtex_file():
with open(GTEx_FILE, "r", newline="") as gtex_file:
    gtex_file
    #for line in gtex_file:


print("====")
list(db.get_snp_batch([("chr2", 43)]))
print("====")
list(db.get_snp_batch([("chr10", 10190)]))

print("""\n\n\n""")
for each in db.get_rsidnum_batch([1,2,3]):
    print("=======")
    print(each)

