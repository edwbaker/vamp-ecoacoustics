"""Writes reference.csv for the portability checks: for each of the 91 pilot recordings of the Acoustic Index User's
Guide and each plugin, the R function's value and the plugin's value in the 0.3.0 benchmark (first pass, plugins built
from source and run through ReVAMP on Linux), with the recording's md5 so that a check can confirm it has the same file.
    python make_reference.py <aiguide.csv of https://doi.org/10.5519/ijg223e8> <folder of the 91 recordings>"""
import csv
import hashlib
import os
import sys

# benchmark id -> plugin, output
PLUGINS = {"aci_sw": ("aci", "aci"), "aci_se": ("aci-acc", "aci"), "adi": ("adi-acc", "adi"), "aei": ("aei-acc", "aei"),
           "bi": ("bi-acc", "bi"), "sh": ("sh", "sh"), "th": ("th", "th"), "h": ("h", "h"), "ndsi": ("ndsi", "ndsi")}

data, folder = sys.argv[1:3]
v = {}
for r in csv.DictReader(open(data, newline="")):
    if r["run"] == "full_length" and r["rep"] == "1":
        v[(r["file"], r["id"], r["kind"])] = r["value"]
files = sorted({f for f, _, _ in v}, key=str.lower)
md5 = {f: hashlib.md5(open(os.path.join(folder, f), "rb").read()).hexdigest() for f in files}
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "reference.csv")
with open(out, "w", newline="") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(["file", "md5", "id", "plugin", "output", "r_value", "benchmark_value"])
    for f in files:
        for i, (p, o) in PLUGINS.items():
            w.writerow([f, md5[f], i, p, o, v[(f, i, "R")], v[(f, i, "vamp")]])
print(f"{len(files)} recordings, {len(files) * len(PLUGINS)} rows -> {out}")
