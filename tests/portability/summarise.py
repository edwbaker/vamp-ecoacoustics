"""Combines the portability checks' results (results-*.csv, one per host) into one table, and prints a Markdown summary:
for each host, the values outside the tolerance of R and the largest difference from the benchmark build.
    python summarise.py <folder of results-*.csv> <combined csv>
Exits 1 if any value is outside the tolerance."""
import csv
import glob
import os
import sys

COLS = ["host", "file", "plugin", "output", "value", "r_value", "benchmark_value", "diff_r", "diff_benchmark", "ok"]
rows = []
for f in sorted(glob.glob(os.path.join(sys.argv[1], "**", "results-*.csv"), recursive=True)):
    rows += list(csv.DictReader(open(f, newline="")))
with open(sys.argv[2], "w", newline="") as fh:
    w = csv.DictWriter(fh, COLS, extrasaction="ignore", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)

num = lambda s: float(s) if s not in ("", "NA", None) else float("inf")
passed = lambda r: r["ok"].upper() == "TRUE"
print("## Release portability\n")
print("| host | values | outside the tolerance of R | largest difference from the benchmark build | largest relative difference from it |")
print("|---|---|---|---|---|")
for h in dict.fromkeys(r["host"] for r in rows):
    x = [r for r in rows if r["host"] == h]
    bench = max(abs(num(r["diff_benchmark"])) for r in x)
    rel = max(abs(num(r["diff_benchmark"])) / abs(float(r["benchmark_value"])) for r in x if float(r["benchmark_value"]))
    print(f"| {h} | {len(x)} | {sum(not passed(r) for r in x)} | {bench:.3g} | {rel:.3g} |")
sys.exit(0 if rows and all(passed(r) for r in rows) else 1)
