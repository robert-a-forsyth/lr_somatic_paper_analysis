#!/usr/bin/env python3
"""
Nextflow pipeline disk usage report generator.

Inputs (resolved relative to this script):
  - ../Run/data.csv                              — samplesheet
  - ../Run/out/pipeline_info/execution_trace_*.txt — nextflow trace (latest if multiple)
  - ../Run/work/<hash>/                          — nextflow work directories
  - ../Run/out/<sample>/                         — nextflow output directories

Output:
  - disk_usage_table.tsv   — 5-column human-readable table
"""

import csv
import glob
import re
import subprocess
import sys
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────────
SCRIPT_DIR   = Path(__file__).resolve().parent
RUN_DIR      = SCRIPT_DIR.parent / "Run"
SAMPLESHEET  = RUN_DIR / "data.csv"
WORK_DIR     = RUN_DIR / "work"
OUT_DIR      = RUN_DIR / "out"
FINAL_REPORT = SCRIPT_DIR / "disk_usage_table.tsv"

trace_matches = sorted(glob.glob(str(OUT_DIR / "pipeline_info" / "execution_trace_*.txt")))
if not trace_matches:
    sys.exit(f"ERROR: no execution trace found in {OUT_DIR / 'pipeline_info'}")
TRACE_FILE = Path(trace_matches[-1])   # use most recent if multiple
print(f"Using trace: {TRACE_FILE}")


# ── Helpers ────────────────────────────────────────────────────────────────────
def human_readable(n_bytes):
    if n_bytes is None:
        return "N/A"
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n_bytes < 1024:
            return f"{n_bytes:.2f} {unit}"
        n_bytes /= 1024
    return f"{n_bytes:.2f} PB"


def du_bytes(path):
    try:
        result = subprocess.run(
            ["du", "-sb", str(path)],
            capture_output=True, text=True, timeout=120,
        )
        if result.returncode == 0:
            return int(result.stdout.split()[0])
    except Exception as e:
        print(f"  WARNING: du failed on {path}: {e}", file=sys.stderr)
    return None


def file_size_bytes(path):
    try:
        return Path(path).stat().st_size
    except Exception as e:
        print(f"  WARNING: cannot stat {path}: {e}", file=sys.stderr)
        return None


def resolve_work_dir(hash_str):
    parts = hash_str.split("/")
    if len(parts) != 2:
        return None
    prefix, suffix = parts
    matches = glob.glob(str(WORK_DIR / prefix / (suffix + "*")))
    return matches[0] if matches else None


def display_name(sample: str) -> str:
    tumor_only = sample.endswith("-TO")
    base = sample[:-3] if tumor_only else sample
    if base.endswith("-PB"):
        base = base[:-3] + "-PacBio"
    return f"{base} ({'tumor-only' if tumor_only else 'paired'})"


# ── Step 1: Parse samplesheet ──────────────────────────────────────────────────
print("Reading samplesheet...")
samples = {}
sample_order = []
with open(SAMPLESHEET, newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        name = row["sample"].strip()
        samples[name] = {
            "bam_tumor":  row["bam_tumor"].strip()  or None,
            "bam_normal": row["bam_normal"].strip() or None,
        }
        sample_order.append(name)

known_samples = set(samples)
print(f"  Found {len(samples)} samples: {sorted(known_samples)}")


# ── Step 2: Parse trace and accumulate work-dir bytes per sample ───────────────
print(f"\nParsing execution trace...")
SAMPLE_RE = re.compile(r"\(([^)]+)\)$")
sample_work_bytes = {s: 0 for s in known_samples}

tasks = []
with open(TRACE_FILE, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        if row["status"] not in ("COMPLETED", "CACHED"):
            continue
        tasks.append(row)

print(f"  {len(tasks)} COMPLETED/CACHED tasks found.")

for i, row in enumerate(tasks, 1):
    hash_str = row["hash"]
    name     = row["name"]

    m = SAMPLE_RE.search(name)
    tag = m.group(1) if m else ""
    sample = tag if tag in known_samples else None
    if sample is None:
        continue

    work_path = resolve_work_dir(hash_str)
    if work_path and Path(work_path).is_dir():
        size = du_bytes(work_path)
        if size is not None:
            sample_work_bytes[sample] += size
    elif work_path is None:
        print(f"  WARNING: cannot resolve work dir for hash {hash_str!r}", file=sys.stderr)

    if i % 50 == 0 or i == len(tasks):
        print(f"  [{i}/{len(tasks)}] processed tasks...", flush=True)

if not any(b > 0 for b in sample_work_bytes.values()):
    sys.exit(
        "ERROR: execution trace has no completed tasks matching any sample "
        "in data.csv — refusing to emit an empty table."
    )


# ── Step 3: Get input BAM file sizes ──────────────────────────────────────────
print("\nMeasuring input BAM file sizes...")
for sname, info in samples.items():
    info["tumor_bam_bytes"]  = file_size_bytes(info["bam_tumor"])  if info["bam_tumor"]  else None
    info["normal_bam_bytes"] = file_size_bytes(info["bam_normal"]) if info["bam_normal"] else None
    t = human_readable(info["tumor_bam_bytes"])
    n = human_readable(info["normal_bam_bytes"]) if info["bam_normal"] else "N/A"
    print(f"  {sname}: tumor={t}, normal={n}")


# ── Step 4: Get output directory sizes ────────────────────────────────────────
print("\nMeasuring output directory sizes...")
for sname in samples:
    out_path = OUT_DIR / sname
    if out_path.is_dir():
        print(f"  Measuring {out_path} ...", flush=True)
        samples[sname]["out_dir_bytes"] = du_bytes(out_path)
    else:
        print(f"  WARNING: output dir not found: {out_path}", file=sys.stderr)
        samples[sname]["out_dir_bytes"] = None
    print(f"    → {human_readable(samples[sname]['out_dir_bytes'])}")


# ── Step 5: Determine output row order ────────────────────────────────────────
emitted = set()
order = []
for name in sample_order:
    if name in emitted:
        continue
    order.append(name)
    emitted.add(name)
    to_name = name + "-TO"
    if not name.endswith("-TO") and to_name in samples and to_name not in emitted:
        order.append(to_name)
        emitted.add(to_name)


# ── Step 6: Write table ────────────────────────────────────────────────────────
print(f"\nWriting report to {FINAL_REPORT}...")
OUT_COLS = [
    "Sample Name",
    "Tumor input BAM size",
    "Normal input BAM size",
    "Work Directory Disk Usage",
    "Output Directory Disk Usage",
]
with open(FINAL_REPORT, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=OUT_COLS, delimiter="\t")
    writer.writeheader()
    for sname in order:
        info = samples[sname]
        wb = sample_work_bytes.get(sname, 0)
        writer.writerow({
            "Sample Name":                 display_name(sname),
            "Tumor input BAM size":        human_readable(info["tumor_bam_bytes"]),
            "Normal input BAM size":       human_readable(info["normal_bam_bytes"]) if info["bam_normal"] else "N/A",
            "Work Directory Disk Usage":   human_readable(wb) if wb else "N/A",
            "Output Directory Disk Usage": human_readable(info["out_dir_bytes"]),
        })

print(f"  {len(order)} rows written.")
print("Done.")
