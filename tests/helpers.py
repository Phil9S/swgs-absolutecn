"""Shared helpers for the swgs-absolutecn pytest suite.

Every workflow run uses its own working directory: the Snakefiles read
config/config.yaml and write logs/ relative to the working directory, and the
R helper scripts read config/config.yaml from it as well.
"""
import csv
import os
import re
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
import yaml

REPO = Path(__file__).resolve().parents[1]
WORKFLOW = REPO / "workflow"
ARGO = WORKFLOW / "argo-auto-workflow.yaml"
SHEET_COLUMNS = ["PATIENT_ID", "SAMPLE_ID", "TP53freq", "smooth", "file", "precPloidy", "precPurity"]
CORES = os.environ.get("SWGS_TEST_CORES", "4")
DROP = object()  # override value that removes a config key

# Config shared by the integration runs (100 kb bins keep runs to about a minute)
BASE = {"bins": 100, "project_name": "test-run_1", "out_dir": "results/"}
# 5 Mb homozygous deletion carried by the DEL_1 test BAM (GRCh37 contig names)
DELETION = ("4", 100_000_000, 105_000_000)


def make_workdir(workdir, samples, drop_columns=(), **overrides):
    """Write config/config.yaml (default config plus overrides) and config/sample_sheet.tsv.

    samples: list of dicts with at least SAMPLE_ID and file. Missing columns default to
    PATIENT_ID P1, TP53freq/precPloidy/precPurity NA and smooth FALSE.
    """
    workdir = Path(workdir)
    (workdir / "config").mkdir(parents=True, exist_ok=True)
    config = yaml.safe_load((REPO / "config" / "default_config.yaml").read_text())
    config["samplesheet"] = "config/sample_sheet.tsv"
    for key, value in overrides.items():
        if value is DROP:
            config.pop(key, None)
        else:
            config[key] = value
    (workdir / "config" / "config.yaml").write_text(yaml.safe_dump(config, sort_keys=False))

    rows = []
    for sample in samples:
        row = {"PATIENT_ID": "P1", "TP53freq": "NA", "smooth": "FALSE", "precPloidy": "NA", "precPurity": "NA"}
        row.update(sample)
        rows.append(row)
    columns = [c for c in SHEET_COLUMNS if c not in drop_columns]
    pd.DataFrame(rows, columns=columns).to_csv(workdir / "config" / "sample_sheet.tsv", sep="\t", index=False)
    return workdir


def read_config(workdir):
    return yaml.safe_load((Path(workdir) / "config" / "config.yaml").read_text())


def run_dir(workdir):
    """Pipeline output directory for the project/bin set in the workdir config."""
    config = read_config(workdir)
    return Path(workdir, config["out_dir"], "sWGS_fitting", f"{config['project_name']}_{config['bins']}kb")


def sample_paths(workdir, sample):
    """Per-sample pipeline files, mirroring the rule outputs."""
    config = read_config(workdir)
    project, bins = config["project_name"], config["bins"]
    base = run_dir(workdir)
    pre = base / "absolute_PRE_down_sampling" / sample
    post = base / "absolute_POST_down_sampling" / sample
    return {
        "link": base / "bams" / f"{sample}.bam",
        "ok": base / "bams" / f"{sample}_bam.ok",
        "clonality": pre / "clonality_results" / f"{project}_{sample}_clonality.tsv",
        "pdf": pre / "clonality_results" / f"{project}_{sample}_clonality.pdf",
        "filtered": pre / f"{project}_{sample}_filtered_results.tsv",
        "fit": pre / f"{project}_{sample}_fit_QC_predownsample.tsv",
        "plot": pre / "plots" / f"{sample}.png",
        "pre_rds": pre / "relative_cn_rds" / f"{project}_{sample}_{bins}kb_relSmoothedCN.rds",
        "fitted": post / f"{project}_{sample}_fitted_QC.tsv",
        "ds_bam": base / "absolute_POST_down_sampling" / "downsampled_bams" / f"{sample}.bam",
        "post_rds": post / "relative_cn_rds" / f"{project}_{sample}_{bins}kb_relSmoothedCN.rds",
        "abs_tsv": post / "abs_cn_rds" / f"{project}_{sample}_{bins}kb_ds_abs_fits.tsv",
        "abs_rds": post / "abs_cn_rds" / f"{project}_{sample}_{bins}kb_ds_absCopyNumber.rds",
        "seg": post / "abs_cn_rds" / f"{project}_{sample}_{bins}kb_ds_absCopyNumber_segTable.tsv",
        "abs_plot": post / "abs_cn_rds" / "plots" / f"{sample}.png",
    }


def relpath(workdir, path):
    """Path relative to the workdir, as Snakemake targets are written."""
    return os.path.relpath(path, workdir)


def run_snakemake(workdir, workflow, *args, ok=True, cores=None):
    """Run a Snakefile in workdir. ok=True/False asserts the outcome, ok=None doesn't check."""
    cmd = ["snakemake", "-s", str(WORKFLOW / workflow), "--cores", str(cores or CORES),
           "--show-failed-logs", *map(str, args)]
    res = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True)
    res.output = res.stdout + res.stderr
    if ok is not None and (res.returncode == 0) != ok:
        tail = "\n".join(res.output.splitlines()[-60:])
        raise AssertionError(f"snakemake {workflow} exit {res.returncode}, expected {'success' if ok else 'failure'}\n"
                             f"cmd: {' '.join(cmd)}\n--- last output ---\n{tail}")
    return res


def job_counts(res):
    """Parse the first 'Job stats' table of a snakemake run into {rule: count}."""
    lines = [line.strip() for line in res.output.splitlines()]
    start = lines.index("Job stats:")
    counts = {}
    for line in lines[start + 3:]:
        parts = line.split()
        if len(parts) != 2 or not parts[1].isdigit():
            break
        counts[parts[0]] = int(parts[1])
        if parts[0] == "total":
            break
    return counts


def read_tsv(path):
    """Read a pipeline TSV the way the R scripts do (no quote or comment characters)."""
    return pd.read_csv(path, sep="\t", quoting=csv.QUOTE_NONE,
                       dtype={"SAMPLE_ID": str, "PATIENT_ID": str, "sample": str, "chromosome": str})


def r_eval(code):
    """Run R code from a temporary file (avoids shell escaping) and return stdout."""
    with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False) as fh:
        fh.write(code)
    try:
        res = subprocess.run(["Rscript", fh.name], capture_output=True, text=True)
    finally:
        os.unlink(fh.name)
    if res.returncode != 0:
        raise AssertionError(f"Rscript failed\n{res.stderr[-3000:]}")
    return res.stdout


def run(cmd, cwd=None, ok=True):
    """Run a command (list) and return the CompletedProcess with .output."""
    res = subprocess.run([str(c) for c in cmd], cwd=cwd, capture_output=True, text=True)
    res.output = res.stdout + res.stderr
    if ok is not None and (res.returncode == 0) != ok:
        raise AssertionError(f"{' '.join(map(str, cmd))} exit {res.returncode}\n{res.output[-3000:]}")
    return res


def samtools_count(bam):
    return int(run(["samtools", "view", "-c", bam]).stdout.strip())


def placeholder_bams(root, names, index="real", ext="bam"):
    """Empty alignment files for dry runs, each a symlink data/<name> -> archive/<name>.

    index: where the index goes - "real" (next to the real file), "stem" (real file with
    .bam replaced by .bai), "link" (next to the symlink only) or "none".
    """
    root = Path(root)
    (root / "archive").mkdir(parents=True, exist_ok=True)
    (root / "data").mkdir(parents=True, exist_ok=True)
    files = {}
    for name in names:
        real = root / "archive" / f"{name}.{ext}"
        link = root / "data" / f"{name}.{ext}"
        real.touch()
        link.symlink_to(real)
        if index == "real":
            Path(f"{real}.bai").touch()
        elif index == "stem":
            (root / "archive" / f"{name}.bai").touch()
        elif index == "link":
            Path(f"{link}.bai").touch()
        files[name] = str(link)
    return files


def argo_template(template):
    """Return (script path, [(flag, value template), ...]) for an Argo template's Rscript call."""
    spec = yaml.safe_load(ARGO.read_text())
    tpl = next(t for t in spec["spec"]["templates"] if t["name"] == template)
    script = tpl["container"]["args"][0]
    rscript = re.search(r"Rscript (\S+)", script).group(1)
    return rscript, re.findall(r'--(\w+) "([^"]*)"', script)


def argo_args(template, values):
    """Build the Rscript argument list for an Argo template, filling {{inputs/workflow.parameters.X}}."""
    _, pairs = argo_template(template)
    args = []
    for flag, value in pairs:
        filled = re.sub(r"\{\{(?:inputs|workflow)\.parameters\.(\w+)\}\}", lambda m: str(values[m.group(1)]), value)
        assert "{{" not in filled, f"unresolved Argo parameter in --{flag}: {value}"
        args += [f"--{flag}", filled]
    return args
