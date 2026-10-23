"""Group A: config, sample sheet and DAG checks as dry runs (no test data, no R jobs run).

Placeholder alignment files are empty files behind symlinks; only Snakemake's parsing,
the checks in workflow/rules/common.smk and the job graph are exercised.
"""
import os
import shutil
import subprocess

import pytest

from helpers import DROP, job_counts, make_workdir, placeholder_bams, run_snakemake

pytestmark = pytest.mark.dryrun

SAMPLES = ("S1", "007")
STAGE_1 = {"bam_check": 2, "sym_link": 2, "relRDS": 2, "gridsearch": 2, "combine_output": 1, "all": 1, "total": 10}
AUTO = {"bam_check": 2, "sym_link": 2, "relRDS": 2, "gridsearch": 2, "autofit": 2, "downsample": 2,
        "ds_relRDS": 2, "rel_to_abs": 2, "all": 1, "total": 17}
STAGE_2 = {"ds_symlink": 2, "split_input": 2, "downsample": 2, "ds_relRDS": 2, "rel_to_abs": 2,
           "combine_results": 1, "all": 1, "total": 12}


def setup(tmp_path, names=SAMPLES, index="real", ext="bam", sample_values=None, drop_columns=(), **overrides):
    files = placeholder_bams(tmp_path, names, index=index, ext=ext)
    samples = [{"SAMPLE_ID": n, "file": files[n], **(sample_values or {}).get(n, {})} for n in names]
    if ext == "cram":
        ref = tmp_path / "archive" / "ref.fa"
        ref.touch()
        overrides = {"filetype": "CRAM", "reference": str(ref), **overrides}
    make_workdir(tmp_path, samples, drop_columns=drop_columns, **overrides)
    return files


def fake_stage1(tmp_path, names=SAMPLES, out_dir="results/", qc_rows=None, project="example_run", bins=30):
    """Write the curated stage_1 QC table and per-sample pre-downsample RDS placeholders."""
    pre = tmp_path / out_dir / "sWGS_fitting" / f"{project}_{bins}kb" / "absolute_PRE_down_sampling"
    for n in names:
        rds = pre / n / "relative_cn_rds"
        rds.mkdir(parents=True, exist_ok=True)
        (rds / f"{project}_{n}_{bins}kb_relSmoothedCN.rds").touch()
    if qc_rows is None:
        qc_rows = [f"{n}\tP1\t2\t1\tTRUE" for n in names]
    if qc_rows is not False:
        (pre / f"{project}_{bins}kb_combined_QC_predownsample.tsv").write_text(
            "SAMPLE_ID\tPATIENT_ID\tploidy\tpurity\tuse\n" + "".join(r + "\n" for r in qc_rows))


def r_has(package):
    res = subprocess.run(["Rscript", "-e", f'quit(status = !requireNamespace("{package}", quietly = TRUE))'])
    return res.returncode == 0


# ---- workflows build the expected job graph --------------------------------

def test_bam_stage1_and_auto_job_graph(tmp_path):
    setup(tmp_path)
    assert job_counts(run_snakemake(tmp_path, "stage_1", "-n")) == STAGE_1
    assert job_counts(run_snakemake(tmp_path, "auto", "-n")) == AUTO


def test_cram_job_graph_declares_index(tmp_path):
    setup(tmp_path, ext="cram")
    cram_stage1 = {k.replace("sym_link", "cram_to_bam"): v for k, v in STAGE_1.items()}
    cram_auto = {k.replace("sym_link", "cram_to_bam"): v for k, v in AUTO.items()}
    assert job_counts(run_snakemake(tmp_path, "stage_1", "-n")) == cram_stage1
    res = run_snakemake(tmp_path, "auto", "-n", "-p")
    assert job_counts(res) == cram_auto
    assert "bams/S1.bam, " in res.output and "bams/S1.bam.bai" in res.output
    assert "samtools index results/sWGS_fitting/example_run_30kb/bams/S1.bam " \
           "results/sWGS_fitting/example_run_30kb/bams/S1.bam.bai" in res.output


def test_stage2_finds_qc_file_without_trailing_slash(tmp_path):
    setup(tmp_path, out_dir="results")
    fake_stage1(tmp_path, out_dir="results")
    assert job_counts(run_snakemake(tmp_path, "stage_2", "-n")) == STAGE_2


def test_signatures_add_compute_jobs(tmp_path):
    setup(tmp_path, computeSigs="TRUE", cinsigimage="phil9s/cinsig:test")
    fake_stage1(tmp_path)
    assert job_counts(run_snakemake(tmp_path, "auto", "-n"))["compute_signatures"] == 2
    stage2 = job_counts(run_snakemake(tmp_path, "stage_2", "-n"))
    assert stage2["compute_signatures"] == 1 and stage2["total"] == 13


@pytest.mark.parametrize("image", ["phil9s/swgs-absolutecn:latest", "docker://phil9s/swgs-absolutecn:latest"])
def test_workflowimage_with_or_without_prefix(tmp_path, image):
    setup(tmp_path, workflowimage=image)
    run_snakemake(tmp_path, "auto", "-n")
    if shutil.which("apptainer") or shutil.which("singularity"):
        run_snakemake(tmp_path, "auto", "-n", "--sdm", "apptainer")


# ---- BAM index check (index must sit next to the real file) -----------------

@pytest.mark.parametrize("index", ["real", "stem"])
def test_index_next_to_real_file_passes(tmp_path, index):
    setup(tmp_path, index=index)
    res = run_snakemake(tmp_path, "stage_1", "-n", "-p")
    # sym_link links straight to the real file
    assert f"ln -s {os.path.realpath(tmp_path)}/archive/S1.bam " in res.output


@pytest.mark.parametrize("index", ["link", "none"])
def test_missing_index_stops_early(tmp_path, index):
    setup(tmp_path, index=index)
    for workflow in ("stage_1", "auto"):
        res = run_snakemake(tmp_path, workflow, "-n", ok=False)
        assert "BAM index not found for sample(s): S1, 007" in res.output


# ---- sample sheet and config errors -----------------------------------------

def test_duplicate_sample_ids(tmp_path):
    setup(tmp_path, names=("S1",))
    sheet = tmp_path / "config" / "sample_sheet.tsv"
    lines = sheet.read_text().splitlines()
    sheet.write_text("\n".join(lines + lines[1:]) + "\n")
    res = run_snakemake(tmp_path, "stage_1", "-n", ok=False)
    assert "Sample sheet contains duplicated sample ids" in res.output


def test_missing_precomputed_columns(tmp_path):
    setup(tmp_path, drop_columns=("precPloidy",))
    res = run_snakemake(tmp_path, "stage_1", "-n", ok=False)
    assert "missing required column(s): precPloidy" in res.output


def test_cram_requires_reference(tmp_path):
    setup(tmp_path, ext="cram", reference="")
    res = run_snakemake(tmp_path, "stage_1", "-n", ok=False)
    assert "CRAM input requires 'reference'" in res.output


@pytest.mark.parametrize("overrides,message", [
    ({"ploidy_min": 5.0, "ploidy_max": 2.0}, "Minimum ploidy exceeds maximum ploidy"),
    ({"purity_min": 0.9, "purity_max": 0.5}, "Minimum purity exceeds maximum purity"),
])
def test_range_limits(tmp_path, overrides, message):
    setup(tmp_path, **overrides)
    assert message in run_snakemake(tmp_path, "stage_1", "-n", ok=False).output


@pytest.mark.parametrize("overrides", [{"bins": 20}, {"fitMethod": "svm"}, {"workflowimage": DROP}])
def test_config_schema_rejects(tmp_path, overrides):
    setup(tmp_path, **overrides)
    res = run_snakemake(tmp_path, "auto", "-n", ok=False)
    assert "validat" in res.output.lower()


@pytest.mark.parametrize("workflow", ["stage_1", "auto"])
def test_all_precomputed_refuses_gridsearch(tmp_path, workflow):
    values = {n: {"precPloidy": 2.0, "precPurity": 1.0} for n in SAMPLES}
    setup(tmp_path, sample_values=values)
    res = run_snakemake(tmp_path, workflow, "-n", ok=False)
    assert "all ploidy and purity values are precomputed" in res.output


@pytest.mark.parametrize("qc_rows,message", [
    (False, "QC file from stage 1 is missing"),
    ([], "QC file is empty"),
    (["S1\tP1\t2\t1\tTRUE", "007\tP1\t2\t1\tNA"], "'use' column contains non-boolean values"),
    # pandas would read True/False as boolean, but the R scripts don't
    (["S1\tP1\t2\t1\tTrue", "007\tP1\t2\t1\tFalse"], "'use' column contains non-boolean values"),
    (["S1\tP1\t2\t1\tTRUE", "S1\tP1\t3\t0.8\tTRUE", "007\tP1\t2\t1\tTRUE"], "duplicated sample ids marked 'TRUE'"),
])
def test_stage2_qc_file_checks(tmp_path, qc_rows, message):
    setup(tmp_path)
    fake_stage1(tmp_path, qc_rows=qc_rows)
    assert message in run_snakemake(tmp_path, "stage_2", "-n", ok=False).output


def test_signatures_without_package_or_image(tmp_path):
    if r_has("CINSignatureQuantification"):
        pytest.skip("CINSignatureQuantification is installed")
    setup(tmp_path, computeSigs="TRUE", cinsigimage="")
    res = run_snakemake(tmp_path, "auto", "-n", ok=False)
    assert "CINSignatureQuantification is not installed" in res.output
