"""B12: CIN signature computation (computeSigs: 'TRUE').

Opt-in: CINSignatureQuantification is GAP-licensed and not bundled. Runs in the active
environment when the package is installed, or in a container when SWGS_TEST_CINSIG_IMAGE
is set (needs apptainer/singularity).
"""
import os
import shutil
import subprocess

import pytest

from helpers import BASE, make_workdir, read_tsv, run_dir, run_snakemake

IMAGE = os.environ.get("SWGS_TEST_CINSIG_IMAGE", "")
RUNTIME = shutil.which("apptainer") or shutil.which("singularity")
INSTALLED = subprocess.run(["Rscript", "-e", 'quit(status = !requireNamespace("CINSignatureQuantification", quietly = TRUE))'],
                           capture_output=True).returncode == 0

pytestmark = [pytest.mark.cinsig,
              pytest.mark.skipif(not (INSTALLED or (IMAGE and RUNTIME)),
                                 reason="needs CINSignatureQuantification or SWGS_TEST_CINSIG_IMAGE with apptainer")]


def test_auto_computes_signatures(tmp_path, bams):
    samples = [{"SAMPLE_ID": s, "file": bams[s]} for s in ("NA12878", "007")]
    make_workdir(tmp_path, samples, **{**BASE, "computeSigs": "TRUE", "cinsigimage": "" if INSTALLED else IMAGE})
    args = [] if INSTALLED else ["--sdm", "apptainer", "--apptainer-args", f"--bind {tmp_path},{os.path.expanduser('~')}"]
    run_snakemake(tmp_path, "auto", *args, cores=8)
    for sample in ("NA12878", "007"):
        folder = run_dir(tmp_path) / "cin_signatures" / sample
        for method in ("DREWS", "MAC"):
            acts = read_tsv(folder / f"test-run_1_{sample}_100kb_{method}_CIN_signature_acts.tsv")
            assert len(acts) == 1
