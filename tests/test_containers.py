"""B11: auto workflow inside the published container (--sdm apptainer).

Opt-in: skipped unless apptainer or singularity is installed. The image defaults to the
config's workflowimage and can be overridden with SWGS_TEST_IMAGE.
"""
import os
import shutil

import pytest

from helpers import BASE, REPO, make_workdir, read_config, run, run_snakemake, sample_paths

RUNTIME = shutil.which("apptainer") or shutil.which("singularity")

pytestmark = [pytest.mark.container,
              pytest.mark.skipif(not RUNTIME, reason="apptainer/singularity not installed")]


def image_uri():
    image = os.environ.get("SWGS_TEST_IMAGE") or read_config_default()["workflowimage"]
    return image if "://" in image else f"docker://{image}"


def read_config_default():
    import yaml
    return yaml.safe_load((REPO / "config" / "default_config.yaml").read_text())


def test_image_has_current_package():
    out = run([RUNTIME, "exec", image_uri(), "Rscript", "-e",
               'cat(as.character(packageVersion("rswgsabsolutecn")))']).stdout.strip()
    assert tuple(int(x) for x in out.split(".")[:3]) >= (1, 1, 1), out


def test_auto_in_container(tmp_path, bams):
    make_workdir(tmp_path, [{"SAMPLE_ID": "NA12878", "file": bams["NA12878"]}],
                 **{**BASE, "workflowimage": image_uri()})
    real = os.path.dirname(os.path.realpath(bams["NA12878"]))
    binds = ",".join(sorted({str(tmp_path), real, os.path.dirname(bams["NA12878"])}))
    run_snakemake(tmp_path, "auto", "--sdm", "apptainer", "--apptainer-args", f"--bind {binds}", cores=8)
    log = (tmp_path / "logs" / "autofit_test-run_1_100kb_NA12878.log").read_text()
    assert "namespace" not in log and "not available" not in log
    assert sample_paths(tmp_path, "NA12878")["abs_tsv"].exists()
    assert read_config(tmp_path)["workflowimage"] == image_uri()
