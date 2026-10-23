"""B8-B9: precomputed ploidy/purity - processPrecomputed.R + stage_2, and partial values."""
import os

import pytest

from helpers import BASE, REPO, make_workdir, read_tsv, relpath, run, run_dir, run_snakemake, sample_paths

pytestmark = pytest.mark.integration

PRECOMPUTED = REPO / "workflow" / "scripts" / "processPrecomputed.R"


@pytest.fixture(scope="module")
def precomputed(tmp_path_factory, bams):
    workdir = tmp_path_factory.mktemp("precomputed")
    samples = [{"SAMPLE_ID": s, "file": bams[s], "precPloidy": 2.0, "precPurity": 1.0} for s in ("NA12878", "007")]
    make_workdir(workdir, samples, **BASE)
    run(["Rscript", PRECOMPUTED], cwd=workdir)
    run_snakemake(workdir, "stage_2", cores=8)
    return {"workdir": workdir, "samples": samples}


def test_precomputed_markers_written(precomputed):
    pre = run_dir(precomputed["workdir"]) / "absolute_PRE_down_sampling"
    assert (pre / "GENERATED_USING_PRECOMPUTED_PL_PU_RDS_FILES_ARE_EMPTY").exists()
    assert (pre / "precomp.ok").exists()


def test_stage2_uses_precomputed_values(precomputed):
    combined = read_tsv(run_dir(precomputed["workdir"]) / "absolute_POST_down_sampling" / "test-run_1_100kb_ds_abs_fits.tsv")
    assert sorted(combined["SAMPLE_ID"]) == ["007", "NA12878"]
    assert (combined["ploidy"] == 2.0).all() and (combined["purity"] == 1.0).all()


def test_precomputed_links_point_to_real_files(precomputed):
    for sample in precomputed["samples"]:
        link = sample_paths(precomputed["workdir"], sample["SAMPLE_ID"])["link"]
        assert os.readlink(link) == os.path.realpath(sample["file"])


def test_missing_bam_is_named(tmp_path, bams):
    samples = [{"SAMPLE_ID": "NA12878", "file": bams["NA12878"], "precPloidy": 2.0, "precPurity": 1.0},
               {"SAMPLE_ID": "GHOST", "file": str(tmp_path / "missing.bam"), "precPloidy": 2.0, "precPurity": 1.0}]
    make_workdir(tmp_path, samples, **BASE)
    res = run(["Rscript", PRECOMPUTED], cwd=tmp_path, ok=False)
    assert "BAM file(s) not found for sample(s): GHOST" in res.output


def test_partial_precomputed_fixes_one_dimension(tmp_path, bams):
    samples = [{"SAMPLE_ID": "NA12878", "file": bams["NA12878"], "precPloidy": 2.0},
               {"SAMPLE_ID": "007", "file": bams["007"], "precPurity": 0.9}]
    make_workdir(tmp_path, samples, **BASE)
    targets = [relpath(tmp_path, sample_paths(tmp_path, s)["clonality"]) for s in ("NA12878", "007")]
    run_snakemake(tmp_path, "auto", *targets, cores=8)
    fixed_ploidy = read_tsv(sample_paths(tmp_path, "NA12878")["clonality"])
    fixed_purity = read_tsv(sample_paths(tmp_path, "007")["clonality"])
    assert set(fixed_ploidy["ploidy"].round(2)) == {2.0} and fixed_ploidy["purity"].nunique() > 1
    assert set(fixed_purity["purity"].round(2)) == {0.9} and fixed_purity["ploidy"].nunique() > 1
