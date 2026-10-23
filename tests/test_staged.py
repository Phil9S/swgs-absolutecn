"""B7: staged workflow - stage_1, curation of the QC table, stage_2.

Curation reuses resources/generate_QC_columns.R (marks the first fit per sample) and then
adds free-text notes with apostrophes and '#', which the R scripts must read literally.
"""
import csv

import pytest

from helpers import BASE, REPO, make_workdir, read_tsv, run, run_dir, run_snakemake, sample_paths

pytestmark = pytest.mark.integration

SAMPLES = ("NA12878", "007", "DEL_1")


def note(i):
    return f"doesn't fit chr17 - it's row {i} #{i}"


@pytest.fixture(scope="module")
def staged(tmp_path_factory, bams):
    workdir = tmp_path_factory.mktemp("staged")
    samples = [{"SAMPLE_ID": "NA12878", "file": bams["NA12878"]},
               {"SAMPLE_ID": "007", "TP53freq": 0.95, "file": bams["007"]},
               {"SAMPLE_ID": "DEL_1", "PATIENT_ID": "P2", "file": bams["DEL_1"]}]
    # out_dir without a trailing slash (get_ds_bams builds the QC path with Path())
    make_workdir(workdir, samples, **{**BASE, "out_dir": "results"})
    run_snakemake(workdir, "stage_1", cores=8)

    qc_path = run_dir(workdir) / "absolute_PRE_down_sampling" / "test-run_1_100kb_combined_QC_predownsample.tsv"
    before = read_tsv(qc_path)
    run(["Rscript", REPO / "resources" / "generate_QC_columns.R"], cwd=workdir)
    curated = read_tsv(qc_path)
    # DEL_1 keeps empty notes, so combining mixes all-NA and text notes across samples
    curated["notes"] = [None if s == "DEL_1" else note(i) for i, s in enumerate(curated["SAMPLE_ID"])]
    # write logicals as TRUE/FALSE - R's read.table doesn't read pandas' True/False as logical
    saved = curated.copy()
    for col in saved.columns:
        if saved[col].dtype == bool:
            saved[col] = saved[col].map({True: "TRUE", False: "FALSE"})
    saved.to_csv(qc_path, sep="\t", index=False, quoting=csv.QUOTE_NONE, na_rep="NA")

    run_snakemake(workdir, "stage_2", cores=8)
    return {"workdir": workdir, "before": before, "curated": curated}


def test_stage1_qc_table_needs_curation(staged):
    before = staged["before"]
    assert set(before["SAMPLE_ID"]) == set(SAMPLES)
    assert before["use"].isna().all()


def test_stage2_uses_the_curated_fits(staged):
    combined = read_tsv(run_dir(staged["workdir"]) / "absolute_POST_down_sampling" / "test-run_1_100kb_ds_abs_fits.tsv")
    assert sorted(combined["SAMPLE_ID"]) == sorted(SAMPLES)
    chosen = staged["curated"][staged["curated"]["use"] == True]  # noqa: E712
    for _, row in chosen.iterrows():
        out = combined[combined["SAMPLE_ID"] == row["SAMPLE_ID"]].iloc[0]
        assert (out["ploidy"], out["purity"]) == pytest.approx((row["ploidy"], row["purity"]))


def test_combined_outputs_cover_all_samples(staged):
    post = run_dir(staged["workdir"]) / "absolute_POST_down_sampling"
    seg = read_tsv(post / "test-run_1_100kb_ds_absCopyNumber_segTable.tsv")
    assert set(seg["sample"]) == set(SAMPLES)
    assert (post / "test-run_1_100kb_ds_absCopyNumber.rds").exists()


def test_notes_with_apostrophes_and_hash_survive(staged):
    curated = staged["curated"]
    for sample in SAMPLES:
        split = read_tsv(sample_paths(staged["workdir"], sample)["fitted"])
        expected = curated.loc[curated["SAMPLE_ID"] == sample, "notes"].fillna("NA").tolist()
        assert split["notes"].fillna("NA").tolist() == expected


def test_report_seg_counts(staged):
    out = run(["Rscript", REPO / "resources" / "report_seg_counts.R", "all"], cwd=staged["workdir"]).output
    assert "pre-downsampled segments" in out and "post-downsampled segments" in out
    assert "NA" not in out.replace("NA12878", "")
