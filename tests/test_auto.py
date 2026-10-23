"""B1-B6: auto workflow runs on the Zenodo NA12878 BAM and BAMs derived from it.

NA12878 is a normal diploid sample, so the expected fit is ploidy 2 / purity 1 with
one segment per autosome. DEL_1 adds a 5 Mb homozygous deletion and MOSAIC_1 adds
2 Mb half-depth blocks (many segments, to exercise smoothing).
"""
import os
import shutil

import numpy as np
import pytest

from helpers import (BASE, DELETION, make_workdir, read_config, read_tsv, r_eval, relpath, run,
                     run_snakemake, sample_paths, samtools_count)

pytestmark = pytest.mark.integration

AUTOFIT_NOTES = "autofit|fitmethod=errorOnly|flagThreshold=0.84|ErrorMetric=clonality"


def selected(fitted):
    rows = fitted[fitted["use"] == True]  # noqa: E712 - pandas boolean column
    assert len(rows) == 1, f"expected one selected fit, found {len(rows)}"
    return rows.iloc[0]


def weighted_median_segval(seg):
    order = np.argsort(seg["segVal"].to_numpy())
    widths = (seg["end"] - seg["start"] + 1).to_numpy()[order]
    cum = np.cumsum(widths) / widths.sum()
    return seg["segVal"].to_numpy()[order][np.searchsorted(cum, 0.5)]


def r_value(code):
    return r_eval(f"suppressMessages(library(Biobase))\n{code}\n").strip()


# ---- B1: default auto run ---------------------------------------------------

@pytest.fixture(scope="module")
def b1(auto_default):
    return auto_default


def test_all_outputs_exist(b1):
    for sample in b1["samples"]:
        for name, path in sample_paths(b1["workdir"], sample["SAMPLE_ID"]).items():
            assert path.exists(), f"{sample['SAMPLE_ID']}: missing {name} ({path})"


def test_one_fit_selected_with_autofit_notes(b1):
    for sample in b1["samples"]:
        fitted = read_tsv(sample_paths(b1["workdir"], sample["SAMPLE_ID"])["fitted"])
        selected(fitted)
        assert fitted["notes"].str.startswith(AUTOFIT_NOTES).all()


def test_normal_genome_is_fitted_diploid(b1):
    for sample in b1["samples"]:
        paths = sample_paths(b1["workdir"], sample["SAMPLE_ID"])
        fit = selected(read_tsv(paths["fitted"]))
        assert 1.8 <= fit["ploidy"] <= 2.2 and fit["purity"] >= 0.9, (sample["SAMPLE_ID"], fit["ploidy"], fit["purity"])
        assert 1.8 <= weighted_median_segval(read_tsv(paths["seg"])) <= 2.2


def test_homozygous_deletion_detected(b1):
    paths = sample_paths(b1["workdir"], "DEL_1")
    chrom, start, end = DELETION
    seg = read_tsv(paths["seg"])
    inside = seg[(seg["chromosome"] == chrom) & (seg["start"] < end) & (seg["end"] > start)]
    assert (inside["segVal"] < 0.5).any(), inside
    loss = selected(read_tsv(paths["fitted"]))["homozygousLoss"]
    assert 4.5e6 <= loss <= 5.5e6, loss


def test_homozygous_loss_counted_in_bins(b1):
    bin_size = read_config(b1["workdir"])["bins"] * 1000
    for sample in b1["samples"]:
        grid = read_tsv(sample_paths(b1["workdir"], sample["SAMPLE_ID"])["clonality"])
        assert (grid["homozygousLoss"] % bin_size == 0).all()


def test_ids_and_tp53_carried_through(b1):
    paths = sample_paths(b1["workdir"], "007")
    assert (read_tsv(paths["fitted"])["SAMPLE_ID"] == "007").all()
    abs_fits = read_tsv(paths["abs_tsv"])
    assert (abs_fits["SAMPLE_ID"] == "007").all()
    assert abs_fits["TP53freq"].iloc[0] == pytest.approx(0.95)
    assert (read_tsv(paths["seg"])["sample"] == "007").all()


def test_segment_counts_are_consistent(b1):
    for sample in b1["samples"]:
        paths = sample_paths(b1["workdir"], sample["SAMPLE_ID"])
        pre_count = int(r_value(f'cat(rswgsabsolutecn::getSegCount(readRDS("{paths["pre_rds"]}")))'))
        assert (read_tsv(paths["fit"])["segments"] == pre_count).all()
        abs_fits = read_tsv(paths["abs_tsv"])
        assert abs_fits["segments.pre"].iloc[0] == pre_count
        assert abs_fits["segments.post"].iloc[0] == len(read_tsv(paths["seg"]))


def test_downsampling_follows_the_selected_fit(b1):
    for sample in b1["samples"]:
        paths = sample_paths(b1["workdir"], sample["SAMPLE_ID"])
        depth = selected(read_tsv(paths["fitted"]))["downsample_depth"]
        total = float(r_value(f'cat(pData(readRDS("{paths["pre_rds"]}"))$total.reads)'))
        ratio = min(round(depth / total, 4), 1)
        if ratio <= 0.96:
            assert not os.path.islink(paths["ds_bam"])
            kept = samtools_count(paths["ds_bam"]) / samtools_count(os.path.realpath(paths["link"]))
            assert kept == pytest.approx(ratio, abs=0.01), (sample["SAMPLE_ID"], kept, ratio)
        else:
            assert os.path.islink(paths["ds_bam"])
        assert os.path.exists(f"{paths['ds_bam']}.bai")


def test_bam_links_point_to_real_files(b1):
    for sample in b1["samples"]:
        source = sample["file"] if os.path.isabs(sample["file"]) else os.path.join(b1["workdir"], sample["file"])
        link = sample_paths(b1["workdir"], sample["SAMPLE_ID"])["link"]
        assert os.readlink(link) == os.path.realpath(source)


def test_rule_logs(b1):
    logs = b1["workdir"] / "logs"
    for log in logs.glob("*.log"):
        text = log.read_text()
        assert "Execution halted" not in text and "Error in" not in text, log.name
    for rule in ("relRDS", "gridsearch_fitting", "downsample", "ds_relRDS", "rel_to_abs"):
        for sample in b1["samples"]:
            log = logs / f"{rule}_test-run_1_100kb_{sample['SAMPLE_ID']}.log"
            assert log.exists() and log.stat().st_size > 0, log.name


def test_grid_spans_config_limits(b1):
    config = read_config(b1["workdir"])
    grid = read_tsv(sample_paths(b1["workdir"], "NA12878")["clonality"])
    assert grid["ploidy"].min() == pytest.approx(config["ploidy_min"])
    assert grid["ploidy"].max() == pytest.approx(config["ploidy_max"])
    assert grid["purity"].min() == pytest.approx(config["purity_min"])
    assert grid["purity"].max() == pytest.approx(config["purity_max"])


# ---- B1a: autofit re-runs with other error metrics --------------------------

@pytest.fixture(scope="module")
def b1_copy(b1, tmp_path_factory):
    copy = tmp_path_factory.mktemp("auto_rerun") / "work"
    shutil.copytree(b1["workdir"], copy, symlinks=True)
    return copy


@pytest.mark.parametrize("metric", ["segvariance", "rmse"])
def test_error_metric_selects_minimum(b1_copy, metric):
    fitted_path = sample_paths(b1_copy, "NA12878")["fitted"]
    run_snakemake(b1_copy, "auto", relpath(b1_copy, fitted_path), "--config", f"errorMetric={metric}",
                  "--forcerun", "autofit")
    fitted = read_tsv(fitted_path)
    assert selected(fitted)[metric] == pytest.approx(fitted[metric].min())
    assert fitted["notes"].str.endswith(f"ErrorMetric={metric}").all()


# ---- B2: random forest autofit at the model's 30 kb calibration -------------

@pytest.fixture(scope="module")
def b2(tmp_path_factory, bams):
    workdir = tmp_path_factory.mktemp("randforest")
    make_workdir(workdir, [{"SAMPLE_ID": "NA12878", "file": bams["NA12878"]}],
                 **{**BASE, "bins": 30, "fitMethod": "randforest"})
    run_snakemake(workdir, "auto", cores=8)
    return workdir


def flags(fitted):
    return fitted["flag"].fillna("").astype(str)


def test_randforest_selects_one_diploid_fit(b2):
    fitted = read_tsv(sample_paths(b2, "NA12878")["fitted"])
    fit = selected(fitted)
    assert (fit["ploidy"], fit["purity"]) == (2.0, 1.0)
    predicted = fitted["pred_class"].astype(str).str.upper() == "TRUE"
    if predicted.any():
        assert fit["clonality"] == pytest.approx(fitted.loc[predicted, "clonality"].min())
    else:
        assert "NOFIT" in str(fit["flag"]) and str(fit["notes"]).startswith("NO_FIT_FOUND|")


def test_randforest_probabilities_and_flags(b2):
    fitted = read_tsv(sample_paths(b2, "NA12878")["fitted"])
    scored = fitted.dropna(subset=["pred_TRUE"])
    assert len(scored) > 0
    assert (scored["pred_TRUE"] + scored["pred_FALSE"]).to_numpy() == pytest.approx(1, abs=0.002)
    low = scored["triageValue"] < 0.84
    assert (flags(scored).str.contains("LOWPROBFLAG") == low).all()
    log = (b2 / "logs" / "autofit_test-run_1_30kb_NA12878.log").read_text()
    assert "not available" not in log and "namespace" not in log


def test_flag_threshold_zero_removes_low_probability_flags(b2):
    fitted_path = sample_paths(b2, "NA12878")["fitted"]
    run_snakemake(b2, "auto", relpath(b2, fitted_path), "--config", "flagThreshold=0", "--forcerun", "autofit")
    fitted = read_tsv(fitted_path)
    assert not flags(fitted).str.contains("LOWPROBFLAG").any()
    assert fitted["notes"].str.contains("flagThreshold=0").all()


# ---- B3: gridsearch filters ---------------------------------------------------

FILTERS = {"filter_underpowered": "TRUE", "filter_homozygous": "TRUE", "homozygous_prop": 6_000_000,
           "homozygous_threshold": 0.6, "ploidy_min": 1.8, "ploidy_max": 3.5, "purity_min": 0.4, "purity_max": 1.0}


@pytest.fixture(scope="module")
def b3(tmp_path_factory, bams):
    workdir = tmp_path_factory.mktemp("filters")
    samples = [{"SAMPLE_ID": s, "file": bams[s]} for s in ("NA12878", "DEL_1")]
    make_workdir(workdir, samples, **{**BASE, **FILTERS})
    run_snakemake(workdir, "auto", cores=8)
    return workdir


def test_filters_keep_only_powered_fits_within_homozygous_limit(b3):
    for sample in ("NA12878", "DEL_1"):
        fit = read_tsv(sample_paths(b3, sample)["fit"])
        assert len(fit) > 0
        assert (fit["powered"] == 1).all()
        assert (fit["homozygousLoss"] <= FILTERS["homozygous_prop"]).all()
    assert (read_tsv(sample_paths(b3, "DEL_1")["fit"])["homozygousLoss"] > 0).all()


def test_grid_uses_configured_ranges(b3):
    grid = read_tsv(sample_paths(b3, "NA12878")["clonality"])
    assert (grid["ploidy"].min(), grid["ploidy"].max()) == pytest.approx((1.8, 3.5))
    assert (grid["purity"].min(), grid["purity"].max()) == pytest.approx((0.4, 1.0))


def test_af_cutoff_filters_on_tp53(tmp_path, bams):
    make_workdir(tmp_path, [{"SAMPLE_ID": "007", "TP53freq": 0.5, "file": bams["007"]}], **BASE)
    filtered_path = sample_paths(tmp_path, "007")["filtered"]
    target = relpath(tmp_path, filtered_path)
    run_snakemake(tmp_path, "auto", target, cores=8)
    diff = (read_tsv(filtered_path)["expected_TP53_AF"] - 0.5).abs()
    assert len(diff) > 0 and (diff <= 0.15 + 1e-9).all()
    run_snakemake(tmp_path, "auto", target, "--config", "af_cutoff=1", "--forcerun", "gridsearch")
    diff = (read_tsv(filtered_path)["expected_TP53_AF"] - 0.5).abs()
    assert (diff > 0.15).any()


def test_homozygous_filter_failure(tmp_path, bams):
    make_workdir(tmp_path, [{"SAMPLE_ID": "DEL_1", "file": bams["DEL_1"]}],
                 **{**BASE, "homozygous_prop": 1_000_000})
    target = relpath(tmp_path, sample_paths(tmp_path, "DEL_1")["fit"])
    res = run_snakemake(tmp_path, "auto", target, ok=False)
    assert "No fits meet filtering criteria - filter_homozygous" in res.output


def test_underpowered_filter_failure(tmp_path, bams):
    make_workdir(tmp_path, [{"SAMPLE_ID": "LOW_1", "file": bams["LOW_1"]}],
                 **{**BASE, "filter_underpowered": "TRUE"})
    target = relpath(tmp_path, sample_paths(tmp_path, "LOW_1")["fit"])
    res = run_snakemake(tmp_path, "auto", target, ok=False)
    assert "No fits meet filtering criteria - filter_underpowered" in res.output


# ---- B4: smoothing, seeds and pairedEnd (relRDS only) -----------------------

@pytest.fixture(scope="module")
def b4(tmp_path_factory, bams):
    """relRDS runs of MOSAIC_1 under different settings; returns per-variant stats."""
    workdir = tmp_path_factory.mktemp("smoothing")
    make_workdir(workdir, [{"SAMPLE_ID": "MOSAIC_1", "file": bams["MOSAIC_1"]}], **BASE)
    paths = sample_paths(workdir, "MOSAIC_1")
    target = relpath(workdir, paths["pre_rds"])
    log = workdir / "logs" / "relRDS_test-run_1_100kb_MOSAIC_1.log"
    results = {}

    def variant(name, *config):
        args = ["--config", *config] if config else []
        run_snakemake(workdir, "auto", target, *args, "--forcerun", "relRDS")
        stats = r_value(f'x <- readRDS("{paths["pre_rds"]}"); '
                        f'cat(rswgsabsolutecn::getSegCount(x), pData(x)$expected.variance)').split()
        keep = workdir / f"{name}.rds"
        shutil.copy(paths["pre_rds"], keep)
        results[name] = {"segments": int(stats[0]), "variance": float(stats[1]), "rds": keep,
                         "segmentations": log.read_text().count("Performing segmentation")}

    variant("base")
    variant("repeat")
    variant("autosmooth", "autoSmooth='TRUE'", "smoothThreshold=44")
    variant("single_end", "pairedEnd='FALSE'")
    variant("no_seed", "use_seed='FALSE'")
    sheet = workdir / "config" / "sample_sheet.tsv"
    sheet.write_text(sheet.read_text().replace("\tFALSE\t", "\tTRUE\t", 1))
    variant("smooth_sample")
    return results


def test_mosaic_has_many_segments(b4):
    assert b4["base"]["segments"] > 44


def test_seeded_segmentation_is_reproducible(b4):
    out = r_value(f'a <- readRDS("{b4["base"]["rds"]}"); b <- readRDS("{b4["repeat"]["rds"]}"); '
                  'cat(identical(assayDataElement(a, "segmented"), assayDataElement(b, "segmented")))')
    assert out == "TRUE"


def test_autosmooth_resegments_above_threshold(b4):
    assert b4["autosmooth"]["segmentations"] > 1
    assert b4["autosmooth"]["segments"] <= b4["base"]["segments"]


def test_sample_sheet_smooth_flag(b4):
    if b4["base"]["segments"] > 300:
        assert b4["smooth_sample"]["segmentations"] > 1
    assert b4["smooth_sample"]["segments"] <= b4["base"]["segments"]


def test_paired_end_changes_expected_variance(b4):
    assert b4["single_end"]["variance"] != pytest.approx(b4["base"]["variance"])


def test_unseeded_run_completes(b4):
    assert b4["no_seed"]["segments"] > 0


# ---- B6: corrupt BAM ------------------------------------------------------

def test_corrupt_bam_fails_bam_check(tmp_path, bams):
    make_workdir(tmp_path, [{"SAMPLE_ID": "CORRUPT_1", "file": bams["CORRUPT_1"]}], **BASE)
    paths = sample_paths(tmp_path, "CORRUPT_1")
    res = run_snakemake(tmp_path, "auto", relpath(tmp_path, paths["ok"]), ok=False)
    assert "BAM check failed" in res.output
    assert not paths["ok"].exists()
    assert paths["ok"].with_name("CORRUPT_1_bam.invalid").exists()
