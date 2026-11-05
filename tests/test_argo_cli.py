"""B13: the Argo workflow's command-line contract.

Each Argo template's Rscript call is rebuilt from workflow/argo-auto-workflow.yaml and run
against the B1 (auto_default) inputs, so a renamed or removed option fails here. Outputs must
match what the Snakemake run produced from the same inputs.
"""
from pathlib import Path

import pandas as pd
import pytest

from helpers import REPO, argo_args, argo_template, read_config, read_tsv, r_eval, run, sample_paths, samtools_count

pytestmark = pytest.mark.integration

SAMPLE = "NA12878"


@pytest.fixture(scope="module")
def ctx(auto_default, tmp_path_factory):
    workdir = auto_default["workdir"]
    config = read_config(workdir)
    out = tmp_path_factory.mktemp("argo_cli")
    values = {**{k: v for k, v in config.items()}, "prplpu": "FALSE",
              "samplesheet": workdir / "config" / "sample_sheet.tsv", "sampleId": SAMPLE}
    return {"workdir": workdir, "paths": sample_paths(workdir, SAMPLE), "out": out, "values": values}


def run_template(ctx, template, **inputs):
    script, _ = argo_template(template)
    local = REPO / "workflow" / "scripts" / Path(script).name
    args = argo_args(template, {**ctx["values"], **inputs})
    return run(["Rscript", local, *args], cwd=ctx["out"])


def same_assays(a, b):
    out = r_eval(f'suppressMessages(library(Biobase)); a <- readRDS("{a}"); b <- readRDS("{b}")\n'
                 'cat(isTRUE(all.equal(assayDataElement(a, "segmented"), assayDataElement(b, "segmented"))), '
                 'isTRUE(all.equal(assayDataElement(a, "copynumber"), assayDataElement(b, "copynumber"))))')
    return out.split() == ["TRUE", "TRUE"]


def test_bam_check(ctx):
    ok = ctx["out"] / f"{SAMPLE}_bam.ok"
    run_template(ctx, "bam-check", bam=ctx["paths"]["link"], fileType="BAM", okOut=ok)
    assert ok.exists()


def test_qdnaseq_mod(ctx):
    rds = ctx["out"] / "qdnaseq.rds"
    run_template(ctx, "qdnaseq-mod", bam=ctx["paths"]["link"], meta=ctx["values"]["samplesheet"], rdsOut=rds)
    assert same_assays(rds, ctx["paths"]["pre_rds"])


def test_gridsearch(ctx):
    outs = {"tsvOut": ctx["out"] / "grid_clonality.tsv", "pdfOut": ctx["out"] / "grid.pdf",
            "plotOut": ctx["out"] / "grid.png", "filtOut": ctx["out"] / "grid_filtered.tsv",
            "fitOut": ctx["out"] / "grid_fit.tsv"}
    run_template(ctx, "gridsearch", rds=ctx["paths"]["pre_rds"], meta=ctx["values"]["samplesheet"], **outs)
    pd.testing.assert_frame_equal(read_tsv(outs["tsvOut"]), read_tsv(ctx["paths"]["clonality"]))
    pd.testing.assert_frame_equal(read_tsv(outs["fitOut"]), read_tsv(ctx["paths"]["fit"]))


def test_autofit(ctx):
    tsv = ctx["out"] / "fitted.tsv"
    run_template(ctx, "autofit", tsv=ctx["paths"]["fit"], tsvOut=tsv)
    pd.testing.assert_frame_equal(read_tsv(tsv), read_tsv(ctx["paths"]["fitted"]))


def test_downsample_bams(ctx):
    bam = ctx["out"] / "downsampled.bam"
    run_template(ctx, "downsample-bams", bam=ctx["paths"]["link"], meta=ctx["paths"]["fitted"],
                 rds=ctx["paths"]["pre_rds"], bamOut=bam)
    assert samtools_count(bam) == samtools_count(ctx["paths"]["ds_bam"])


def test_rel_to_abs(ctx):
    outs = {"tsvOut": ctx["out"] / "abs.tsv", "rdsOut": ctx["out"] / "abs.rds",
            "tabOut": ctx["out"] / "seg.tsv", "plotOut": ctx["out"] / "abs.png"}
    run_template(ctx, "qdnaseq-rel-to-abs", rds=ctx["paths"]["post_rds"], meta=ctx["paths"]["fitted"], **outs)
    pd.testing.assert_frame_equal(read_tsv(outs["tsvOut"]), read_tsv(ctx["paths"]["abs_tsv"]))
    pd.testing.assert_frame_equal(read_tsv(outs["tabOut"]), read_tsv(ctx["paths"]["seg"]))


def test_combine_results_takes_comma_separated_lists(ctx, auto_default):
    samples = [s["SAMPLE_ID"] for s in auto_default["samples"]]
    paths = [sample_paths(ctx["workdir"], s) for s in samples]
    out = {k: ctx["out"] / f"combined.{k}" for k in ("tsv", "rds", "tab")}
    run(["Rscript", REPO / "workflow" / "scripts" / "combineResults.R",
         "--tsv", ",".join(str(p["abs_tsv"]) for p in paths),
         "--rds", ",".join(str(p["abs_rds"]) for p in paths),
         "--tab", ",".join(str(p["seg"]) for p in paths),
         "--tsvOut", out["tsv"], "--rdsOut", out["rds"], "--tabOut", out["tab"]], cwd=ctx["out"])
    expected = pd.concat([read_tsv(p["abs_tsv"]) for p in paths], ignore_index=True)
    pd.testing.assert_frame_equal(read_tsv(out["tsv"]), expected)
    assert sorted(read_tsv(out["tab"])["sample"].unique()) == sorted(samples)
    combined = r_eval(f'cat(Biobase::sampleNames(readRDS("{out["rds"]}")), sep = ",")').strip()
    assert sorted(combined.split(",")) == sorted(samples)
