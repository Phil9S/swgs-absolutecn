"""Checks for the shared test helpers."""
import pytest

from helpers import argo_template, make_workdir, read_config, read_tsv

pytestmark = pytest.mark.dryrun


def test_make_workdir_keeps_leading_zero_ids_and_overrides(tmp_path):
    make_workdir(tmp_path, [{"SAMPLE_ID": "007", "PATIENT_ID": "001", "file": "x.bam"}],
                 bins=100, filter_underpowered="TRUE")
    sheet = read_tsv(tmp_path / "config" / "sample_sheet.tsv")
    assert sheet.loc[0, "SAMPLE_ID"] == "007"
    assert sheet.loc[0, "PATIENT_ID"] == "001"
    config = read_config(tmp_path)
    assert config["bins"] == 100
    assert config["filter_underpowered"] == "TRUE"
    assert config["samplesheet"] == "config/sample_sheet.tsv"


def test_argo_template_parsing():
    script, pairs = argo_template("qdnaseq-rel-to-abs")
    assert script == "/app/workflow/scripts/qdnaseqRelToAbs.R"
    assert [flag for flag, _ in pairs] == ["rds", "meta", "genome", "sampleName",
                                          "tsvOut", "rdsOut", "tabOut", "plotOut"]
