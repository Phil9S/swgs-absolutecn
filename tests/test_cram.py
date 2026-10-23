"""B10: CRAM input - conversion to BAM with the declared index, and file type mismatch.

Uses a chr21-only CRAM of the Zenodo BAM and the GRCh37 chr21 reference. A full CRAM run
needs the complete hs37d5 reference and only runs when SWGS_TEST_REFERENCE points at it.
"""
import os
import shutil

import pytest

from helpers import BASE, make_workdir, relpath, run, run_snakemake, sample_paths, samtools_count

pytestmark = [pytest.mark.cram, pytest.mark.integration]


def test_cram_converted_to_indexed_bam(tmp_path, chr21):
    make_workdir(tmp_path, [{"SAMPLE_ID": "C21", "file": chr21["cram"]}],
                 **{**BASE, "filetype": "CRAM", "reference": chr21["fasta"]})
    paths = sample_paths(tmp_path, "C21")
    run_snakemake(tmp_path, "stage_1", relpath(tmp_path, paths["ok"]), relpath(tmp_path, paths["link"]), "--notemp")
    cram_reads = int(run(["samtools", "view", "-c", "-T", chr21["fasta"], chr21["cram"]]).stdout)
    assert samtools_count(paths["link"]) == cram_reads > 0
    assert os.path.exists(f"{paths['link']}.bai")
    assert not os.path.islink(paths["link"])


def test_cram_given_as_bam_is_rejected(tmp_path, chr21):
    (tmp_path / "inputs").mkdir()
    cram = tmp_path / "inputs" / "MIS.cram"
    shutil.copy(chr21["cram"], cram)
    # placeholder index so the BAM index check passes and bam_check sees the CRAM
    (tmp_path / "inputs" / "MIS.cram.bai").touch()
    make_workdir(tmp_path, [{"SAMPLE_ID": "MIS", "file": str(cram)}], **BASE)
    paths = sample_paths(tmp_path, "MIS")
    run_snakemake(tmp_path, "stage_1", relpath(tmp_path, paths["ok"]), ok=False)
    invalid = paths["ok"].with_name("MIS_bam.invalid")
    assert "Filetype specified as BAM, detected CRAM" in invalid.read_text()


@pytest.mark.skipif(not os.environ.get("SWGS_TEST_REFERENCE"), reason="set SWGS_TEST_REFERENCE to an hs37d5 FASTA")
def test_full_cram_stage1(tmp_path, source_bam):
    reference = os.environ["SWGS_TEST_REFERENCE"]
    cram = tmp_path / "NA12878.cram"
    run(["samtools", "view", "-C", "-T", reference, "-o", cram, source_bam])
    run(["samtools", "index", cram])
    make_workdir(tmp_path, [{"SAMPLE_ID": "NA12878", "file": str(cram)}],
                 **{**BASE, "filetype": "CRAM", "reference": reference})
    run_snakemake(tmp_path, "stage_1", cores=8)
