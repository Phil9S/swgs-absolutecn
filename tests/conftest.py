"""Session fixtures: test data (download, md5 check, cache), derived BAMs and shared workflow runs.

Test data is cached in $SWGS_TEST_DATA (default ~/.cache/swgs-absolutecn-tests) so it is
downloaded once. Derived BAMs are built from the Zenodo NA12878 BAM with samtools and cached
next to it, keyed on the source md5.
"""
import gzip
import hashlib
import json
import os
import shutil
import urllib.request
from pathlib import Path

import pytest

from helpers import BASE, DELETION, make_workdir, run, run_snakemake

DATA = json.loads((Path(__file__).parent / "test_data.json").read_text())
CACHE = Path(os.environ.get("SWGS_TEST_DATA", Path.home() / ".cache" / "swgs-absolutecn-tests"))
ZENODO = DATA["zenodo"]
SOURCE_MD5 = ZENODO["files"][ZENODO["bam"]]["md5"]

# 2 Mb half-depth blocks every 20 Mb on the autosomes (MOSAIC_1 - many segments for smoothing)
MOSAIC_START, MOSAIC_STEP, MOSAIC_BLOCK = 10_000_000, 20_000_000, 2_000_000
AUTOSOMES = [str(i) for i in range(1, 23)]


def md5sum(path):
    digest = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch(url, dest, md5=None):
    """Download url to dest unless a copy with the right md5 is already cached."""
    dest = Path(dest)
    if dest.exists() and (not md5 or md5sum(dest) == md5):
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    try:
        with urllib.request.urlopen(url, timeout=120) as response, open(part, "wb") as fh:
            shutil.copyfileobj(response, fh, 1 << 20)
    except OSError as err:
        part.unlink(missing_ok=True)
        pytest.skip(f"test data unavailable ({url}): {err}")
    if md5 and md5sum(part) != md5:
        part.unlink()
        raise AssertionError(f"md5 mismatch for {url} - expected {md5}")
    if dest.is_symlink() or dest.exists():
        dest.unlink()
    part.replace(dest)
    return dest


def contig_lengths(bam):
    """{contig: length} from samtools idxstats (excludes the '*' unmapped line)."""
    lengths = {}
    for line in run(["samtools", "idxstats", bam]).stdout.splitlines():
        name, length = line.split("\t")[:2]
        if name != "*":
            lengths[name] = int(length)
    return lengths


def write_bed(path, intervals):
    Path(path).write_text("".join(f"{c}\t{s}\t{e}\n" for c, s, e in intervals if e > s))
    return path


def derive(path, make, index=True):
    """Build a derived BAM once (via a temporary name) and index it."""
    path = Path(path)
    if path.exists() and (not index or Path(f"{path}.bai").exists()):
        return path
    tmp = path.with_name(path.stem + ".tmp" + path.suffix)
    make(tmp)
    if index:
        run(["samtools", "index", tmp])
        Path(f"{tmp}.bai").replace(f"{path}.bai")
    tmp.replace(path)
    return path


@pytest.fixture(scope="session")
def source_bam():
    """The Zenodo NA12878 0.1X BAM and its index, verified by md5."""
    folder = CACHE / f"zenodo-{ZENODO['record']}"
    for name, meta in ZENODO["files"].items():
        fetch(ZENODO["url"].format(name=name), folder / name, meta["md5"])
    return folder / ZENODO["bam"]


@pytest.fixture(scope="session")
def derived_dir(source_bam):
    folder = CACHE / f"derived-{SOURCE_MD5[:8]}"
    folder.mkdir(parents=True, exist_ok=True)
    return folder


@pytest.fixture(scope="session")
def bams(source_bam, derived_dir):
    """Test BAMs, all indexed next to the real file."""
    src = str(source_bam)
    lengths = contig_lengths(src)

    def subsample(fraction, seed):
        return lambda out: run(["samtools", "view", "-b", "--subsample", fraction,
                                "--subsample-seed", seed, "-o", out, src])

    def deletion(out):
        chrom, start, end = DELETION
        keep = [(c, 0, n) for c, n in lengths.items() if c != chrom] + [(chrom, 0, start), (chrom, end, lengths[chrom])]
        bed = write_bed(derived_dir / "deletion_keep.bed", keep)
        run(["samtools", "view", "-b", "-L", bed, "-o", out, src])

    def mosaic(out):
        half, full = [], []
        for c, n in lengths.items():
            pos = 0
            blocks = range(MOSAIC_START, n - MOSAIC_BLOCK, MOSAIC_STEP) if c in AUTOSOMES else []
            for start in blocks:
                full.append((c, pos, start))
                half.append((c, start, start + MOSAIC_BLOCK))
                pos = start + MOSAIC_BLOCK
            full.append((c, pos, n))
        full_bed = write_bed(derived_dir / "mosaic_full.bed", full)
        half_bed = write_bed(derived_dir / "mosaic_half.bed", half)
        a, b = Path(f"{out}.a.bam"), Path(f"{out}.b.bam")
        run(["samtools", "view", "-b", "-L", full_bed, "-o", a, src])
        run(["samtools", "view", "-b", "-L", half_bed, "--subsample", "0.5", "--subsample-seed", "3", "-o", b, src])
        run(["samtools", "merge", "-f", "-o", out, a, b])
        a.unlink()
        b.unlink()

    corrupt = derived_dir / "CORRUPT_1.bam"
    if not corrupt.exists():
        with open(src, "rb") as fh:
            corrupt.write_bytes(fh.read(1 << 20))
        # placeholder index so the corrupt file reaches bam_check
        Path(f"{corrupt}.bai").write_bytes(b"")

    return {
        "NA12878": src,
        "007": str(derive(derived_dir / "SUB50.bam", subsample("0.5", "7"))),
        "LOW_1": str(derive(derived_dir / "LOW5.bam", subsample("0.05", "11"))),
        "DEL_1": str(derive(derived_dir / "DEL_1.bam", deletion)),
        "MOSAIC_1": str(derive(derived_dir / "MOSAIC_1.bam", mosaic)),
        "CORRUPT_1": str(corrupt),
    }


@pytest.fixture(scope="session")
def chr21(source_bam, derived_dir):
    """GRCh37 chr21 reference (uncompressed, faidx) and a chr21-only CRAM of the source BAM."""
    ref = DATA["chr21_reference"]
    gz = fetch(ref["url"], CACHE / "reference" / ref["name"], ref["md5"] or None)
    fasta = gz.with_suffix("")
    if not Path(f"{fasta}.fai").exists():
        with gzip.open(gz, "rb") as src, open(fasta, "wb") as dst:
            shutil.copyfileobj(src, dst, 1 << 20)
        run(["samtools", "faidx", fasta])
    cram = derived_dir / "chr21.cram"
    if not Path(f"{cram}.crai").exists():
        # keep only the chr21 @SQ line so the chr21-only reference covers the header
        header = [line for line in run(["samtools", "view", "-H", source_bam]).stdout.splitlines()
                  if not line.startswith("@SQ") or "\tSN:21\t" in line + "\t"]
        reads = run(["samtools", "view", source_bam, "21"]).stdout
        sam = derived_dir / "chr21.sam"
        sam.write_text("\n".join(header) + "\n" + reads)
        run(["samtools", "view", "-C", "-T", fasta, "-o", cram, sam])
        run(["samtools", "index", cram])
        sam.unlink()
    return {"fasta": str(fasta), "cram": str(cram)}


# ---- shared workflow runs -------------------------------------------------

@pytest.fixture(scope="session")
def auto_default(tmp_path_factory, bams):
    """B1: auto workflow with default settings at 100 kb on three samples.

    007 uses a relative sample-sheet path that is a symlink to the cached BAM,
    DEL_1 carries the 5 Mb homozygous deletion.
    """
    workdir = tmp_path_factory.mktemp("auto_default")
    (workdir / "inputs").mkdir()
    (workdir / "inputs" / "007.bam").symlink_to(bams["007"])
    samples = [
        {"SAMPLE_ID": "NA12878", "PATIENT_ID": "P1", "file": bams["NA12878"]},
        {"SAMPLE_ID": "007", "PATIENT_ID": "P1", "TP53freq": 0.95, "file": "inputs/007.bam"},
        {"SAMPLE_ID": "DEL_1", "PATIENT_ID": "P2", "file": bams["DEL_1"]},
    ]
    make_workdir(workdir, samples, **BASE)
    res = run_snakemake(workdir, "auto", cores=8)
    return {"workdir": workdir, "samples": samples, "result": res, "bams": bams}
