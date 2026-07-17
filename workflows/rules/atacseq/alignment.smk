import hashlib
import os


# --- CONFIGURATION ---
REF_GENOME = config["ref_genome"]
BOWTIE2_INDEX_DIR = config.get("bowtie2_index_dir")
CONFIGURED_BOWTIE2_INDEX_PREFIX = config.get("bowtie2_index_prefix")
REFERENCE_DIR = config["reference_dir"]
TRIMMED_DIR = config["trimmed_dir"]
ALIGNMENT_DIR = config["alignment_dir"]
BOWTIE2_CONFIG = config.get("bowtie2", {})
BOWTIE2_EXTRA_ARGS = BOWTIE2_CONFIG.get("extra_args", "")
ACTIVE_ATAC_PIPELINE = config.get("pipeline") == "atacseq"

STANDARD_INDEX_SUFFIXES = [
    ".1.bt2",
    ".2.bt2",
    ".3.bt2",
    ".4.bt2",
    ".rev.1.bt2",
    ".rev.2.bt2",
]
LARGE_INDEX_SUFFIXES = [
    ".1.bt2l",
    ".2.bt2l",
    ".3.bt2l",
    ".4.bt2l",
    ".rev.1.bt2l",
    ".rev.2.bt2l",
]


def _bowtie2_index_files(prefix, suffixes):
    return [f"{prefix}{suffix}" for suffix in suffixes]


PREBUILT_BOWTIE2_INDEX = CONFIGURED_BOWTIE2_INDEX_PREFIX is not None

if PREBUILT_BOWTIE2_INDEX:
    BOWTIE2_INDEX_PREFIX = CONFIGURED_BOWTIE2_INDEX_PREFIX
    known_suffixes = STANDARD_INDEX_SUFFIXES + LARGE_INDEX_SUFFIXES
    if ACTIVE_ATAC_PIPELINE and any(
        BOWTIE2_INDEX_PREFIX.endswith(suffix) for suffix in known_suffixes
    ):
        raise ValueError(
            "ATAC-seq bowtie2_index_prefix must not include a .bt2 or .bt2l file suffix: "
            f"{BOWTIE2_INDEX_PREFIX}"
        )

    standard_files = _bowtie2_index_files(
        BOWTIE2_INDEX_PREFIX, STANDARD_INDEX_SUFFIXES
    )
    large_files = _bowtie2_index_files(BOWTIE2_INDEX_PREFIX, LARGE_INDEX_SUFFIXES)
    standard_present = any(os.path.exists(path) for path in standard_files)
    large_present = any(os.path.exists(path) for path in large_files)
    standard_complete = all(os.path.isfile(path) for path in standard_files)
    large_complete = all(os.path.isfile(path) for path in large_files)

    if ACTIVE_ATAC_PIPELINE and standard_present and large_present:
        raise ValueError(
            "ATAC-seq prebuilt Bowtie2 index is mixed or ambiguous because .bt2 and "
            f".bt2l files share the prefix: {BOWTIE2_INDEX_PREFIX}"
        )
    if standard_complete:
        BOWTIE2_INDEX_FORMAT = "standard"
        BOWTIE2_INDEX_FILES = standard_files
    elif large_complete:
        BOWTIE2_INDEX_FORMAT = "large"
        BOWTIE2_INDEX_FILES = large_files
    elif ACTIVE_ATAC_PIPELINE:
        missing_standard = [path for path in standard_files if not os.path.isfile(path)]
        missing_large = [path for path in large_files if not os.path.isfile(path)]
        raise ValueError(
            "ATAC-seq bowtie2_index_prefix does not identify a complete Bowtie2 index. "
            "Provide all six .bt2 files or all six .bt2l files. "
            f"Missing standard files: {', '.join(missing_standard)}. "
            f"Missing large files: {', '.join(missing_large)}."
        )
    else:
        # The ATAC module is parsed even when another pipeline is selected.
        BOWTIE2_INDEX_FORMAT = "standard"
        BOWTIE2_INDEX_FILES = standard_files
else:
    large_index = BOWTIE2_CONFIG.get("large_index", False)
    if not isinstance(large_index, bool):
        if ACTIVE_ATAC_PIPELINE:
            raise ValueError("ATAC-seq bowtie2.large_index must be true or false.")
        large_index = False

    ref_basename = os.path.basename(REF_GENOME)
    for suffix in [".fa.gz", ".fasta.gz", ".fna.gz", ".fa", ".fasta", ".fna"]:
        if ref_basename.endswith(suffix):
            ref_basename = ref_basename[: -len(suffix)]
            break
    if not ref_basename:
        ref_basename = "genome"

    BOWTIE2_INDEX_PREFIX = os.path.join(BOWTIE2_INDEX_DIR, ref_basename)
    if large_index:
        BOWTIE2_INDEX_FORMAT = "large"
        BOWTIE2_INDEX_FILES = _bowtie2_index_files(
            BOWTIE2_INDEX_PREFIX, LARGE_INDEX_SUFFIXES
        )
    else:
        BOWTIE2_INDEX_FORMAT = "standard"
        BOWTIE2_INDEX_FILES = _bowtie2_index_files(
            BOWTIE2_INDEX_PREFIX, STANDARD_INDEX_SUFFIXES
        )

index_identity = hashlib.sha256(
    f"{os.path.abspath(BOWTIE2_INDEX_PREFIX)}|{BOWTIE2_INDEX_FORMAT}".encode()
).hexdigest()[:12]
BOWTIE2_VALIDATION_MARKER = os.path.join(
    REFERENCE_DIR, f"bowtie2_index_{index_identity}.validated.OK"
)


if not PREBUILT_BOWTIE2_INDEX:
    rule bowtie2_build:
        """
        Builds and tracks all six Bowtie2 index components.
        """
        input:
            ref = REF_GENOME
        output:
            index = BOWTIE2_INDEX_FILES
        params:
            prefix = BOWTIE2_INDEX_PREFIX,
            large_arg = "--large-index" if BOWTIE2_INDEX_FORMAT == "large" else ""
        threads: 1
        log:
            os.path.join("logs", config["pipeline"], "bowtie2_build", "bowtie2_build.log")
        shell:
            """
            set -euo pipefail
            pixi run bowtie2-build {params.large_arg} "{input.ref}" "{params.prefix}" \
                > {log}.out 2> {log}.err
            """


rule bowtie2_validate:
    """
    Validates the complete generated or prebuilt Bowtie2 index once.
    """
    input:
        index = BOWTIE2_INDEX_FILES
    output:
        marker = BOWTIE2_VALIDATION_MARKER
    params:
        prefix = BOWTIE2_INDEX_PREFIX
    threads: 1
    log:
        os.path.join("logs", config["pipeline"], "bowtie2_validate", "bowtie2_validate.log")
    shell:
        """
        set -euo pipefail
        pixi run bowtie2-inspect -n "{params.prefix}" > {log}.out 2> {log}.err
        touch "{output.marker}"
        """


rule bowtie2_align_pe:
    """
    Aligns trimmed paired-end reads with Bowtie2 and sorts the resulting BAM.
    """
    input:
        r1 = os.path.join(TRIMMED_DIR, "{sample}_R1.trimmed.fq.gz"),
        r2 = os.path.join(TRIMMED_DIR, "{sample}_R2.trimmed.fq.gz"),
        index = BOWTIE2_INDEX_FILES,
        index_validated = BOWTIE2_VALIDATION_MARKER
    output:
        bam = os.path.join(ALIGNMENT_DIR, "{sample}_pe.sorted.bam")
    params:
        extra = BOWTIE2_EXTRA_ARGS,
        index_prefix = BOWTIE2_INDEX_PREFIX
    threads: 4
    log:
        os.path.join("logs", config["pipeline"], "bowtie2_align", "{sample}_pe.log")
    shell:
        """
        set -euo pipefail
        (pixi run bowtie2 -p {threads} {params.extra} -x "{params.index_prefix}" \
            -1 "{input.r1}" -2 "{input.r2}" | \
        pixi run samtools view -bS - | \
        pixi run samtools sort -@ {threads} -o "{output.bam}" -) \
            > {log}.out 2> {log}.err
        """
