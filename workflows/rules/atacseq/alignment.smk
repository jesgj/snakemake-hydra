import hashlib
import os
import shlex


# --- CONFIGURATION ---
REF_GENOME = config["ref_genome"]
BOWTIE2_INDEX_DIR = config.get("bowtie2_index_dir")
CONFIGURED_BOWTIE2_INDEX_PREFIX = config.get("bowtie2_index_prefix")
REFERENCE_DIR = config["reference_dir"]
TRIMMED_DIR = config["trimmed_dir"]
ALIGNMENT_DIR = config["alignment_dir"]
READ_GROUP_QC_DIR = config["read_group_qc_dir"]
SAMPLES_INFO = config["samples_info"]
BOWTIE2_CONFIG = config.get("bowtie2", {})
BOWTIE2_EXTRA_ARGS = BOWTIE2_CONFIG.get("extra_args", "")
ACTIVE_ATAC_PIPELINE = config.get("pipeline") == "atacseq"
BOWTIE2_MANIFEST_HELPER = os.path.join("src", "atac_bowtie2_manifest.py")


def _atac_alignment_run_info(wildcards):
    return SAMPLES_INFO[wildcards.sample]["runs"][wildcards.run]


def _atac_sample_run_bams(wildcards):
    return [
        os.path.join(ALIGNMENT_DIR, "runs", f"{wildcards.sample}__{run}.sorted.bam")
        for run in SAMPLES_INFO[wildcards.sample]["runs"]
    ]


def _atac_sample_run_validations(wildcards):
    return [
        os.path.join(
            READ_GROUP_QC_DIR,
            f"{wildcards.sample}__{run}.read_group_validation.tsv",
        )
        for run in SAMPLES_INFO[wildcards.sample]["runs"]
    ]


def _atac_merged_read_group_args(wildcards):
    arguments = []
    for run_info in SAMPLES_INFO[wildcards.sample]["runs"].values():
        read_group = run_info["read_group"]
        values = [
            read_group["id"],
            read_group["sample"],
            read_group["library"],
            read_group["platform"],
            read_group["platform_unit"],
        ]
        arguments.append("--expected " + " ".join(shlex.quote(value) for value in values))
    return " ".join(arguments)

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
PREBUILT_BOWTIE2_INDEX_MANIFEST = (
    f"{CONFIGURED_BOWTIE2_INDEX_PREFIX}.manifest.json"
    if PREBUILT_BOWTIE2_INDEX
    else None
)

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

if PREBUILT_BOWTIE2_INDEX and ACTIVE_ATAC_PIPELINE:
    if not os.path.isfile(PREBUILT_BOWTIE2_INDEX_MANIFEST):
        raise ValueError(
            "ATAC-seq prebuilt bowtie2_index_prefix requires a provenance manifest. "
            f"Missing manifest: {PREBUILT_BOWTIE2_INDEX_MANIFEST}. Create it from "
            "known index provenance with: pixi run --locked python3 "
            f"{BOWTIE2_MANIFEST_HELPER} create --prefix "
            f"{BOWTIE2_INDEX_PREFIX!r} --reference {REF_GENOME!r} "
            f"--format {BOWTIE2_INDEX_FORMAT} --builder-version <bowtie2 version>"
        )

index_identity = hashlib.sha256(
    f"{os.path.abspath(BOWTIE2_INDEX_PREFIX)}|{BOWTIE2_INDEX_FORMAT}".encode()
).hexdigest()[:12]
validation_scope = "_manifest" if PREBUILT_BOWTIE2_INDEX else ""
BOWTIE2_VALIDATION_MARKER = os.path.join(
    REFERENCE_DIR,
    f"bowtie2_index_{index_identity}{validation_scope}.validated.OK",
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
        index = BOWTIE2_INDEX_FILES,
        manifest = PREBUILT_BOWTIE2_INDEX_MANIFEST if PREBUILT_BOWTIE2_INDEX else [],
        ref_genome = REF_GENOME if PREBUILT_BOWTIE2_INDEX else [],
        helper = BOWTIE2_MANIFEST_HELPER if PREBUILT_BOWTIE2_INDEX else []
    output:
        marker = BOWTIE2_VALIDATION_MARKER
    params:
        prefix = BOWTIE2_INDEX_PREFIX,
        prebuilt = "true" if PREBUILT_BOWTIE2_INDEX else "false",
        index_format = BOWTIE2_INDEX_FORMAT
    threads: 1
    log:
        os.path.join("logs", config["pipeline"], "bowtie2_validate", "bowtie2_validate.log")
    shell:
        """
        set -euo pipefail
        if {params.prebuilt}; then
            pixi run python3 {input.helper:q} validate \
                --manifest {input.manifest:q} \
                --prefix {params.prefix:q} \
                --reference {input.ref_genome:q} \
                --format {params.index_format:q} --check-dictionary \
                > {log}.out 2> {log}.err
        else
            : > {log}.out
            : > {log}.err
        fi
        pixi run bowtie2-inspect -n {params.prefix:q} >> {log}.out 2>> {log}.err
        touch {output.marker:q}
        """


rule bowtie2_align_run_pe:
    """
    Aligns trimmed paired-end reads with Bowtie2 and sorts the resulting BAM.
    """
    input:
        r1 = os.path.join(TRIMMED_DIR, "{sample}__{run}_R1.trimmed.fq.gz"),
        r2 = os.path.join(TRIMMED_DIR, "{sample}__{run}_R2.trimmed.fq.gz"),
        index = BOWTIE2_INDEX_FILES,
        index_validated = BOWTIE2_VALIDATION_MARKER
    output:
        bam = os.path.join(ALIGNMENT_DIR, "runs", "{sample}__{run}.sorted.bam")
    params:
        extra = BOWTIE2_EXTRA_ARGS,
        index_prefix = BOWTIE2_INDEX_PREFIX,
        rg_id = lambda wildcards: _atac_alignment_run_info(wildcards)["read_group"]["id"],
        rg_sm = lambda wildcards: _atac_alignment_run_info(wildcards)["read_group"]["sample"],
        rg_lb = lambda wildcards: _atac_alignment_run_info(wildcards)["read_group"]["library"],
        rg_pl = lambda wildcards: _atac_alignment_run_info(wildcards)["read_group"]["platform"],
        rg_pu = lambda wildcards: _atac_alignment_run_info(wildcards)["read_group"]["platform_unit"],
        rg_sm_arg = lambda wildcards: f'SM:{_atac_alignment_run_info(wildcards)["read_group"]["sample"]}',
        rg_lb_arg = lambda wildcards: f'LB:{_atac_alignment_run_info(wildcards)["read_group"]["library"]}',
        rg_pl_arg = lambda wildcards: f'PL:{_atac_alignment_run_info(wildcards)["read_group"]["platform"]}',
        rg_pu_arg = lambda wildcards: f'PU:{_atac_alignment_run_info(wildcards)["read_group"]["platform_unit"]}',
        bowtie_threads = lambda wildcards, threads: max(1, threads // 2),
        sort_extra_threads = lambda wildcards, threads: max(0, threads - max(1, threads // 2) - 1)
    threads: 4
    log:
        out=os.path.join("logs", config["pipeline"], "bowtie2_align", "{sample}__{run}.log.out"),
        err=os.path.join("logs", config["pipeline"], "bowtie2_align", "{sample}__{run}.log.err")
    shell:
        """
        set -euo pipefail
        (
          BOWTIE_CMD=(pixi run bowtie2 -p {params.bowtie_threads} {params.extra}
            -x {params.index_prefix:q} --rg-id {params.rg_id:q}
            --rg {params.rg_sm_arg:q} --rg {params.rg_lb_arg:q}
            --rg {params.rg_pl_arg:q} --rg {params.rg_pu_arg:q}
            -1 {input.r1:q} -2 {input.r2:q})
          if [ {threads} -eq 1 ]; then
            TMP_SAM={output.bam:q}.alignment.tmp.sam
            trap 'rm -f "$TMP_SAM"' EXIT
            "${{BOWTIE_CMD[@]}}" -S "$TMP_SAM"
            pixi run samtools sort -@ 0 -o {output.bam:q} "$TMP_SAM"
          else
            "${{BOWTIE_CMD[@]}}" | \
            pixi run samtools sort -@ {params.sort_extra_threads} -o {output.bam:q} -
          fi
        ) > {log.out:q} 2> {log.err:q}
        """


rule validate_bowtie2_run_read_group:
    """
    Validates the ATAC BAM read-group header and every alignment RG tag.
    """
    input:
        bam = os.path.join(ALIGNMENT_DIR, "runs", "{sample}__{run}.sorted.bam"),
        validator = "src/validate_bam_read_group.py"
    output:
        report = os.path.join(
            READ_GROUP_QC_DIR, "{sample}__{run}.read_group_validation.tsv"
        )
    params:
        rg_id = lambda wildcards: _atac_alignment_run_info(wildcards)["read_group"]["id"],
        rg_sm = lambda wildcards: _atac_alignment_run_info(wildcards)["read_group"]["sample"],
        rg_lb = lambda wildcards: _atac_alignment_run_info(wildcards)["read_group"]["library"],
        rg_pl = lambda wildcards: _atac_alignment_run_info(wildcards)["read_group"]["platform"],
        rg_pu = lambda wildcards: _atac_alignment_run_info(wildcards)["read_group"]["platform_unit"]
    log:
        os.path.join(
            "logs", config["pipeline"], "read_group_validation", "{sample}__{run}.log"
        )
    shell:
        """
        set -euo pipefail
        pixi run samtools view -h {input.bam:q} | \
        pixi run python3 src/validate_bam_read_group.py \
            --expected-id {params.rg_id:q} \
            --expected-sample {params.rg_sm:q} \
            --expected-library {params.rg_lb:q} \
            --expected-platform {params.rg_pl:q} \
            --expected-platform-unit {params.rg_pu:q} \
            > {output.report:q} 2> {log:q}.err
        """


rule samtools_merge_runs:
    """
    Merges validated coordinate-sorted run BAMs for one biological sample.
    """
    input:
        bams=_atac_sample_run_bams,
        read_group_validations=_atac_sample_run_validations,
        validator="src/validate_merged_bam_read_groups.py"
    output:
        bam=os.path.join(ALIGNMENT_DIR, "{sample}_pe.sorted.bam"),
        read_group_validation=os.path.join(
            READ_GROUP_QC_DIR, "{sample}_merged.read_group_validation.tsv"
        )
    params:
        expected_read_groups=_atac_merged_read_group_args
    threads: 4
    log:
        os.path.join("logs", config["pipeline"], "samtools_merge_runs", "{sample}.log")
    shell:
        """
        set -euo pipefail
        pixi run samtools merge -f -@ {threads} -o {output.bam:q} {input.bams:q} \
            > {log:q}.out 2> {log:q}.err
        pixi run samtools quickcheck -v {output.bam:q} >> {log:q}.out 2>> {log:q}.err
        pixi run samtools view -h {output.bam:q} | \
        pixi run python3 src/validate_merged_bam_read_groups.py \
            {params.expected_read_groups} > {output.read_group_validation:q} \
            2>> {log:q}.err
        """
