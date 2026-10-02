import os
import re
import shutil
import subprocess
from collections import defaultdict


def validate_required_path(path, config_name):
    """
    Ensures a required config path value is present.
    """
    if not path:
        raise ValueError(f"Missing '{config_name}' in the config.")


def validate_existing_file(path, config_name):
    """
    Ensures a config path points to an existing file.
    """
    validate_required_path(path, config_name)
    if not os.path.isfile(path):
        raise ValueError(f"Config '{config_name}' must point to an existing file: {path}")


def validate_existing_parent_dir(path, config_name):
    """
    Ensures a config path uses an existing parent directory.
    """
    validate_required_path(path, config_name)
    parent_dir = os.path.dirname(path) or "."
    if not os.path.isdir(parent_dir):
        raise ValueError(
            f"Config '{config_name}' must use an existing parent directory: {parent_dir}"
        )


def validate_sambamba_filter(expression, config_name):
    """
    Validates a Sambamba filter expression with the installed Sambamba parser.
    """
    if not isinstance(expression, str):
        raise ValueError(
            f"Config '{config_name}' must be a Sambamba filter expression string."
        )

    expression = expression.strip()
    if not expression:
        return expression

    sambamba = shutil.which("sambamba")
    if sambamba is None:
        raise RuntimeError(
            f"Sambamba is required to validate config '{config_name}'. "
            "Run Snakemake through 'pixi run'."
        )

    try:
        result = subprocess.run(
            [
                sambamba,
                "view",
                "-S",
                "-c",
                "-F",
                expression,
                os.devnull,
            ],
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=10,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(
            f"Timed out while validating config '{config_name}' with Sambamba."
        ) from exc

    if result.returncode != 0:
        details = [line.strip() for line in result.stderr.splitlines() if line.strip()]
        detail = details[-1] if details else "Sambamba rejected the expression"
        raise ValueError(
            f"Invalid Sambamba filter in config '{config_name}': {detail}. "
            f"Expression: {expression!r}"
        )

    return expression


def discover_samples(raw_fastqs_dir):
    """
    Discovers samples from a directory of raw FASTQ files.
    It scans the raw_fastqs_dir for FASTQ files and groups them by sample name.
    
    Assumes file naming convention like: 
    - {sample_name}_R1.fastq.gz, {sample_name}_R2.fastq.gz
    - {sample_name}_1.fastq.gz, {sample_name}_2.fastq.gz
    - {sample_name}.R1.fastq.gz, {sample_name}.R2.fastq.gz
    - Single-end files: {sample_name}.fastq.gz
    """
    samples = defaultdict(lambda: {'R1': None, 'R2': None})
    
    # Regex to capture sample name, read number, and extension
    # It can handle _R1, _R2, _1, _2, and other common separators.
    pattern = re.compile(r"(.+?)[_.-]?(R[12]|[12])(_001)?\.(fastq|fq)\.gz")

    files = os.listdir(raw_fastqs_dir)
    
    # First pass for paired-end files
    for filename in files:
        match = pattern.match(filename)
        if match:
            sample_name = match.group(1)
            read_identifier = match.group(2)
            read = 'R1' if read_identifier in ['R1', '1'] else 'R2'
            
            if samples[sample_name][read] is None:
                samples[sample_name][read] = os.path.join(raw_fastqs_dir, filename)

    # Second pass for single-end files (those that didn't match PE pattern)
    matched_files = {f for s in samples.values() for f in s.values() if f}
    for filename in files:
        full_path = os.path.join(raw_fastqs_dir, filename)
        if full_path not in matched_files and filename.endswith(('.fastq.gz', '.fq.gz')):
            sample_name = re.sub(r'(\.fastq\.gz|\.fq\.gz)$', '', filename)
            if sample_name not in samples:
                samples[sample_name]['R1'] = full_path

    # Finalize types and convert to dict
    final_samples = {}
    for sample_name, info in samples.items():
        clean_name = sample_name.rstrip('-_.')
        if info['R1'] and info['R2']:
            info['type'] = 'PE'
        elif info['R1']:
            info['type'] = 'SE'
        else:
            continue # Skip if no R1 file was found
        
        # remove None values
        final_info = {k:v for k,v in info.items() if v is not None}
        final_samples[clean_name] = final_info

    return final_samples

def prepare_sample_data(config):
    """
    Retrieves sample information from the config, or discovers it from the filesystem.
    Returns a tuple of (samples_info, samples_list).
    """
    samples_info_from_config = config.get("samples_info", {})
    
    # If samples_info is provided in config, use it directly.
    if samples_info_from_config:
        samples_info = samples_info_from_config
    else:
        # Otherwise, try to discover from raw_fastqs_dir
        raw_fastqs_dir = config.get("raw_fastqs_dir")
        if raw_fastqs_dir and os.path.isdir(raw_fastqs_dir):
            samples_info = discover_samples(raw_fastqs_dir)
        else:
            samples_info = {}
        
    samples = list(samples_info.keys())
    config['samples_info'] = samples_info
    return samples_info, samples


def normalize_read_group_metadata(samples_info, config_name):
    """
    Adds validated Bowtie2 read-group metadata to each sample.

    Existing manifests remain valid by using the sample name for ID, SM, LB,
    and PU, and ILLUMINA for PL. Run-aware manifests should override these
    defaults with a nested ``read_group`` mapping.
    """
    normalized_samples = {}
    read_group_ids = {}
    field_defaults = {
        "id": lambda sample: sample,
        "sample": lambda sample: sample,
        "library": lambda sample: sample,
        "platform": lambda sample: "ILLUMINA",
        "platform_unit": lambda sample: sample,
    }

    for sample, info in samples_info.items():
        normalized_info = dict(info)
        configured_read_group = normalized_info.get("read_group", {})
        if configured_read_group is None:
            configured_read_group = {}
        if not isinstance(configured_read_group, dict):
            raise ValueError(
                f"Config '{config_name}.{sample}.read_group' must be a mapping."
            )

        read_group = {}
        for field, default_factory in field_defaults.items():
            value = configured_read_group.get(field, default_factory(sample))
            if not isinstance(value, str) or not value.strip():
                raise ValueError(
                    f"Config '{config_name}.{sample}.read_group.{field}' "
                    "must be a non-empty string."
                )
            if any(character in value for character in ["\t", "\n", "\r"]):
                raise ValueError(
                    f"Config '{config_name}.{sample}.read_group.{field}' "
                    "must not contain tabs or newlines."
                )
            read_group[field] = value

        read_group_id = read_group["id"]
        if read_group_id in read_group_ids:
            raise ValueError(
                f"Config '{config_name}' assigns read-group ID {read_group_id!r} "
                f"to both {read_group_ids[read_group_id]!r} and {sample!r}. "
                "Use one unique ID per run or lane."
            )
        read_group_ids[read_group_id] = sample

        normalized_info["read_group"] = read_group
        normalized_samples[sample] = normalized_info

    return normalized_samples
