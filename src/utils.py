import os
import re
from collections import defaultdict

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