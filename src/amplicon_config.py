"""Validate machine-independent amplicon configuration without filesystem writes."""

import copy
import csv
import os
import re
from pathlib import Path


def integer(value, name, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError(f"amplicon.{name} must be an integer >= {minimum}.")
    return value


def pair(value, name, minimum=0):
    if not isinstance(value, list) or len(value) != 2:
        raise ValueError(f"amplicon.{name} must contain [forward, reverse].")
    return [integer(v, name, minimum) for v in value]


def discover(directory):
    samples = {}
    pattern = re.compile(r"(.+?)[_.-](?:R?([12]))(?:_001)?\.(?:fastq|fq)(?:\.gz)?$")
    for path in sorted(Path(directory).iterdir()):
        if not path.is_file() or not re.search(r"\.(fastq|fq)(\.gz)?$", path.name):
            continue
        match = pattern.fullmatch(path.name)
        sample, mate = (match[1], f"R{match[2]}") if match else (re.sub(r"\.(fastq|fq)(\.gz)?$", "", path.name), "R1")
        info = samples.setdefault(sample, {})
        if mate in info:
            raise ValueError(f"Ambiguous amplicon FASTQs for {sample}/{mate}; use samples_info or sample_sheet.")
        info[mate] = str(path)
    return samples


def normalize(config):
    config = copy.deepcopy(config)
    config.setdefault("output_dir", "results/amplicon")
    config.setdefault("log_dir", "logs/amplicon")
    for key in ("output_dir", "log_dir"):
        if not isinstance(config[key], str) or not config[key].strip():
            raise ValueError(f"amplicon.{key} must be a nonempty path.")
    config["threads"] = integer(config.get("threads", 4), "threads", 1)
    raw = config.get("samples_info") or {}
    sheet = config.get("sample_sheet")
    if sheet and raw:
        raise ValueError("Use either amplicon.sample_sheet or samples_info, not both.")
    if sheet:
        with open(sheet, newline="") as handle:
            rows = csv.DictReader(handle, delimiter="\t")
            if not {"sample", "R1"}.issubset(rows.fieldnames or []):
                raise ValueError("amplicon.sample_sheet needs tab-separated sample, R1, and optional R2, mode, batch columns.")
            for row in rows:
                name = row.pop("sample")
                if name in raw:
                    raise ValueError(f"Duplicate amplicon sample in sample_sheet: {name}")
                raw[name] = {k: v for k, v in row.items() if v}
    if not raw:
        directory = config.get("raw_fastqs_dir")
        if directory and os.path.isdir(directory):
            raw = discover(directory)
    if not isinstance(raw, dict) or not raw:
        raise ValueError("No amplicon samples found. Set samples_info, sample_sheet, or raw_fastqs_dir.")
    defaults = {
        "trunc_len": [0, 0], "trim_left": [0, 0], "max_ee": [2, 2],
        "trunc_q": 2, "min_len": 50, "min_overlap": 12, "max_mismatch": 0,
        "chimera_method": "consensus", "pool": False, "nbases": 100000000,
        "seed": 1,
    }
    supplied = config.get("dada2", {})
    if not isinstance(supplied, dict) or set(supplied) - set(defaults):
        raise ValueError("Unknown amplicon.dada2 options; see config/config.yaml.")
    config["dada2"] = {**defaults, **supplied}
    d = config["dada2"]
    for key in ("trunc_len", "trim_left"):
        d[key] = pair(d[key], f"dada2.{key}")
    if not isinstance(d["max_ee"], list) or len(d["max_ee"]) != 2 or any(type(v) not in (int, float) or not 0 <= v < float("inf") for v in d["max_ee"]):
        raise ValueError("amplicon.dada2.max_ee needs two finite nonnegative numbers.")
    for key in ("trunc_q", "min_len", "min_overlap", "max_mismatch", "nbases", "seed"):
        d[key] = integer(d[key], f"dada2.{key}", 1 if key in ("min_len", "min_overlap", "nbases") else 0)
    if d["chimera_method"] not in ("consensus", "pooled", "per-sample"):
        raise ValueError("amplicon.dada2.chimera_method must be consensus, pooled, or per-sample.")
    if d["pool"] not in (False, True, "pseudo") or type(d["pool"]) not in (bool, str):
        raise ValueError("amplicon.dada2.pool must be a boolean or 'pseudo'.")
    samples = {}
    for sample, info in raw.items():
        if not isinstance(sample, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sample):
            raise ValueError(f"Unsafe amplicon sample ID: {sample!r}.")
        if not isinstance(info, dict) or not info.get("R1"):
            raise ValueError(f"amplicon sample {sample} needs R1 (or has an orphan R2).")
        info = copy.deepcopy(info)
        for mate in ("R1", "R2"):
            if info.get(mate) and (not isinstance(info[mate], str) or not os.path.isfile(info[mate])):
                raise ValueError(f"amplicon sample {sample}: {mate} is not an existing FASTQ file.")
        has_pair = bool(info.get("R2"))
        if info.get("type") not in (None, "PE" if has_pair else "SE"):
            raise ValueError(f"amplicon sample {sample}: declared type disagrees with R1/R2 layout.")
        info["mode"] = info.get("mode", config.get("mode", "auto"))
        if info["mode"] == "auto":
            info["mode"] = "merge" if has_pair else "forward"
        if info["mode"] not in ("merge", "forward", "reverse") or (not has_pair and info["mode"] != "forward"):
            raise ValueError(f"amplicon sample {sample}: merge/reverse need paired inputs; mode must be auto, merge, forward, or reverse.")
        info["batch"] = info.get("batch", "default")
        if not isinstance(info["batch"], str) or not info["batch"].strip():
            raise ValueError(f"amplicon sample {sample} needs a nonempty batch identifier.")
        if not isinstance(config.get("primers", {}), dict) or not isinstance(info.get("primers", {}), dict):
            raise ValueError(f"amplicon sample {sample}: primers must be a mapping.")
        primers = {**config.get("primers", {}), **info.get("primers", {})}
        if set(primers) - {"mode", "forward", "reverse", "error_rate", "min_overlap"}:
            raise ValueError(f"Unknown primer options for amplicon sample {sample}.")
        if primers.get("mode") not in ("remove", "already_removed"):
            raise ValueError(f"amplicon sample {sample}: explicitly set primers.mode to remove or already_removed.")
        if primers["mode"] == "remove":
            for mate in ("forward", "reverse") if has_pair or primers.get("reverse") else ("forward",):
                if not isinstance(primers.get(mate), str) or not re.fullmatch(r"[ACGTRYSWKMBDHVN]+", primers[mate], re.I):
                    raise ValueError(f"amplicon sample {sample} needs a {mate} IUPAC primer in read orientation.")
        rate = primers.get("error_rate", 0.1)
        if type(rate) not in (int, float) or not 0 <= rate < 1:
            raise ValueError("amplicon.primers.error_rate must be in [0, 1).")
        primers["error_rate"] = rate
        primers["min_overlap"] = integer(primers.get("min_overlap", 10), "primers.min_overlap", 1)
        if primers["mode"] == "remove" and any(primers["min_overlap"] > len(primers[m]) for m in ("forward", "reverse") if primers.get(m)):
            raise ValueError(f"amplicon sample {sample}: primers.min_overlap exceeds primer length.")
        info["primers"] = primers
        overrides = info.get("filter", {})
        filter_keys = {"trunc_len", "trim_left", "max_ee", "min_len", "trunc_q"}
        if not isinstance(overrides, dict) or set(overrides) - filter_keys:
            raise ValueError(f"Unsupported per-sample filter for {sample}.")
        # Validate overrides using the same contract as global options.
        filters = {k: copy.deepcopy(d[k]) for k in filter_keys}
        filters.update(overrides)
        for key in ("trunc_len", "trim_left"):
            filters[key] = pair(filters[key], f"samples_info.{sample}.filter.{key}")
        if not isinstance(filters["max_ee"], list) or len(filters["max_ee"]) != 2 or any(type(v) not in (int, float) or not 0 <= v < float("inf") for v in filters["max_ee"]):
            raise ValueError(f"Invalid max_ee override for {sample}.")
        for key in ("min_len", "trunc_q"):
            filters[key] = integer(filters[key], f"samples_info.{sample}.filter.{key}", 1 if key == "min_len" else 0)
        info["filter"] = filters
        samples[sample] = info
    config["samples_info"] = dict(sorted(samples.items()))
    if len({s["mode"] for s in samples.values()}) > 1:
        raise ValueError("Amplicon ASVs must cover a comparable region: do not mix merge, forward, and reverse modes in one cohort. Select forward for all samples or use separate runs.")
    if not isinstance(config.get("taxonomy", {}), dict):
        raise ValueError("amplicon.taxonomy must be a mapping.")
    taxonomy = {"enabled": False, "training_set": None, "reference_version": None, "min_boot": 50, "try_rc": True, **config.get("taxonomy", {})}
    if set(taxonomy) - {"enabled", "training_set", "reference_version", "min_boot", "try_rc"}:
        raise ValueError("Unknown amplicon.taxonomy options.")
    for key in ("enabled", "try_rc"):
        if type(taxonomy[key]) is not bool:
            raise ValueError(f"amplicon.taxonomy.{key} must be boolean.")
    taxonomy["min_boot"] = integer(taxonomy["min_boot"], "taxonomy.min_boot")
    if taxonomy["min_boot"] > 100:
        raise ValueError("amplicon.taxonomy.min_boot must be <= 100.")
    if taxonomy["enabled"]:
        if not taxonomy["training_set"] or not taxonomy["reference_version"]:
            raise ValueError("Enabled amplicon taxonomy needs training_set and reference_version.")
    config["taxonomy"] = taxonomy
    return config
