import os


# Prefer an ignored machine-local config without recursively merging placeholders.
CONFIG_FILE = os.environ.get(
    "HYDRA_CONFIG",
    "config/config.local.yaml"
    if os.path.isfile("config/config.local.yaml")
    else "config/config.yaml",
)
configfile: CONFIG_FILE

if config.get("pipeline") not in ("rnaseq", "wgbs", "chip_cr", "atacseq", "amplicon"):
    raise ValueError("config.pipeline must be rnaseq, wgbs, chip_cr, atacseq, or amplicon.")

# Import modules for different pipelines
rnaseq_config = config["rnaseq"].copy()
rnaseq_config["pipeline"] = config["pipeline"]
module rnaseq:
    snakefile: "workflows/rnaseq.smk"
    config: rnaseq_config
    #prefix: "rnaseq_"

wgbs_config = config["wgbs"].copy()
wgbs_config["pipeline"] = config["pipeline"]
module wgbs:
    snakefile: "workflows/wgbs.smk"
    config: wgbs_config
    #prefix: "wgbs_"

chip_cr_config = config["chip_cr"].copy()
chip_cr_config["pipeline"] = config["pipeline"]
module chip_cr:
    snakefile: "workflows/chip_cr.smk"
    config: chip_cr_config

atacseq_config = config["atacseq"].copy()
atacseq_config["pipeline"] = config["pipeline"]
module atacseq:
    snakefile: "workflows/atacseq.smk"
    config: atacseq_config

amplicon_config = config.get("amplicon", {}).copy()
amplicon_config["pipeline"] = config["pipeline"]
module amplicon:
    snakefile: "workflows/amplicon.smk"
    config: amplicon_config

# Conditional rule import and rule_all definition
if config["pipeline"] == "rnaseq":
    use rule * from rnaseq as rnaseq_*
    use rule all from rnaseq as rnaseq_all

    rule all:
        default_target: True
        input: rules.rnaseq_all.input
elif config["pipeline"] == "wgbs":
    use rule * from wgbs as wgbs_*
    use rule all from wgbs as wgbs_all

    # Make it the default
    rule all:
        default_target: True
        input: rules.wgbs_all.input
elif config["pipeline"] == "chip_cr":
    use rule * from chip_cr as chip_cr_*
    use rule all from chip_cr as chip_cr_all
    
    rule all:
        default_target: True
        input: rules.chip_cr_all.input
elif config["pipeline"] == "atacseq":
    use rule * from atacseq as atacseq_*
    use rule all from atacseq as atacseq_all

    rule all:
        default_target: True
        input: rules.atacseq_all.input
elif config["pipeline"] == "amplicon":
    use rule * from amplicon as amplicon_*

    rule all:
        default_target: True
        input: rules.amplicon_all.input
