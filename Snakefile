# Load config file
configfile: "config/config.yaml"

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

# Conditional rule import and rule_all definition
if config["pipeline"] == "rnaseq":
    use rule * from rnaseq as rnaseq_*
    use rule all from rnaseq as rnaseq_all

    rule all:
        input: rules.rnaseq_all.input
elif config["pipeline"] == "wgbs":
    use rule * from wgbs as wgbs_*
    use rule all from wgbs as wgbs_all

    # Make it the default
    rule all:
        input: rules.wgbs_all.input
elif config["pipeline"] == "chip_cr":
    use rule * from chip_cr as chip_cr_*
    use rule all from chip_cr as chip_cr_all
    
    rule all:
        input: rules.chip_cr_all.input
else:
    rule all:
        input: "ERROR: Write 'rnaseq', 'wgbs', or 'chip_cr' in config.yaml"
