# workflows/rules/mbias.smk
import os

# --- CONFIGURATION ---
SORTED_FILTERED_BAM_DIR = config["sorted_filtered_bam_dir"]
REF_GENOME = config["ref_genome"]
MBIAS_DIR = config["mbias_dir"]
SAMPLES = list(config['samples_info'].keys())
METHYLDACKEL_MBIAS_EXTRA_ARGS = config.get("methyldackel_mbias", {}).get("extra_args", "")

# --- RULES ---

rule methyldackel_mbias:
    """
    Runs MethylDackel mbias on a sorted and filtered BAM file to determine methylation bias.
    The suggested options for MethylDackel extract are saved to a file.
    """
    input:
        bam=os.path.join(SORTED_FILTERED_BAM_DIR, "{sample}_pe.filtered.sorted.bam"),
        ref=REF_GENOME
    output:
        mbias=os.path.join(MBIAS_DIR, "{sample}.mbias.txt"),
        options=os.path.join(MBIAS_DIR, "{sample}.options.txt")
    params:
        extra=METHYLDACKEL_MBIAS_EXTRA_ARGS,
        prefix=os.path.join(MBIAS_DIR, "{sample}")
    threads: 16
    log:
        os.path.join("logs", config["pipeline"], "methyldackel_mbias", "{sample}.log")
    shell:
        """
        pixi run MethylDackel mbias -@ {threads} {params.extra} {input.ref} {input.bam} {params.prefix}  > {output.mbias} 2> {log}
        # The options are needed for the next step, so we extract them from the log.
        # If grep fails, it will have a non-zero exit code, so we add `|| true` to prevent the script from failing
        options=$(grep "Suggested inclusion options:" {log} | sed 's/Suggested inclusion options: //') || true
        echo -n "${{options}}" > {output.options}
        """

rule aggregate_mbias_options:
    """
    Aggregates the mbias options from all samples into a single TSV file.
    """
    input:
        expand(os.path.join(MBIAS_DIR, "{sample}.options.txt"), sample=SAMPLES)
    output:
        tsv=os.path.join(MBIAS_DIR, "all_samples_mbias_options.tsv")
    shell:
        """
        echo -e 'sample\toptions' > {output.tsv}
        for f in {input}; do
            sample=$(basename $f .options.txt)
            options=$(cat $f)
            echo -e "$sample\t$options" >> {output.tsv}
        done
        """
