rule multiBamSummary_generic:
    """
    Generic rule for deepTools multiBamSummary. Override I/O in parent workflow.
    """
    input:
        bams=["path/to/input.bam"],
        bais=["path/to/input.bam.bai"]
    output:
        npz="path/to/output.npz"
    params:
        extra=""
    threads: 8
    log:
        "logs/multiBamSummary.log"
    shell:
        "pixi run multiBamSummary bins -b {input.bams} -o {output.npz} -p {threads} {params.extra} > {log}.out 2> {log}.err"


rule plotCorrelation_generic:
    """
    Generic rule for deepTools plotCorrelation. Override I/O in parent workflow.
    """
    input:
        npz="path/to/input.npz"
    output:
        heatmap="path/to/output.png",
        matrix="path/to/output.tab"
    params:
        extra=""
    threads: 1
    log:
        "logs/plotCorrelation.log"
    shell:
        "pixi run plotCorrelation -in {input.npz} -o {output.heatmap} --outFileCorMatrix {output.matrix} {params.extra} > {log}.out 2> {log}.err"


rule plotFingerprint_generic:
    """
    Generic rule for deepTools plotFingerprint. Override I/O in parent workflow.
    """
    input:
        bams=["path/to/input.bam"],
        bais=["path/to/input.bam.bai"]
    output:
        plot="path/to/output.png",
        metrics="path/to/output.tab"
    params:
        extra=""
    threads: 8
    log:
        "logs/plotFingerprint.log"
    shell:
        "pixi run plotFingerprint -b {input.bams} -o {output.plot} --outRawCounts {output.metrics} -p {threads} {params.extra} > {log}.out 2> {log}.err"
