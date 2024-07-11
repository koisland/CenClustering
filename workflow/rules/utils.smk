rule wget:
    output:
        "",
    params:
        url="",
    resources:
        mem=1,
    log:
        "",
    shell:
        """
        wget --no-verbose {params.url} -O {output} 2> {log}
        """


rule convert_bbed_to_bed:
    input:
        "",
    output:
        "",
    conda:
        "../env/tools.yaml"
    resources:
        mem=4,
    log:
        "",
    shell:
        """
        bigBedToBed {input} {output} 2> {log}
        """
