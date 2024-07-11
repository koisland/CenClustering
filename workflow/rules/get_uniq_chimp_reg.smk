
rule fmt_chimp_sd_bed:
    input:
        rules.cvt_chimp_sd_bbed.output,
    output:
        "results/tracks/mPanTro3.sd.filtered.bed",
    resources:
        mem=2,
    shell:
        """
        cut -f 1,2,3,4 {input} > {output}
        """


rule fmt_chimp_rm_bed:
    input:
        script="workflow/scripts/filter_rm.py",
        bed=rules.cvt_chimp_rm_bbed.output,
    output:
        "results/tracks/mPanTro3.rm.filtered.bed",
    conda:
        "../env/py.yaml"
    resources:
        mem=2,
    shell:
        """
        python {input.script} -i <(cut -f 14 {input.bed} | sed 's/,/\\n/g') > {output}
        """


rule fmt_chimp_censat_bed:
    input:
        bed=rules.cvt_chimp_censat_bbed.output,
    output:
        "results/tracks/mPanTro3.censat.filtered.bed",
    params:
        allowed_censat_prefixes="|".join(["mon", "ct", "dhor"]),
    resources:
        mem=2,
    shell:
        """
        grep -Pv "{params.allowed_censat_prefixes}" {input} | cut -f 1,2,3,4 > {output}
        """

rule index_chimp_genome:
    input:
        rules.get_chimp_asm.output,
    output:
        "data/asm/mPanTro3.fa.gz.fai"
    log:
        "logs/index_chimp_genome.log"
    conda:
        "../env/tools.yaml"
    shell:
        """
        samtools faidx {input}
        """

rule find_uniq_regions_chimp:
    input:
        chimp_genome_size=rules.index_chimp_genome.output,
        rm_bed=rules.fmt_chimp_rm_bed.output,
        segdup_bed=rules.fmt_chimp_sd_bed.output,
        censat_bed=rules.fmt_chimp_censat_bed.output,
    output:
        uniq_bed="results/uniq/mPanTro3_regions.bed",
    conda:
        "../env/tools.yaml"
    resources:
        mem=4,
    log:
        "logs/find_uniq_regions_chimp.log",
    shell:
        """
        bedtools complement \
        -i <(cat {input.rm_bed} {input.segdup_bed} {input.censat_bed} | sort -k 1,1 -k2,2n) \
        -g <(cut -f 1,2  {input.chimp_genome_size} | sort -k1,1) > {output} 2> {log}
        """
