
rule get_ref_regions:
    input:
        ref=rules.get_chm13_asm.output,
        bed=join(config["chr_arm_bed_dir"], "{chrom}_{arm}.bed"),
    output:
        fa="results/fastas/CHM13/{chrom}_{arm}.fa",
    conda:
        "../env/tools.yaml"
    log:
        "logs/get_ref_{chrom}_{arm}_regions.log",
    shell:
        """
        seqtk subseq {input.ref} {input.bed} > {output} 2> {log}
        """


rule get_chimp_regions:
    input:
        bam=expand(rules.alignment_idx.output.bam, ref="CHM13", sm="mPanTro3"),
        bam_idx=expand(rules.alignment_idx.output.bam_idx, ref="CHM13", sm="mPanTro3"),
        bed=join(config["chr_arm_bed_dir"], "{chrom}_{arm}.bed"),
    output:
        fa="results/fastas/mPanTro3/{chrom}_{arm}.fa",
    singularity:
        "docker://eichlerlab/subseqfa:1.0"
    log:
        "logs/get_chimp_{chrom}_{arm}_regions.log",
    shell:
        """
        while IFS='' read -r line; do
            rgn=$(echo $line | awk '{{print $1":"$2"-"$3}}')
            subseqfa -b -v -r "${{rgn}}" {input.bam} >> {output} 2> {log}
        done < {input.bed}
        """


use rule get_chimp_regions as get_sm_regions with:
    input:
        bam=lambda wc: ALIGNMENTS[str(wc.sm)][str(wc.asm)],
        bam_idx=lambda wc: f"{ALIGNMENTS[str(wc.sm)][str(wc.asm)]}.bai",
        bed=join(config["chr_arm_bed_dir"], "{chrom}_{arm}.bed"),
    output:
        fa="results/fastas/{sm}/{chrom}_{arm}_{asm}.fa",
    params:
        seq_prefix="{sm}_",
    log:
        "logs/get_{sm}_{chrom}_{arm}_{asm}_regions.log",


rule add_prefix_to_seq_chimp:
    input:
        rules.get_chimp_regions.output,
    output:
        fa="results/fastas/mPanTro3/{chrom}_{arm}_renamed.fa",
    params:
        seq_prefix="mPanTro3_",
        extra="",
    conda:
        "../env/tools.yaml"
    log:
        "logs/add_prefix_to_seq_chimp_{chrom}_{arm}.log",
    shell:
        """
        seqkit replace -p ^ -r "{params.seq_prefix}" {input} {params.extra} > {output}
        """


use rule add_prefix_to_seq_chimp as add_prefix_to_seq_sm with:
    input:
        rules.get_sm_regions.output,
    output:
        fa="results/fastas/{sm}/{chrom}_{arm}_{asm}_renamed.fa",
    log:
        "logs/add_prefix_to_seq_{sm}_{chrom}_{arm}_{asm}.log",
    params:
        seq_prefix="{sm}_",
        # Strip unnecessary information.
        # Only {sample}_{hap} or {sample}_unassigned
        extra=lambda wc: (
            "| "
            + "sed -e 's/_haplotype/_h/g' -e 's/_h/_hap/g' "
            + "| "
            + r"seqkit replace -p '^(\S+)_(hap\d|unassigned)(.+)$' -r '${1}_${2}"
            + f"_{wc.asm}' "
            + "| "
            + "seqkit rmdup"
        ),


rule merge_regions_by_arm:
    input:
        fa_sm=lambda wc: chain(*[
            expand(rules.add_prefix_to_seq_sm.output, sm=sm, arm=[wc.arm], chrom=[wc.chrom], asm=[ASSEMBLIES[sm]])
            for sm in SAMPLES
        ]),
        fa_chimp=rules.add_prefix_to_seq_chimp.output,
    output:
        fa="results/msa/{chrom}_{arm}.fa",
        faidx="results/msa/{chrom}_{arm}.fa.fai",
    conda:
        "../env/tools.yaml"
    log:
        "logs/merge_regions_by_arm_{chrom}_{arm}.log",
    shell:
        """
        cat {input.fa_sm} {input.fa_chimp} > {output.fa}
        samtools faidx {output.fa} 2> {log}
        """


rule remove_dupes:
    input:
        script="workflow/scripts/dedupe_fai.py",
        faidx=rules.merge_regions_by_arm.output.faidx,
    output:
        faidx="results/msa/{chrom}_{arm}.fa.fai",
    params:
        region_len=20_000,
        heuristic="min",
    conda:
        "../env/py.yaml"
    log:
        "logs/remove_dupes_{chrom}_{arm}.log",
    shell:
        """
        python {input.script} -i {input.faidx} --heuristic {params.heuristic} --len {params.region_len} > {output} 2> {log}
        """


rule subset_fa:
    input:
        fa=rules.merge_regions_by_arm.output.fa,
        faidx=rules.remove_dupes.output.faidx,
    output:
        fa="results/msa/{chrom}_{arm}_dedup.fa",
        faidx="results/msa/{chrom}_{arm}_dedup.fa.fai",
    params:
        omit="mPanTro3_{chrom}_hap2",
    conda:
        "../env/tools.yaml"
    log:
        "logs/subset_fa_{chrom}_{arm}.log",
    shell:
        """
        {{ seqtk subseq {input.fa} <(cut -f 1 {input.faidx} | grep -v "{params.omit}" ) | \
        seqkit rmdup ;}} > {output.fa} 2> {log}
        samtools faidx {output.fa}
        """


rule msa:
    input:
        fa=rules.subset_fa.output.fa,
    output:
        msa="results/msa/{chrom}_{arm}_msa.fa",
        faidx="results/msa/{chrom}_{arm}_msa.fa.fai",
    conda:
        "../env/tools.yaml"
    threads: 20
    log:
        "logs/msa_{chrom}_{arm}.log",
    shell:
        """
        mafft --auto --thread {threads} {input.fa} > {output.msa} 2> {log}
        samtools faidx {output.msa}
        """


rule generate_tree:
    input:
        msa=rules.msa.output.msa,
        faidx=rules.msa.output.faidx,
    output:
        touch("results/msa/{chrom}_{arm}.done"),
    threads: 20
    params:
        outgroup_pattern="mPanTro3",
        mode="MFP",
        bootstrap=1000,
    conda:
        "../env/tools.yaml"
    log:
        "logs/generate_tree_{chrom}_{arm}.log",
    shell:
        """
        iqtree2 -redo -T {threads} -s {input.msa} \
        -m {params.mode} -B {params.bootstrap} \
        -o $(grep "{params.outgroup_pattern}" {input.faidx} | head -n1 | cut -f 1 | sed "s/:/_/g") > {log}
        """


rule phylogeny_all:
    input:
        [
            [
                expand(rules.merge_regions_by_arm.output, chrom=chrom, arm=CHR_ARMS[chrom]),
                expand(rules.msa.output, chrom=chrom, arm=CHR_ARMS[chrom]),
                expand(rules.generate_tree.output, chrom=chrom, arm=CHR_ARMS[chrom])
            ]
            for chrom in CHR_ARMS
        ]
    default_target: True
