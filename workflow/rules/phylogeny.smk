ARMS = ["parm", "qarm"]
HEURISTIC = ["1", "2"]
with open("hgsvc3_samples.list", "rt") as fh:
    SAMPLES = [l.strip() for l in fh.readlines()]

wildcard_constraints:
    arm="|".join(ARMS),
    sm="|".join(SAMPLES),
    heuristic="|".join(HEURISTIC)


rule get_ref_regions:
    input:
        ref="/project/logsdon_shared/data/reference/T2T-CHM13v2.fasta",
        bed="chr2_{arm}.bed",
    output:
        fa="results/fastas/CHM13/chr2_{arm}.fa",
    conda:
        "../env/tools.yaml"
    log:
        "logs/get_ref_{arm}_regions.log"
    shell:
        """
        seqtk subseq {input.ref} {input.bed} > {output} 2> {log}
        """


rule get_chimp_regions:
    input:
        bam="results/CHM13/bam/mPanTro3.bam",
        bam_idx="results/CHM13/bam/mPanTro3.bam.bai",
        bed="chr2_{arm}.bed",
    output:
        fa="results/fastas/mPanTro3/{arm}.fa",
    singularity:
        "docker://eichlerlab/subseqfa:1.0"
    log:
        "logs/get_chimp_{arm}_regions.log"
    shell:
        """
        while IFS='' read -r line; do
            rgn=$(echo $line | awk '{{print $1":"$2"-"$3}}')
            subseqfa -b -v -r "${{rgn}}" {input.bam} >> {output} 2> {log}
        done < {input.bed}
        """


use rule get_chimp_regions as get_sm_regions with:
    input:
        bam="/project/logsdon_shared/data/hgsvc3/alignments_t2tv2_batch3/{sm}.bam",
        bam_idx="/project/logsdon_shared/data/hgsvc3/alignments_t2tv2_batch3/{sm}.bam.bai",
        bed="chr2_{arm}.bed",
    output:
        fa="results/fastas/{sm}/{arm}.fa",
    params:
        seq_prefix="{sm}_",
    log:
        "logs/get_{sm}_{arm}_regions.log"


rule add_prefix_to_seq_chimp:
    input:
        rules.get_chimp_regions.output
    output:
        fa="results/fastas/mPanTro3/{arm}_renamed.fa",
    params:
        seq_prefix="mPanTro3_",
        extra="",
    conda:
        "../env/tools.yaml"
    log:
        "logs/add_prefix_to_seq_chimp_{arm}.log"
    shell:
        """
        seqkit replace -p ^ -r "{params.seq_prefix}" {input} {params.extra} > {output}
        """


use rule add_prefix_to_seq_chimp as add_prefix_to_seq_sm with:
    input:
        rules.get_sm_regions.output
    output:
        fa="results/fastas/{sm}/{arm}_renamed.fa",
    log:
        "logs/add_prefix_to_seq_{sm}_{arm}.log"
    params:
        seq_prefix="{sm}_",
        # Strip unnecessary information.
        # Only {sample}_{hap}
        extra=lambda wc: "| seqkit replace -p '^(\S+)_hap\w+(\d)-(.+)$' -r '${1}_hap${2}' | seqkit rmdup"


rule filter_samples_by_list:
    input:
        filter_file="mark_samples.list",
        fa_sm=lambda wc: expand(
            rules.add_prefix_to_seq_sm.output,
            sm=SAMPLES,
            arm=[wc.arm]
        ),
    output:
        fa="results/msa/{arm}_filtered.fa"
    conda:
        "../env/tools.yaml"
    log:
        "logs/filter_samples_by_list_{arm}.log"
    shell:
        """
        seqtk subseq <(cat {input.fa_sm}) {input.filter_file} > {output} 2> {log}
        """


rule merge_regions_by_arm:
    input:
        fa_sm=rules.filter_samples_by_list.output,
        # fa_sm=lambda wc: expand(
        #     rules.add_prefix_to_seq_sm.output,
        #     sm=SAMPLES,
        #     arm=[wc.arm]
        # ),
        fa_chimp=rules.add_prefix_to_seq_chimp.output,
    output:
        fa="results/msa/{arm}.fa",
        faidx="results/msa/{arm}.fa.fai",
    conda:
        "../env/tools.yaml"
    log:
        "logs/merge_regions_by_arm_{arm}.log"
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
        faidx="results/msa/{arm}_{heuristic}.fa.fai",
    params:
        region_len=20_000,
        heuristic=lambda wc: "min" if wc.heuristic == "1" else "max"
    conda:
        "../env/py.yaml"
    log:
        "logs/remove_dupes_{arm}_{heuristic}.log"
    shell:
        """
        python {input.script} -i {input.faidx} --heuristic {params.heuristic} --len {params.region_len} > {output} 2> {log}
        """

rule subset_fa:
    input:
        fa=rules.merge_regions_by_arm.output.fa,
        faidx=rules.remove_dupes.output.faidx
    output:
        fa="results/msa/{arm}_{heuristic}_dedup.fa",
        faidx="results/msa/{arm}_{heuristic}_dedup.fa.fai",
    params:
        omit="mPanTro3_chr12_hap2"
    conda:
        "../env/tools.yaml"
    log:
        "logs/subset_fa_{arm}_{heuristic}.log"
    shell:
        """
        seqtk subseq {input.fa} <(cut -f 1 {input.faidx} | grep -v "{params.omit}" ) > {output.fa} 2> {log}
        samtools faidx {output.fa}
        """

rule msa:
    input:
        fa=rules.subset_fa.output.fa
    output:
        msa="results/msa/{arm}_{heuristic}_msa.fa",
        faidx="results/msa/{arm}_{heuristic}_msa.fa.fai",
    conda:
        "../env/tools.yaml"
    threads:
        20
    log:
        "logs/msa_{arm}_{heuristic}.log"
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
        touch("results/msa/{arm}_{heuristic}.done")
    threads:
        20
    params:
        outgroup_pattern="mPanTro3",
        mode="MFP",
        bootstrap=1000,
    conda:
        "../env/tools.yaml"
    log:
        "logs/generate_tree_{arm}_{heuristic}.log"
    shell:
        """
        iqtree2 -redo -T {threads} -s {input.msa} \
        -m {params.mode} -B {params.bootstrap} \
        -o $(grep "{params.outgroup_pattern}" {input.faidx} | head -n1 | cut -f 1 | sed "s/:/_/g") > {log}
        """


rule phylogeny_all:
    input:
        expand(rules.merge_regions_by_arm.output, arm=ARMS),
        expand(rules.msa.output, arm=ARMS, heuristic=HEURISTIC),
        expand(rules.generate_tree.output, arm=ARMS, heuristic=HEURISTIC),
    default_target: True
