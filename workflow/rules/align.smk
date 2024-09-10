include: "data.smk"


ALN_CFG = {
    "ref": {"CHM13": rules.get_chm13_asm.output},
    "sm": {"mPanTro3": rules.get_chimp_asm.output},
    "aln_threads": 24,
    # https://github.com/koisland/asm-to-reference-alignment/blob/remove_sm_num_index/config/clint.yaml
    "mm2_opts": "-x asm20 --secondary=no -s 25000 -K 8G",
}


# Align assemblies to reference.
# Pull alignment workflow from github.
# See output files here. Most of everything gets dumped in subdirs from results/{ref}/{ftype}/{asm}.{ftype}:
# https://github.com/mrvollger/asm-to-reference-alignment/blob/main/workflow/rules/reference_alignment.smk
module align_asm_to_ref:
    snakefile:
        github(
            "koisland/asm-to-reference-alignment",
            path="workflow/Snakefile",
            branch="minimal",
        )
    config:
        ALN_CFG


use rule * from align_asm_to_ref as asm_ref_*


rule alignment_idx:
    input:
        bam=rules.asm_ref_alignment.output,
    output:
        bam="results/{ref}/bam/{sm}.bam",
        bam_idx="results/{ref}/bam/{sm}.bam.bai",
    conda:
        "../env/tools.yaml"
    log:
        "logs/asm_{sm}_{ref}_alignment_idx.log",
    shell:
        """
        samtools sort {input.bam} -o {output.bam} 2> {log}
        samtools index {output.bam} 2>> {log}
        """


rule align_all:
    input:
        expand(
            rules.asm_ref_all.input, ref=ALN_CFG["ref"].keys(), sm=ALN_CFG["sm"].keys()
        ),
        expand(
            rules.alignment_idx.output,
            ref=ALN_CFG["ref"].keys(),
            sm=ALN_CFG["sm"].keys(),
        ),
