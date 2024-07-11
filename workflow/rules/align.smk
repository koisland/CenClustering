# Generate table
ASM_TABLE = f"/tmp/table_{hash(workflow.basedir)}.asm.tbl"
with open(ASM_TABLE, "wt") as tbl_fh:
    tbl_fh.write("sample\tasm\n")
    tbl_fh.write(f"mPanTro3\t{rules.get_chimp_asm.output}\n")


ALN_CFG = {
    "ref": {"CHM13": rules.get_chm13_asm.output},
    "tbl": ASM_TABLE,
    "aln_threads": 4,
    # https://github.com/koisland/asm-to-reference-alignment/blob/remove_sm_num_index/config/clint.yaml
    "mm2_opts": "-x asm20 --secondary=no -s 25000 -K 8G",
    "second_aln": "no",
    ## THE FOLLOWING OPTIONS ARE ONLY USED FOR GENE CONVERSION ANALYSIS AND NOT ALIGNMENT AND CAN BE IGNORRED
    # subset the gene conversion analysis to just these regions on the reference
    # bed: /net/eichler/vol26/projects/chm13_t2t/nobackups/Assembly_analysis/SEDEF/chm13_v1.1_plus38Y.SDs.bed
    "break_paf": 10000,
    "bed": "config/SDs.and.lowid.bed",
    "window": 1000,
    "slide": 100,
    "min_aln_len": 1000000,
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
            branch="remove_sm_num_index",
        )
    config:
        ALN_CFG


use rule * from align_asm_to_ref exclude dipcall, ideogram as asm_ref_*


# Override rules from ^ because input getter functions run before concat_asm so cannot detect file.
use rule alignment from align_asm_to_ref as asm_ref_alignment with:
    input:
        ref=rules.get_chm13_asm.output,
        query=rules.get_chimp_asm.output,
    resources:
        mem=30,


use rule alignment2 from align_asm_to_ref as asm_ref_alignment2 with:
    input:
        ref_fasta=rules.get_chm13_asm.output,
        query=rules.get_chimp_asm.output,
        # Weird. Blank if passing reference rule output from above. Use string instead.
        aln="temp/{ref}/{sm}.bam",


use rule all from align_asm_to_ref as align_asm_ref_all with:
    input:
        ra=rules.asm_ref_reference_alignment.input,
        SafFire=expand(
            rules.asm_ref_SafFire.output,
            sm=["mPanTro3"],
            ref=["CHM13"],
        ),