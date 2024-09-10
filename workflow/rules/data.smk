include: "utils.smk"


use rule wget as get_chm13_asm with:
    output:
        "data/asm/T2T-CHM13.fa.gz",
    log:
        "logs/get_asm_chm13v2.0.log",
    params:
        url="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz",


use rule wget as get_chimp_asm with:
    output:
        "data/asm/mPanTro3.fa.gz",
    log:
        "logs/get_asm_mPanTro3.log",
    params:
        url="https://genomeark.s3.amazonaws.com/species/Pan_troglodytes/mPanTro3/assembly_curated/mPanTro3.analysis-dip.20231122.fasta.gz",
