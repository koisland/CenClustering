include: "utils.smk"


use rule wget as get_chimp_rm with:
    output:
        "data/tracks/mPanTro3.rm.bb",
    log:
        "logs/get_rm_mPanTro3.log",
    params:
        url="https://genomeark.s3.amazonaws.com/species/Pan_troglodytes/mPanTro3/assembly_curated/repeats/mPanTro3_v2.0.RepeatMasker_v1.2.bb",


use rule wget as get_chimp_censat with:
    output:
        "data/tracks/mPanTro3.censat.bb",
    log:
        "logs/get_censat_mPanTro3.log",
    params:
        url="https://genomeark.s3.amazonaws.com/species/Pan_troglodytes/mPanTro3/assembly_curated/repeats/mPanTro3_v2.0_CenSat_v1.2.bb",


use rule wget as get_chimp_sd with:
    output:
        "data/tracks/mPanTro3.sd.bb",
    log:
        "logs/get_sd_mPanTro3.log",
    params:
        url="https://genomeark.s3.amazonaws.com/species/Pan_troglodytes/mPanTro3/assembly_curated/repeats/mPanTro3_v2.0.SD_v1.0.bb",


# use rule wget as get_chimp_cen with:
#     output:
#         "data/tracks/mPanTro3.hap1hap2.ALR.500kbp.bed",
#     log:
#         "logs/get_cen_mPanTro3.log",
#     params:
#         url="https://eichlerlab.gs.washington.edu/help/glogsdon/Shared_with_T2T_primates_CenSat/centromeric_regions/mPanTro3.hap1hap2.ALR.500kbp.bed",


use rule convert_bbed_to_bed as cvt_chimp_rm_bbed with:
    input:
        rules.get_chimp_rm.output,
    output:
        "data/tracks/mPanTro3.rm.bed",
    log:
        "logs/cvt_chimp_rm_bbed.log",


use rule convert_bbed_to_bed as cvt_chimp_sd_bbed with:
    input:
        rules.get_chimp_sd.output,
    output:
        "data/tracks/mPanTro3.sd.bed",
    log:
        "logs/cvt_chimp_sd_bbed.log",


use rule convert_bbed_to_bed as cvt_chimp_censat_bbed with:
    input:
        rules.get_chimp_censat.output,
    output:
        "data/tracks/mPanTro3.censat.bed",
    log:
        "logs/cvt_chimp_censat_bbed.log",


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
