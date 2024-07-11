
# filter aln bed
rule filter_format_aln_bed:
    input:
        script="workflow/scripts/filter_aln.py",
        aln_bed="results/CHM13/SafFire/mPanTro3.bed"
    output:
        "results/aln/chm13_to_chimp_filtered.bed"
    conda:
        "../env/py.yaml"
    log:
        "logs/filter_format_aln_bed.log"
    shell:
        """
        python {input.script} -i {input.aln_bed} > {output} 2> {log}
        """

# hap1 chr2
# chr2	91700730	91982439

rule get_chm13_hor_array_coords:
    input:
        "data/chm13/chm13_cen_coords.bed"
    output:
        "results/intersect/chm13_cen_coords_minmax.bed"
    conda:
        "../env/tools.yaml"
    log:
        "logs/get_chm13_hor_array_coords.log"
    shell:
        """
        cat {input} | bedtools groupby -g 1 -c 2,3 -o min,max > {output} 2> {log}
        """

rule get_chm13_hor_array_edge_coords:
    input:
        script="workflow/scripts/get_pqarm_edge.py",
        cen_coords=rules.get_chm13_hor_array_coords.output
    output:
        "results/intersect/chm13_cen_coords_edge.bed"
    params:
        bp_edge=1_000_000
    conda:
        "../env/py.yaml"
    log:
        "logs/get_chm13_hor_array_edge_coords.log"
    shell:
        """
        python {input.script} -i {input.cen_coords} -e {params.bp_edge} > {output} 2> {log}
        """

 # intersect uniq regions with ^
rule intersect_uniq_chimp_regions_with_chm13_aln:
    input:
        chm13_edge_bed=rules.get_chm13_hor_array_edge_coords.output,
        chm13_aln_bed=rules.filter_format_aln_bed.output,
        uniq_chimp_bed=rules.find_uniq_regions_chimp.output,
    output:
        "results/intersect/overlapping_uniq_regions_chm13_chimp.bed"
    params:
        uniq_chimp_region_len_thr=10000
    log:
        "logs/intersect_uniq_chimp_regions_with_chm13_aln.log"
    conda:
        "../env/tools.yaml"
    shell:
        """
        {{ bedtools intersect \
            -a {input.chm13_aln_bed} \
            -b <(awk -v OFS="\\t" '{{
                len=$3 - $2;
                if (len > {params.uniq_chimp_region_len_thr}) {{ print }}
            }}' {input.uniq_chimp_bed}) \
            -wa -wb -F 1.0 | \
            sort -k 1 -k 2n | \
            awk -v OFS="\\t" '{{ print $5,$6,$7,$1,$2,$3,$19 }}' | \
            bedtools intersect -wa -wb -a - -b {input.chm13_edge_bed} | \
            sort | \
            uniq \
        ;}} > {output} 2> {log}
        """
# TODO: group by chr and arm. minmax for rows that map correctly to chimp