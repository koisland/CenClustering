# Snakemake-CenPhylogeny
Workflow to generate phylogenetic tree for centromeres with chimp (`mPanTro3`) as an outgroup.

It will:
* Align CHM13 to chimp.
* Extract query sequence mapping to both chimp and CHM13 from a provided alignment of an assembly to CHM13
* Run a multiple-sequence alignment of extracted sequences.
* Generate a phylogeny tree.

Only chr2 and chr20 supported currently.

### Usage
```bash
snakemake -p -c20 --use-conda --use-singularity --singularity-args "--bind /project/logsdon_shared/data" --rerun-incomplete
```
