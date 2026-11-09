# GABAinducesPeroxisome

R scripts for the RNA-seq part of a study on how leaf development changes peroxisome activity, and on the
role of GABA in that change.

## Study summary

Published RNA-seq data from successive developmental stages of maize (*Zea mays*) and rice (*Oryza sativa*)
were used to follow the expression of genes in 21 peroxisome metabolic pathways.

- Adult leaves were enriched for photorespiration genes compared with juvenile stages. Genes for ROS/RNS
  metabolism, NADPH regeneration and polyamine catabolism were also enriched at later stages of leaf
  differentiation.
- The gene most often upregulated in differentiated leaves, across all datasets of both species, was
  *BETAINE ALANINE DEHYDROGENASE* (*BADH*). BADH converts 4-aminobutyraldehyde to GABA during polyamine
  catabolism.
- qRT-PCR in developing Einkorn wheat (*Triticum monococcum*) seedlings confirmed that *BADH* and
  *CATALASE3* (*CAT3*) are upregulated in older seedlings.
- Adding GABA to the seedlings increased the number of peroxisomes and upregulated *BADH* and *PEX11C*,
  a peroxin needed for peroxisome fission.
- Conclusion: GABA may help regulate the peroxisome fission machinery during leaf differentiation.

This repository holds only the RNA-seq analysis (maize and rice). The Einkorn qRT-PCR and GABA-treatment
experiments are not part of it.

## Data

| Species | Dataset | Samples |
|---|---|---|
| Rice | GSE155932 | Five developmental stages (S1–S5), two replicates each |
| Rice | GSE57950 (HHZ) | Control vs drought (used for WGCNA) |
| Maize | SRP014652 | Leaf position and stage (8th, 11th, 13th and immature leaves at V9; 13th leaf at VT and R2; topmost, base and tip leaf sections at V3–V7) and kernel time course (0, 6, 12, 30 days after pollination) |
| Maize | GSE151455 | Kernels at 10 vs 65 days after pollination |

Inputs are featureCounts tables named `<dataset>_raw_counts_hisat-sorted.txt` in a `counts/` folder next to
each analysis folder (the count tables are not in this repository).

## Layout

Both `OryzaSativa/` and `ZeaMays/` use the same folders:

| Folder / script | What it does |
|---|---|
| `deseq2/` | DESeq2 differential expression. Each script defines one contrast (for example `GSE155932_deseq_s2_vs_s1.R`, `SRP014652_v_stage_8vs11.R`). Genes with fewer than 100 reads in total are removed; the scripts give the DE gene list (`dn_<dataset>`) and variance-stabilised counts (`vstc_dds_<dataset>`). `ZeaMays/deseq2/time_series_deg.R` tests the leaf series as a time course |
| `peroxisome/`, `per_heatmap/` | Z-score heatmaps of peroxisome genes (rice–maize homolog table) across stages, clustered and non-clustered, some annotated by peroxisome function |
| `autophagy/` | The same heatmaps for autophagy (ATG) homologs |
| `Bargraph_SPGUN.R` (rice), `Bargraph_VPGU.R` (maize) | For each sample, counts the peroxisome genes with above-average expression (z > 0) in each metabolic pathway and draws stacked bar graphs (the pathway representation per stage) |
| `TOPGO/`, `TopGO/` | GO enrichment (topGO, classic Fisher test) of differentially expressed genes |
| `kegg/` | KEGG pathway enrichment with goseq, corrected for gene-length bias |
| `OryzaSativa/wgcna/wgcna.R` | WGCNA co-expression network: keeps the most variable genes after the variance-stabilising transform, finds modules and relates them to condition |

## Running

Requirements: R with DESeq2, pheatmap, ggplot2, ggsci, plyr/dplyr, topGO, goseq, WGCNA.

1. Run the `deseq2/` script for the contrast you want.
2. In the same R session, run the heatmap, bar graph, topGO or KEGG script. These scripts use the objects
   from step 1.

Annotation files (Phytozome annotation tables, MSU7 GO annotation, peroxisome and autophagy homolog
tables) are read from absolute local paths. Change those paths to your own copies before running.
