# Population Genomics of Miltefosine Resistance in *Leishmania infantum*

## Project Status
🚧 **Active Development** - Living research project demonstrating ongoing bioinformatics capabilities

## Overview

*Leishmania infantum* causes visceral leishmaniasis, affecting thousands globally. This project investigates the evolutionary dynamics of miltefosine resistance through population genomic analysis of the Miltefosine Sensitivity Locus (MSL) on chromosome 31.

**Central Research Question:** Is the spread of MSL deletions (conferring drug resistance) driven by positive selection or neutral demographic processes (drift/founder effects)?

## Key Findings (Preliminary)

- Clear Old World/New World phylogenetic split
- 12-13 distinct populations identified via ADMIXTURE
- MSL deletions predominantly in Brazilian populations
- Ploidy variation (diploid vs tetraploid) correlates with geography

## Analytical Approach

### 1. Hierarchical Population Structure (ADMIXTURE)
- **All samples + outgroup (K=12):** Species-level structure
- **L. infantum only (K=13):** Fine-scale cryptic populations
- **Americas only (K=13):** High-resolution New World structure

Cross-validation error minimization determined optimal K for each dataset.

### 2. Multi-Scale Phylogenetic Analysis
Comparing phylogenetic signal at three scales to detect selection:

- **Whole-genome phylogenies:** Neutral evolutionary baseline
- **100kb windows around MSL:** Regional selection signal
- **20kb windows around MSL:** Fine-scale locus-specific evolution

Phylogenies stratified by ploidy level (diploid vs tetraploid strains).

**Rationale:** If MSL is under selection, chr31 phylogeny will differ from genome-wide patterns.

### 3. Geographic Distribution Analysis
*Integration in progress*

### 4. Hardy-Weinberg Equilibrium Testing
*Integration in progress* - Tests whether MSL allele frequencies deviate from neutral expectations within populations.

## Clinical Relevance

MSL deletions confer resistance to **miltefosine**, the only oral treatment for leishmaniasis. Understanding whether resistant strains spread via:
- **Selection** (faster than expected, requires intervention)
- **Drift** (random spread, reflects demographic history)

...has direct implications for treatment policy and surveillance strategies.

## Methods Summary

**Phylogenetics:**
- IQ-TREE2 maximum likelihood with 1000 ultrafast bootstrap replicates
- GTR+F+ASC model (ascertainment bias correction for SNP data)
- ggtree (R) for visualization with geographic/ploidy annotation

**Population Structure:**
- ADMIXTURE v1.3 with K=1-20 tested
- Cross-validation error minimization for optimal K selection
- High-confidence assignments: ≥90% ancestry threshold

**Data:**
- 671 *L. infantum* isolates, global sampling
- Whole-genome sequencing data (VCF format)
- Geographic metadata, MSL copy number (coverage-based)

## Repository Structure

```
├── data/
│   ├── metadata/              # Sample information, population assignments
│   ├── population_structure/  # ADMIXTURE input files (.Q files)
│   └── phylogenetic/          # Tree files (.treefile, .contree)
├── scripts/
│   ├── 01_population_structure_admixture.R
│   ├── 02_phylogenetic_whole_genome.R
│   ├── 03_phylogenetic_100kb_windows.R
│   └── 04_phylogenetic_20kb_windows.R
├── results/
│   ├── population_structure/  # Population assignments, summaries
│   ├── phylogenies/           # Tree files
│   └── tables/                # Summary statistics
├── figures/
│   ├── admixture/             # Population structure plots
│   └── phylogenetics/         # Annotated phylogenetic trees
└── docs/
    ├── methods.md             # Detailed methodology
    └── data_sources.md        # Data provenance
```

## Software & Dependencies

**R (v4.4.2):**
- tidyverse, ggtree, treeio, ape, phytools
- RColorBrewer, ggtext

**External:**
- ADMIXTURE v1.3
- IQ-TREE2
- vcftools (preprocessing)

## Data Sources

Genomic data: 671 *L. infantum* isolates from global sampling (Jeffares lab, University of York, unpublished).

## Project Timeline

- **Oct 2024:** Initial phylogenetic analyses (MSc thesis)
- **Nov 2025:** ADMIXTURE population structure completed
- **Dec 2025:** GitHub portfolio publication (this repository)
- **Q1 2026:** Integration of HWE/geographical analyses
- **2026:** Manuscript submission (collaborative, multi-author)

## Author

**Mathew Johansson**  
MSc Bioinformatics, University of York  
Voluntary Research Associate, Jeffares Lab  
GitHub: [@MathewJohansson](https://github.com/MathewJohansson)

## Acknowledgments

- Dr. Daniel Jeffares (University of York) - Supervision, data generation
- Zeynep [Surname] - Geographical and HWE analyses collaboration

## License

Code: MIT License (upon manuscript publication)  
Data: Available upon reasonable request pending publication

---

**Note:** This is an active research project. Analyses are being refined for manuscript preparation. The repository demonstrates bioinformatics workflow, reproducible research practices, and clinical/evolutionary genomics expertise suitable for junior bioinformatics positions.
