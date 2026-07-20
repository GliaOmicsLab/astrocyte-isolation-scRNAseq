**Welcome!**
---
These scripts guide you through the analysis of single-cell (scRNA-seq) data generated from astrocyte enrichment experiments using three different isolation protocols.
These data support the publication titled: **Tissue dissociation protocols dictate reactive astrocyte recovery in single-cell transcriptomics** (DOI: ).

The raw data are publicly available on [NCBI GEO](https://). The processed single-cell RNA-sequencing data, curated as an annotated Seurat object, are available on [Zenodo](https://zenodo.org/records/).


## Experimental Overview
Capturing dynamic, morphologically complex cell states such as reactive astrocytes via standard single-cell workflows is heavily restricted by tissue-dissociation biases. To address this challenge, we developed and benchmarked an optimized, 
low-shear tissue processing workflow (**WorthGentle**) (a step-by-step protocol is available on [Protocol.io](https://www.protocols.io/private/)).

The results were compared between protocols:
-	**WorthGentle** – an optimized in-house protocol utilizing papain-based tissue dissociation, wide-bore pipette homogenization, and low-speed (60 x g) ovomucoid separation to isolate intact cells. 
-	**Miltenyi** – a standardized commercial protocol utilizing automated papain-based dissociation, and high-speed (3000 x g)density-gradient centrifugation for debris removal. 
-	**WorthMech** – an in-house high-shear control variant identical to WorthGentle, but utilizing standard-bore pipette tips to isolate the specific impact of physical shear stress. 
 
### Samples - cortical tissue was collected from female and male mouse models across three experimental groups:
- **control** - C57Bl/6J mice
- **permanent middle cerebral artery occlusion** - C57Bl/6J mice with tissues harvestd 7 days post-stroke 
- **tauopathy** - transgenic mice expressing human P301S tau [Tg(Thy1-MAPT*P301S)2541Godt]
---

## Data processing 
To maintain a reproducible and modular workflow, the project directory is organized into three distinct functional layers:

- **`OF` — Functions** — scripts containing functions/color palettes/orders/labels of variables to maintain consistency across notebooks
- **`1DP` — Data Processing** —  markdown files dedicated to raw data ingestion, quality filtering, downstream analysis, (e.g., creating Seurat objects, DEG analysis...).
- **`2V` — Visualization**  — markdown files focused on generating plots (e.g., UMAPs, dot plots, violin plots...)
---

### Description of markdown scripts 

1. **1DP_00_raw_data_to_seurat.Rmd** - processing of raw data to Seurat object 
2. **1DP_01_all_cells.Rmd** - all isolated cells - quality filtering, normalization, SCTransform clustering, doublet exclusion, cell-type annotation; visalisation in **2V_01_all_cells.Rmd**
3. **1DP_02_astro_subset.Rmd** - astrocytes - normalization, SCTransform clustering, cell-type annotation; visalisation in **2V_02_astro_subset.Rmd** and **2V_06_technical_parameters.Rmd**
4. **1DP_03_astro_DEA_GO.Rmd** - astrocytes - targeted differential gene expression analysis of reactive astrocytes isolated via the WorthGentle versus Miltenyi protocols, followed by gene 
ontology biological process enrichment analysis of the significant protocol-specific genes; visualisation in **2V_03_astro_DEA_GO.Rmd**
5. **1DP_04_validation_experiment.Rmd** - raw 384-well plate data processing, and target population enumeration to quantify the recovery efficiency of 
*Actb*<sup>+</sup>/*Aldh1l1*<sup>+</sup>/*Gfap*<sup>+</sup> reactive astrocytes relative to the *Actb*<sup>+</sup>/*Aldh1l1*<sup>+</sup> pan-astrocytic baseline across the Miltenyi and WorthGentle protocols; visalisation in **2V_04_validation_experiment.Rmd**
6. **1DP_05_visium.Rmd** - spatial transcriptomics data visium - data downloaded from the public repository to track *Gfap* expression dynamics following permanent middle cerebral artery occlusion; visualisation in **2V_05_visium.Rmd**
7. **2V_07_public_data.Rmd** - visualisation of meta-analysis of astrocyte recovery yields across publicly available scRNA-seq datasets from mouse stroke models
---

### Figures in article 
This section details the workflow for generating figures in the publication.

**Fig. 1: Reactive astrocytes are selectively underrepresented in scRNA-seq datasets ** 
- **A** - BioRender
- **B-E, I** - 2V_02_astro_subset.Rmd
- **F-G** - 2V_05_visium.Rmd
- **H** - 2V_04_validation_experiment.Rmd
- **J** - 2V_07_public_data.Rmd

**Fig. 2: Protocol-dependent gene expression profiles** 
- **A-C** - 2V_03_astro_DEA_GO.Rmd
- **D** - microscopic data

**SFig. 1: Characterization of cells obtained in scRNA-seq dataset** 
- **A-C** - 2V_01_all_cells.Rmd

**SFig. 2: The proportion of reactive astrocytes in scRNA-seq datasets is affected by isolation protocols** 
- **A, C** - 2V_02_astro_subset.Rmd
- **B** - 2V_06_technical_parameters.Rmd

**SFig. 3: Gating strategy for cell sorting** 
- **A-B** - FlowJo - FACS data