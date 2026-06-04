# LEAP

**L**UG-**E**ncoded **A**ncestral **P**rojection — a phylogeny-based computational framework for inferring the transcriptional states of unobserved internal nodes from cell phylogenies via lineage-specific upregulated genes (LUGs).

## Overview

LEAP integrates cell lineage tracing with single-cell RNA sequencing (scRNA-seq) to reconstruct the transcriptomes of ancestral internal nodes in a developmental cell lineage tree. It identifies a set of lineage-committed upregulated genes (LUGs) from terminal cell transcriptomes, maps them back onto internal nodes of the phylogeny, and enables lineage-informed longitudinal analysis of cell state dynamics throughout development.

The method facilitates the discovery of subtle fate specializations and hidden developmental states that lack prominent transcriptional differentiation, supporting the construction of lineage-resolved cell atlases in complex organisms.

<img width="1924" height="849" alt="image" src="https://github.com/user-attachments/assets/f922c00b-b1d6-4bf4-b22e-2dd76ca3b3a2" />

## Installation

```r
devtools::install_github("shadowdeng1994/LEAP")
```

## Dependencies

- R (>= 4.0.0)
- [TarCA.beta](https://github.com/shadowdeng1994/TarCA.beta)
- tidyverse
- parallel
- castor
- Seurat

## Workflow

### 1. Load packages

```r
library(LEAP)
library(TarCA.beta)
library(parallel)
library(Seurat)
```

### 2. Load demo data

Demo data can be downloaded from [LEAP_sourcedata](https://github.com/shadowdeng1994/LEAP_sourcedata). It includes:

- **Query data**: a reconstructed cell phylogeny (`Var.Tree`), tip annotations (`Var.Ann`), and a gene expression matrix (`Var.Expression`)
- **Reference data**: a reference expression matrix (`Data.Ref_Expression`) and metadata (`Data.Ref_metadata`) for label transfer

```r
load("DemoData/Query.RData")
load("DemoData/Reference.RData")
```

### 3. Create TarCA object

Convert the phylogenetic tree and tip annotations into a TarCA object, which serves as the central data structure for all downstream analyses.

```r
Data.TarCA <- CreateTarCAObject(Var.Tree, Var.Ann)
```

### 4. Reconstruct internal node transcriptome

This is the core step of LEAP. It performs three sub-steps automatically:

1. **LUG detection** — identifies lineage-specific upregulated genes using a permutation test (1000 shuffles by default)
2. **Phylogenetic depth conversion** — maps LUG activation events onto the phylogenetic depth of internal nodes
3. **Internal state inference** — reconstructs the binarized transcriptome of each internal node

```r
Var.InternalMatrix <- GetInternalMatrix(
  Data.TarCA,
  Var.Expression,
  Threads = 90,
  SaveObject = TRUE
)
```

### 5. Binarize reference data

Convert the reference scRNA-seq dataset into a binarized matrix using the same LUG genes, grouped by developmental time bins or other metadata categories.

```r
tmp.overlapLUG <- intersect(rownames(Var.InternalMatrix), rownames(Data.Ref_Expression))

Var.ReferenceMatrix <- GetReferenceMatrix(
  Data.Ref_Expression[tmp.overlapLUG, ],
  Data.Ref_metadata,
  GroupBy = "embryo.time.bin",
  Threads = 90,
  SaveObject = TRUE
)
```

### 6. Define transfer labels

Create combined labels from reference metadata (e.g., time bin and lineage) for transferring to internal nodes.

```r
tmp.Label <- setNames(
  paste0(Data.Ref_metadata$embryo.time.bin, "|", Data.Ref_metadata$lineage),
  rownames(Data.Ref_metadata)
)
tmp.Label <- tmp.Label[colnames(Var.ReferenceMatrix)]
```

### 7. Run label transfer

Use Seurat's label transfer framework to assign developmental states to the reconstructed internal nodes.

```r
Var.Predictions <- RunLabelTransfer(
  Var.InternalMatrix,
  Var.ReferenceMatrix,
  tmp.Label,
  SaveObject = TRUE
)
```

### 8. Extract predicted states

Map the transferred labels onto the cell phylogeny, propagating annotations to unsolved descendant nodes.

```r
Var.PredictedState <- GetPredictedState(
  Data.TarCA,
  Var.Predictions,
  RemoveRoot = "Root",
  SaveObject = TRUE
)
```

## Key Functions

| Function | Description |
|---|---|
| `CreateTarCAObject` | Create a TarCA object from a phylogenetic tree and tip annotations |
| `GetInternalMatrix` | Reconstruct LUG-transcriptome of internal nodes (main workflow) |
| `GetReferenceMatrix` | Binarize reference dataset by expression percentiles within groups |
| `RunLabelTransfer` | Transfer labels from reference to query data via Seurat |
| `GetPredictedState` | Map inferred states to nodes and their descendants |
| `fun.ExtractLUG` | Infer LUGs based on permutation test |
| `fun.CalculatePotency` | Calculate Shannon entropy-based potency for each internal node |
| `fun.AssignLaterState` | Propagate predicted states through the internal node hierarchy |

## Output

All results are saved to the `LEAPOutput/` directory when `SaveObject = TRUE`:

- `Var.LUGs.rds` — detected lineage-specific upregulated genes with empirical P-values
- `Var.PhyloDepthOfLUGs.rds` — phylogenetic depth of LUG activation events
- `Var.InternalState.rds` — reconstructed binarized transcriptome of internal nodes
- `Var.RefMatrix.rds` — binarized reference matrix
- `Var.QuerySeurat.rds` / `Var.RefSeurat.rds` — Seurat objects for query and reference
- `Var.Anchors.rds` — transfer anchors
- `Var.Predictions.rds` — label transfer predictions
- `Var.PredictedState.rds` — final predicted states mapped onto the phylogeny

## License

MIT
