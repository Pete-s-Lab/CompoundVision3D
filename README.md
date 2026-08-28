# CV3D

**CV3D** is an R package for extracting and analysing three-dimensional compound-eye surface geometry. It provides the numerical methods used by the companion CV3D UI and can also be used directly from R.

- R package repository: https://github.com/Pete-s-Lab/CV3D
- Graphical workflow: https://github.com/Pete-s-Lab/CV3D-UI

## Main functionality

The package includes tools to:

- import and represent corneal surface meshes;
- calculate local surface heights and optional local-height normalization;
- condense thresholded surface points into facet-position candidates;
- identify eye-boundary facets and construct an edge-aware first-ring neighbour graph;
- estimate facet size from retained neighbour-centre distances;
- estimate facet surface normals, including the regularised facet-centre envelope method;
- calculate local inter-facet angles, Snyder eye parameter, sampling frequency, and anatomical acuity estimates;
- align point clouds to specimen-level anatomical landmarks;
- project facet viewing directions to a spherical field of view; and
- produce standardized eye visualizations and QC outputs.

The end-to-end Fiji/ImageJ + Blender + R workflow is managed by **CV3D UI**. The UI is therefore the recommended entry point for complete datasets, while this repository contains the R analysis package itself.

## Installation

Install the development version from GitHub with:

```r
remotes::install_github("Pete-s-Lab/CV3D", build_vignettes = TRUE)
```

or:

```r
devtools::install_github("Pete-s-Lab/CV3D", build_vignettes = TRUE)
```

Then load the package with:

```r
library(CV3D)
```

## Current geometric workflow

In the current CV3D UI workflow, manually checked facet centres are followed by a dedicated neighbour-selection step. `find_neighbours_edge_aware()` detects likely boundary facets from local angular coverage and constructs a conservative reciprocal neighbour graph. This stored graph is subsequently used for facet-size and optical calculations.

Facet surface normals can be estimated with `get_facet_normals_envelope()`, the current UI default. The method reconstructs a regularised triangular envelope directly from facet centres, independently of the optical neighbour graph, and derives a locally weighted surface normal for each original facet centre. The original `get_facet_normals()` estimator remains available for backward compatibility and method comparison.

See `vignettes/CV3D.Rmd` for a worked overview of the current package workflow and the relevant function documentation for algorithmic details.

## Units

CV3D currently assumes that spatial input coordinates are in **micrometres (µm)**. Consequently, length-derived outputs are reported in µm and the Snyder eye parameter is reported in µm·rad. Functions that operate on arbitrary geometric radii expect those radii in the same coordinate units as their input point clouds.

## Example

The package contains reduced example data for demonstrating the R workflow:

```r
library(CV3D)

data(cv3d_example_facets)

facets <- find_neighbours_edge_aware(
  cv3d_example_facets,
  edge_gap_threshold_deg = 90,
  verbose = FALSE
)

normals <- get_facet_normals_envelope(
  cv3d_example_facets,
  envelope_factor = 1.25,
  verbose = FALSE
)
```

The example values are intended to demonstrate the API; dataset-specific parameters should be selected and validated for the data being analysed.

## Documentation

Use R's help system for individual functions, for example:

```r
?find_neighbours_edge_aware
?get_facet_normals_envelope
?calculate_local_heights
?get_optic_properties
```

The package vignette gives a compact workflow overview:

```r
vignette("CV3D")
```

For the complete interactive workflow, file naming, QC steps, and export structure, see the tutorial in the CV3D-UI repository.

## Issues

Please report package problems at:

https://github.com/Pete-s-Lab/CV3D/issues
