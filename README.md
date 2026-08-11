# CV3D

CV3D is an R package for extracting and analysing compound-eye surface geometry from three-dimensional data. It implements the R-based analysis components of the CV3D workflow.

## What CV3D does

CV3D provides functions to:

- import triangle centres and surface normals from ASCII STL meshes;
- calculate local surface heights;
- condense thresholded surface-height points into facet candidates;
- identify neighbouring facets and estimate facet size;
- estimate facet surface normals;
- calculate local inter-facet angles, eye parameter, sampling frequency, and anatomical acuity estimates;
- align eye data to anatomical landmarks;
- create standardised face-on 2D or 3D eye views for visualisation and QC; and
- project facet viewing directions onto a spherical field of view.


## Spatial units

CV3D currently assumes that all input mesh and point-cloud coordinates are expressed in micrometres (µm). Facet-size estimates are therefore reported in µm, and Snyder's eye parameter in µm·rad.

## CV3D workflow

The R package can be used directly in R, but it is designed to work together with the [CV3D UI](https://github.com/Pete-s-Lab/CV3D-UI), which guides the complete workflow from 3D image data to analysis-ready compound-eye measurements.

## Installation

```r
# install.packages("remotes")
remotes::install_github("Pete-s-Lab/CV3D")
```

## Basic example

```r
data(cv3d_example_facets)

facets <- find_neighbours(cv3d_example_facets)
facet_sizes <- calculate_facet_size(facets)

head(facet_sizes)
```

For standardised face-on QC and comparison plots, `view_eye_face_on()` can display centred eye clouds in base R or, when installed, `rgl`. Added comparison clouds can be positioned in units of the reference eye's width or height, for example `translate_x = 1.1`.

## Documentation

Individual functions are documented in R:

```r
?CV3D
?find_neighbours
?calculate_facet_size
?view_eye_face_on
```

The package vignette provides a worked example of the R workflow.

## Citation

A formal citation for CV3D will be added following publication of the associated methods paper.

## Issues

Please report bugs and feature requests through the [GitHub issue tracker](https://github.com/Pete-s-Lab/CV3D/issues).
