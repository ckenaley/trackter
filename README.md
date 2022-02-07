
[![R-CMD-check](https://github.com/ckenaley/trackter/workflows/R-CMD-check/badge.svg)](https://github.com/ckenaley/trackter/actions)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/ckenaley/trackter?branch=master&svg=true)](https://ci.appveyor.com/project/ckenaley/trackter)
[![](https://www.r-pkg.org/badges/version/trackter?color=orange)](https://cran.r-project.org/package=trackter)[![Codecov
test
coverage](https://codecov.io/gh/ckenaley/trackter/branch/master/graph/badge.svg)](https://codecov.io/gh/ckenaley/trackter?branch=master)

<!-- README.md is generated from README.Rmd. Please edit that file -->

    ## Warning in register(): Can't find generic `scale_type` in package ggplot2 to
    ## register S3 method.

# *trackter*

## Description

*trackter* is an R package for semiautomated tracking and analysis of 2D
kinematics from video and image data. The core functions of *trackter*
automatically detect a region of interest (ROI) and compute important
kinematic and shape parameters based on the ROI’s contour. These
functions use thresholding and segmentation to identify ROIs and, thus,
moderately contrasted images are required. Presented below is some of
*trackter*’s functionality in the context of fish locomotion.

Please report any bugs or performance issues.

## Getting started

For tutorial and tips on usage, please check out *trackter*’s [github
page](https://ckenaley.github.io/trackter/)

## Installation

The release version of *trackter* can be installed with:

``` r
install.packages("trackter")
```

The current development version can be installed with:

``` r
    require(devtools)
    install_github("ckenaley/trackter")
    require(trackter)
```

## External dependencies

The core functions of *trackter* that extract shape and contour data
from images ( `kin.simple` and `kin.search`) depend upon *EBImage*,
available on the Bioconductor repository. The current build and
development versions of *trackter* install this dependency. If it does
not install, it can be done so easily with just a few lines of code:

``` r
  if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
   BiocManager::install("EBImage")
```

*trackter* also contains several functions for image and video
processing. These functions depend on the popular `FFmpeg` package and
it must be installed if the user intends to use them. Installation is
platform-dependent. I found the [`FFmpeg` wiki installation and
compilation guide to be quite
useful](https://trac.ffmpeg.org/wiki/CompilationGuide).

## Features

**Automated kinematic analysis**

-   Fast and accurate contour and shape analysis of ROIs.
-   ROI detection with search parameters including position and size.
-   Relevant functions: `kin.search`, `kin.simple`, `kin.free`

**Tools for kinematic analysis of swimming animals**

-   Calculate midline (propulsive) wavelength, trailing-edge frequency,
    paired-fin position.
-   Relevant functions: `amp.freq`, `halfwave`, `wave`, and `fin.kin`

\_\_Tools for image and video processing

-   Access `FFmpeg` functionality, including filters and codecs, to
    extract frames, stitch videos, and edit images and videos.

-   Relevant functions: `images.to.videos`, `vid.to.images`

**Other miscellaneous, low-level tools for kinematic analysis**

-   Compute distances in 2d space, angles, heading/bearing, convert
    radians to degrees and vice versa.
-   Relevant functions: `dist.2d`, `cosine.ang`, `bearing.xy`, `deg`,
    `rad`

## Bugs and feedback

Please report issues [here](https://github.com/ckenaley/trackter/issues)
