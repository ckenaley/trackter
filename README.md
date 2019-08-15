
[![Build
Status](https://travis-ci.com/ckenaley/trackter.svg?branch=master)](https://travis-ci.com/ckenaley/trackter)

# *trackter*

## Description

*trackter* is an R package for semiautomated tracking and analysis of 2D
kinematics from video and image data. The core functions of *trackter*
automatically detect a region of interest (ROI) and compute important
kinematic and shape parameters based on the ROI’s contour. Below, some
of *trackter*’s functionality is detailed and examples of its usage in
the context of fish locomotion are presented.

Please report any bugs or performance issues—this page is currently
under
    development.

## Installation

``` r
    require(devtools)
```

    ## Loading required package: devtools

``` r
    install_github("ckenaley/trackter")
```

    ## Skipping install of 'trackter' from a github remote, the SHA1 (d1472697) has not changed since last install.
    ##   Use `force = TRUE` to force installation

``` r
    require(trackter)
```

    ## Loading required package: trackter

## External dependencies

The core functions of *trackter* that extract shape and contour data
from images ( `kin.simple` and `kin.search`) depend upon *EBImage* It
can be installed easily with just a few lines of code:

``` r
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("EBImage")
```

    ## Bioconductor version 3.9 (BiocManager 1.30.4), R 3.6.0 (2019-04-26)

    ## Installing package(s) 'EBImage'

    ## Update old packages: 'boot', 'clipr', 'cluster', 'devtools', 'foreign',
    ##   'git2r', 'httr', 'nlme', 'openssl', 'pkgbuild', 'RcppArmadillo',
    ##   'remotes', 'testthat', 'tinytex', 'usethis', 'xml2'

*trackter* also contains several functions for image and video
processing. These functions depend on the popular `FFmpeg` package.
Installation is platform-dependent. I found the [`FFmpeg` wiki
installation and compilation guide to be quite
useful](https://trac.ffmpeg.org/wiki/CompilationGuide).

## Features

## Examples

Load *trackter*

``` r
#library(trackter)
#download example avi video
f <- "https://github.com/ckenaley/exampledata/blob/master/Pseed5_BCF_exp.avi?raw=true"
download.file(f,"sunfish_test.avi")

print(file.exists("Pseed5_BCF_exp.avi"))
```

Next, a filter is designated that rescales the images extracted from the
video with `video.to.images2`which is then called to perform the
extraction. This filter rescales the images to 600 pixels wide,
maintaining the aspect ratio. %Users may want to consult `FFmpeg`
documentation for filtering video data at
<https://ffmpeg.org/ffmpeg-filters.html>.

``` r
#extract images, reduce to 600 px wide
# filter
filt.red <- " -vf scale=600:-1 " 
vid.to.images2(vid.path="sunfish_test.avi",filt = filt.red)
```
