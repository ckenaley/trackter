---
output:
  md_document:
    variant: markdown_github
---

[![Build Status](https://travis-ci.com/ckenaley/trackter.svg?branch=master)](https://travis-ci.com/ckenaley/trackter)


```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "figure/",
  fig.height = 1
)
```
# *trackter*

## Description

*trackter* is an R package for semiautomated tracking and analysis of 2D kinematics from video and image data. The core functions of *trackter* automatically detect a region of interest (ROI) and compute important kinematic and shape parameters based on the ROI's contour. Below, some of *trackter*'s functionality is detailed and examples of its usage in the context of fish locomotion are presented.

Please report any bugs or performance issues---this page is currently
under development.

## Installation

```r{}
    require(devtools)
    install_github("ckenaley/trackter")
    require(trackter)  
```

## External dependencies

The core functions of *trackter* that extract shape and contour data from images ( `kin.simple` and `kin.search`) depend upon *EBImage* It can be installed easily with just a few lines of code: 

```r{}
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("EBImage")
```

*trackter* also contains several functions for image and video processing. These functions depend on the popular `FFmpeg` package.
Installation is platform-dependent. I found the [`FFmpeg` wiki installation and compilation guide to be quite useful](https://trac.ffmpeg.org/wiki/CompilationGuide).

## Workflow 

![](workflow.pdf)

The following generalized workflow is summarized above. Users can begin analyses with an image sequence or a video. Supported image formats include .jpg, .tif, and .png. If the starting point of analysis is a video file, `video.to.images` or `video.to.images2` to extract images to a ``processed_images" subdirectory. Supported video formats for `video.to.images` and `video.to.images2` include .avi and .mp4 files. `video.to.images2` supports added `FFmpeg` functionality, allowing users to specify filters. 

Two `kin` functions in *trackter* can be used to extract midline and position data from ROIs in each image of a sequence. `kin.simple` and `kin.search` create binary images and use imported functions from *EBImage* to threshold and segment each image. The user can specify a threshold value (0--1) with ``thresh" or determine the value with Otsu's method (``thresh='otsu'"). Each `kin` function allows users to filter the potential ROIs by specifying a minimum size (in \% of the 2D pixel field) with ``size.min" and whether potential ROIs are free of the images edges (``edges=TRUE"). `kin.simple` is more straightforward, finding and analyzing ROIs that pass size and edge criteria. 


For image fields that contain other large contrasted objects (e.g, a wall or structure intended to alter flow), there is `kin.search`. This function allows the user to specify one of four search criteria which find the ROI according to: 


	* A centroid that is the shortest linear distance to the center of the field (``search.for=offset'').
	* A centroid x position that is closest to the x position of the center of the field (``search.for=offset.x").
	* A centroid y position that is closest to the y position of the center of the field (``search.for=offset.y").
	* The largest ROI (``search.for=largest")


	Both `kin` functions return a list of data tables that include frame-specific position, midline, and contour data of detected ROIs. Data returned in the ``kin.dat" data table include frame number, x and y positions of the head and tail, and amplitude of the tail relative to a theoretical midline established by a linear prediction from the parameter ``ant.per" (the proportion of the body represented by a stiff head). Midline data returned in the ``midline" data table include x and y positions of the ROI midline, y position of the theoretical linear midline, and a smoothed y position of the midline. The midline of the body contour (``y.m") is calculated as the vertical midpoint between minimum and maximum y values at each x position along the body. ``y.m" is passed through a user-specified smoothing operation (either a LOESS or smooth-spline fit) to calculate the smoothed midline (``y.pred"). The theoretical linear midline is subtracted from smoothed midline values to produce "wave.y". 
	
	Users can pass the data returned by the \textit{kin} routines to a number of other functions for analysis. `amp.freq` calculates the trailing-edge frequency given amplitude values and a sample rate. Functions `wave` and `halfwave` compute wavelength and half-wavelength from a wave-like data table. `fin.kin` allows the user to specify an apriori fin position (relative position along the body) to extract left and right fin contours and amplitude data. `fin.kin` also returns a composite contour of the full-body minus the fins and a smoothed midline based on these data. 

Lastly, if the user specifies ``save=TRUE" in the `kin` routines, processed images, complete with predicted, smoothed midlines, are saved as .jpg files in the ``processed_images" subdirectory. Using `images.to.video` or `images.to.video2`, these images can be compiled in a video file. `images.to.video2` permits use of FFmpegs flexible filters and video codecs. 

