---
title: "Automated Kinematic Analysis in R"
author: "Christopher Kenaley"
date: "1/24/2018"
output:
  html_document:
    toc: true
    toc_float:
      collapsed: false
      smooth_scroll: false
      theme: united
    highlight: tango
---
#Introduction

trackter extracts 2D kinemtics from video and data. Although trackter was developed for analysis of swimming fishes, this package should prove useful forthose interested in waveforms of any undulating or oscillating body in a pixel field. The intent is to find body midline positions and, from that, compute tail-beat frequency and amplitude. The basic premise here is that an image of back-lit body can be transformed into a binary image rather easily and those negative (i.e., black pixels) can be used to infer the position of the midline. The midline data are then plotted onto the original images from the video and those images are stiched together to produce an MPEG complete with undulating fish and midline data. So you go from . . .

<img src="trout1_64_2018-01-23-130029-0000001.jpg" alt="trout swimming", width = "20%"> $\rightarrow$ <img src="trout1_64_2018-01-23-130029-0000.avi_001_bin.jpg" alt="trout swimming", width = "20%"> $\rightarrow$ <img src="trout1_64_2018-01-23-130029-0000_001_data.jpg" alt="trout swimming", width = "20%">

 . . . and, with a little more effort, this . . . 
 <blockquote class="twitter-video" data-lang="en"><p lang="en" dir="ltr">RT <a href="https://twitter.com/stemdevevo?ref_src=twsrc%5Etfw">@stemdevevo</a>: RT <a href="https://twitter.com/kenaley?ref_src=twsrc%5Etfw">@kenaley</a>: Know why I love working in R? &#39;Cause in about an hour I learned how to do this.... <a href="https://twitter.com/hashtag/RStats?src=hash&amp;ref_src=twsrc%5Etfw">#RStats</a> <a href="https://twitter.com/hashtag/6bl?src=hash&amp;ref_src=twsrc%5Etfw">#6bl</a>/s <a href="https://t.co/1UQyhF259Q">pic.twitter.com/1UQyhF259Q</a></p>&mdash; Pranay Roy ☃️ (@pranayroy01) <a href="https://twitter.com/pranayroy01/status/956082489918881792?ref_src=twsrc%5Etfw">January 24, 2018</a></blockquote>
<script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script>

# What you need

The script assumes the following:

 1. you have an image sequence or video of a back-lit body or at least highly constrasted images or videos
 2. you have [ffmpeg](https://www.ffmpeg.org/) installed

The back lighting---in this example, accomplished with an [LED light box](http://a.co/hIbIxNq) placed under the swim tunnel---creates some pretty good contrast. However, the scrips in trackter have been tested on some dark fish swimming on light back grounds and the results are comparable.

If you have an image sequence, you don't need ffmpeg installed. To stich outputted images (those with overlayed mideline track) you could make an AVI file from an image stack with some simple opertions in something like [imageJ](https://imagej.nih.gov/ij/)). By using trackter's \code[kin.img], near total automation from image extracted from the original vidoe to stitching overlayed images together after the binary images have been worked up.


# Install and load

Let's use `pacman` to install the required packages that are available through the normal mean (i.e. CRAN), that's if you don't have these already.

```{r, echo=F}
if (!require("pacman")) install.packages("pacman")
pacman::p_load(jpeg, png, ggplot2, grid,plyr,data.table,BiocInstaller)
```

Now we'll also need to install the [EBImage package](https://github.com/aoles/EBImage), on which the script depends for thresholding and object identification. That can be done with `biocLite` from `BiocInstaller`.

```{r}
#source("http://bioconductor.org/biocLite.R")
#biocLite("EBImage")
library("EBImage")

```

# Preliminaries
Let's set the working directory. 

```{r}
setwd("~/Google Drive/trout")
```

Now let's create a few important functions. First, base R doesn't contain a simple standard error function, odly enough, so let's create one: 


```{r}
se <- function(x){sd(x,na.rm = T)/sqrt(length(x))}

```

Because once  goal of this project was to determine tail beat amplitude, we need a point of reference. For this we'll compute a midline determined by a prediction from an linear model of the some anterior part of the body. That is, we assume the head is stiff and thus is will serve as the point of reference to determine the amplitude. The midline is projected along the body length and, once we have the tail position, we can determine the distance between the tail and that line. To make this easier, we  need a function that determines the distance between an oscilating point (i.e., the tip of the tail) and the midline so that we can determin the amplitude.  This requires just a little geometry: 

```{r}
#function to find distance between a point and a line
line.point.2d <- function(pt1,pt2,pt0){
  ((pt2[2]-pt1[2])*pt0[1]-(pt2[1]-pt1[1])*pt0[2]+pt2[1]*pt1[2]-pt2[2]*pt1[1])/sqrt((pt2[2]-pt1[2])^2+(pt2[1]-pt1[1])^2)}
```

To find amplitude and the frequency, we find the position of peaks within a series of position data. Here's a function to do just that inspired by [this exchange](https://stackoverflow.com/questions/34205515/finding-local-maxima-and-minima-in-r).

```{r}
find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}
```

#Getting down to business
Now that we have all this out of the way we can load and analyze a video. What's presented below if written in a way that could easily be wrapped in a function that could be applied to many videos in directory.

First, identify the video you'd like to work up and drop the "avi" extension. You can download the video [from this link] if you'd like.

```{r}
trial <- "trout1_32_test.avi"
trial <- gsub(".avi","",trial)
```

Now, store the working directory info so it's easier to pass that to `ffmpeg`. For this, we'll need to get rid of any spaces in the directory name, as is the case for Google drive. We'll also store where the video images should be placed.

```{r}
curdir <- gsub("Google Drive","\"Google Drive\"",getwd())
image.dir <- paste0(curdir,"/images/")
```

Next, we'll extract the images from AVI file and store them in the images directory using `ffmpeg` with a system call. The quality is reduced by half to speed processing. See the [`ffmpeg` documentation]( https://trac.ffmpeg.org/wiki/Encode/MPEG-4) for more info.

```{r}
  system(paste0("ffmpeg -i ", curdir,"/",trial, ".avi -q:v 15 ", image.dir,"/",trial,"%3d.jpg")) #see https://trac.ffmpeg.org/wiki/Encode/MPEG-4
```
Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
