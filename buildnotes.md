
# Notes for build

  - use precompiled .Rmd for vignettes ([a
    la](https://ropensci.org/blog/2019/12/08/precompute-vignettes/))
    because of build fails on win-builder.
      - knit tutorial with `knitr::knit("vignettes/tutorial.Rmd.orig",
        output = "vignettes/tutorial.Rmd")`
      - add \*.Rmd.orig to Rbuildignore
      - update and reknit on new submission/commit/release -set knit to
        `fig.path=""`
      - manually move figures from trackter root to vignette with code
        below
      - add r code for users with
        `knitr::purl("vignettes/tutorial.Rmd.orig", output =
        "vignettes/tutorial.R")`

<!-- end list -->

``` r
knitr::knit("vignettes/tutorial.Rmd.orig", output = "vignettes/tutorial.Rmd")

knitr::purl("vignettes/tutorial.Rmd.orig", output = "vignettes/tutorial.R")
```

``` r
f = list.files(full.names = TRUE,pattern = ".png")
f2 <- gsub("./","./vignettes/",f)
file.copy(f,f2)
sapply(f,unlink)
```

# Notes for testing