
# Notes for build

  - use precompiled .Rmd for vignettes ([a
    la](https://ropensci.org/blog/2019/12/08/precompute-vignettes/))
    because of build fails on win-builder.
      - knit tutorial with `knitr::knit("vignettes/tutorial.Rmd.orig",
        output = "vignettes/tutorial.Rmd")`
      - add \*.Rmd.orig to Rbuildignore
      - update and reknit on new submission/commit/release 
      -set knit to
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

# Notes for testing/CI

  - follow these pointers for `pkgdown` site build with CI on Travis
  - check travis.yaml with <https://config.travis-ci.com/explore>
  - imagemagick requirements related to ubuntu/linux
    <https://srvanderplas.netlify.app/post/2019-05-07-travis-tricks/>

# Notes for gitting/site

  - `usethis::use_pkgdown()` for `pkgdown`
  - `pkgdown::build_articles` and `pkgdown::build_site`
  - see:
    <https://www.datacamp.com/community/tutorials/cd-package-docs-pkgdown-travis>
  - see:
    <https://cran.r-project.org/web/packages/pkgdown/vignettes/pkgdown.html>
  - abandoned deploy, now going with static, but often builds of site
    locally with push to git [a
    la](https://sahirbhatnagar.com/blog/2020/03/03/creating-a-website-for-your-r-package/)
      - `use_travis_deploy()` broken with:

<!-- end list -->

    Error: `github_token()` was deprecated in usethis 2.0.0 and is now defunct.
    Call `gh::gh_token()` to retrieve a GitHub personal access token
    Call `gh_token_help()` if you need help getting or configuring your token
