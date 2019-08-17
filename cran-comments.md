## Test environments
* local OS X install, R 3.6.0
* OS X (on travis-ci) R 3.6.1
* win-builder (on appveyor-ci, devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking dependencies in R code ... NOTE
  Namespace in Imports field not imported from: 'R6'

  R6 is a build-time dependency.

## Downstream dependencies