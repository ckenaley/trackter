## Test environments
* local OS X install, R 3.6.0
* OS X (on travis-ci, release)
* win-builder (on appveyor-ci, devel and release)
* win-builder (devel and release)

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs. 

## Downstream dependencies
There are currently no downstream dependencies for this package.

## Comments
This is the third submission of this package. Upon second submission, the following was noted:

* References describing the methods should be included in the description. ACTION: no references are included in descriptions. References to methods include only those to functions in dependencies and these are linked with 'link{function}'.
* Ensure that your functions do not write by default or in your 
examples/vignettes/tests in the user's home filespace. ACTION: Functions now require user input to specify where any output is written and 'tempdir()'
used throughout tests and examples.
* '\\dontrun' should be only used if the example really cannot be executed. ACTION: '\\dontrun' is reserved for examples of functions for which the installation of additional software (ffmpeg)  can't be assumed and for examples that run >5s. For the latter, shorter additional examples are now provided.
* Please always write TRUE and FALSE instead of T and F. ACTION: All instances of 'T' or 'F' spelled out.


win-builder throws one NOTE: 'Possibly mis-spelled words in DESCRIPTION:
  ROIs (5:69)
  midline (5:109)'. 
  
These are words common to the image analysis lexicon.