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
This is the second submission of this package. Upon first submission, errors related to either tests reliant on FFmpeg external dependencies or  'flush(stde)); flush(stdout()) Error: unexpected ')' in "flush(stde))" '. The former fixed by running test only if FFmpeg is installed; the latter by removing 'unlink()' from examples. FFmpeg has been added to SystemRequirements in the DESCRIPTION.

win-builder throws one NOTE: 'Possibly mis-spelled words in DESCRIPTION:
  ROIs (5:69)
  midline (5:109)'. 
  
These are words common to the image analysis lexicon.