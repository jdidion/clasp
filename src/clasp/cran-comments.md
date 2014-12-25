## Test environments

* local OS X install, R 3.1.1
* windows, R unstable (via win-builder)
* tried to build on two different linux systems where I don't have admin access; cannot install dependency Rmpfr due to out of date libmpfr on both systems

## R CMD check results

There were no ERRORs, WARNINGs or NOTEs in my local build. win-builder had one warning about no code that exercises the package. There is indeed a vignette on in the github repo, but it takes too long for it to run for it to be included in the package.
