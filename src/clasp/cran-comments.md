## Test environments

* local OS X install, R 3.1.2 (via devtools)
* windows, R unstable (via win-builder)
* tried to build on two different linux systems where I don't have admin access; cannot install dependency Rmpfr due to out of date libmpfr on both systems

## R CMD check results

There were no ERRORs. There was one WARNING about no code that exercises the package. There is indeed a vignette on in the github repo, but it takes too long for it to run for it to be included in the package.

There were two NOTEs:

>Maintainer: ‘John Didion <john.didion@nih.gov>’
>New submission
>Non-FOSS package license (CC BY-NC-SA 4.0)

and

>No repository set, so cyclic dependency check skipped

The latter appears to be due to a bug in devtools: https://github.com/hadley/devtools/issues/602