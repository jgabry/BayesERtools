## Test environments

* local Mac OSX install, R 4.4.1
* Windows, macOS, Ubuntu 22.04.1 on Github actions (devel and release)
* win-builder (devel and release)


## R CMD check results

0 errors | 0 warnings | 1 note

Days since last update: 1

## Note

Fixed the CRAN error in the previous version:

- Update package dependency that caused error in unit tests
- Use Suggests packages conditionally in tests and vignettes 
