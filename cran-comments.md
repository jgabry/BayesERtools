## Test environments

* local Mac OSX install, R 4.4.1
* Windows, macOS, Ubuntu 24.04.1 on Github actions (devel and release)
* win-builder (devel and release)


## R CMD check results

0 errors | 0 warnings | 1 note

Maintainer: 'Kenta Yoshida <yoshida.kenta.6@gmail.com>'

New submission

Possibly misspelled words in DESCRIPTION:
  Emax (16:62, 17:45)

## Note

Thank you for the review. CRAN cookbook was very helpful.

Fixed the issues raised by the CRAN team in the previous submission:

- Expanded DESCRIPTION
- Removed LICENSE file
- Added return value to the function documentation
- Removed or replaced \dontrun{}
