## Test environments
* local R installation, R 4.3.1
* local ubuntu 20.04.6
* https://github.com/doccstat/fastcpd/blob/main/.github/workflows/check-standard.yaml

## R CMD check results

❯ checking CRAN incoming feasibility ... [16s] NOTE
  Maintainer: 'Xingchi Li <anthony.li@stat.tamu.edu>'

  New submission

  Possibly misspelled words in DESCRIPTION:
    Xianyang (16:18)
    Zhang (16:27)

❯ checking examples ... [39s] NOTE
  Examples with CPU (user + system) or elapsed time > 5s
           user system elapsed
  fastcpd 32.17   2.06   35.08

❯ checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    ''NULL''

❯ checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'

0 errors ✔ | 0 warnings ✔ | 4 notes ✖
