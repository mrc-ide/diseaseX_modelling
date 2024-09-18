# diseaseX_modelling

🦠🌍😷Modelling the deployment of a broad-spectrum vaccine against a
hypothetical "Disease X", caused by novel sarbecovirus with pandemic
potential 🦠🌍😷

## Overview

This repository contains the code used to evaluate and explore the
potential impact of a broadly protective sarbecovirus vaccine (BPSV)
during a hypothetical future pandemic caused by a novel sarbecovirus
pathogen ("SARS-X"). This work is not yet published, but has been
pre-printed and is available on medRxiv:
<https://www.medrxiv.org/content/10.1101/2024.08.12.24311730>

## Repo Contents

-   [analysis](./analysis): Code running the analyses and generating the
    figures featured in the paper.
-   [data](./data): Contains publicly available COVID-19 case and
    vaccination rate data from Our World in Data.
-   [figures](./figures): Containing .PDF and .ai versions of paper
    figures.
-   [functions](./functions): Custom functions required for the analyses
    presented in the paper.
-   [outputs](./outputs): Containing .rds outputs from the analyses.

## Software Requirements & Installation Guide

The results in the associated preprint were generated using R version
`4.4.0` "Puppy Cup". Other than the required R packages (specified in
each script and [main.R](./main.R)), running the code contained in this
repository requires that you clone and install the R package available
at <https://github.com/cwhittaker1000/squire.page.sarsX>. This
repository contains the codebase implementing the dynamical transmission
model used in this study, wrapped into a custom R package. In order to
this install this package, users are advised to clone this repository,
navigate to the directory in an R session and then install the package
using `devtools::install(".")` (installation typically takes order of
single-digit minutes) Further details of this mathematical modelling
framework are described in
[https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(22)00320-6/fulltext](https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(22)00320-6/fulltext){.uri}.

## Instructions for Use

Each of the files contained in [analysis](./analysis) implements a set
of analyses comprising a figure from the paper - assuming relevant R
packages are installed as described above and the working directory is
set to the home directory of this repo, users will be able to run these
scripts and reproduce the outputs. Exact runtime varies depending on the
exact set of analyses and number of associated simulations being carried
out. Single simulations from the branching process based framework
typically take order of 20 seconds for a population of 10,000
individuals; single simulations from the compartmental transmission
modelling framework take order of 1 second.

Any issues, please post an issue on this Github repository or feel free
to reach out at
[charles.whittaker16\@imperial.ac.uk](mailto:charles.whittaker16@imperial.ac.uk){.email}
😄
