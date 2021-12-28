
<img src="figs/Fig1.png" alt="mcfly logo" width="200px" align="right"/>

# There and Back to the Present: An Evolutionary Tale on Biological Diversity

[![DOI:10.1101/2021.12.11.472171](http://img.shields.io/badge/DOI-10.1101/2021.12.11.472171-B31B1B.svg)](https://www.biorxiv.org/content/10.1101/2021.12.11.472171v1)

## Authors

**Leandro Duarte**, Gabriel Nakamura, Vanderlei Debastiani, Renan
Maestri, Maria João Veloso da Costa Ramos Pereira, Marcus Cianciaruso,
José Alexandre F. Diniz-Filho

## Repository overview

This repository contains data and R scripts that allows to fully
reproduce the analyses presented in the manuscript entitled “There and
Back to the Present: An Evolutionary Tale on Biological Diversity”, in
which we introduce an individual-based simulation approach coupled with
Approximate Bayesian Computation (ABC) that allows to parameterize
adaptation rates of species niche positions along the evolution of a
monophyletic lineage, and the intensity of dispersal limitation between
local assemblages potentially connected by dispersal (metacommunity).
The analytic tool was implemented in an R package called
[mcfly](https://github.com/GabrielNakamura/mcfly). A preprint version of
this work can be found
[here](https://www.biorxiv.org/content/10.1101/2021.12.11.472171v1).

If you want to use this repository you can clone or download to your
computer machine using the following code:

To clone using github from command line:

``` bash
git clone https://github.com/GabrielNakamura/Data_MS-mcfly your_repo_title
```

To download using R:

``` r
download.file(url = "https://github.com/GabrielNakamura/Data_MS-mcfly/archive/main.zip", destfile = "Data_MS-mcfly.zip")
```

## Repository layout

The script and data are organized in folders containing files that are
described in the next section. Names of R scripts present the following
structure: **Number\_UpperCaseLetter\_ScriptName.R**

  - The number indicates the sequence in which the script must be
    opened, therefore, the script that innitiate with 01 are the first
    scripts that has to be opened;
  - The upper case letter indicates what the script does. C indicate
    scripts with codes to read libraries and data; scripts with D
    contain code for running the analysis;
  - The last part is the name of the script, that provide a glimpse of
    what the code contained in script aim to do.

For example: script 01\_C\_dataPhylostomidae.R must be the first script
to be opened since it reads data and libraries related to empirical
communities (Phylostomidae).

### `R` folder

The `R` folder contains all the scripts needed to read data and
libraries and running the analysis contained in the manuscript.

### `Data` folder

The `Data` folder contains all data necessary to run the simulations
(sensibility analysis) and empirical analysis (with Phylostomidae
communities)

### `Output` folder

This folder contains data from all analysis (sensibility and empirical).
Those files can be used to produce the Figures present in the
manuscript.

### `figs` folder

This folder contains all Figures presented in the manuscript. Figure
names corresponds are the same as the latest [preprint version of this
manuscript](https://www.biorxiv.org/content/10.1101/2021.12.11.472171v2).
All the Figures were made using
[DataGraph](https://apps.apple.com/us/app/datagraph/id407412840?mt=12)
software, however all data needed to construct the same figures can be
find in `Data` and `Outpub` folders.

## Additional software

As an additional tool we developed the R package `mcfly`. The latest
development version of the package can be downloaded using:

``` r
devtools::install_github(repo = "GabrielNakamura/mcfly", ref = "master")
```

**Important**: All analyses described in the preprint was done using
previous versions of `mcfly`

Additionaly to R,
[DataGraph](https://apps.apple.com/us/app/datagraph/id407412840?mt=12)
software was used to produce all the Figures present in the manuscript.
Files used in additional software can be found in `DataGraph` folder.

## Getting help

If you find a problem or a bug, please open an [Issue at
GitHub](https://github.com/GabrielNakamura/Data_MS-mcfly/issues) with a
minimal reproducible example.
