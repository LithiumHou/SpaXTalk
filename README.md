## About
![Figure](/1.jpg)
This directory contains the code and resources of the following paper: 
Screening communication crosstalk in space

## What is SpaXTalk?
SpaXTalk is a statistical-learning approach that use spatial transcriptomics (ST) data to:
- model the cell-cell communication and intracellular signal transduction at single-cell resolution with spatial constraints
- identify the regulatory pathways that are activated by CCC signals
- quantify the level of crosstalk that vary across the space

## How to install SpaXTalk?
R >= 4.3.1 is required to correctly install SpaXTalk and other dependencies. We strongly suggest that you use RStudio to install and run SpaXTalk!

To install the SpaXTalk R package, you may either install from remote or from local. The installation takes ~2 minutes on a normal PC
<details>
  <summary>OPTION 1: remote installation</summary>

Run the following command in R:

```R
devtools::install_github("LithiumHou/SpaXTalk", dependencies = T, upgrade = "always")
```

Note: using `devtools::install_github` in Rstudio sometimes causes a github's token issue. In this case, you may need to generate a token. Please see [here](https://usethis.r-lib.org/articles/git-credentials.html). Alternatively, you may try local installation (see below).

</details>

<details>
  <summary>OPTION 2: install from local</summary>
You may download or clone the SpaXTalk repository to your device and run:
  
```R
install.packages(package_file,repos = NULL,type = "source")   # Replace package file with the path where you store the SpaXTalk repository
```

</details>

## How to use SpaXTalk?

We provide a step-by-step tutorial to show the functionality of SpaXTalk [here](https://github.com/LithiumHou/SpaXTalk/blob/main/vignettes/demo.md).
