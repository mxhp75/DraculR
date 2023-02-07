# DraculR
## A public repository for the DraculR Shiny Application.

DraculR is a new lightweight tool for haemolysis detection in human plasma miR-Seq datasets.

This repository contains all code and data required to run DraculR, including an example dataset taken from our previous publication.[^fn1]

The test dataset can be found <a href="https://github.com/mxhp75/DraculR/tree/main/dataExample">here</a>. If you choose to test DraculR using the example data, please remember to select "10" as the number of your smallest group for filtering, and select the pregnancy associated miRNA:

- hsa-miR-186-5p
- hsa-miR-425-5p
- hsa-miR-25-3p
- hsa-miR-363-3p
- hsa-miR-183-5p
- hsa-miR-451a
- hsa-miR-182-5p
- hsa-miR-191-5p
- hsa-miR-194-5p
- hsa-miR-20b-5p

for removal from the calculation before calculating the Haemolysis Metric for these data.

You may access the DraculR Shiny.io page <a href="https://mxhp75.shinyapps.io/DraculR">here</a>.

### How to clone this repository

If you plan on using DraculR regularly, we suggest you clone this repository and run locally (you will need RStudio for this).

The command `$ git clone` is used to make an exact copy of an existing Git repository on your local machine (or wherever you choose to save the cloned repository). Your new copy will be saved into a local directory with the same name and directory structure as the GitHub repository you are cloning. By cloning this repository you also create remote tracked branches, and checkout an initial branch locally. By default, Git clone will create a reference to the remote repository called _origin_.

To create a local copy of the DraculR tool simply open your terminal and type the following:

```bash
  $ cd folder/to/clone-into/
  $ git clone https://github.com/mxhp75/DraculR
```

Please ensure the following packages are installed on your local machine:

```r
library(dplyr)
library(plyr)
library(ggplot2)
library(patchwork)
library(shiny)
library(ggrepel)
library(scales)
library(tidyr)
library(magrittr)
library(reshape)
library(edgeR)
library(readr)
library(DT)
library(psych)
library(tools)
library(shinyhelper)
```

You can now open the `draculR.R` script and Run from RStudio.

[^fn1]: Smith, Melanie D., Shalem Y. Leemaqz, Tanja Jankovic-Karasoulos, Dale McAninch, Dylan McCullough, James Breen, Claire T. Roberts, and Katherine A. Pillman. 2022. “Haemolysis Detection in MicroRNA-Seq from Clinical Plasma Samples.” Genes 13 (7). https://doi.org/10.3390/genes13071288.
