coveffectsplot
========
![CRAN](http://www.r-pkg.org/badges/version-last-release/coveffectsplot)
[![Travis-CI Build Status](https://travis-ci.org/smouksassi/coveffectsplot.svg?branch=master)](https://travis-ci.org/smouksassi/coveffectsplot)
![DOWNLOADS](http://cranlogs.r-pkg.org/badges/grand-total/coveffectsplot)

A function and a Shiny App that Produce Forest Plots to Visualize Covariate Effects as commonly used in pharmacometrics population PK/PD reports.

### Installation and Running information
```
# Install from CRAN:
install.packages("coveffectsplot")

# Or the development version from GitHub:
#install.packages("devtools")
devtools::install_github('smouksassi/coveffectsplot')
coveffectsplot::run_interactiveforestplot()
```
### Example
Several example data are provided. if you want to bring your own data it should have at a minimum the following columns with the exact names:

paramname: Parameter on which the effects are shown e.g. CL, Cmax, AUC etc.

covname: Covariate name that the effects belong to e.g. Weight, SEX, Dose etc.

label: Covariate value that the effects of which is shown e.g. 50 kg, 50 kg\90 kg (here the reference value is contained in the label)

mid: Middle value for the effects usually the median from the uncertainty distribution

lower: Lower value for the effects usually the 2.5% or 5% from the uncertainty distribution

upper: Upper value for the effects usually the 97.5% or 95% from the uncertainty distribution

You might also choose to have a covname with value All (or other appropirate value) to illustrate and show the uncertainty on the reference value in a separate facet. Additionally, you might  want to have a covname with value BSV to illustrate and show the the between subject variability (BSV) spread. The example data show where does 90 and 50% of the patients will be given the model BSV estimate for the shown paramname. The vignette will walk you through how to compute and build the required data that the shiny app or the function `forest_plot`expects. There is some data management steps that the app does automatically chossing to call the function in script will require you to take responsibilty to build the table LABEL and to control the ordering of the variables

When you load the app and press the button to load the example data you get:


![example plot with the included dataset](./inst/shiny/img/snapshotforest.png)
