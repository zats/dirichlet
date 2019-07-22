# A new statistical method to analyze Morris Water Maze data using Dirichlet distribution

This is a python package forked from Eric Sue's [dirichlet](https://github.com/ericsuh/dirichlet) package which adds a test of uniformity for a single data set (*i.e.* equal fractions in each category) and plotting capability.

This was intended for the particular case of Morris Water Maze data as presented in Maugard and Doux, 2018. The test includes an approximate Bartlett-type correction for small samples as it is based on a likelihood-ratio test, which is only accurate for large samples.

## Use 
The package contains a python module called `dirichlet` which can be used with python online or locally, or with R (see below).

### With python

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/xuod/dirichlet/master)

We provide a jupyter notebook [example.ipynb](https://nbviewer.jupyter.org/github/xuod/dirichlet/blob/master/example/example.ipynb) to show how to use the `dirichlet` module to perform the uniformity test and produce plots with test data and/or your data.

You can run it with python either:
- online on [Binder](https://mybinder.org/v2/gh/xuod/dirichlet/master) (no installation required, you can even edit or copy/paste your data);
- locally by downloading the repository as a zip file and unzipping it. The jupyter notebook [example.ipynb](https://nbviewer.jupyter.org/github/xuod/dirichlet/blob/master/example/example.ipynb) in the `example` subfolder with test data and plots. If you have never used jupyter, see the [jupyter website](http://jupyter.org/).

<img width="20%" height="20%" src="https://github.com/xuod/dirichlet/blob/master/example/3Tg.png"> <img width="20%" height="20%" src="https://github.com/xuod/dirichlet/blob/master/example/wt.png">

### With R

You can import the python `dirichlet` module using the `reticulate` package for R. Here's an example from the R console:
```
library(reticulate)
setwd('/path/to/dirichlet/example')
data_3tg<-data.matrix(read.csv('3Tg.csv'))
dirichlet <- import_from_path("dirichlet", path="../")
dirichlet$test_uniform(data_3tg, label='3Tg', do_MWM_correction=TRUE, verbose=TRUE)
```
The data should be a 2D matrix with samples in rows. The `plot` function can also be used but the plot is not as good.

